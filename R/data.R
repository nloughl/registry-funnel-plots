# Load extraction outputs and derive the inputs for funnel plots -------------------------------
#
# Inputs are the CSVs written by the regextract pipeline (regitry-outlier-detection repo):
#   UKA_long.csv          one row per device x time point
#   UKA_casemix_long.csv  one row per case-mix stratum x sex x time point

load_extraction <- function(data_dir, procedure = "UKA") {
  read <- function(f) readr::read_csv(file.path(data_dir, f), show_col_types = FALSE, guess_max = 1e5)
  # drop trademark signs (e.g. EPRD "SIGMA(TM)"): they break string handling in non-UTF-8 R sessions
  clean <- function(d) dplyr::mutate(d, dplyr::across(dplyr::where(is.character),
                                                      ~ gsub("\u2122", "", .x, fixed = TRUE, useBytes = TRUE)))
  read <- function(f) clean(readr::read_csv(file.path(data_dir, f), show_col_types = FALSE, guess_max = 1e5))
  list(
    device  = read(paste0(procedure, "_long.csv")),
    casemix = if (file.exists(file.path(data_dir, paste0(procedure, "_casemix_long.csv"))))
                read(paste0(procedure, "_casemix_long.csv")) else NULL
  )
}

has_estimate <- function(d) !is.na(d$estimate) & d$value_status %in% c("ok", "ci_not_estimable")

# Short, readable point labels ---------------------------------------------------------------
short_device_label <- function(d) {
  lab <- gsub("\\[[^]]*\\]", " ", d$device_label)           # NJR [M.Fem:M.Tib] tags
  lab <- trimws(gsub("\\s+", " ", lab))
  both <- !is.na(d$femoral) & !is.na(d$tibial)
  same <- both & tolower(d$femoral) == tolower(d$tibial)
  lab[both] <- paste(d$femoral[both], "/", d$tibial[both])
  lab[same] <- d$femoral[same]
  lab <- ifelse(!is.na(d$compartment) & d$compartment == "lateral", paste0(lab, " (lat)"), lab)
  # EPRD lists one design in several fixation/bearing sections: add a short tag
  sec <- tolower(d$section %||% NA_character_)
  tag <- paste0(ifelse(grepl("mobile", sec), "MB", ifelse(grepl("fixed", sec), "FB", "")),
                ifelse(grepl("uncemented", sec), " uncem", ifelse(grepl("hybrid", sec), " hyb",
                       ifelse(grepl("cemented", sec), " cem", ""))))
  ifelse(!is.na(sec) & nzchar(trimws(tag)), paste0(lab, " (", trimws(tag), ")"), lab)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Device families (for "one design across registries" plots) --------------------------------
load_families <- function(path) readr::read_csv(path, show_col_types = FALSE, comment = "#")

assign_family <- function(d, families) {
  fam <- rep(NA_character_, nrow(d))
  for (i in seq_len(nrow(families))) {
    hit <- is.na(fam) & grepl(families$pattern[i], d$device_label, ignore.case = TRUE)
    fam[hit] <- families$family[i]
  }
  dplyr::mutate(d, family = fam)
}

# Registry means ("standardising the means") -------------------------------------------------
#
# Each registry's reference rate at time t, used as the funnel centre. Priority:
#   1. manual override (config/registry_means_override.csv)
#   2. the registry's own summary row(s): total / average rows with an estimate
#      (several, e.g. LROI cemented + uncemented, are pooled weighted by n)
#   3. case-mix subtotals at t (AOANJRR sex subtotals), pooled weighted by n
#   4. fallback: n-weighted mean of the listed device rows (flagged - not a registry-reported value)

pool <- function(p, n) sum(p * n) / sum(n)

# CI of an n-weighted pool of published estimates (all in percent): each part's SE is backed out of its
# 95% CI width, SE_i = (ucl - lcl) / 3.92, and SE_pool = sqrt(sum((w_i * SE_i)^2)) with w_i = n_i / sum(n).
# Approximate (KM CIs are not symmetric, parts treated as independent) - the same spirit as the pooled mean.
pool_ci <- function(est, lcl, ucl, n) {
  if (any(is.na(c(lcl, ucl, n)))) return(c(NA_real_, NA_real_))
  if (length(est) == 1) return(c(lcl, ucl) / 100)
  w <- n / sum(n)
  se <- sqrt(sum((w * (ucl - lcl) / 3.92)^2))
  (pool(est, n) + c(-1, 1) * 1.96 * se) / 100
}

mean_row <- function(reg, p, n, source, ci = c(NA_real_, NA_real_), ci_source = NA_character_) {
  tibble::tibble(registry = reg, p_mean = p, p_lcl = ci[1], p_ucl = ci[2], n_total = n,
                 mean_source = source, ci_source = if (all(is.na(ci))) "none" else ci_source)
}

#' Registry means with their 95% CI where one can be had:
#'   override        lcl_pct / ucl_pct columns if filled
#'   summary row(s)  published CI; several rows (LROI cemented + uncemented) -> pool_ci()
#'   case-mix        AOANJRR male + female subtotal CIs -> pool_ci()
#'   fallback        none (EPRD: mean of listed devices, no published total)
#' `p_ref` (the funnel reference) is the mean; set_reference() switches it to a CI bound.
registry_means <- function(device, casemix = NULL, time_yr, overrides = NULL) {
  regs <- sort(unique(device$registry))
  out <- purrr::map_dfr(regs, function(reg) {
    dv <- dplyr::filter(device, registry == reg, time_yr == !!time_yr)
    listed_n <- sum(dv$n_total[dv$row_type %in% c("device", "other")], na.rm = TRUE)

    if (!is.null(overrides)) {
      ov <- dplyr::filter(overrides, registry == reg, time_yr == !!time_yr)
      if (nrow(ov)) {
        ci <- if (all(c("lcl_pct", "ucl_pct") %in% names(ov))) c(ov$lcl_pct[1], ov$ucl_pct[1]) / 100 else c(NA, NA)
        return(mean_row(reg, ov$mean_pct[1] / 100, ov$n_total[1], paste("override:", ov$source[1]),
                        as.numeric(ci), "override"))
      }
    }
    summ <- dplyr::filter(dv, row_type %in% c("total", "average"), !is.na(estimate))
    if (nrow(summ)) {
      n <- dplyr::coalesce(summ$n_total, if (nrow(summ) == 1) listed_n else NA_real_)
      return(mean_row(reg, pool(summ$estimate / 100, n), sum(n),
                      paste0("registry summary row: ", paste(summ$device_label, collapse = " + ")),
                      pool_ci(summ$estimate, summ$lcl, summ$ucl, n),
                      if (nrow(summ) == 1) "published" else "pooled from published CIs (approx.)"))
    }
    if (!is.null(casemix)) {
      cm <- dplyr::filter(casemix, registry == reg, time_yr == !!time_yr, row_type == "subtotal", !is.na(estimate))
      if (nrow(cm)) return(mean_row(reg, pool(cm$estimate / 100, cm$n_total), sum(cm$n_total),
                                    paste0("case-mix subtotals pooled (", cm$table_id[1], ")"),
                                    pool_ci(cm$estimate, cm$lcl, cm$ucl, cm$n_total),
                                    "pooled from published CIs (approx.)"))
    }
    dd <- dplyr::filter(dv, row_type == "device", !is.na(estimate), !is.na(n_total))
    mean_row(reg, pool(dd$estimate / 100, dd$n_total), sum(dd$n_total),
             "FALLBACK: n-weighted mean of listed devices")
  })
  # carry report metadata for tables/captions
  meta <- device |>
    dplyr::distinct(registry, report_year, metric_type) |>
    dplyr::group_by(registry) |>
    dplyr::summarise(report_year = paste(unique(report_year), collapse = "/"),
                     metric_type = paste(unique(metric_type), collapse = "/"), .groups = "drop")
  dplyr::left_join(out, meta, by = "registry") |> set_reference("mean")
}

REFERENCES <- c("Registry mean" = "mean", "Lower 95% CI of mean" = "lcl", "Upper 95% CI of mean" = "ucl")
reference_label <- function(ref) c(mean = "registry mean", lcl = "lower CI of registry mean",
                                   ucl = "upper CI of registry mean")[ref]

#' Choose the funnel reference: "mean", "lcl" or "ucl". Registries without a CI keep their mean
#' (reference_used says what each registry actually got).
set_reference <- function(means, reference = c("mean", "lcl", "ucl")) {
  reference <- match.arg(reference)
  bound <- switch(reference, mean = means$p_mean, lcl = means$p_lcl, ucl = means$p_ucl)
  dplyr::mutate(means, reference = .env$reference,
                p_ref = dplyr::coalesce(.env$bound, p_mean),
                reference_used = ifelse(is.na(.env$bound), "mean", .env$reference))
}

#' One row per registry x reference (mean and any CI bounds that exist), ready for build_limit_curves():
#' limits are centred on p_ref (the bound) and standardised against it too, so every funnel is symmetric about 0.
#' fill_missing = TRUE keeps registries without a CI (EPRD) in the CI rows, using their mean (as set_reference does).
reference_refs <- function(means, references = c("mean", "lcl", "ucl"), n_max, phi = NULL, fill_missing = FALSE) {
  purrr::map_dfr(references, function(r) {
    m <- set_reference(means, r)
    if (!fill_missing) m <- dplyr::filter(m, reference_used == .env$r)
    if (!is.null(phi)) m <- dplyr::left_join(m, phi, by = "registry")
    dplyr::transmute(m, registry, reference = .env$r, p_ref, p_centre = p_ref, n_max = n_max,
                     phi = if ("phi" %in% names(m)) dplyr::coalesce(phi, 1) else 1)
  })
}

# Device points ---------------------------------------------------------------------------------
device_points <- function(device, means, time_yr, families = NULL) {
  d <- device |>
    dplyr::filter(row_type == "device", time_yr == !!time_yr, !is.na(n_total)) |>
    dplyr::filter(has_estimate(dplyr::pick(dplyr::everything()))) |>
    dplyr::mutate(label = short_device_label(dplyr::pick(dplyr::everything())), p = estimate / 100) |>
    dplyr::left_join(dplyr::select(means, registry, p_ref, p_mean, reference_used), by = "registry") |>
    # standardised against the reference in use (registry mean, or a CI bound of it) - the same value
    # the funnel is centred on, so reference = 0 (difference) or 1 (ratio)
    dplyr::mutate(dev = p - p_ref, ratio = p / p_ref)
  if (!is.null(families)) d <- assign_family(d, families)
  d
}

# Case-mix strata ------------------------------------------------------------------------------
#
# AOANJRR KP22: sex x age.  NJR 3.K6: design block (unicondylar cemented / uncemented-hybrid)
# x sex x age; `njr_level = "group"` uses the design blocks, "stratum" the medial/lateral rows.
casemix_strata <- function(casemix, time_yr, njr_level = c("group", "stratum")) {
  njr_level <- match.arg(njr_level)
  casemix |>
    dplyr::filter(time_yr == !!time_yr, !is.na(n_total), !is.na(age_group), !is.na(sex)) |>
    dplyr::filter(has_estimate(dplyr::pick(dplyr::everything()))) |>
    dplyr::filter(dplyr::case_when(
      is.na(design_group) ~ row_type == "stratum",
      TRUE ~ row_type == njr_level)) |>
    dplyr::mutate(
      p = estimate / 100,
      design = dplyr::case_when(
        is.na(design_group) ~ "all UKA",
        njr_level == "group" ~ sub("^All unicondylar, ", "", design_group),
        TRUE ~ paste(sub("^All unicondylar, ", "", design_group), design_subgroup, sep = ": ")),
      age_group = factor(age_group, levels = c("<55", "55-64", "65-74", ">=75")),
      stratum = paste(design, sex, age_group, sep = " | ")
    )
}

# Procedure-period strata (LROI K052B) ---------------------------------------------------------------
#
# One row per 2-year period of primary UKA (2009-2010 ... 2021-2022) with a KM estimate at time_yr.
# NOTE the metric: LROI reports MAJOR revision (first revision of the femur or tibia) here, while its
# device tables report any revision - period rates run a little below the device-table rates.
period_strata <- function(casemix, time_yr) {
  if (is.null(casemix) || !"procedure_period" %in% names(casemix)) return(NULL)
  d <- casemix |>
    dplyr::filter(!is.na(procedure_period), time_yr == !!time_yr, !is.na(n_total)) |>
    dplyr::filter(has_estimate(dplyr::pick(dplyr::everything())))
  if (!nrow(d)) return(NULL)
  d |>
    dplyr::arrange(period_start) |>
    dplyr::mutate(p = estimate / 100,
                  period = factor(procedure_period, levels = unique(procedure_period)))
}

#' Sequential colours for procedure periods (oldest = dark, newest = light); registry colours stay reserved.
period_colours <- function(periods) {
  periods <- as.character(unique(periods))
  stats::setNames(grDevices::hcl.colors(length(periods), "viridis"), periods)
}

#' Periods whose funnel a point is outside of: for each device and period, "above" / "within" / "below"
#' the period's limit at the device's n (limits built around the period's own rate).
period_status <- function(points, periods, level = 0.998, method = "normal") {
  tidyr::crossing(dplyr::select(points, registry, label, n_total, p),
                  dplyr::select(periods, period, p_period = p)) |>
    dplyr::mutate(lim = purrr::map2(p_period, n_total, ~ funnel_limits(.x, .y, level, method)),
                  upper = purrr::map_dbl(lim, "upper"), lower = purrr::map_dbl(lim, "lower"),
                  status = dplyr::case_when(p > upper ~ "above", p < lower ~ "below", TRUE ~ "within")) |>
    dplyr::select(-lim)
}
