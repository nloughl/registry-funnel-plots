# Load extraction outputs and derive the inputs for funnel plots -------------------------------
#
# Inputs are the CSVs written by the regextract pipeline (regitry-outlier-detection repo):
#   UKA_long.csv          one row per device x time point
#   UKA_casemix_long.csv  one row per case-mix stratum x sex x time point

load_extraction <- function(data_dir, procedure = "UKA") {
  read <- function(f) readr::read_csv(file.path(data_dir, f), show_col_types = FALSE, guess_max = 1e5)
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

registry_means <- function(device, casemix = NULL, time_yr, overrides = NULL) {
  regs <- sort(unique(device$registry))
  out <- purrr::map_dfr(regs, function(reg) {
    dv <- dplyr::filter(device, registry == reg, time_yr == !!time_yr)
    listed_n <- sum(dv$n_total[dv$row_type %in% c("device", "other")], na.rm = TRUE)

    if (!is.null(overrides)) {
      ov <- dplyr::filter(overrides, registry == reg, time_yr == !!time_yr)
      if (nrow(ov)) return(tibble::tibble(registry = reg, p_ref = ov$mean_pct[1] / 100,
                                          n_total = ov$n_total[1], mean_source = paste("override:", ov$source[1])))
    }
    summ <- dplyr::filter(dv, row_type %in% c("total", "average"), !is.na(estimate))
    if (nrow(summ)) {
      n <- dplyr::coalesce(summ$n_total, if (nrow(summ) == 1) listed_n else NA_real_)
      return(tibble::tibble(registry = reg, p_ref = pool(summ$estimate / 100, n), n_total = sum(n),
                            mean_source = paste0("registry summary row: ", paste(summ$device_label, collapse = " + "))))
    }
    if (!is.null(casemix)) {
      cm <- dplyr::filter(casemix, registry == reg, time_yr == !!time_yr, row_type == "subtotal", !is.na(estimate))
      if (nrow(cm)) return(tibble::tibble(registry = reg, p_ref = pool(cm$estimate / 100, cm$n_total),
                                          n_total = sum(cm$n_total),
                                          mean_source = paste0("case-mix subtotals pooled (", cm$table_id[1], ")")))
    }
    dd <- dplyr::filter(dv, row_type == "device", !is.na(estimate), !is.na(n_total))
    tibble::tibble(registry = reg, p_ref = pool(dd$estimate / 100, dd$n_total), n_total = sum(dd$n_total),
                   mean_source = "FALLBACK: n-weighted mean of listed devices")
  })
  # carry report metadata for tables/captions
  meta <- device |>
    dplyr::distinct(registry, report_year, metric_type) |>
    dplyr::group_by(registry) |>
    dplyr::summarise(report_year = paste(unique(report_year), collapse = "/"),
                     metric_type = paste(unique(metric_type), collapse = "/"), .groups = "drop")
  dplyr::left_join(out, meta, by = "registry")
}

# Device points ---------------------------------------------------------------------------------
device_points <- function(device, means, time_yr, families = NULL) {
  d <- device |>
    dplyr::filter(row_type == "device", time_yr == !!time_yr, !is.na(n_total)) |>
    dplyr::filter(has_estimate(dplyr::pick(dplyr::everything()))) |>
    dplyr::mutate(label = short_device_label(dplyr::pick(dplyr::everything())), p = estimate / 100) |>
    dplyr::left_join(dplyr::select(means, registry, p_ref), by = "registry") |>
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
