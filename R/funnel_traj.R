# Funnel trajectories: the standard funnel at several follow-up years ---------------------------------------
#
# Used by funnel_trajectories.Rmd. Nothing new in the funnel itself: each registry's funnel is centred on its
# overall registry mean at year t (registry_means()), the limits come from funnel_limits(), and the points are
# the models at detail level 0 (rows as listed) or 1 (brand, n-weighted). What is added:
#   * the funnel's test statistic as a z-score, so a model can be followed across years;
#   * a simple trajectory label per model (persistent / early / late / isolated, high or low);
#   * a sensitivity analysis replacing the total number of procedures n by the effective sample size
#     n_eff = p (1 - p) / SE^2, with SE from the published KM 95% CI (accounts for censoring).

#' Model points for one follow-up year, with the registry mean as reference and both sample sizes.
#'   level 0: one point per registry row; level 1: one point per brand per registry (n-weighted, as detail
#'   level 1), with Var(p) = sum(w_i^2 Var(p_i)) for the n_eff sensitivity analysis.
#' Var(p_i) = ((ucl - lcl) / 3.92)^2 from the published CI (proportion scale).
tr_points <- function(device, casemix, t, families, models, level = 1, overrides = NULL) {
  means <- registry_means(device, casemix, t, overrides) |> dplyr::filter(is.finite(p_mean), n_total > 0)
  pts <- device_points(device, means, t, families) |>
    dplyr::filter(family %in% models) |>
    dplyr::mutate(var_p = ifelse(is.na(lcl) | is.na(ucl), NA_real_, ((ucl - lcl) / 100 / 3.92)^2))
  if (!nrow(pts)) return(NULL)
  if (level == 1) {
    pts <- pts |>
      dplyr::group_by(registry, family) |>
      dplyr::summarise(n_rows = dplyr::n(), label = dplyr::first(family),
                       p = sum(n_total * p) / sum(n_total),
                       var_p = sum((n_total / sum(n_total))^2 * var_p),
                       n_total = sum(n_total),
                       p_ref = dplyr::first(p_ref), p_mean = dplyr::first(p_mean), .groups = "drop")
  } else {
    pts <- dplyr::mutate(pts, n_rows = 1L)
  }
  pts |>
    dplyr::mutate(time_yr = t, n_eff = p * (1 - p) / var_p,
                  dev = p - p_ref, ratio = p / p_ref) |>
    dplyr::left_join(dplyr::select(means, registry, n_registry = n_total, mean_source), by = "registry")
}

#' Funnel z-score for a point at volume n: the statistic the funnel limits are drawn from, so
#' |z| > 1.96 <=> outside the 95% funnel and |z| > 3.09 <=> outside the 99.8% funnel.
#'   normal: z = (p - p0) / sqrt(p0 (1 - p0) / n)
#'   exact : z = qnorm(mid-p binomial tail of x = round(n p) under Bin(n, p0)), matching exact limits
tr_z <- function(p, p0, n, method = "normal") {
  if (method == "exact") {
    x <- round(n * p)
    pr <- stats::pbinom(x - 1, round(n), p0) + 0.5 * stats::dbinom(x, round(n), p0)
    return(stats::qnorm(pmin(pmax(pr, 1e-12), 1 - 1e-12)))
  }
  (p - p0) / sqrt(p0 * (1 - p0) / n)
}

#' Limit curves for every registry x year (centred on the registry mean), over an n grid.
tr_curves <- function(pts, levels, method, n_col = "n_total") {
  refs <- pts |>
    dplyr::group_by(registry, time_yr) |>
    dplyr::summarise(p_ref = dplyr::first(p_ref), n_max = max(.data[[n_col]], na.rm = TRUE) * 1.3, .groups = "drop")
  n_min <- max(10, floor(min(pts[[n_col]], na.rm = TRUE) * 0.5))
  build_limit_curves(refs, levels, method, n_min = n_min)
}

#' Trajectory label from the z-scores at the core years (in order). High side first; low side mirrors it.
#'   Persistent : beyond the 95% limit at >= 2 consecutive years
#'   Early      : beyond it only at the first two years
#'   Late       : beyond it only from the third year on (e.g. 5 and/or 10 years), not consecutive
#'   Isolated   : any other pattern of non-consecutive exceedances
#'   "(99.8%)"  : beyond the 99.8% limit at least once in that direction
tr_classify <- function(z, years, z95 = stats::qnorm(0.975), z998 = stats::qnorm(0.999)) {
  side <- function(s, dir) {
    ex <- !is.na(z) & s * z > z95
    if (!any(ex)) return(NA_character_)
    consec <- length(ex) > 1 && any(ex[-1] & ex[-length(ex)])
    lab <- if (consec) "Persistent" else if (all(which(ex) <= 2)) "Early" else if (all(which(ex) >= 3)) "Late" else "Isolated"
    paste0(lab, " ", dir, if (any(!is.na(z) & s * z > z998)) " (99.8%)" else "")
  }
  out <- stats::na.omit(c(side(1, "high"), side(-1, "low")))
  if (length(out)) paste(out, collapse = "; ") else "Within"
}

#' Status text per point: "above 99.8%", "above 95%", "within", "below 95%", "below 99.8%".
tr_status <- function(z, z95 = stats::qnorm(0.975), z998 = stats::qnorm(0.999)) {
  dplyr::case_when(is.na(z) ~ NA_character_, z > z998 ~ "above 99.8%", z > z95 ~ "above 95%",
                   z < -z998 ~ "below 99.8%", z < -z95 ~ "below 95%", TRUE ~ "within")
}
