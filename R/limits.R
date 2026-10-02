# Funnel-plot control limits ------------------------------------------------------
#
# For a reference (registry mean) proportion p0 and volume n, the limits give the range a device
# with true rate p0 would fall in with probability `level`. The width depends on p0 and n only;
# devices are then compared against them.
#
# method = "normal": p0 +/- z * sqrt(p0 (1 - p0) / n)   (the approximation used in v1)
# method = "overdispersed": normal limits widened by sqrt(phi), where phi is the registry's
#                    overdispersion factor (Spiegelhalter 2005, multiplicative model; see estimate_phi())
# method = "exact" : binomial quantiles with the interpolation of Spiegelhalter (2005, Stat Med
#                    24:1185), which stays accurate for small n / small p0 where the normal
#                    approximation under-covers (see normal_binom_dist_comparison.R).
# Both treat a KM cumulative revision % as a binomial proportion at time t, ignoring censoring:
# the limits are approximate, and the smaller the number still at risk at t, the more so.

exact_bound <- function(p0, n, prob) {
  rp <- stats::qbinom(prob, n, p0)
  a  <- (stats::pbinom(rp, n, p0) - prob) / stats::dbinom(rp, n, p0)
  (rp - a) / n
}

#' Limits for one reference rate over a vector of n. Returns absolute limits and the two
#' standardised versions: difference (limit - p0) and ratio (limit / p0).
funnel_limits <- function(p0, n, levels = c(0.95, 0.998), method = c("normal", "exact", "overdispersed"),
                          phi = 1, centre = p0) {
  method <- match.arg(method)
  phi_used <- if (method == "overdispersed") max(phi, 1) else 1
  infl <- sqrt(phi_used)
  tidyr::crossing(n = n, level = levels) |>
    dplyr::mutate(
      alpha = 1 - level,
      lower = if (method == "exact") exact_bound(p0, n, alpha / 2) else
        p0 - infl * stats::qnorm(1 - alpha / 2) * sqrt(p0 * (1 - p0) / n),
      upper = if (method == "exact") exact_bound(p0, n, 1 - alpha / 2) else
        p0 + infl * stats::qnorm(1 - alpha / 2) * sqrt(p0 * (1 - p0) / n),
      p0 = p0,
      # standardised scales keep the unclipped lower limit, so both sides of the funnel show;
      # on the absolute scale a revision rate cannot go below 0
      # `centre` is the value the standardised scales are relative to; by default the funnel's own
      # reference p0, so every funnel is centred on 0 (difference) / 1 (ratio)
      dev_lower = lower - centre, dev_upper = upper - centre,
      ratio_lower = lower / centre, ratio_upper = upper / centre,
      lower = pmax(lower, 0),
      level_lab = level_label(level),
      method = .env$method, phi = .env$phi_used
    ) |>
    dplyr::select(-alpha)
}

#' Log-spaced n grid (dense where the funnel bends, sparse where it is flat).
n_grid <- function(n_min = 20, n_max = 1e5, length_out = 400) {
  unique(round(10^seq(log10(n_min), log10(n_max), length.out = length_out)))
}

#' Limit curves for many references at once.
#' refs: one row per reference, with `p_ref`, `n_max` and any key columns (registry, stratum...).
build_limit_curves <- function(refs, levels = c(0.95, 0.998), method = "normal", n_min = 20) {
  keys <- setdiff(names(refs), c("p_ref", "n_max", "phi", "p_centre"))
  purrr::pmap_dfr(refs, function(...) {
    r <- list(...)
    funnel_limits(r$p_ref, n_grid(n_min, r$n_max), levels, method, phi = r$phi %||% 1,
                  centre = r$p_centre %||% r$p_ref) |>
      dplyr::bind_cols(tibble::as_tibble(r[keys]))
  })
}

#' Flag each point against its own reference limits, e.g. "above 99.8%", "above 95%", "within 95%".
classify_points <- function(points, levels = c(0.95, 0.998), method = "normal") {
  lv <- sort(levels, decreasing = TRUE)
  phis <- if ("phi" %in% names(points)) points$phi else rep(1, nrow(points))
  status <- purrr::pmap_chr(list(points$p_ref, points$n_total, points$p, phis), function(p0, n, p, ph) {
    lim <- dplyr::arrange(funnel_limits(p0, n, lv, method, ph), dplyr::desc(level))
    for (i in seq_len(nrow(lim))) {
      if (p > lim$upper[i]) return(paste("above", lim$level_lab[i]))
      if (p < lim$lower[i]) return(paste("below", lim$level_lab[i]))
    }
    paste("within", level_label(min(lv)))
  })
  dplyr::mutate(points, funnel_status = status)
}

# Outlier reporting ----------------------------------------------------------------------------------
#
# A device is a (high) outlier when its rate is above the UPPER control limit of its own registry at
# its own volume n. For each outlier:
#   delta_p        = p - p_ref                (device rate minus the reference: registry mean, or a CI bound of it)
#   excess         = n * delta_p              (revisions above what the registry mean would give)
# The cumulative excess revision is the sum of `excess` over all outliers.

outlier_report <- function(points, level = 0.998, method = "normal") {
  if (!nrow(points)) return(points[0, ])
  phis <- if ("phi" %in% names(points)) points$phi else rep(1, nrow(points))
  up <- purrr::pmap_dbl(list(points$p_ref, points$n_total, phis),
                        function(p0, n, ph) funnel_limits(p0, n, level, method, ph)$upper)
  out <- points |>
    dplyr::mutate(upper_limit = up) |>
    dplyr::filter(p > upper_limit) |>          # ONLY devices outside (above) the control limit
    dplyr::mutate(delta_p = p - p_ref,
                  excess = n_total * delta_p) |>
    dplyr::arrange(dplyr::desc(excess))
  # guard: every row counted in the outlier metrics must be above its registry's upper limit
  stopifnot(all(out$p > out$upper_limit))
  out
}

#' Per-registry summary (+ overall row) of an outlier report.
#' `points` is every device on the plot, so registries without outliers still get a 0 row.
#' n_implants and both excess sums are over OUTLIERS ONLY (devices above the upper control limit);
#' devices inside the funnel contribute nothing, even if their rate is above the registry mean.
outlier_summary <- function(outliers, points) {
  per_reg <- points |>
    dplyr::distinct(registry) |>
    dplyr::left_join(
      outliers |>
        dplyr::group_by(registry) |>
        dplyr::summarise(n_outliers = dplyr::n(), n_implants = sum(n_total),
                         cumulative_excess = sum(excess),
                         .groups = "drop"),
      by = "registry") |>
    dplyr::left_join(dplyr::count(points, registry, name = "n_devices"), by = "registry") |>
    dplyr::mutate(dplyr::across(c(n_outliers, n_implants, cumulative_excess),
                                ~ dplyr::coalesce(.x, 0))) |>
    dplyr::arrange(registry)
  dplyr::bind_rows(
    per_reg,
    dplyr::summarise(per_reg, registry = "All registries", n_devices = sum(n_devices),
                     n_outliers = sum(n_outliers), n_implants = sum(n_implants),
                     cumulative_excess = sum(cumulative_excess))
  ) |>
    dplyr::select(registry, n_devices, n_outliers, n_implants, cumulative_excess)
}


# Overdispersion ------------------------------------------------------------------------------------
#
# Registry data usually vary more between devices than binomial sampling allows. Spiegelhalter's
# (2005) multiplicative model estimates phi = mean(z^2) over a registry's devices, with z-scores
# winsorised at the 10th/90th percentiles so the outliers being screened for do not inflate it.
# phi <= 1 means no overdispersion (limits unchanged).
estimate_phi <- function(points, winsor = 0.1) {
  points |>
    dplyr::mutate(p0 = if ("p_mean" %in% names(points)) p_mean else p_ref) |>   # dispersion around the mean
    dplyr::filter(!is.na(p), !is.na(n_total), !is.na(p0)) |>
    dplyr::group_by(registry) |>
    dplyr::summarise(
      n_devices_phi = dplyr::n(),
      phi = {
        z <- (p - p0) / sqrt(p0 * (1 - p0) / n_total)
        if (length(z) >= 5) {
          q <- stats::quantile(z, c(winsor, 1 - winsor))
          z <- pmin(pmax(z, q[1]), q[2])
        }
        max(mean(z^2), 1)
      }, .groups = "drop")
}
