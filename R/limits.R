# Funnel-plot control limits ------------------------------------------------------
#
# For a reference (registry mean) proportion p0 and volume n, the limits give the range a device
# with true rate p0 would fall in with probability `level`. The width depends on p0 and n only;
# devices are then compared against them.
#
# method = "normal": p0 +/- z * sqrt(p0 (1 - p0) / n)   (the approximation used in v1)
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
funnel_limits <- function(p0, n, levels = c(0.95, 0.998), method = c("normal", "exact")) {
  method <- match.arg(method)
  tidyr::crossing(n = n, level = levels) |>
    dplyr::mutate(
      alpha = 1 - level,
      lower = if (method == "normal") {
        p0 - stats::qnorm(1 - alpha / 2) * sqrt(p0 * (1 - p0) / n)
      } else exact_bound(p0, n, alpha / 2),
      upper = if (method == "normal") {
        p0 + stats::qnorm(1 - alpha / 2) * sqrt(p0 * (1 - p0) / n)
      } else exact_bound(p0, n, 1 - alpha / 2),
      p0 = p0,
      # standardised scales keep the unclipped lower limit, so both sides of the funnel show;
      # on the absolute scale a revision rate cannot go below 0
      dev_lower = lower - p0, dev_upper = upper - p0,
      ratio_lower = lower / p0, ratio_upper = upper / p0,
      lower = pmax(lower, 0),
      level_lab = level_label(level),
      method = method
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
  keys <- setdiff(names(refs), c("p_ref", "n_max"))
  purrr::pmap_dfr(refs, function(...) {
    r <- list(...)
    funnel_limits(r$p_ref, n_grid(n_min, r$n_max), levels, method) |>
      dplyr::bind_cols(tibble::as_tibble(r[keys]))
  })
}

#' Flag each point against its own reference limits, e.g. "above 99.8%", "above 95%", "within 95%".
classify_points <- function(points, levels = c(0.95, 0.998), method = "normal") {
  lv <- sort(levels, decreasing = TRUE)
  status <- purrr::pmap_chr(list(points$p_ref, points$n_total, points$p), function(p0, n, p) {
    lim <- dplyr::arrange(funnel_limits(p0, n, lv, method), dplyr::desc(level))
    for (i in seq_len(nrow(lim))) {
      if (p > lim$upper[i]) return(paste("above", lim$level_lab[i]))
      if (p < lim$lower[i]) return(paste("below", lim$level_lab[i]))
    }
    paste("within", level_label(min(lv)))
  })
  dplyr::mutate(points, funnel_status = status)
}
