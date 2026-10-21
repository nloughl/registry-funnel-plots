# Funnels over follow-up time ----------------------------------------------------------------------------
#
# Used by funnel_over_time.Rmd. Keeps the funnel's within-registry comparison (each model against its own
# registry) and adds the time and cross-registry machinery of Combescure et al. (2025):
#   1. per registry row and year: log cumulative hazard ratio vs the REST of the registry (cloglog scale)
#      and its z-score - the quantity a funnel tests, written so it can be followed over years;
#   2. interval (conditional) hazard ratios between follow-up years, to separate early from late failure;
#   3. pooling of the log HRs across registries per year (variant nested in registry, shared comparator),
#      and one model across years (log HR linear in log t, random intercepts and slopes);
#   4. trajectory classification and its false-alarm rate under the null (simulation).
# Needs crr_meta.R (cll_from_crr, crr_inputs, crr_row_cov) and data.R (registry_means).

# Comparator: the rest of the registry ---------------------------------------------------------------------
#' Registry mean minus the model's rows, per registry x year.
#'   p_rest = (N p_reg - sum n_o p_o) / (N - sum n_o)           (KM rates treated as additive: approximate)
#'   Var(p_rest) ~= Var(p_reg) * N / (N - sum n_o)               (variance scales with 1 / size)
#' Var(p_reg) comes from the registry mean's published (or pooled) 95% CI. Only rows with an estimate at t
#' are subtracted; rows without one stay in the "rest" (noted in `unsubtracted_n`).
ft_comparators <- function(device, casemix, inp, years, overrides = NULL, max_unsubtracted = 0.10) {
  purrr::map_dfr(years, function(t) {
    m <- registry_means(device, casemix, t, overrides)
    m <- dplyr::filter(m, is.finite(p_mean), n_total > 0)
    if (!nrow(m)) return(NULL)
    purrr::map_dfr(seq_len(nrow(m)), function(i) {
      r <- m[i, ]
      rows_all <- dplyr::distinct(dplyr::filter(inp, registry == r$registry), row_id, n_total)
      rows_t <- dplyr::filter(inp, registry == r$registry, time_yr == t)
      n_o <- sum(rows_t$n_total); N <- r$n_total
      if (n_o >= N || !nrow(rows_t)) return(NULL)
      p_rest <- (N * r$p_mean - sum(rows_t$n_total * rows_t$estimate / 100)) / (N - n_o)
      v_reg <- if (is.na(r$p_lcl)) NA_real_ else ((r$p_ucl - r$p_lcl) / 3.92)^2
      v_rest <- v_reg * N / (N - n_o)
      deriv <- 1 / ((1 - p_rest) * (-log(1 - p_rest)))         # d cloglog / dp
      tibble::tibble(registry = r$registry, time_yr = t, n_registry = N, p_registry = r$p_mean,
                     p_reg_lcl = r$p_lcl, p_reg_ucl = r$p_ucl, mean_source = r$mean_source,
                     n_model = n_o, model_share = n_o / N, rows_subtracted = nrow(rows_t),
                     unsubtracted_n = sum(rows_all$n_total) - n_o,
                     n_rest = N - n_o, p_rest = p_rest, crr_rest = 100 * p_rest,
                     y_rest = cll_from_crr(100 * p_rest), v_rest_y = v_rest * deriv^2)
    })
  }) |> dplyr::filter(is.finite(y_rest), p_rest > 0) |>
    # when model rows without an estimate at t make up > max_unsubtracted of the "rest", the comparator
    # still contains much of the model itself (e.g. NJR 20 y, LROI 16 y): drop that registry-year
    dplyr::mutate(comparator_ok = unsubtracted_n / n_rest <= max_unsubtracted)
}

# Cumulative contrasts -------------------------------------------------------------------------------------
#' One row per registry row (variant) x year: log cumulative HR vs the rest of the registry, SE, z, HR.
#' The variant and the rest of the registry are disjoint, so their sampling errors are independent.
ft_contrasts <- function(inp, comp) {
  inp |>
    dplyr::inner_join(dplyr::select(comp, registry, time_yr, crr_rest, y_rest, v_rest_y, p_registry), by = c("registry", "time_yr")) |>
    dplyr::mutate(lhr = yi - y_rest, se = sqrt(vi + v_rest_y), z = lhr / se,
                  hr = exp(lhr), hr_lcl = exp(lhr - 1.96 * se), hr_ucl = exp(lhr + 1.96 * se), level = "variant")
}

#' Brand-level contrast per registry x year: the registry's rows collapsed, n-weighted as in detail level 1
#' (crr_nw_collapse()), vs the rest of the registry.
ft_brand_contrasts <- function(inp, comp) {
  inp |>
    dplyr::group_by(model, registry, time_yr) |>
    dplyr::summarise(cl = list(crr_nw_collapse(estimate, vi, n_total)), n_total = sum(n_total),
                     n_rows = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(yi = purrr::map_dbl(cl, "y"), vi = purrr::map_dbl(cl, "v")) |>
    dplyr::select(-cl) |>
    dplyr::mutate(estimate = crr_from_cll(yi), row_id = paste0(registry, ": ", model, " (all rows)"),
                  label = paste0(model, " (all rows)")) |>
    dplyr::inner_join(dplyr::select(comp, registry, time_yr, crr_rest, y_rest, v_rest_y, p_registry), by = c("registry", "time_yr")) |>
    dplyr::mutate(lhr = yi - y_rest, se = sqrt(vi + v_rest_y), z = lhr / se,
                  hr = exp(lhr), hr_lcl = exp(lhr - 1.96 * se), hr_ucl = exp(lhr + 1.96 * se), level = "brand")
}

# Interval (conditional) contrasts ---------------------------------------------------------------------------
#' Revision between consecutive follow-up years t1 < t2 among knees unrevised at t1.
#' Interval cumulative hazard dH = H(t2) - H(t1), H = -log(1 - CRR/100), with H(0) = 0.
#' Var(dH) ~= Var(H2) - Var(H1) (independent increments), Var(H) = v * H^2 from the cloglog variance v.
#' When published-CI rounding makes that difference <= 0, Var(dH) ~= Var(H2) * dH / H2 is used (flagged).
#' Interval log HR = log(dH_model / dH_rest); variance by the delta method.
ft_interval_one <- function(est_pct, v, years, cuts) {
  keep <- years %in% cuts
  est_pct <- est_pct[keep]; v <- v[keep]; years <- years[keep]
  o <- order(years); est_pct <- est_pct[o]; v <- v[o]; years <- years[o]
  H <- -log(1 - est_pct / 100); VH <- v * H^2
  H0 <- c(0, H); V0 <- c(0, VH); y0 <- c(0, years)
  dH <- diff(H0); dV <- diff(V0)
  approx <- dV <= 0
  dV[approx] <- VH[approx] * dH[approx] / H[approx]
  tibble::tibble(t1 = utils::head(y0, -1), t2 = years, dH = dH, var_dH = dV, var_approx = approx)
}

ft_intervals <- function(contr, cuts = c(1, 3, 5, 10, 15, 20)) {
  contr |>
    dplyr::group_by(level, model, registry, row_id, label) |>
    dplyr::group_modify(function(d, k) {
      m <- ft_interval_one(d$estimate, d$vi, d$time_yr, cuts)
      r <- ft_interval_one(d$crr_rest, d$v_rest_y, d$time_yr, cuts)
      dplyr::inner_join(m, r, by = c("t1", "t2"), suffix = c("_m", "_r"))
    }) |>
    dplyr::ungroup() |>
    dplyr::filter(dH_m > 0, dH_r > 0) |>
    dplyr::mutate(interval = paste0(t1, "–", t2, " y"),
                  lhr = log(dH_m / dH_r), se = sqrt(var_dH_m / dH_m^2 + var_dH_r / dH_r^2), z = lhr / se,
                  hr = exp(lhr), hr_lcl = exp(lhr - 1.96 * se), hr_ucl = exp(lhr + 1.96 * se),
                  crr_cond_m = 100 * (1 - exp(-dH_m)), crr_cond_r = 100 * (1 - exp(-dH_r)))
}

# Pooling across registries ----------------------------------------------------------------------------------
#' Sampling covariance of the variant contrasts at one year: rows of the same registry share the
#' comparator, so Cov(lhr_i, lhr_j) = Var(y_rest) for i != j in the same registry.
ft_V_year <- function(d) {
  V <- diag(d$se^2, nrow(d))
  for (i in seq_len(nrow(d))) for (j in seq_len(nrow(d)))
    if (i != j && d$registry[i] == d$registry[j]) V[i, j] <- d$v_rest_y[i]
  V
}

#' Pooled log HR at one year (or one interval): REML, random effects registry / row (dropped when not
#' estimable), t-based CI with containment df when >= k_reg_min registries, z otherwise.
ft_pool <- function(d, k_reg_min = 3, V = NULL) {
  k_reg <- dplyr::n_distinct(d$registry)
  if (k_reg < 2) return(NULL)
  multi <- any(table(d$registry) > 1)
  random <- if (multi) list(~ 1 | registry, ~ 1 | row_id) else list(~ 1 | registry)
  use_t <- k_reg >= k_reg_min
  args <- list(yi = d$lhr, V = if (is.null(V)) d$se^2 else V, random = random, data = d, method = "REML")
  if (use_t) args <- c(args, test = "t", dfs = "contain")
  f <- tryCatch(suppressWarnings(do.call(metafor::rma.mv, args)), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  pr <- stats::predict(f)
  tibble::tibble(k = f$k, k_reg = k_reg, registries = paste(sort(unique(d$registry)), collapse = ", "),
                 ci_method = if (use_t) "t (contain df)" else "z", exploratory = !use_t,
                 lhr = f$b[1], hr = exp(f$b[1]), hr_lcl = exp(f$ci.lb), hr_ucl = exp(f$ci.ub),
                 pi_lcl = exp(pr$pi.lb), pi_ucl = exp(pr$pi.ub), pval = f$pval,
                 Q = f$QE, Q_p = f$QEp, sigma2 = paste(sprintf("%s %.3f", f$s.names, f$sigma2), collapse = "; "))
}

#' Joint model across years for the variant contrasts.
#'   lhr(t) = a + b log t + (u0_r + u1_r log t) + (w0_j + w1_j log t) + e,  e ~ N(0, V)
#' V combines, within a registry: each row's own covariance across years (crr_row_cov) and the shared
#' comparator's covariance across years (all rows and years of the registry).
#' a = log HR at 1 year; b = change in log HR per unit log t (b > 0: the excess grows with follow-up).
ft_joint <- function(contr, slopes = TRUE) {
  d <- contr |> dplyr::arrange(registry, row_id, time_yr) |> dplyr::group_by(row_id) |>
    dplyr::filter(!(estimate == dplyr::lag(estimate, default = -1) & lcl == dplyr::lag(lcl, default = -1))) |>
    dplyr::ungroup()
  blocks <- lapply(split(d, d$registry)[unique(d$registry)], function(b) {
    yrs <- sort(unique(b$time_yr))
    rest <- dplyr::distinct(b, time_yr, crr_rest, v_rest_y) |> dplyr::arrange(time_yr)
    Cr <- crr_row_cov(rest$crr_rest, rest$v_rest_y); dimnames(Cr) <- list(rest$time_yr, rest$time_yr)
    M <- Cr[as.character(b$time_yr), as.character(b$time_yr), drop = FALSE]
    for (rid in unique(b$row_id)) {
      ix <- which(b$row_id == rid)
      M[ix, ix] <- M[ix, ix] + crr_row_cov(b$estimate[ix], b$vi[ix])
    }
    M
  })
  V <- as.matrix(Matrix::bdiag(blocks))
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 1e-10 * max(ev)) V <- as.matrix(Matrix::nearPD(V, keepDiag = TRUE)$mat)
  d$lt <- log(d$time_yr)
  rand <- if (slopes) list(~ lt | registry, ~ lt | row_id) else list(~ 1 | registry, ~ 1 | row_id)
  args <- list(yi = d$lhr, V = V, mods = ~ lt, random = rand, data = d, method = "REML", test = "t", dfs = "contain")
  if (slopes) args$struct <- c("GEN", "GEN")
  f <- tryCatch(suppressWarnings(do.call(metafor::rma.mv, args)), error = function(e) NULL)
  if (is.null(f) && slopes) return(ft_joint(contr, slopes = FALSE))
  list(fit = f, data = d, V = V, slopes = slopes)
}

ft_joint_predict <- function(jm, t) {
  p <- stats::predict(jm$fit, newmods = log(t))
  tibble::tibble(time_yr = t, hr = exp(p$pred), hr_lcl = exp(p$ci.lb), hr_ucl = exp(p$ci.ub))
}

# Trajectory classification -----------------------------------------------------------------------------------
#' For one line (a variant or a registry's brand-level contrast): flags on the core years and the late
#' interval. High side shown; the low side mirrors it.
#'   Persistent high : z > z95 at >= 2 consecutive core years
#'   Late high       : interval z > z95 for the last core interval (e.g. 5-10 y) without persistent high
#'   Early high      : z > z95 only at the first core years (1 and/or 3 y), not later
#'   Isolated high   : any other single exceedance
#'   (and the same for low; "Within" otherwise). "Strong" = any z beyond z998.
ft_classify <- function(z_by_year, z_late, core, z95 = stats::qnorm(0.975), z998 = stats::qnorm(0.999)) {
  z <- z_by_year[as.character(core)]
  side <- function(s) {
    ex <- !is.na(z) & s * z > z95
    consec <- any(ex[-1] & ex[-length(ex)])
    late <- !is.na(z_late) && s * z_late > z95
    early_only <- any(ex) && all(core[ex] <= core[min(2, length(core))])
    c(persistent = consec, late = late && !consec, early = early_only && !consec && !late,
      isolated = any(ex) && !consec && !late && !early_only,
      strong = any(!is.na(z) & s * z > z998))
  }
  h <- side(1); l <- side(-1)
  name <- function(x, dir) {
    if (x["persistent"]) paste("Persistent", dir) else if (x["late"]) paste("Late", dir) else
      if (x["early"]) paste("Early", dir) else if (x["isolated"]) paste("Isolated", dir) else NA_character_
  }
  hl <- name(h, "high"); ll <- name(l, "low")
  hl <- if (!is.na(hl) && h["strong"]) paste0(hl, " (beyond 99.8%)") else hl
  ll <- if (!is.na(ll) && l["strong"]) paste0(ll, " (beyond 99.8%)") else ll
  out <- paste(stats::na.omit(c(hl, ll)), collapse = "; ")
  if (out == "") "Within" else out
}

#' False-alarm rate of each flagging rule under the null (true HR = 1 at every year), by simulation from
#' the line's own correlation of z across the core years (from the joint covariance) and of the late
#' interval z with the cumulative z's (approximated as independent of earlier years, which holds for
#' independent increments). Returns the share of simulated null trajectories flagged high by each rule.
ft_false_alarm <- function(R, core, n_sim = 20000, z95 = stats::qnorm(0.975), seed = 1) {
  set.seed(seed)
  Z <- MASS::mvrnorm(n_sim, rep(0, ncol(R)), R)
  zl <- stats::rnorm(n_sim)
  ex <- Z > z95
  consec <- if (ncol(ex) > 1) apply(ex[, -1, drop = FALSE] & ex[, -ncol(ex), drop = FALSE], 1, any) else rep(FALSE, n_sim)
  c(single_year = mean(ex[, min(3, ncol(ex))]),            # e.g. 5-year funnel alone
    any_year = mean(apply(ex, 1, any)),
    persistent = mean(consec),
    late_interval = mean(zl > z95 & !consec),
    any_rule = mean(apply(ex, 1, any) | zl > z95))
}
