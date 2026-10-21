# Leave-one-out and influence analysis of registry means ------------------------------------------
#
# Used by loo_influence.Rmd. With few models (rare joints), one model can drive the registry mean and
# camouflage outliers. Two ways to recompute the mean without a model:
#
#   subtraction  remove model i's procedures from the PUBLISHED registry total:
#                p_(-i) = (N p_pub - n_i p_i) / (N - n_i)
#                patient-weighted (same estimand as the published mean), keeps unlisted "other" models
#   meta-analytic  random-effects meta-analysis of the model rates (metafor::rma on the complementary log-log
#                scale by default, or logit; SE from each model's published Kaplan-Meier CI), then
#                metafor::leave1out() / influence(); every rate is back-transformed to the proportion scale
#                model-weighted (the expected rate of a typical MODEL, not of a typical patient)
#
# Reference options per model (all evaluated with the same control limits):
#   pub      published registry mean
#   loo_sub  the model's own subtraction leave-one-out mean
#   re       random-effects pooled mean (all models)
#   loo_re   the model's own random-effects leave-one-out mean
#   min_loo  the lowest subtraction leave-one-out mean in the registry
#   min_loo_re the lowest random-effects (metafor) leave-one-out mean in the registry: the "least
#            upward-biased" mean recommended in supervisor guidance (2026-10) for registries with few models

LOO_REFS <- c(pub = "Published mean", loo_sub = "Own leave-one-out (subtraction)",
              re = "Random-effects mean", loo_re = "Own leave-one-out (random effects)",
              min_loo = "Lowest leave-one-out mean (subtraction)",
              min_loo_re = "Lowest leave-one-out mean (random effects)")
LOO_REF_COLOURS <- c(pub = "#222222", loo_sub = "#1B9E77", re = "#2166AC", loo_re = "#7570B3", min_loo = "#D95F02",
                     min_loo_re = "#E7298A")
LOO_STATUS_COLOURS <- c(above = "#B2182B", within = "#9E9E9E", below = "#2166AC")

# Transformation of the revision proportion p used for the meta-analysis --------------------------------
#   cloglog  log(-log(1 - p)) = log(-log S): the scale registries build Kaplan-Meier CIs on (Kalbfleisch &
#            Prentice; the published CIs are symmetric on it) and the log cumulative hazard; used for pooled
#            UKA survival by Heaps et al. 2019 (J Arthroplasty, metafor + REML). Default.
#   logit    log(p / (1 - p)): metafor's standard for proportions (measure "PLO"); sensitivity analysis.
# f = transform, inv = back-transform to p, dfdp = derivative (delta method).
MA_SCALES <- list(
  cloglog = list(label = "complementary log-log", f = function(p) log(-log(1 - p)),
                 inv = function(x) 1 - exp(-exp(x)), dfdp = function(p) 1 / (-(1 - p) * log(1 - p))),
  logit   = list(label = "logit", f = stats::qlogis, inv = stats::plogis, dfdp = function(p) 1 / (p * (1 - p))))
ma_scale <- function(scale) {
  if (!scale %in% names(MA_SCALES)) stop("scale must be one of: ", paste(names(MA_SCALES), collapse = ", "))
  MA_SCALES[[scale]]
}

#' Model rows for one joint at time_yr, with estimates and variances on the chosen scale.
#' SE_i = [f(ucl_i) - f(lcl_i)] / (2 * 1.96), from the published Kaplan-Meier 95% CI.
#' AOANJRR elbow stems (printed per diagnosis) are pooled per stem: n-weighted rate; each row's SE is taken
#' from its own CI on the chosen scale, moved to the proportion scale (delta method), combined n-weighted,
#' and moved back to the chosen scale at the pooled rate.
loo_models <- function(ext, cfg, time_yr, z = stats::qnorm(0.975), scale = "cloglog") {
  sc <- ma_scale(scale)
  m <- rare_means(ext, cfg, time_yr)
  if (!nrow(m)) return(NULL)
  d <- rare_points(ext, m, cfg, time_yr, elbow_by_diagnosis = TRUE) |>
    dplyr::mutate(lcl = lcl / 100, ucl = ucl / 100)
  by_dx <- !is.na(d$diagnosis) & d$registry == "AOANJRR" & ext$joint == "ELBOW"
  single <- d[!by_dx, ] |>
    dplyr::mutate(
      ok = !is.na(lcl) & !is.na(ucl) & lcl > 0 & ucl < 1 & p > 0 & p < 1,
      yi = ifelse(ok, sc$f(p), NA_real_),
      sei = ifelse(ok, (sc$f(ucl) - sc$f(lcl)) / (2 * z), NA_real_),
      se_source = ifelse(ok, paste0("published CI (", sc$label, " scale)"), NA_character_),
      model = family, n_rows = 1L)
  pooled <- d[by_dx, ] |>
    dplyr::mutate(se_p = ((sc$f(ucl) - sc$f(lcl)) / (2 * z)) / sc$dfdp(p)) |>   # row SE, proportion scale
    dplyr::group_by(registry, family) |>
    dplyr::summarise(p = sum(p * n_total) / sum(n_total),
                     se_p = sqrt(sum((n_total / sum(n_total))^2 * se_p^2)),
                     n_total = sum(n_total), n_rows = dplyr::n(),
                     diagnosis = paste(unique(diagnosis), collapse = " + "), .groups = "drop") |>
    dplyr::mutate(lcl = NA_real_, ucl = NA_real_, ok = p > 0 & p < 1,
                  yi = sc$f(p), sei = se_p * sc$dfdp(p),
                  se_source = paste0("pooled across diagnoses (delta method, ", sc$label, ")"), model = family) |>
    dplyr::select(-se_p)
  out <- dplyr::bind_rows(single, pooled)
  # fallback (no usable CI, or a 0% estimate): events = p * n, 0.5 continuity correction, binomial SE (delta)
  fb <- is.na(out$yi) | is.na(out$sei)
  if (any(fb)) {
    e <- out$p[fb] * out$n_total[fb]; nn <- out$n_total[fb]
    pc <- (e + 0.5) / (nn + 1)
    out$yi[fb]  <- sc$f(pc)
    out$sei[fb] <- sqrt(pc * (1 - pc) / nn) * sc$dfdp(pc)
    out$se_source[fb] <- "binomial approximation (0.5 correction)"
  }
  # duplicate model names within a registry (e.g. two Infinity talar combinations): keep the row label
  out <- out |>
    dplyr::group_by(registry, model) |>
    dplyr::mutate(model = if (dplyr::n() > 1) label else model) |>
    dplyr::ungroup() |>
    dplyr::mutate(model = ifelse(n_rows > 1, paste0(model, " [", n_rows, " dx]"), model),
                  vi = sei^2, joint = ext$joint, time_yr = time_yr, scale = scale) |>
    dplyr::select(joint, time_yr, registry, model, n_total, p, lcl, ucl, yi, vi, sei, se_source, scale,
                  n_rows, diagnosis, low_at_risk)
  list(models = out, means = m)
}

#' Meta-analysis and leave-one-out for one registry (k >= 3 models).
#' REML random effects with the Knapp-Hartung adjustment; a common-effect (inverse-variance) fit for comparison.
loo_fit <- function(md, N, p_pub, test = "knha", scale = "cloglog") {
  inv <- ma_scale(scale)$inv
  md <- dplyr::arrange(md, dplyr::desc(n_total))
  k <- nrow(md)
  re <- metafor::rma(yi, vi, data = md, slab = model, method = "REML", test = test)   # slab: model names for metafor's own plots
  ce <- metafor::rma(yi, vi, data = md, method = "EE")
  l1 <- as.data.frame(metafor::leave1out(re, transf = inv))                     # back to the proportion scale
  inf <- stats::influence(re)   # S3 method influence.rma.uni from metafor
  infd <- as.data.frame(inf$inf)
  w_re <- stats::weights(re) / 100; w_ce <- stats::weights(ce) / 100
  models <- md |>
    dplyr::mutate(
      share_n = n_total / N, w_ce = w_ce, w_re = w_re,
      # subtraction leave-one-out from the published total
      loo_sub = (N * p_pub - n_total * p) / (N - n_total),
      loo_sub_se = sqrt(loo_sub * (1 - loo_sub) / (N - n_total)),
      d_sub = loo_sub - p_pub,
      # meta-analytic leave-one-out (back-transformed)
      loo_re = l1$estimate, loo_re_lb = l1$ci.lb, loo_re_ub = l1$ci.ub,
      d_re = loo_re - inv(re$b[1]),
      tau2_without = l1$tau2, I2_without = l1$I2,
      rstudent = infd$rstudent, dffits = infd$dffits, cook_d = infd$cook.d,
      cov_r = infd$cov.r, hat = infd$hat, tau2_del = infd$tau2.del,
      influential = as.character(infd$inf) %in% "*")
  summ <- tibble::tibble(
    k = k, N = N, p_pub = p_pub,
    scale = scale,
    p_re = inv(re$b[1]), p_re_lb = inv(re$ci.lb), p_re_ub = inv(re$ci.ub),
    p_ce = inv(ce$b[1]),
    tau2 = re$tau2, I2 = re$I2, QE = re$QE, QEp = re$QEp,
    pi_lb = inv(stats::predict(re)$pi.lb), pi_ub = inv(stats::predict(re)$pi.ub),
    loo_sub_min = min(models$loo_sub), loo_sub_max = max(models$loo_sub),
    loo_re_min = min(models$loo_re), loo_re_max = max(models$loo_re),
    drop_for_min = models$model[which.min(models$loo_sub)],
    drop_for_min_re = models$model[which.min(models$loo_re)],
    top_model = models$model[which.max(models$share_n)], top_share = max(models$share_n),
    top_w_re = models$w_re[which.max(models$share_n)])
  list(re = re, ce = ce, models = models, summary = summ)
}

#' Run loo_models() + loo_fit() for every registry with >= k_min models.
loo_run <- function(ext, cfg, time_yr, k_min = 3, test = "knha", scale = "cloglog") {
  lm <- loo_models(ext, cfg, time_yr, scale = scale)
  if (is.null(lm)) return(NULL)
  regs <- lm$models |> dplyr::count(registry) |> dplyr::filter(n >= k_min) |> dplyr::pull(registry)
  fits <- purrr::map(stats::setNames(regs, regs), function(r) {
    mm <- dplyr::filter(lm$means, registry == r)
    loo_fit(dplyr::filter(lm$models, registry == r), mm$n_total, mm$p_mean, test, scale)
  })
  list(
    joint = ext$joint, time_yr = time_yr, scale = scale, fits = fits, means = lm$means,
    skipped = setdiff(unique(c(lm$means$registry, lm$models$registry)), regs),
    models = purrr::imap_dfr(fits, function(f, r) dplyr::mutate(f$models, registry = r)),
    summary = purrr::imap_dfr(fits, function(f, r) dplyr::mutate(f$summary, registry = r, .before = 1)) |>
      dplyr::left_join(dplyr::select(lm$means, registry, p_lcl, p_ucl), by = "registry") |>
      dplyr::mutate(joint = ext$joint, time_yr = time_yr, .before = 1))
}

#' Each model's reference under every option, with its limits and status.
#' method: funnel limit method ("exact" recommended for small volumes); level: flagging level.
loo_status <- function(run, level = 0.998, method = "exact") {
  mods <- run$models |>
    dplyr::left_join(dplyr::select(run$summary, registry, p_pub, p_re, loo_sub_min, loo_re_min), by = "registry")
  long <- dplyr::bind_rows(
    dplyr::mutate(mods, ref = "pub",     p0 = p_pub),
    dplyr::mutate(mods, ref = "loo_sub", p0 = loo_sub),
    dplyr::mutate(mods, ref = "re",      p0 = p_re),
    dplyr::mutate(mods, ref = "loo_re",  p0 = loo_re),
    dplyr::mutate(mods, ref = "min_loo", p0 = loo_sub_min),
    dplyr::mutate(mods, ref = "min_loo_re", p0 = loo_re_min))
  lim <- purrr::map2_dfr(long$p0, long$n_total, function(p0, n) {
    l <- funnel_limits(p0, n, level, method); tibble::tibble(lower = l$lower, upper = l$upper,
                                                             dev_lower = l$dev_lower, dev_upper = l$dev_upper)
  })
  dplyr::bind_cols(long, lim) |>
    dplyr::mutate(ref = factor(ref, levels = names(LOO_REFS)),
                  dev = p - p0,
                  status = dplyr::case_when(p > upper ~ "above", p < lower ~ "below", TRUE ~ "within"),
                  status = factor(status, levels = c("above", "within", "below")),
                  level = level, method = method)
}

#' Meta-analytic outlier test: the externally standardised (leave-one-out) residual from influence().
#' z_i = (y_i - mu_(-i)) / sqrt(v_i + tau2_(-i) + Var(mu_(-i)))  -> flagged if |z| > qnorm(1 - (1 - level)/2)
loo_rstudent_status <- function(run, level = 0.998) {
  zc <- stats::qnorm(1 - (1 - level) / 2)
  run$models |>
    dplyr::mutate(z_crit = zc,
                  status = dplyr::case_when(rstudent > zc ~ "above", rstudent < -zc ~ "below", TRUE ~ "within"),
                  status = factor(status, levels = c("above", "within", "below")))
}

#' Real tail probability of the normal and exact upper limits under the binomial (calibration check).
limit_calibration <- function(p0, n, level = 0.998) {
  target <- (1 - level) / 2
  purrr::map_dfr(n, function(nn) {
    lims <- dplyr::bind_rows(
      dplyr::mutate(funnel_limits(p0, nn, level, "normal"), method = "normal"),
      dplyr::mutate(funnel_limits(p0, nn, level, "exact"), method = "exact"))
    lims |> dplyr::transmute(
      n = nn, p0 = p0, method,
      upper, lower,
      # P(observed proportion > upper) and P(< lower) for a model whose true rate is p0
      p_above = 1 - stats::pbinom(floor(upper * nn + 1e-9), nn, p0),
      p_below = ifelse(lower <= 0, 0, stats::pbinom(ceiling(lower * nn - 1e-9) - 1, nn, p0)),
      target = target)
  })
}

#' Sensitivity of the meta-analysis to the transformation scale.
#' Fits every joint x follow-up year on each scale and returns one row per joint, year, registry and scale:
#' pooled mean (CI), prediction interval, tau2, I2, the lowest leave-one-out mean and the model left out,
#' and the models flagged against the lowest leave-one-out mean (funnel limit `method` at `level`).
loo_scale_compare <- function(exts, cfg, years, scales = c("cloglog", "logit"), k_min = 3, test = "knha",
                              level = 0.90, method = "normal") {
  grid <- tidyr::expand_grid(j = names(exts), t = years, sc = scales)
  purrr::pmap_dfr(grid, function(j, t, sc) {
    r <- loo_run(exts[[j]], cfg, t, k_min = k_min, test = test, scale = sc)
    if (is.null(r) || !nrow(r$summary)) return(NULL)
    fl <- loo_status(r, level, method) |>
      dplyr::filter(ref == "min_loo_re", status != "within") |>
      dplyr::group_by(registry) |>
      dplyr::summarise(n_flag = dplyr::n(), flagged = paste(sort(paste0(model, " (", status, ")")), collapse = "; "),
                       .groups = "drop")
    r$summary |>
      dplyr::select(joint, time_yr, registry, scale, k, p_re, p_re_lb, p_re_ub, pi_lb, pi_ub, tau2, I2,
                    loo_re_min, drop_for_min_re) |>
      dplyr::left_join(fl, by = "registry") |>
      dplyr::mutate(n_flag = dplyr::coalesce(n_flag, 0L), flagged = dplyr::coalesce(flagged, ""))
  })
}

#' One row per joint, year and registry: the two scales side by side, with differences in percentage points
#' and whether the lowest-LOO model and the flagged models agree.
loo_scale_wide <- function(cmp, a = "cloglog", b = "logit") {
  x <- dplyr::filter(cmp, scale == a); y <- dplyr::filter(cmp, scale == b)
  dplyr::inner_join(x, y, by = c("joint", "time_yr", "registry", "k"), suffix = c("_a", "_b")) |>
    dplyr::mutate(d_re_pp = 100 * (p_re_a - p_re_b),
                  d_loo_pp = 100 * (loo_re_min_a - loo_re_min_b),
                  same_drop = drop_for_min_re_a == drop_for_min_re_b,
                  same_flags = flagged_a == flagged_b)
}

# plots -------------------------------------------------------------------------------------------

#' Forest plot: model rates (published CI) with the published mean (black) and the random-effects mean
#' (blue band = its CI, dotted = prediction interval). Point size ~ n.
plot_loo_forest <- function(run, ...) {
  inv <- ma_scale(if (is.null(run$scale)) "cloglog" else run$scale)$inv
  md <- run$models |>
    dplyr::group_by(registry) |>
    dplyr::mutate(model = factor(model, levels = rev(model[order(-n_total)]))) |>
    dplyr::ungroup() |>
    dplyr::mutate(lo = dplyr::coalesce(lcl, inv(yi - 1.96 * sei)),
                  hi = dplyr::coalesce(ucl, inv(yi + 1.96 * sei)))
  s <- run$summary
  ggplot2::ggplot(md, ggplot2::aes(y = model)) +
    ggplot2::geom_rect(data = s, ggplot2::aes(xmin = p_re_lb, xmax = p_re_ub, ymin = -Inf, ymax = Inf),
                       inherit.aes = FALSE, fill = LOO_REF_COLOURS[["re"]], alpha = 0.12) +
    ggplot2::geom_vline(data = s, ggplot2::aes(xintercept = pi_lb), colour = LOO_REF_COLOURS[["re"]], linetype = "13") +
    ggplot2::geom_vline(data = s, ggplot2::aes(xintercept = pi_ub), colour = LOO_REF_COLOURS[["re"]], linetype = "13") +
    ggplot2::geom_vline(data = s, ggplot2::aes(xintercept = p_re), colour = LOO_REF_COLOURS[["re"]], linewidth = 0.7) +
    ggplot2::geom_vline(data = s, ggplot2::aes(xintercept = p_pub), colour = LOO_REF_COLOURS[["pub"]], linewidth = 0.7) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0.25, colour = "grey40") +
    ggplot2::geom_point(ggplot2::aes(x = p, size = n_total, fill = influential), shape = 22, colour = "grey15") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#D95F02", `FALSE` = "grey70"),
                               labels = c(`TRUE` = "influential (metafor rule)", `FALSE` = "not influential"), name = NULL) +
    ggplot2::scale_size_area(max_size = 7, labels = scales::comma, name = "n") +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::facet_wrap(ggplot2::vars(registry), scales = "free") +
    ggplot2::labs(x = paste0(run$time_yr, "-year cumulative revision (95% CI)"), y = NULL, ...) +
    theme_funnel()
}

#' Leave-one-out plot: the registry mean recomputed without each model, by both methods.
#' methods: which leave-one-out means to draw, "loo_sub" (subtraction) and/or "loo_re" (metafor).
plot_loo_means <- function(run, methods = c("loo_sub", "loo_re"), ...) {
  md <- run$models |>
    dplyr::group_by(registry) |>
    dplyr::mutate(model = factor(model, levels = rev(model[order(-n_total)]))) |>
    dplyr::ungroup()
  s <- run$summary
  ggplot2::ggplot(md, ggplot2::aes(y = model)) +
    ggplot2::geom_rect(data = s, ggplot2::aes(xmin = p_lcl, xmax = p_ucl, ymin = -Inf, ymax = Inf),
                       inherit.aes = FALSE, fill = "grey50", alpha = 0.15) +
    ggplot2::geom_vline(data = s, ggplot2::aes(xintercept = p_pub), colour = LOO_REF_COLOURS[["pub"]], linewidth = 0.7) +
    ggplot2::geom_vline(data = s, ggplot2::aes(xintercept = p_re), colour = LOO_REF_COLOURS[["re"]], linewidth = 0.7, linetype = "22") +
    (if ("loo_re" %in% methods) list(
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = loo_re_lb, xmax = loo_re_ub), height = 0.2,
                              colour = LOO_REF_COLOURS[["loo_re"]], alpha = 0.6),
      ggplot2::geom_point(ggplot2::aes(x = loo_re, colour = "loo_re"), size = 2.6, shape = 17))) +
    (if ("loo_sub" %in% methods) ggplot2::geom_point(ggplot2::aes(x = loo_sub, colour = "loo_sub"), size = 2.6)) +
    ggplot2::scale_colour_manual(values = LOO_REF_COLOURS[methods], labels = LOO_REFS[methods],
                                 name = "Mean without this model") +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    ggplot2::facet_wrap(ggplot2::vars(registry), scales = "free") +
    ggplot2::labs(x = "Registry mean recomputed without the model on the y axis", y = "Model left out", ...) +
    theme_funnel()
}

#' Influence: Cook's distance per model, with metafor's flag.
plot_loo_cook <- function(run, ...) {
  md <- run$models |>
    dplyr::group_by(registry) |>
    dplyr::mutate(model = factor(model, levels = rev(model[order(-n_total)]))) |>
    dplyr::ungroup()
  ggplot2::ggplot(md, ggplot2::aes(x = cook_d, y = model, fill = influential)) +
    ggplot2::geom_col(width = 0.7, colour = "grey20", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = stats::qchisq(0.5, 1), linetype = "22") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#D95F02", `FALSE` = "grey70"),
                               labels = c(`TRUE` = "influential", `FALSE` = "not influential"), name = NULL) +
    ggplot2::facet_wrap(ggplot2::vars(registry), scales = "free") +
    ggplot2::labs(x = "Cook's distance (dashed: metafor cut-off, median of chi-square 1 df = 0.45)", y = NULL, ...) +
    theme_funnel()
}

#' Status of every model under every reference: one panel per reference x registry.
#' Points at p - p0 (own reference); the bracket at each point is that model's own control limit
#' (exact by default), so a point outside its bracket is flagged. Panels with one shared reference
#' (published, random effects, lowest leave-one-out) also show the full funnel.
plot_loo_status <- function(st, run, x_log = TRUE, label = TRUE, ...) {
  shared <- st |>
    dplyr::filter(ref %in% c("pub", "re", "min_loo", "min_loo_re")) |>
    dplyr::distinct(registry, ref, p0, level, method)
  nmax <- max(st$n_total) * 1.4
  curves <- purrr::pmap_dfr(shared, function(registry, ref, p0, level, method) {
    funnel_limits(p0, n_grid(10, nmax), level, method) |>
      dplyr::mutate(registry = registry, ref = ref)
  }) |>
    dplyr::mutate(ref = factor(ref, levels = names(LOO_REFS)))
  xr <- range(st$n_total) * c(0.6, 1.5)
  curves <- dplyr::filter(curves, n >= xr[1])
  lab <- if (label) dplyr::filter(st, status != "within") else st[0, ]
  p <- ggplot2::ggplot(st, ggplot2::aes(x = n_total)) +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = n, y = dev_upper), colour = "grey55") +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = n, y = dev_lower), colour = "grey55") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = dev_lower, ymax = dev_upper), width = 0.06, colour = "grey30", linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(y = dev, fill = status), shape = 21, size = 2.5, colour = "grey15", stroke = 0.3) +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(y = dev, label = model), size = 2.6, max.overlaps = Inf,
                             seed = 1, min.segment.length = 0) +
    ggplot2::scale_fill_manual(values = LOO_STATUS_COLOURS, drop = FALSE, name = "Status") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    (if (x_log) ggplot2::scale_x_log10(labels = scales::comma) else ggplot2::scale_x_continuous(labels = scales::comma)) +
    ggplot2::coord_cartesian(xlim = xr) +
    ggplot2::facet_grid(ggplot2::vars(ref), ggplot2::vars(registry), scales = "free",
                        labeller = ggplot2::labeller(ref = ggplot2::as_labeller(stats::setNames(
                          stringr::str_wrap(LOO_REFS, 18), names(LOO_REFS))))) +
    ggplot2::labs(x = "Procedure volume (n)", y = "Model rate minus its reference", ...) +
    theme_funnel()
  p
}

#' Small-volume limits: published-mean funnel with normal vs exact limits, points at p - p_pub.
plot_limit_methods <- function(run, level = 0.998, x_log = TRUE, ...) {
  s <- run$summary
  nmax <- max(run$models$n_total) * 1.4
  xr <- range(run$models$n_total) * c(0.6, 1.5)
  curves <- purrr::pmap_dfr(dplyr::select(s, registry, p_pub), function(registry, p_pub) {
    dplyr::bind_rows(
      dplyr::mutate(funnel_limits(p_pub, n_grid(10, nmax), level, "normal"), lim = "normal"),
      dplyr::mutate(funnel_limits(p_pub, n_grid(10, nmax), level, "exact"),  lim = "exact")) |>
      dplyr::mutate(registry = registry)
  }) |> dplyr::filter(n >= xr[1])
  pts <- run$models |> dplyr::mutate(dev = p - p_pub_of(registry, s))
  ggplot2::ggplot() +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = n, y = dev_upper, linetype = lim, colour = lim)) +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = n, y = dev_lower, linetype = lim, colour = lim)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = n_total, y = dev), size = 2.2, colour = "grey20") +
    ggrepel::geom_text_repel(data = pts, ggplot2::aes(x = n_total, y = dev, label = model), size = 2.5,
                             max.overlaps = Inf, seed = 1, min.segment.length = 0, colour = "grey30") +
    ggplot2::scale_linetype_manual(values = c(normal = "22", exact = "solid"), name = "Limits") +
    ggplot2::scale_colour_manual(values = c(normal = "#D95F02", exact = "#222222"), name = "Limits") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    (if (x_log) ggplot2::scale_x_log10(labels = scales::comma) else ggplot2::scale_x_continuous(labels = scales::comma)) +
    ggplot2::coord_cartesian(xlim = xr) +
    ggplot2::facet_wrap(ggplot2::vars(registry), scales = "free") +
    ggplot2::labs(x = "Procedure volume (n)", y = "Model rate minus published mean", ...) +
    theme_funnel()
}
p_pub_of <- function(registry, s) s$p_pub[match(registry, s$registry)]

#' Meta-analytic test as a flat "z funnel": externally standardised residual vs n.
plot_rstudent <- function(rs, level = 0.998, x_log = TRUE, ...) {
  zc <- stats::qnorm(1 - (1 - level) / 2); z95 <- stats::qnorm(0.975)
  ggplot2::ggplot(rs, ggplot2::aes(x = n_total, y = rstudent)) +
    ggplot2::geom_hline(yintercept = c(-zc, zc), linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = c(-z95, z95), linewidth = 0.4, linetype = "22") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(fill = status), shape = 21, size = 2.6, colour = "grey15", stroke = 0.3) +
    ggrepel::geom_text_repel(ggplot2::aes(label = model), size = 2.5, max.overlaps = Inf, seed = 1,
                             min.segment.length = 0, colour = "grey30") +
    ggplot2::scale_fill_manual(values = LOO_STATUS_COLOURS, drop = FALSE, name = "Status") +
    (if (x_log) ggplot2::scale_x_log10(labels = scales::comma) else ggplot2::scale_x_continuous(labels = scales::comma)) +
    ggplot2::facet_wrap(ggplot2::vars(registry), scales = "free_x") +
    ggplot2::labs(x = "Procedure volume (n)",
                  y = "Leave-one-out standardised residual (z)", ...) +
    theme_funnel()
}

# lower control limits for screening ------------------------------------------------------------------

#' Power: probability that a model whose true rate is `ratio` x the reference is above the upper limit,
#' for each level and n (binomial; the limits use the reference p0).
limit_power <- function(p0, ratio, n, levels, method = "exact") {
  tidyr::crossing(n = n, level = levels, ratio = ratio) |>
    dplyr::mutate(upper = purrr::pmap_dbl(list(n, level), function(nn, lv) funnel_limits(p0, nn, lv, method)$upper),
                  p_true = pmin(p0 * ratio, 0.999),
                  power = 1 - stats::pbinom(floor(upper * n + 1e-9), n, p_true),
                  false_alarm = 1 - stats::pbinom(floor(upper * n + 1e-9), n, p0))
}

#' Highest level whose upper (or lower) limit each model is outside, for several references.
screen_status <- function(run, levels, refs = c("pub", "min_loo_re"), method = "exact") {
  lv <- sort(levels, decreasing = TRUE)
  st <- purrr::map_dfr(lv, function(l) loo_status(run, l, method)) |>
    dplyr::filter(ref %in% refs)
  st |>
    dplyr::group_by(joint, registry, model, ref) |>
    dplyr::summarise(n_total = dplyr::first(n_total), p = dplyr::first(p), p0 = dplyr::first(p0), dev = dplyr::first(dev),
                     above = if (any(status == "above")) max(level[status == "above"]) else NA_real_,
                     below = if (any(status == "below")) max(level[status == "below"]) else NA_real_,
                     .groups = "drop") |>
    dplyr::mutate(screen = dplyr::case_when(
                    !is.na(above) ~ paste("above", level_label(above)),
                    !is.na(below) ~ paste("below", level_label(below)),
                    TRUE ~ paste("within", level_label(min(levels)))),
                  ref = factor(ref, levels = names(LOO_REFS)))
}

#' Funnel plots with several (lower) control limits, one row per reference and one column per registry.
#' Points are coloured by the highest level whose limit they are outside.
plot_screen_funnels <- function(ss, run, levels, method = "exact", x_log = TRUE, ...) {
  refs <- ss |> dplyr::distinct(registry, ref, p0)
  # shared references only: one p0 per registry x ref
  stopifnot(!any(duplicated(dplyr::select(refs, registry, ref))))
  nmax <- max(ss$n_total) * 1.4
  xr <- range(ss$n_total) * c(0.6, 1.5)
  curves <- purrr::pmap_dfr(refs, function(registry, ref, p0) {
    funnel_limits(p0, n_grid(10, nmax), levels, method) |> dplyr::mutate(registry = registry, ref = ref)
  }) |>
    dplyr::filter(n >= xr[1]) |>
    dplyr::mutate(ref = factor(ref, levels = names(LOO_REFS)), lvl = factor(level_label(level), levels = level_label(sort(levels))))
  lv_sorted <- sort(levels, decreasing = TRUE)
  cats <- c(paste("above", level_label(lv_sorted)), paste("within", level_label(min(levels))),
            paste("below", level_label(rev(rev(lv_sorted)))))
  reds <- grDevices::colorRampPalette(c("#B2182B", "#F4A582"))(length(levels))
  blues <- grDevices::colorRampPalette(c("#2166AC", "#92C5DE"))(length(levels))
  cols <- stats::setNames(c(reds, "#9E9E9E", blues), cats)
  ltys <- stats::setNames(c("solid", "42", "13", "1111")[seq_along(levels)], level_label(sort(levels)))
  ss <- dplyr::mutate(ss, screen = factor(screen, levels = cats))
  lab <- dplyr::filter(ss, !startsWith(as.character(screen), "within"))
  ggplot2::ggplot() +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = n, y = dev_upper, linetype = lvl, group = interaction(registry, ref, lvl))) +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = n, y = dev_lower, linetype = lvl, group = interaction(registry, ref, lvl))) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    ggplot2::geom_point(data = ss, ggplot2::aes(x = n_total, y = dev, fill = screen), shape = 21, size = 2.7,
                        colour = "grey15", stroke = 0.3) +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(x = n_total, y = dev, label = model), size = 2.6,
                             max.overlaps = Inf, seed = 1, min.segment.length = 0) +
    ggplot2::scale_linetype_manual(values = ltys, name = paste0("Control limit (", method, ")")) +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE, name = "Model") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    (if (x_log) ggplot2::scale_x_log10(labels = scales::comma) else ggplot2::scale_x_continuous(labels = scales::comma)) +
    ggplot2::coord_cartesian(xlim = xr) +
    ggplot2::facet_grid(ggplot2::vars(ref), ggplot2::vars(registry), scales = "free",
                        labeller = ggplot2::labeller(ref = ggplot2::as_labeller(stats::setNames(
                          stringr::str_wrap(LOO_REFS, 20), names(LOO_REFS))))) +
    ggplot2::labs(x = "Procedure volume (n)", y = "Model rate minus reference", ...) +
    theme_funnel()
}
