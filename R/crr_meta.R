# Cumulative revision rate (CRR) meta-analysis across registries -------------------------------------
#
# Used by crr_meta_analysis.Rmd. Adapts the approach of Combescure et al. (2025, EFORT Open Rev
# 10:277-285, doi:10.1530/EOR-2024-0020) from hip cup/stem combinations to knee models:
#   * each registry row of a model (a "variant", e.g. NJR single peg vs twin peg Oxford) is one estimate;
#   * revision-free survival is put on the complementary log-log scale and the SE is backed out of the
#     published 95% CI;
#   * estimates are pooled per follow-up year with a random-effects model (REML), variant nested in registry;
#   * heterogeneity: Cochran's Q; sources: meta-regression on fixation + registry (exp(coef) = HR).
# Added here (not in the paper): t-based CIs with containment df, a single model across all follow-up
# years (addresses the paper's own stated limitation), and sensitivity analyses.

# Transformations --------------------------------------------------------------------------------------
#' CRR in percent -> complementary log-log of revision-free survival: log(-log(S)), S = 1 - CRR/100.
#' Increasing in CRR; equals log of the cumulative hazard H(t) = -log S(t).
cll_from_crr <- function(crr_pct) log(-log(1 - crr_pct / 100))
#' Inverse: cloglog value -> CRR in percent.
crr_from_cll <- function(y) 100 * (1 - exp(-exp(y)))

#' Model rows for the meta-analysis.
#' `device` = load_extraction()$device; returns one row per registry row x follow-up year with an estimate
#' and a usable 95% CI, with yi / vi on the cloglog scale and attributes from config/device_attributes.csv.
crr_inputs <- function(device, families, rules, models, attribute_source = "design") {
  d <- assign_family(device, families)
  d$label <- short_device_label(d)
  d <- dplyr::filter(d, family %in% models, row_type == "device")
  bad <- with(d, is.na(estimate) | is.na(lcl) | is.na(ucl) | lcl <= 0 | ucl >= 100 |
                 !(lcl <= estimate & estimate <= ucl))
  dropped <- d[bad, ]
  d <- d[!bad, ]
  d <- derive_attributes(d, rules, attribute_source)
  z <- stats::qnorm(0.975)
  out <- d |>
    dplyr::mutate(
      model = family,
      lateral = grepl("lateral", device_label, ignore.case = TRUE) | compartment %in% "lateral",
      # SIRIS reports cemented and hybrid together; kept as its own level in the meta-regression and
      # counted with "cemented" in the fixation subgroups (most Oxford "cemented/hybrid" are fully cemented)
      fixation_group = dplyr::case_when(fixation %in% c("cemented", "cemented/hybrid") ~ "Cemented",
                                        fixation == "cementless" ~ "Cementless",
                                        fixation == "hybrid" ~ "Hybrid", TRUE ~ NA_character_),
      yi = cll_from_crr(estimate),
      y_lcl = cll_from_crr(lcl), y_ucl = cll_from_crr(ucl),
      sei = (y_ucl - y_lcl) / (2 * z),
      vi = sei^2,
      # symmetry of the published CI on the cloglog scale: upper half-width / lower half-width
      # (1 = symmetric, as assumed when the SE is backed out of the CI)
      ci_asym = (y_ucl - yi) / (yi - y_lcl),
      row_id = paste0(registry, ": ", label),
      registry_model = paste(registry, model, sep = " | "))
  attr(out, "dropped") <- dropped
  out
}

# Brand-level collapse (detail level 1) -------------------------------------------------------------------
#' Collapse rows to one estimate, weighting by procedures as detail level 1 does:
#'   p = sum(w_i p_i), w_i = n_i / sum(n),  Var(p) = sum(w_i^2 Var(p_i))   (rows are disjoint patient groups)
#' Var(p_i) is the row's cloglog variance mapped back to the CRR scale (delta method). Returned on the
#' cloglog scale (y, v) and as CRR (%). Inverse-variance weighting on the cloglog scale is NOT used: there an
#' estimate's precision tracks its number of revisions, so high-revision variants would get more weight
#' and pull the brand rate up.
crr_nw_collapse <- function(estimate, vi, n) {
  p <- estimate / 100
  d <- 1 / ((1 - p) * (-log(1 - p)))                  # d cloglog / dp
  var_p <- vi / d^2
  w <- n / sum(n)
  pp <- sum(w * p); vp <- sum(w^2 * var_p)
  dd <- 1 / ((1 - pp) * (-log(1 - pp)))
  list(y = log(-log(1 - pp)), v = vp * dd^2, crr = 100 * pp)
}

# Per-follow-up-year random-effects models ----------------------------------------------------------------
#' Fit one pooled estimate. Random effects: registry, and variant (row) within registry; a level is dropped
#' when it cannot be estimated (one registry, or one row per registry), as in Combescure et al.
#' With several models (an "all models" estimate) a model level is added on top.
#' test = "t": t-distribution with containment df (dfs = "contain"), recommended with few registries;
#' test = "z": normal CI, as in the paper.
crr_fit <- function(dd, test = c("t", "z")) {
  test <- match.arg(test)
  k_reg <- dplyr::n_distinct(dd$registry); k <- nrow(dd); k_mod <- dplyr::n_distinct(dd$model)
  if (k == 0) return(NULL)
  if (k == 1) return(list(fit = NULL, k = 1, k_reg = 1, k_mod = 1, b = dd$yi, ci = dd$yi + c(-1, 1) * 1.96 * dd$sei,
                          structure = "single estimate"))
  multi_row <- any(table(dd$registry_model) > 1)
  random <- list()
  if (k_mod > 1) random <- c(random, ~ 1 | model)
  if (k_reg > 1) random <- c(random, if (k_mod > 1) ~ 1 | registry_model else ~ 1 | registry)
  if (multi_row || length(random) == 0) random <- c(random, ~ 1 | row_id)
  args <- list(yi = dd$yi, V = dd$vi, random = random, data = dd, method = "REML")
  if (test == "t") args <- c(args, test = "t", dfs = "contain")
  fit <- tryCatch(do.call(metafor::rma.mv, args), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  structure <- paste(vapply(random, function(f) sub(".*\\| ", "", deparse(f)), ""), collapse = " / ")
  pr <- stats::predict(fit)
  list(fit = fit, k = k, k_reg = k_reg, k_mod = k_mod, b = as.numeric(fit$b), ci = c(fit$ci.lb, fit$ci.ub),
       pi = c(pr$pi.lb, pr$pi.ub), structure = structure)
}

#' Heterogeneity summary of a crr_fit(): Cochran's Q, and a multilevel I^2 (share of total variance
#' due to each random-effect level; typical within-study variance as in Higgins & Thompson / Nakagawa).
crr_het <- function(cf) {
  f <- cf$fit
  if (is.null(f)) return(tibble::tibble(Q = NA_real_, Q_df = NA_real_, Q_p = NA_real_, tau_total = NA_real_,
                                        I2_total = NA_real_, sigma2_txt = NA_character_))
  W <- diag(1 / f$vi); X <- stats::model.matrix(f)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  v_typ <- (f$k - f$p) / sum(diag(P))
  s2 <- f$sigma2
  tibble::tibble(Q = f$QE, Q_df = f$k - f$p, Q_p = f$QEp, tau_total = sqrt(sum(s2)),
                 I2_total = 100 * sum(s2) / (sum(s2) + v_typ),
                 sigma2_txt = paste(sprintf("%s %.3f", f$s.names, s2), collapse = "; "))
}

#' Pooled CRR per group x follow-up year. `groups` is a named list of filters (e.g. All, Cemented ...).
#' A cross-registry pooled estimate needs >= 2 registries; with fewer than `k_reg_min` registries the
#' z-based CI is used and the row is flagged (a t CI with 1 df is uninformative).
crr_pool_table <- function(inp, groups, time_yrs, k_reg_min = 3, test = "t") {
  purrr::map_dfr(names(groups), function(g) {
    dg <- dplyr::filter(inp, !!groups[[g]])
    purrr::map_dfr(time_yrs, function(t) {
      dd <- dplyr::filter(dg, time_yr == t)
      k_reg <- dplyr::n_distinct(dd$registry)
      if (k_reg < 2) return(NULL)
      use_test <- if (k_reg >= k_reg_min) test else "z"
      cf <- crr_fit(dd, use_test)
      if (is.null(cf)) return(NULL)
      dplyr::bind_cols(
        tibble::tibble(group = g, time_yr = t, k = cf$k, k_reg = k_reg, n_total = sum(dplyr::distinct(dd, row_id, n_total)$n_total),
                       registries = paste(sort(unique(dd$registry)), collapse = ", "),
                       structure = cf$structure, ci_method = if (use_test == "t") "t (contain df)" else "z",
                       exploratory = k_reg < k_reg_min,
                       crr = crr_from_cll(cf$b), crr_lcl = crr_from_cll(cf$ci[1]), crr_ucl = crr_from_cll(cf$ci[2]),
                       pi_lcl = crr_from_cll(cf$pi[1]), pi_ucl = crr_from_cll(cf$pi[2])),
        crr_het(cf))
    })
  })
}

#' Leave-one-registry-out pooled CRR (group "All") at each time point.
crr_loro <- function(inp, time_yrs, test = "t") {
  purrr::map_dfr(time_yrs, function(t) {
    dd <- dplyr::filter(inp, time_yr == t)
    regs <- sort(unique(dd$registry))
    if (length(regs) < 3) return(NULL)
    purrr::map_dfr(regs, function(r) {
      d2 <- dplyr::filter(dd, registry != r)
      cf <- crr_fit(d2, if (dplyr::n_distinct(d2$registry) >= 3) test else "z")
      tibble::tibble(time_yr = t, left_out = r, crr = crr_from_cll(cf$b),
                     crr_lcl = crr_from_cll(cf$ci[1]), crr_ucl = crr_from_cll(cf$ci[2]))
    })
  })
}

#' Sensitivity: collapse each registry's rows to one n-weighted estimate first (crr_nw_collapse(), detail
#' level 1's one-point-per-brand with a CI), then pool registries with a two-level REML model (Knapp-Hartung CI).
crr_collapsed <- function(inp, time_yrs) {
  purrr::map_dfr(time_yrs, function(t) {
    dd <- dplyr::filter(inp, time_yr == t)
    if (dplyr::n_distinct(dd$registry) < 2) return(NULL)
    per_reg <- dd |> dplyr::group_by(registry) |>
      dplyr::summarise(cl = list(crr_nw_collapse(estimate, vi, n_total)), .groups = "drop") |>
      dplyr::mutate(yi = purrr::map_dbl(cl, "y"), vi = purrr::map_dbl(cl, "v"))
    k <- nrow(per_reg)
    f <- metafor::rma(yi, vi, data = per_reg, method = "REML", test = if (k >= 3) "knha" else "z")
    tibble::tibble(time_yr = t, k_reg = k, crr = crr_from_cll(f$b[1]), crr_lcl = crr_from_cll(f$ci.lb),
                   crr_ucl = crr_from_cll(f$ci.ub), tau = sqrt(f$tau2))
  })
}

# Meta-regression (Table 3/4 analogue) -----------------------------------------------------------------
#' At one follow-up year: yi ~ fixation + registry, random effect for the row. exp(coef) on the cloglog
#' scale = ratio of cumulative hazards = hazard ratio under proportional hazards. Reference levels are
#' the categories with the largest number of procedures (as in the paper).
crr_metareg <- function(inp, t, test = "t") {
  dd <- dplyr::filter(inp, time_yr == t, !is.na(fixation))
  if (dplyr::n_distinct(dd$registry) < 3) return(NULL)
  ref_of <- function(v) names(sort(tapply(dd$n_total, v, sum), decreasing = TRUE))[1]
  dd$fixation <- stats::relevel(factor(dd$fixation), ref_of(dd$fixation))
  dd$registry <- stats::relevel(factor(dd$registry), ref_of(dd$registry))
  mods <- if (nlevels(droplevels(dd$fixation)) > 1) ~ fixation + registry else ~ registry
  args <- list(yi = dd$yi, V = dd$vi, mods = mods, random = ~ 1 | row_id, data = dd, method = "REML")
  if (test == "t") args <- c(args, test = "t")
  f <- tryCatch(do.call(metafor::rma.mv, args), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  cf <- stats::coef(summary(f))
  joint_p <- function(pat) { i <- grep(pat, rownames(cf)); if (length(i)) stats::anova(f, btt = i)$QMp else NA_real_ }
  fix_p <- joint_p("^fixation"); reg_p <- joint_p("^registry")
  tibble::tibble(time_yr = t, term = rownames(cf), est = cf$estimate, ci.lb = cf$ci.lb, ci.ub = cf$ci.ub,
                 pval = cf$pval) |>
    dplyr::filter(term != "intrcpt") |>
    dplyr::mutate(variable = ifelse(grepl("^fixation", term), "Fixation", "Registry"),
                  level = sub("^(fixation|registry)", "", term),
                  reference = ifelse(variable == "Fixation", levels(dd$fixation)[1], levels(dd$registry)[1]),
                  HR = exp(est), HR_lcl = exp(ci.lb), HR_ucl = exp(ci.ub),
                  QM_p = f$QMp, fixation_p = fix_p, registry_p = reg_p, QE = f$QE, QE_p = f$QEp, sigma2_resid = f$sigma2, k = f$k)
}

# Single model across follow-up years --------------------------------------------------------------------
#' Approximate sampling covariance of one row's cloglog estimates across follow-up years.
#' Cumulative-hazard estimates have (approximately) independent increments, so for t_j <= t_k
#'   Cov(H_j, H_k) ~= Var(H_j)  and, by the delta method, Cov(log H_j, log H_k) ~= Var(H_j) / (H_j H_k),
#' with Var(H_j) = v_j H_j^2 (v_j = variance of log H_j from the CI). Var(H) is forced non-decreasing and
#' each block is projected to the nearest positive-definite matrix when ties make it singular.
crr_row_cov <- function(est_pct, vi) {
  H <- -log(1 - est_pct / 100)
  varH <- cummax(vi * H^2)
  n <- length(H)
  M <- outer(seq_len(n), seq_len(n), function(i, j) varH[pmin(i, j)] / (H[i] * H[j]))
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 1e-10 * max(ev)) M <- as.matrix(Matrix::nearPD(M, keepDiag = TRUE, ensureSymmetry = TRUE)$mat)
  M
}

#' Joint model: all rows and all follow-up years in one rma.mv.
#'   yi = b0 + b1 log(t) + b2 log(t)^2 + u_registry(t) + w_row(t) + e,   e ~ N(0, V)
#' log H(t) linear in log(t) is a Weibull cumulative hazard; the quadratic term lets the hazard shape
#' change. Random intercepts and log(t) slopes per registry and per row (struct "GEN").
#' Repeated identical estimates (KM plateaus carried forward) are dropped - they add no information and
#' make V singular.
crr_joint <- function(inp, slopes = TRUE) {
  dd <- inp |> dplyr::arrange(row_id, time_yr) |> dplyr::group_by(row_id) |>
    dplyr::filter(!(estimate == dplyr::lag(estimate, default = -1) & lcl == dplyr::lag(lcl, default = -1))) |>
    dplyr::ungroup()
  V <- as.matrix(Matrix::bdiag(lapply(split(dd, dd$row_id)[unique(dd$row_id)], function(b) crr_row_cov(b$estimate, b$vi))))
  dd$lt <- log(dd$time_yr); dd$lt2 <- dd$lt^2
  rand <- if (slopes) list(~ lt | registry, ~ lt | row_id) else list(~ 1 | registry, ~ 1 | row_id)
  args <- list(yi = dd$yi, V = V, mods = ~ lt + lt2, random = rand, data = dd, method = "REML",
               test = "t", dfs = "contain")
  if (slopes) args$struct <- c("GEN", "GEN")
  f <- tryCatch(suppressWarnings(do.call(metafor::rma.mv, args)), error = function(e) NULL)
  if (is.null(f) && slopes) return(crr_joint(inp, slopes = FALSE))
  list(fit = f, data = dd, slopes = slopes)
}

crr_joint_predict <- function(jm, t) {
  p <- stats::predict(jm$fit, newmods = cbind(log(t), log(t)^2))
  tibble::tibble(time_yr = t, crr = crr_from_cll(p$pred), crr_lcl = crr_from_cll(p$ci.lb),
                 crr_ucl = crr_from_cll(p$ci.ub))
}

# Formatting ---------------------------------------------------------------------------------------------
fmt_crr <- function(est, lo, hi, d = 2) ifelse(is.na(est), "", sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), est, lo, hi))

# Design variants (finer than fixation; the analogue of the paper's "associated implant") --------------------
#' One label per registry row, shared across registries where the reports allow:
#'   lateral rows -> "Lateral"; a stated non-standard femoral surface (e.g. TiNbN) -> "<material>-coated";
#'   "twin peg" in the label -> "Twin peg (cemented)"; otherwise the fixation ("Cemented", "Cementless",
#'   "Hybrid", SIRIS "Cemented/hybrid"). Registries that do not name the peg design are counted as "Cemented".
crr_design_variant <- function(inp) {
  fx <- dplyr::recode(inp$fixation, cemented = "Cemented", cementless = "Cementless", hybrid = "Hybrid",
                      `cemented/hybrid` = "Cemented/hybrid (SIRIS)", .missing = "Fixation unknown")
  mat_stated <- !is.na(inp$material) & inp$material_source %in% "stated" & inp$material != "CoCr"
  dplyr::case_when(
    inp$lateral ~ "Lateral",
    mat_stated ~ paste0(inp$material, "-coated (", tolower(fx), ")"),
    grepl("twin peg", inp$device_label, ignore.case = TRUE) ~ paste0("Twin peg (", tolower(fx), ")"),
    TRUE ~ fx)
}

#' Table 3/4 analogue: multivariable meta-regression at one follow-up year.
#' yi ~ registry + <factor>, random effect per row (REML). Reference = level with the most procedures.
#' Returns one row per level (reference rows included, HR = 1) with procedures, adjusted HR (95% CI), the
#' level's p-value, and on the reference row the omnibus p of that factor (as in Combescure et al. Tables 3-4).
crr_table34 <- function(inp, t, factor_col = "fixation", test = "t") {
  dd <- dplyr::filter(inp, time_yr == t, !is.na(.data[[factor_col]]))
  if (dplyr::n_distinct(dd$registry) < 3) return(NULL)
  dd$f <- dd[[factor_col]]
  ref_of <- function(v) names(sort(tapply(dd$n_total, v, sum), decreasing = TRUE))[1]
  dd$registry <- stats::relevel(factor(dd$registry), ref_of(dd$registry))
  dd$f <- stats::relevel(factor(dd$f), ref_of(dd$f))
  has_f <- nlevels(droplevels(dd$f)) > 1                # e.g. ZUK: every row cemented -> registry only
  args <- list(yi = dd$yi, V = dd$vi, mods = if (has_f) ~ registry + f else ~ registry, random = ~ 1 | row_id,
               data = dd, method = "REML")
  if (test == "t") args <- c(args, test = "t")
  fit <- tryCatch(do.call(metafor::rma.mv, args), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- stats::coef(summary(fit))
  omni <- function(prefix) { i <- grep(paste0("^", prefix), rownames(cf)); if (length(i)) stats::anova(fit, btt = i)$QMp else NA_real_ }
  block <- function(var, prefix, lab) {
    lv <- levels(dd[[var]])
    n <- tapply(dd$n_total, dd[[var]], sum)[lv]
    i <- match(paste0(prefix, lv), rownames(cf))          # exact match (data.frame row names would partial-match)
    pick <- function(col) ifelse(is.na(i), NA_real_, cf[[col]][i])
    tibble::tibble(time_yr = t, variable = lab, level = lv, n = as.numeric(n), reference = lv == lv[1],
                   HR = ifelse(lv == lv[1], 1, exp(pick("estimate"))), HR_lcl = exp(pick("ci.lb")), HR_ucl = exp(pick("ci.ub")),
                   pval = pick("pval"), omnibus_p = ifelse(lv == lv[1], omni(prefix), NA_real_))
  }
  out <- dplyr::bind_rows(block("registry", "registry", "Registry"),
                          if (has_f) block("f", "f", c(fixation = "Fixation", design_variant = "Design variant")[[factor_col]] %||% factor_col))
  attr(out, "fit") <- fit
  dplyr::mutate(out, k = fit$k, df_resid = fit$k - fit$p, QE = fit$QE, QE_p = fit$QEp, sigma2_resid = fit$sigma2)
}

# One estimate per registry ----------------------------------------------------------------------------------
#' Collapse each registry's rows at each follow-up year to one estimate, n-weighted (crr_nw_collapse(); the
#' detail level 1 point, with a CI). The CI ignores between-variant heterogeneity within the registry, so it
#' is too narrow when the registry's variants differ (Q_within measures how much they do).
crr_registry_collapse <- function(inp) {
  inp |>
    dplyr::group_by(model, registry, time_yr) |>
    dplyr::group_modify(function(d, k) {
      cl <- crr_nw_collapse(d$estimate, d$vi, d$n_total)
      tibble::tibble(n_rows = nrow(d), y = cl$y, v = cl$v, n_total = sum(d$n_total),
                     variants = paste(d$label, collapse = " ; "),
                     Q_within = if (nrow(d) > 1) { m <- sum(d$yi / d$vi) / sum(1 / d$vi); sum((d$yi - m)^2 / d$vi) } else NA_real_)
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(crr = crr_from_cll(y), crr_lcl = crr_from_cll(y - 1.96 * sqrt(v)), crr_ucl = crr_from_cll(y + 1.96 * sqrt(v)))
}

#' Share of the total weight each registry gets in the pooled estimate at year t:
#' 3-level model (rows nested in registry) vs treating every row as an independent study vs collapsed.
crr_registry_weights <- function(inp, t) {
  dd <- dplyr::filter(inp, time_yr == t)
  if (dplyr::n_distinct(dd$registry) < 2) return(NULL)
  w_share <- function(fit, reg) tapply(stats::weights(fit), reg, sum)
  f3 <- crr_fit(dd, "t")$fit
  fi <- metafor::rma(yi, vi, data = dd, method = "REML")
  cl <- dd |> dplyr::group_by(registry) |>
    dplyr::summarise(c = list(crr_nw_collapse(estimate, vi, n_total)), .groups = "drop") |>
    dplyr::mutate(yi = purrr::map_dbl(c, "y"), vi = purrr::map_dbl(c, "v"))
  fc <- metafor::rma(yi, vi, data = cl, method = "REML")
  tibble::tibble(time_yr = t, registry = sort(unique(dd$registry)),
                 rows = as.integer(table(dd$registry)[sort(unique(dd$registry))]),
                 w_independent = as.numeric(w_share(fi, dd$registry)[sort(unique(dd$registry))]),
                 w_nested = as.numeric(w_share(f3, dd$registry)[sort(unique(dd$registry))]),
                 w_collapsed = as.numeric(stats::setNames(stats::weights(fc), cl$registry)[sort(unique(dd$registry))]))
}
