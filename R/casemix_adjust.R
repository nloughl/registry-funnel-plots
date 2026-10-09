# Case-mix adjustment: shifted funnels and model point ranges ---------------------------------------
#
# Used by casemix_adjustment.Rmd. Patient case mix = age-sex class (AOANJRR KP22; NJR 3.K6 with its
# cemented and uncemented/hybrid blocks pooled per class, n-weighted, so fixation - a model attribute -
# is not treated as case mix).
#
# For each registry three references are compared at follow-up time t:
#   mean  the registry mean p_mean (as in the main report)
#   low   the rate of the LOWEST-risk age-sex class  p_low   (most favourable patients)
#   high  the rate of the HIGHEST-risk age-sex class p_high  (least favourable patients)
#
# Option A, "adjusted funnels": points stay at p - p_mean; each reference's funnel is centred on
#   (p_ref - p_mean) with the width of p_ref, i.e. limits(p_ref, n) - p_mean.
# Option B, "point ranges": funnels stay centred on 0 (width of p_ref only); each model gets one
#   deviation per reference, p - p_ref, drawn as three markers or as a band around the main point.
# The two are the same comparison drawn differently:  p > upper(p_ref, n)  <=>  p - p_ref > upper(p_ref, n) - p_ref.
# classify_casemix() is shared, so both options flag exactly the same models.

CM_REFS <- c(mean = "Registry mean", low = "Lowest-risk class", high = "Highest-risk class")
CM_REF_COLOURS <- c(mean = "#222222", low = "#2166AC", high = "#B2182B")

CM_STATUS <- c(
  robust_high = "High under every reference",
  sens_high   = "High vs mean, not vs highest-risk class",
  cond_high   = "High only vs lowest-risk class",
  within      = "Within under every reference",
  cond_low    = "Low only vs highest-risk class",
  sens_low    = "Low vs mean, not vs lowest-risk class",
  robust_low  = "Low under every reference")
CM_STATUS_COLOURS <- c(
  "High under every reference"              = "#B2182B",
  "High vs mean, not vs highest-risk class" = "#EF8A62",
  "High only vs lowest-risk class"          = "#F4C7A1",
  "Within under every reference"            = "#9E9E9E",
  "Low only vs highest-risk class"          = "#B8D8EA",
  "Low vs mean, not vs lowest-risk class"   = "#67A9CF",
  "Low under every reference"               = "#2166AC")

#' Age-sex class rates per registry at time_yr (NJR design blocks pooled per class, n-weighted).
casemix_classes <- function(casemix, time_yr, min_n = 0) {
  casemix_strata(casemix, time_yr, njr_level = "group") |>
    dplyr::group_by(registry, sex, age_group) |>
    dplyr::summarise(p = sum(p * n_total) / sum(n_total), n_total = sum(n_total),
                     blocks = paste(unique(design), collapse = " + "),
                     table_id = paste(unique(table_id), collapse = "/"), .groups = "drop") |>
    dplyr::filter(n_total >= min_n) |>
    dplyr::mutate(class = paste(age_group, sex))
}

#' Lowest- and highest-risk class per registry, with the registry mean.
casemix_extremes <- function(classes, means) {
  classes |>
    dplyr::group_by(registry) |>
    dplyr::summarise(
      n_classes = dplyr::n(),
      low_class = class[which.min(p)], p_low = min(p), n_low = n_total[which.min(p)],
      high_class = class[which.max(p)], p_high = max(p), n_high = n_total[which.max(p)],
      p_classes_pooled = sum(p * n_total) / sum(n_total), n_classes_total = sum(n_total),
      table_id = dplyr::first(table_id), .groups = "drop") |>
    dplyr::inner_join(dplyr::select(means, registry, p_mean, mean_n = n_total), by = "registry")
}

# long: one row per registry x reference (mean / low / high)
cm_ref_long <- function(ext) {
  tibble::tibble(
    registry = rep(ext$registry, 3),
    ref = rep(names(CM_REFS), each = nrow(ext)),
    p_ref = c(ext$p_mean, ext$p_low, ext$p_high),
    p_mean = rep(ext$p_mean, 3))
}

#' Limit curves for the three references.
#' shifted = TRUE  (option A): centred on p_ref - p_mean, so dev = limit - p_mean
#' shifted = FALSE (option B): centred on 0, so dev = limit - p_ref (width only)
cm_ref_curves <- function(ext, levels, method = "normal", n_max = 1e5, shifted = TRUE) {
  refs <- cm_ref_long(ext) |>
    dplyr::mutate(p_centre = if (shifted) p_mean else p_ref, n_max = n_max) |>
    dplyr::select(registry, ref, p_ref, p_centre, n_max)
  build_limit_curves(refs, levels, method) |>
    dplyr::mutate(ref = factor(ref, levels = names(CM_REFS)))
}

#' Envelope over ALL age-sex classes (option A: shifted to p_class - p_mean).
cm_class_envelope <- function(classes, means, levels, method = "normal", n_max = 1e5, shifted = TRUE) {
  refs <- classes |>
    dplyr::inner_join(dplyr::select(means, registry, p_mean), by = "registry") |>
    dplyr::transmute(registry, class, p_ref = p, p_centre = if (shifted) p_mean else p, n_max = n_max)
  build_limit_curves(refs, levels, method) |>
    dplyr::group_by(registry, level, n) |>
    dplyr::summarise(up_lo = min(dev_upper), up_hi = max(dev_upper),
                     lo_lo = min(dev_lower), lo_hi = max(dev_lower), .groups = "drop")
}

#' Model points with one deviation per reference and their status under each.
#' d_mean = p - p_mean, d_low = p - p_low (largest), d_high = p - p_high (smallest).
classify_casemix <- function(points, ext, level = 0.998, method = "normal") {
  lim <- function(p0, n) funnel_limits(p0, n, level, method)
  pts <- points |>
    dplyr::select(-dplyr::any_of(c("p_low", "p_high", "low_class", "high_class"))) |>
    dplyr::inner_join(dplyr::select(ext, registry, p_low, p_high, low_class, high_class), by = "registry")
  bounds <- purrr::pmap_dfr(list(pts$p_mean, pts$p_low, pts$p_high, pts$n_total), function(pm, pl, ph, n) {
    a <- lim(pm, n); b <- lim(pl, n); c <- lim(ph, n)
    tibble::tibble(up_mean = a$upper, lo_mean = a$lower, up_low = b$upper, lo_low = b$lower,
                   up_high = c$upper, lo_high = c$lower)
  })
  dplyr::bind_cols(pts, bounds) |>
    dplyr::mutate(
      d_mean = p - p_mean, d_low = p - p_low, d_high = p - p_high,
      above_mean = p > up_mean, above_low = p > up_low, above_high = p > up_high,
      below_mean = p < lo_mean, below_low = p < lo_low, below_high = p < lo_high,
      cm_key = dplyr::case_when(
        above_high ~ "robust_high",
        above_mean ~ "sens_high",
        above_low  ~ "cond_high",
        below_low  ~ "robust_low",
        below_mean ~ "sens_low",
        below_high ~ "cond_low",
        TRUE ~ "within"),
      cm_status = factor(CM_STATUS[cm_key], levels = CM_STATUS),
      excess_mean = ifelse(above_mean, n_total * d_mean, 0),
      excess_high = ifelse(above_high, n_total * d_high, 0),
      excess_low  = ifelse(above_low,  n_total * d_low,  0))
}

#' Per-registry counts by status, and excess revisions of the flagged models under each reference.
casemix_summary <- function(cl) {
  cl |>
    dplyr::group_by(registry) |>
    dplyr::summarise(
      models = dplyr::n(),
      high_vs_lowest = sum(above_low), high_vs_mean = sum(above_mean), high_vs_highest = sum(above_high),
      low_vs_highest = sum(below_high), low_vs_mean = sum(below_mean), low_vs_lowest = sum(below_low),
      excess_vs_lowest = sum(excess_low), excess_vs_mean = sum(excess_mean), excess_vs_highest = sum(excess_high),
      .groups = "drop")
}

# shared plot pieces ------------------------------------------------------------------------------

cm_y_range <- function(cl, curves, n_lo, pad = 0.1) {
  y <- c(cl$d_low, cl$d_high, cl$d_mean)
  at <- dplyr::filter(curves, n >= n_lo) |> dplyr::filter(n == min(n))
  y <- c(y, at$dev_upper, at$dev_lower)
  r <- range(y, na.rm = TRUE)
  r + c(-1, 1) * diff(r) * pad
}

cm_x_range <- function(cl) range(cl$n_total, na.rm = TRUE) * c(0.7, 1.3)

# "flagged": models that are high under any reference, or low under every reference
cm_label_data <- function(cl, label_points) switch(label_points,
  all = cl, none = cl[0, ], cl[cl$cm_key %in% c("robust_high", "sens_high", "cond_high", "robust_low"), ])

cm_scales <- function(x_log) list(
  ggplot2::scale_linetype_manual(values = c(mean = "solid", low = "22", high = "22"), labels = CM_REFS,
                                 name = "Funnel reference"),
  ggplot2::scale_colour_manual(values = CM_REF_COLOURS, labels = CM_REFS, name = "Funnel reference"),
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), breaks = scales::pretty_breaks(8)),
  if (x_log) ggplot2::scale_x_log10(labels = scales::comma) else ggplot2::scale_x_continuous(labels = scales::comma))

#' Option A: shifted funnels (mean / lowest-risk / highest-risk), shaded envelope over all classes,
#' model points at p - p_mean coloured by case-mix status.
plot_cm_adjusted <- function(cl, curves, envelope, level, label_points = "flagged", x_log = TRUE,
                             facet = TRUE, ...) {
  cv <- dplyr::filter(curves, .data$level == !!level, registry %in% cl$registry)
  ev <- dplyr::filter(envelope, .data$level == !!level, registry %in% cl$registry)
  xr <- cm_x_range(cl)
  yr <- cm_y_range(cl, cv, xr[1])
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = ev, ggplot2::aes(x = n, ymin = up_lo, ymax = up_hi), fill = "#B2182B", alpha = 0.10) +
    ggplot2::geom_ribbon(data = ev, ggplot2::aes(x = n, ymin = lo_lo, ymax = lo_hi), fill = "#2166AC", alpha = 0.10) +
    ggplot2::geom_line(data = cv, ggplot2::aes(x = n, y = dev_upper, colour = ref, linetype = ref,
                                               group = interaction(registry, ref)), linewidth = 0.7) +
    ggplot2::geom_line(data = cv, ggplot2::aes(x = n, y = dev_lower, colour = ref, linetype = ref,
                                               group = interaction(registry, ref)), linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    ggplot2::geom_point(data = cl, ggplot2::aes(x = n_total, y = d_mean, fill = cm_status),
                        shape = 21, size = 2.8, colour = "grey15", stroke = 0.4) +
    ggrepel::geom_text_repel(data = cm_label_data(cl, label_points),
                             ggplot2::aes(x = n_total, y = d_mean, label = label), size = 2.7,
                             max.overlaps = Inf, seed = 1, min.segment.length = 0) +
    ggplot2::scale_fill_manual(values = CM_STATUS_COLOURS, drop = FALSE, name = "Model status") +
    cm_scales(x_log) +
    ggplot2::coord_cartesian(xlim = xr, ylim = yr) +
    ggplot2::labs(x = "Procedure volume (n)", y = "Model rate minus registry mean", ...) +
    theme_funnel() +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 3.5), order = 2),
                    colour = ggplot2::guide_legend(order = 1), linetype = ggplot2::guide_legend(order = 1))
  if (facet) p <- p + ggplot2::facet_wrap(ggplot2::vars(registry), scales = "free")
  p
}

#' Option A, upper limits only: the three shifted UPPER limits (mean / lowest-risk / highest-risk),
#' no shaded envelope and no lower limits. One-sided screening for poor performers, so the "low"
#' statuses are not shown: any model not above an upper limit is drawn grey ("Not above any upper limit").
#' The y range is built from the drawn layers only (upper limits + points), not from the lower limits.
CM_STATUS_UPPER <- c(CM_STATUS[c("robust_high", "sens_high", "cond_high")], not_high = "Not above any upper limit")
CM_STATUS_UPPER_COLOURS <- c(CM_STATUS_COLOURS[CM_STATUS[c("robust_high", "sens_high", "cond_high")]],
                             "Not above any upper limit" = "#BDBDBD")

plot_cm_adjusted_upper <- function(cl, curves, level, label_points = "flagged", x_log = TRUE,
                                   facet = TRUE, pad = 0.1, ...) {
  high_keys <- c("robust_high", "sens_high", "cond_high")
  cl <- dplyr::mutate(cl,
    up_key    = ifelse(cm_key %in% high_keys, cm_key, "not_high"),
    up_status = factor(CM_STATUS_UPPER[up_key], levels = CM_STATUS_UPPER))
  cv <- dplyr::filter(curves, .data$level == !!level, registry %in% cl$registry)
  xr <- cm_x_range(cl)
  # y range: points (vs mean), 0, and the upper limits across the visible x range
  vis <- dplyr::filter(cv, n >= xr[1], n <= xr[2])
  yr <- range(c(0, cl$d_mean, vis$dev_upper), na.rm = TRUE)
  yr <- yr + c(-1, 1) * diff(yr) * pad
  lab <- switch(label_points, all = cl, none = cl[0, ], cl[cl$up_key %in% high_keys, ])
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    ggplot2::geom_line(data = cv, ggplot2::aes(x = n, y = dev_upper, colour = ref, linetype = ref,
                                               group = interaction(registry, ref)), linewidth = 0.7) +
    ggplot2::geom_point(data = cl, ggplot2::aes(x = n_total, y = d_mean, fill = up_status),
                        shape = 21, size = 2.8, colour = "grey15", stroke = 0.4) +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(x = n_total, y = d_mean, label = label), size = 2.7,
                             max.overlaps = Inf, seed = 1, min.segment.length = 0) +
    ggplot2::scale_fill_manual(values = CM_STATUS_UPPER_COLOURS, drop = FALSE, name = "Model status") +
    cm_scales(x_log) +
    ggplot2::coord_cartesian(xlim = xr, ylim = yr) +
    ggplot2::labs(x = "Procedure volume (n)", y = "Model rate minus registry mean", ...) +
    theme_funnel() +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 3.5), order = 2),
                    colour = ggplot2::guide_legend(order = 1), linetype = ggplot2::guide_legend(order = 1))
  if (facet) p <- p + ggplot2::facet_wrap(ggplot2::vars(registry), scales = "free")
  p
}

#' Option B: funnels centred on 0 (width of each reference only); each model drawn as a range of
#' deviations. style = "markers": three markers (vs mean / lowest-risk / highest-risk) joined by a line;
#'        style = "band":    the main point (vs mean) with a shaded band spanning the other two.
#'        style = "dots":    the main point (black) with a blue dot above (vs lowest-risk class) and a red dot
#'                           below (vs highest-risk class), joined by a thin line; no status colours, the reader
#'                           compares each dot with the funnel line of the same colour.
plot_cm_ranges <- function(cl, curves, level, style = c("markers", "band", "dots"), label_points = "flagged",
                           x_log = TRUE, facet = TRUE, ...) {
  style <- match.arg(style)
  cv <- dplyr::filter(curves, .data$level == !!level, registry %in% cl$registry)
  xr <- cm_x_range(cl)
  yr <- cm_y_range(cl, cv, xr[1])
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = cv, ggplot2::aes(x = n, y = dev_upper, colour = ref, linetype = ref,
                                               group = interaction(registry, ref)), linewidth = 0.6) +
    ggplot2::geom_line(data = cv, ggplot2::aes(x = n, y = dev_lower, colour = ref, linetype = ref,
                                               group = interaction(registry, ref)), linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40")
  if (style == "dots") {
    dt <- cl |>
      dplyr::select(registry, label, n_total, d_mean, d_low, d_high) |>
      tidyr::pivot_longer(c(d_mean, d_low, d_high), names_to = "ref", values_to = "d") |>
      dplyr::mutate(ref = factor(sub("^d_", "", ref), levels = names(CM_REFS)))
    p <- p +
      ggplot2::geom_linerange(data = cl, ggplot2::aes(x = n_total, ymin = d_high, ymax = d_low),
                              colour = "grey60", linewidth = 0.35) +
      ggplot2::geom_point(data = dt, ggplot2::aes(x = n_total, y = d, colour = ref, shape = ref),
                          size = 2.6) +
      ggplot2::scale_shape_manual(values = c(mean = 16, low = 16, high = 16), labels = CM_REFS,
                                  name = "Funnel reference")
  } else if (style == "markers") {
    mk <- cl |>
      dplyr::select(registry, label, n_total, cm_status, d_mean, d_low, d_high) |>
      tidyr::pivot_longer(c(d_mean, d_low, d_high), names_to = "ref", values_to = "d") |>
      dplyr::mutate(ref = factor(sub("^d_", "", ref), levels = names(CM_REFS)))
    p <- p +
      ggplot2::geom_linerange(data = cl, ggplot2::aes(x = n_total, ymin = d_high, ymax = d_low),
                              colour = "grey55", linewidth = 0.4) +
      ggplot2::geom_point(data = mk, ggplot2::aes(x = n_total, y = d, colour = ref, shape = ref), size = 2.4) +
      ggplot2::scale_shape_manual(values = c(mean = 16, low = 24, high = 25), labels = CM_REFS, name = "Funnel reference")
  } else {
    p <- p +
      ggplot2::geom_linerange(data = cl, ggplot2::aes(x = n_total, ymin = d_high, ymax = d_low, colour = NULL),
                              colour = "grey35", alpha = 0.22, linewidth = 3.2, lineend = "round") +
      ggplot2::geom_point(data = cl, ggplot2::aes(x = n_total, y = d_low), colour = CM_REF_COLOURS[["low"]], shape = 95, size = 4) +
      ggplot2::geom_point(data = cl, ggplot2::aes(x = n_total, y = d_high), colour = CM_REF_COLOURS[["high"]], shape = 95, size = 4) +
      ggplot2::geom_point(data = cl, ggplot2::aes(x = n_total, y = d_mean, fill = cm_status),
                          shape = 21, size = 2.8, colour = "grey15", stroke = 0.4) +
      ggplot2::scale_fill_manual(values = CM_STATUS_COLOURS, drop = FALSE, name = "Model status") +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 3.5), order = 2))
  }
  p <- p +
    ggrepel::geom_text_repel(data = cm_label_data(cl, label_points),
                             ggplot2::aes(x = n_total, y = d_mean, label = label), size = 2.7,
                             max.overlaps = Inf, seed = 1, min.segment.length = 0, nudge_x = 0.04) +
    cm_scales(x_log) +
    ggplot2::coord_cartesian(xlim = xr, ylim = yr) +
    ggplot2::labs(x = "Procedure volume (n)", y = "Model rate minus reference rate", ...) +
    theme_funnel() +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1), linetype = ggplot2::guide_legend(order = 1),
                    shape = ggplot2::guide_legend(order = 1))
  if (facet) p <- p + ggplot2::facet_wrap(ggplot2::vars(registry), scales = "free")
  p
}
