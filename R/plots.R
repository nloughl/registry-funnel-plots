# Funnel plots ------------------------------------------------------------------------------------
#
# One core function, `funnel_plot()`, draws limit curves (+ optional points) for any grouping.
# The named variants below are thin wrappers that pick the data, grouping and defaults, so a new
# variation is usually a few lines: filter the curves/points and call funnel_plot() with other options.
#
# y scale:  "difference" -> device rate minus reference rate (reference = 0)   [default]
#           "ratio"      -> device rate / reference rate (reference = 1)
#           "absolute"   -> revision rate itself (reference = registry mean)

y_cols <- function(y) switch(y,
  difference = list(lower = "dev_lower", upper = "dev_upper", point = "dev", ref = 0,
                    lab = "Deviation from registry mean", fmt = scales::percent_format(accuracy = 0.1)),
  ratio      = list(lower = "ratio_lower", upper = "ratio_upper", point = "ratio", ref = 1,
                    lab = "Ratio to registry mean", fmt = scales::label_number(accuracy = 0.1, suffix = "x")),
  absolute   = list(lower = "lower", upper = "upper", point = "p", ref = NA,
                    lab = "Cumulative revision", fmt = scales::percent_format(accuracy = 0.1)),
  stop("y must be 'difference', 'ratio' or 'absolute'"))

#' Core plot.
#' curves   output of build_limit_curves() (must contain the `colour_by` column)
#' points   optional device/stratum points with n_total, p, dev, ratio, label
#' levels   which control-limit levels to draw, e.g. 0.95, 0.998 or both
#' sides    "both" or "upper" (upper only = one-sided screening for poor performers)
funnel_plot <- function(curves, points = NULL, means = NULL,
                        colour_by = "registry", colour_scale = NULL,
                        levels = c(0.95, 0.998), sides = c("both", "upper"), y = "difference",
                        label_points = TRUE, label_filter = NULL, facet = NULL,
                        xlim = NULL, ylim = NULL, n_floor = NULL, x_log = FALSE, ribbon = FALSE,
                        title = NULL, subtitle = NULL, caption = NULL, x_lab = "Procedure volume (n)") {
  sides <- match.arg(sides)
  yc <- y_cols(y)
  cv <- dplyr::filter(curves, level %in% levels)
  if (!nrow(cv)) stop("no curves for levels ", paste(levels, collapse = ", "))

  p <- ggplot2::ggplot(cv, ggplot2::aes(x = n, colour = .data[[colour_by]]))
  if (ribbon) {
    p <- p + ggplot2::geom_ribbon(
      data = dplyr::filter(cv, level == max(levels)),
      ggplot2::aes(ymin = .data[[yc$lower]], ymax = .data[[yc$upper]], fill = .data[[colour_by]]),
      alpha = 0.08, colour = NA, show.legend = FALSE)
  }
  p <- p + ggplot2::geom_line(ggplot2::aes(y = .data[[yc$upper]], linetype = level_lab,
                                           group = interaction(.data[[colour_by]], level_lab,
                                                               if (!is.null(facet)) .data[[facet]] else 1)),
                              linewidth = 0.6)
  if (sides == "both") {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data[[yc$lower]], linetype = level_lab,
                                             group = interaction(.data[[colour_by]], level_lab,
                                                                 if (!is.null(facet)) .data[[facet]] else 1)),
                                linewidth = 0.6)
  }
  if (y == "absolute") {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = p0), linewidth = 0.4, alpha = 0.7)
  } else {
    p <- p + ggplot2::geom_hline(yintercept = yc$ref, colour = "black", linewidth = 0.4)
  }

  if (!is.null(points) && nrow(points)) {
    has_risk <- "low_at_risk" %in% names(points) && any(points$low_at_risk %in% TRUE)
    pt_map <- ggplot2::aes(x = n_total, y = .data[[yc$point]], colour = .data[[colour_by]])
    if (has_risk) {
      points$low_at_risk <- factor(points$low_at_risk %in% TRUE, levels = c(FALSE, TRUE))
      pt_map$shape <- quote(low_at_risk)
    }
    p <- p + ggplot2::geom_point(data = points, pt_map, size = 2.3, inherit.aes = FALSE)
    if (has_risk) {
      p <- p + ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                                           labels = c(`TRUE` = "<=250 at risk"),
                                           breaks = "TRUE", name = NULL, drop = TRUE)
    }
    if (label_points) {
      lp <- if (is.null(label_filter)) points else dplyr::filter(points, !!label_filter)
      p <- p + ggrepel::geom_text_repel(
        data = lp, ggplot2::aes(x = n_total, y = .data[[yc$point]], label = label, colour = .data[[colour_by]]),
        size = 2.8, show.legend = FALSE, max.overlaps = Inf, seed = 1, min.segment.length = 0,
        segment.size = 0.3, box.padding = 0.35, inherit.aes = FALSE)
    }
  }

  if (!is.null(colour_scale)) {
    p <- p + colour_scale
  } else if (colour_by == "registry" && !is.null(means)) {
    p <- p + scale_colour_registry(means) + scale_fill_registry(means, guide = "none")
  }
  p <- p + scale_linetype_level(levels) +
    ggplot2::scale_y_continuous(labels = yc$fmt, breaks = scales::pretty_breaks(n = 8))
  p <- p + if (x_log) ggplot2::scale_x_log10(labels = scales::comma) else ggplot2::scale_x_continuous(labels = scales::comma)
  if (!is.null(facet)) p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet]]))
  rng <- dynamic_ranges(cv, points, yc, sides, n_floor, x_log)
  p + ggplot2::coord_cartesian(xlim = xlim %||% rng$x, ylim = ylim %||% rng$y) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption, x = x_lab, y = yc$lab) +
    theme_funnel()
}

#' Axis windows worked out from the data, so every point and both sides of every funnel are in view.
#'  x: 0 (or n_floor on a log axis) to the largest device volume (or the curves' n range).
#'  y: symmetric about the reference (0 or 1), wide enough for all points and for the funnel down to
#'     n_floor - the smallest device volume on the plot (funnels diverge as n -> 0, so some floor is needed).
dynamic_ranges <- function(curves, points = NULL, yc, sides = "both", n_floor = NULL, x_log = FALSE,
                           x_pad = 1.05, y_pad = 1.08) {
  has_pts <- !is.null(points) && nrow(points) > 0
  n_floor <- n_floor %||% if (has_pts) min(points$n_total, na.rm = TRUE) else 100
  x_hi <- (if (has_pts) max(points$n_total, na.rm = TRUE) else max(curves$n)) * x_pad
  x_lo <- if (x_log) n_floor / x_pad else 0
  vis <- curves[curves$n >= n_floor & curves$n <= x_hi, , drop = FALSE]
  ys <- c(vis[[yc$upper]], if (sides == "both") vis[[yc$lower]], if (has_pts) points[[yc$point]])
  ys <- ys[is.finite(ys)]
  if (is.na(yc$ref)) {                       # absolute scale: plain data range
    r <- range(ys)
    y <- r + c(-1, 1) * diff(r) * (y_pad - 1)
  } else {
    h <- max(abs(ys - yc$ref)) * y_pad
    y <- if (sides == "both") yc$ref + c(-h, h) else yc$ref + c(-0.08 * h, h)
  }
  list(x = c(x_lo, x_hi), y = y, n_floor = n_floor)
}

# Named variants ----------------------------------------------------------------------------------

#' Per-registry panels: each registry's own funnel (no devices), mean set to 0 (or 1 for ratio) so
#' both sides of the funnel show. Strip labels carry the registry mean and n.
plot_registry_panels <- function(curves, means, levels, y = "difference", n_floor = 100, x_log = TRUE, ...) {
  lab <- registry_legend_labels(means)
  funnel_plot(curves, NULL, means, levels = levels, y = y, n_floor = n_floor, ribbon = TRUE, x_log = x_log, ...) +
    ggplot2::facet_wrap(~registry, labeller = ggplot2::as_labeller(lab)) +
    ggplot2::guides(colour = "none")
}

#' All registries' limits overlaid after setting every registry mean to 0.
plot_standardised_overlay <- function(curves, means, levels, points = NULL, ...) {
  funnel_plot(curves, points, means, levels = levels, ...)
}

#' One device family (e.g. Oxford) across registries, against each registry's own limits.
plot_family <- function(curves, points, means, family, levels, ...) {
  pts <- dplyr::filter(points, .data$family == !!family)
  if (!nrow(pts)) stop("no devices in family '", family, "'")
  funnel_plot(dplyr::filter(curves, registry %in% pts$registry), pts,
              dplyr::filter(means, registry %in% pts$registry), levels = levels, ...)
}

#' Reference grid: one row per reference (registry mean, lower / upper 95% CI of the mean), one column per
#' registry. Every panel is standardised to its own reference (0 = that reference), so all funnels are
#' symmetric about 0 (normal / overdispersed limits) and each point is plotted against the reference of its
#' panel: a point outside a panel's funnel is an outlier under that reference.
#' ref_curves: build_limit_curves(reference_refs(...)) - must have `reference` and `p0` columns.
#' Points outside the funnel of their panel are labelled.
plot_reference_grid <- function(ref_curves, means, level, points = NULL, y = "difference",
                                n_floor = NULL, x_log = TRUE, ...) {
  yc <- y_cols(y)
  ref_lab <- c(mean = "Mean", lcl = "Lower CI", ucl = "Upper CI")
  cv <- dplyr::filter(ref_curves, .data$level == !!level, registry %in% means$registry) |>
    dplyr::mutate(reference = factor(reference, levels = names(ref_lab), labels = ref_lab))
  pts <- NULL
  if (!is.null(points) && nrow(points)) {
    refs <- dplyr::distinct(cv, registry, reference, p0)
    pts <- points |>
      dplyr::select(-dplyr::any_of(c("dev", "ratio", "reference"))) |>
      dplyr::inner_join(refs, by = "registry", relationship = "many-to-many") |>
      dplyr::mutate(dev = p - p0, ratio = p / p0)
    # outside this panel's funnel? (interpolate the panel's limit curve at the point's n)
    pts$outside <- purrr::pmap_lgl(list(pts$registry, pts$reference, pts$n_total, pts[[yc$point]]),
      function(rg, rf, n, yv) {
        c1 <- cv[cv$registry == rg & cv$reference == rf, ]
        yv > stats::approx(c1$n, c1[[yc$upper]], n, rule = 2)$y || yv < stats::approx(c1$n, c1[[yc$lower]], n, rule = 2)$y
      })
  }
  p <- ggplot2::ggplot(cv, ggplot2::aes(x = n, colour = registry)) +
    ggplot2::geom_line(ggplot2::aes(y = .data[[yc$upper]]), linewidth = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = .data[[yc$lower]]), linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = yc$ref, linewidth = 0.4)
  if (!is.null(pts)) {
    p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = n_total, y = .data[[yc$point]], shape = outside), size = 2.2) +
      ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c(`FALSE` = "inside", `TRUE` = "outside"),
                                  name = "Point vs funnel") +
      ggrepel::geom_text_repel(data = dplyr::filter(pts, outside),
                               ggplot2::aes(x = n_total, y = .data[[yc$point]], label = label),
                               size = 2.5, show.legend = FALSE, max.overlaps = Inf, seed = 1, min.segment.length = 0)
  }
  rng <- dynamic_ranges(cv, pts, yc, "both", n_floor, x_log)
  p + ggplot2::facet_grid(reference ~ registry) +
    scale_colour_registry(means, ci = TRUE, name = "Registry (mean, 95% CI, n)") +
    ggplot2::scale_y_continuous(labels = yc$fmt, breaks = scales::pretty_breaks(n = 6)) +
    (if (x_log) ggplot2::scale_x_log10(labels = function(x) ifelse(x >= 1000, paste0(x / 1000, 'k'), x))
     else ggplot2::scale_x_continuous(labels = scales::comma)) +
    ggplot2::coord_cartesian(xlim = rng$x, ylim = rng$y) +
    ggplot2::labs(x = "Procedure volume (n)", y = sub("registry mean", "reference", yc$lab), ...) +
    theme_funnel() + ggplot2::theme(legend.position = "bottom", legend.box = "vertical") +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, title.position = "top"))
}

#' Case-mix: one set of limits per age-sex class (e.g. "<55 Female"), one panel per registry/design.
#' Paired colours: one hue per age group, dark = female, light = male (registry colours are reserved).
casemix_class_levels <- function() {
  as.vector(t(outer(c("<55", "55-64", "65-74", ">=75"), c("Male", "Female"), paste)))
}

plot_casemix_panels <- function(cm_curves, level, sides = "both", y = "difference", n_floor = 100,
                                x_log = TRUE, ...) {
  yc <- y_cols(y)
  cv <- dplyr::filter(cm_curves, .data$level == !!level) |>
    dplyr::mutate(class = factor(paste(age_group, sex), levels = casemix_class_levels()))
  rng <- dynamic_ranges(cv, NULL, yc, sides, n_floor, x_log)
  p <- ggplot2::ggplot(cv, ggplot2::aes(x = n, colour = class, group = interaction(class, design, registry))) +
    ggplot2::geom_line(ggplot2::aes(y = .data[[yc$upper]]), linewidth = 0.6)
  if (sides == "both") p <- p + ggplot2::geom_line(ggplot2::aes(y = .data[[yc$lower]]), linewidth = 0.6)
  p + ggplot2::geom_hline(yintercept = yc$ref, linewidth = 0.4) +
    ggplot2::facet_wrap(ggplot2::vars(registry, design), nrow = 1,
                        labeller = ggplot2::labeller(.multi_line = FALSE)) +
    ggplot2::scale_colour_brewer(palette = "Paired", name = "Age-sex class", drop = FALSE) +
    ggplot2::scale_y_continuous(labels = yc$fmt) +
    (if (x_log) ggplot2::scale_x_log10(labels = scales::comma) else ggplot2::scale_x_continuous(labels = scales::comma)) +
    ggplot2::coord_cartesian(xlim = rng$x, ylim = rng$y) +
    ggplot2::labs(x = "Procedure volume (n)", y = sub("registry", "age-sex class", yc$lab), ...) +
    theme_funnel()
}

#' Case-mix envelope: for each registry, the band spanned by its age-sex class limits
#' (narrowest to widest), in the registry colour, with the registry-mean limit on top.
plot_casemix_envelope <- function(cm_curves, reg_curves, means, level, points = NULL, y = "difference",
                                  n_floor = 100, ...) {
  yc <- y_cols(y)
  env <- cm_curves |>
    dplyr::filter(.data$level == !!level) |>
    dplyr::group_by(registry, n) |>
    dplyr::summarise(up_lo = min(.data[[yc$upper]]), up_hi = max(.data[[yc$upper]]),
                     lo_lo = min(.data[[yc$lower]]), lo_hi = max(.data[[yc$lower]]), .groups = "drop")
  means <- dplyr::filter(means, registry %in% env$registry)
  rc <- dplyr::filter(reg_curves, .data$level == !!level, registry %in% env$registry)
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = env, ggplot2::aes(x = n, ymin = up_lo, ymax = up_hi, fill = registry), alpha = 0.25) +
    ggplot2::geom_ribbon(data = env, ggplot2::aes(x = n, ymin = lo_lo, ymax = lo_hi, fill = registry), alpha = 0.25) +
    ggplot2::geom_line(data = rc, ggplot2::aes(x = n, y = .data[[yc$upper]], colour = registry), linewidth = 0.6) +
    ggplot2::geom_line(data = rc, ggplot2::aes(x = n, y = .data[[yc$lower]], colour = registry), linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = yc$ref, linewidth = 0.4)
  if (!is.null(points) && nrow(points)) {
    pts <- dplyr::filter(points, registry %in% env$registry)
    p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = n_total, y = .data[[yc$point]], colour = registry), size = 2.3) +
      ggrepel::geom_text_repel(data = pts, ggplot2::aes(x = n_total, y = .data[[yc$point]], label = label, colour = registry),
                               size = 2.8, show.legend = FALSE, max.overlaps = Inf, seed = 1, min.segment.length = 0)
  }
  pts <- if (!is.null(points) && nrow(points)) dplyr::filter(points, registry %in% env$registry) else NULL
  rng <- dynamic_ranges(dplyr::filter(cm_curves, .data$level == !!level), pts, yc, "both", n_floor)
  p + scale_colour_registry(means) + scale_fill_registry(means, guide = "none") +
    ggplot2::scale_y_continuous(labels = yc$fmt) +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::coord_cartesian(xlim = rng$x, ylim = rng$y) +
    ggplot2::labs(x = "Procedure volume (n)", y = yc$lab, ...) +
    theme_funnel()
}

#' Procedure-period funnels: one funnel per period (e.g. LROI 2009-2010 ... 2021-2022), each built around
#' that period's own revision rate, with the registry-mean funnel (dashed) and optional device points.
#' y = "absolute": funnels and points at their true rates (the funnels shift with the period's rate).
#' y = "difference" / "ratio": every period funnel centred on 0 / 1 (shows how the width changes);
#'   points are then relative to the registry mean.
plot_period_funnels <- function(period_curves, reg_curves, means, level, points = NULL, y = "absolute",
                                n_floor = NULL, x_log = TRUE, label_points = TRUE, ...) {
  yc <- y_cols(y)
  cv <- dplyr::filter(period_curves, .data$level == !!level)
  rc <- dplyr::filter(reg_curves, .data$level == !!level, registry %in% unique(cv$registry))
  cols <- period_colours(levels(cv$period))
  p <- ggplot2::ggplot(cv, ggplot2::aes(x = n, colour = period, group = period)) +
    ggplot2::geom_line(ggplot2::aes(y = .data[[yc$upper]]), linewidth = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = .data[[yc$lower]]), linewidth = 0.6) +
    ggplot2::geom_line(data = rc, ggplot2::aes(x = n, y = .data[[yc$upper]]), colour = "black",
                       linetype = "22", linewidth = 0.7, inherit.aes = FALSE) +
    ggplot2::geom_line(data = rc, ggplot2::aes(x = n, y = .data[[yc$lower]]), colour = "black",
                       linetype = "22", linewidth = 0.7, inherit.aes = FALSE)
  if (is.na(yc$ref)) {
    p <- p + ggplot2::geom_hline(yintercept = unique(rc$p0), colour = "black", linewidth = 0.4)
  } else p <- p + ggplot2::geom_hline(yintercept = yc$ref, colour = "black", linewidth = 0.4)
  if (!is.null(points) && nrow(points)) {
    p <- p + ggplot2::geom_point(data = points, ggplot2::aes(x = n_total, y = .data[[yc$point]]),
                                 colour = registry_colour(points$registry[1]), size = 2.3, inherit.aes = FALSE)
    if (label_points) p <- p + ggrepel::geom_text_repel(
      data = points, ggplot2::aes(x = n_total, y = .data[[yc$point]], label = label),
      colour = registry_colour(points$registry[1]), size = 2.8, max.overlaps = Inf, seed = 1,
      min.segment.length = 0, inherit.aes = FALSE)
  }
  rng <- dynamic_ranges(dplyr::bind_rows(cv, rc), points, yc, "both", n_floor, x_log)
  if (is.na(yc$ref)) rng$y[1] <- max(0, rng$y[1])
  p + ggplot2::scale_colour_manual(values = cols, name = paste0("Procedure period\n(", level_label(level), " limits)")) +
    ggplot2::scale_y_continuous(labels = yc$fmt, breaks = scales::pretty_breaks(n = 8)) +
    (if (x_log) ggplot2::scale_x_log10(labels = scales::comma) else ggplot2::scale_x_continuous(labels = scales::comma)) +
    ggplot2::coord_cartesian(xlim = rng$x, ylim = rng$y) +
    ggplot2::labs(x = "Procedure volume (n)", y = yc$lab, ...) +
    theme_funnel()
}

#' Save helper: file name from a slug, into out_dir.
save_plot <- function(p, slug, out_dir, width = 10, height = 6, dpi = 300) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(out_dir, paste0(slug, ".png"))
  ggplot2::ggsave(f, p, width = width, height = height, dpi = dpi)
  invisible(f)
}
