# Registry colours ---------------------------------------------------------------
# One fixed colour per registry, used by every plot, so a registry is recognisable across figures.
# Okabe-Ito based (colour-blind safe). Add new registries here; never reuse a colour.

REGISTRY_COLOURS <- c(
  AOANJRR = "#E69F00",  # orange
  NJR     = "#0072B2",  # blue
  LROI    = "#D55E00",  # vermillion
  SIRIS   = "#CC79A7",  # reddish purple
  EPRD    = "#009E73",  # bluish green
  AJRR    = "#56B4E9",  # sky blue
  NZJR    = "#882255",  # wine
  SAR     = "#332288",  # indigo
  NAR     = "#44AA99",  # teal
  CJRR    = "#999999"   # grey
)

registry_colour <- function(registries) {
  missing <- setdiff(unique(registries), names(REGISTRY_COLOURS))
  if (length(missing)) stop("No colour defined for registry: ", paste(missing, collapse = ", "),
                            ". Add it to REGISTRY_COLOURS in R/palette.R")
  REGISTRY_COLOURS[unique(registries)]
}

#' Legend label per registry: "NJR: 3.21%, n = 206,689"
registry_legend_labels <- function(means) {
  stats::setNames(
    sprintf("%s: %s, n = %s", means$registry,
            scales::percent(means$p_ref, accuracy = 0.01),
            scales::comma(means$n_total)),
    means$registry
  )
}

#' Colour scale keyed on the registry code, labelled with mean rate and n.
#' Mapping colour to `registry` (not to the label text) keeps colours fixed whatever the labels say.
scale_colour_registry <- function(means, name = "Registry (mean revision, n)", ...) {
  regs <- means$registry
  ggplot2::scale_colour_manual(values = registry_colour(regs), breaks = regs,
                               labels = registry_legend_labels(means)[regs], name = name, ...)
}

scale_fill_registry <- function(means, name = "Registry (mean revision, n)", ...) {
  regs <- means$registry
  ggplot2::scale_fill_manual(values = registry_colour(regs), breaks = regs,
                             labels = registry_legend_labels(means)[regs], name = name, ...)
}

# Control-limit line types: always the same for the same level
level_label <- function(level) paste0(format(level * 100, trim = TRUE, drop0trailing = TRUE), "%")

# fixed line type per level, so 95% and 99.8% look the same in every figure
LEVEL_LINETYPES <- c("95%" = "dashed", "99.8%" = "solid", "99%" = "dotdash", "90%" = "dotted")

scale_linetype_level <- function(levels, name = "Control limit") {
  labs <- level_label(sort(unique(levels)))
  vals <- LEVEL_LINETYPES[labs]
  vals[is.na(vals)] <- "longdash"
  ggplot2::scale_linetype_manual(values = stats::setNames(unname(vals), labs), breaks = labs, name = name)
}

theme_funnel <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "right",
                   plot.title.position = "plot",
                   plot.caption.position = "plot",
                   plot.caption = ggplot2::element_text(hjust = 0, size = ggplot2::rel(0.75), colour = "grey30"),
                   strip.background = ggplot2::element_rect(fill = "grey95"))
}
