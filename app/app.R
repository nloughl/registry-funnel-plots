# Interactive registry funnel plots ---------------------------------------------------------------
#
# Run from the repo root:   shiny::runApp("app")
# Data: the regextract outputs (UKA_long.csv, UKA_casemix_long.csv). Default location is the sibling
# repo ../regitry-outlier-detection/outputs/UKA; set REGISTRY_DATA_DIR to point elsewhere.
#
# Everything statistical (registry means, limits, outliers, phi) comes from ../R, the same code the
# funnel_plots.Rmd report uses, so the app and the report always agree.

library(shiny)
library(bslib)
library(plotly)
library(dplyr)

APP_DIR  <- normalizePath(if (file.exists("app.R")) "." else "app")
ROOT     <- normalizePath(file.path(APP_DIR, ".."))
for (f in list.files(file.path(ROOT, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

DATA_DIR  <- Sys.getenv("REGISTRY_DATA_DIR",
                        file.path(ROOT, "..", "regitry-outlier-detection", "outputs", "UKA"))
PROCEDURE <- "UKA"

ext       <- load_extraction(DATA_DIR, PROCEDURE)
families  <- load_families(file.path(ROOT, "config", "device_families.csv"))
overrides <- readr::read_csv(file.path(ROOT, "config", "registry_means_override.csv"),
                             comment = "#", show_col_types = FALSE)
attr_rules <- load_attribute_rules(file.path(ROOT, "config", "device_attributes.csv"))

# follow-up years reported by at least two registries (device rows with an estimate)
YEARS <- ext$device |>
  filter(row_type == "device", has_estimate(pick(everything()))) |>
  distinct(registry, time_yr) |>
  count(time_yr) |>
  filter(n >= 2) |>
  pull(time_yr) |>
  sort()
REGISTRIES <- sort(unique(ext$device$registry))
# registries with an age-sex case-mix table (LROI's rows in the case-mix file are procedure periods, not age-sex)
CASEMIX_REGS <- if (is.null(ext$casemix)) character() else sort(unique(ext$casemix$registry[!is.na(ext$casemix$age_group)]))

METHODS <- c("Normal approximation (Wald)" = "normal",
             "Exact binomial (Spiegelhalter)" = "exact")
#             "Overdispersion-adjusted (multiplicative)" = "overdispersed")

# ------------------------------------------------------------------------------------------------ UI
source(file.path(APP_DIR, "rare_tab.R"), local = TRUE)   # Rare joints tab; local = TRUE: it uses ROOT / DATA_DIR defined above

uka_ui <- nav_panel(paste("Knee", paste0("(", PROCEDURE, ")")), layout_sidebar(
  fillable = FALSE,
  sidebar = sidebar(
    width = 320,
    selectInput("year", "Follow-up year", choices = YEARS, selected = if (3 %in% YEARS) 3 else YEARS[1]),
    selectInput("family", "Model", choices = c("All models" = "__all__", families$family), selected = "Oxford"),
    selectInput("detail", "Detail level", choices = DETAIL_LEVELS, selected = 1),
    conditionalPanel("input.detail == '3'",
      selectInput("attr3", "Level-3 attribute", choices = LEVEL3_ATTRS)),
    conditionalPanel("input.detail == '2' || input.detail == '3'",
      radioButtons("fix_filter", "Fixation shown",
                   choices = c("All" = "all", "Cemented only" = "cemented", "Cementless only" = "cementless"),
                   selected = "all", inline = TRUE)),
    conditionalPanel("input.detail != '0'",
      radioButtons("attr_src", "Attributes",
                   choices = c("Stated in report" = "stated", "Stated + design knowledge" = "design"),
                   selected = "stated", inline = TRUE)),
    checkboxGroupInput("registries", "Registries", choices = REGISTRIES, selected = REGISTRIES, inline = TRUE),
    radioButtons("level", "Control limit",
                 choices = c("99.8%" = "0.998", "95%" = "0.95", "Both" = "both"), selected = "0.998", inline = TRUE),
    selectInput("method", "Funnel distribution", choices = METHODS),
    radioButtons("reference", "Funnel reference", choices = REFERENCES, selected = "mean"),
    checkboxInput("ref_overlay", "Also draw funnels around the mean and both CI bounds", FALSE),
    radioButtons("centre", "Funnel",
                 choices = c("Registry mean" = "mean", "Registry mean + case-mix envelope" = "casemix",
                             "Registry mean + LROI funnels by procedure period" = "period"),
                 selected = "mean"),
    radioButtons("layout", "Layout",
                 choices = c("Overlay (standardised)" = "overlay", "Panel per registry (rate)" = "panels"),
                 selected = "overlay"),
    conditionalPanel("input.layout == 'overlay'",
      radioButtons("yscale", "Y scale", choices = c("Difference from mean" = "difference", "Ratio to mean" = "ratio"),
                   selected = "difference", inline = TRUE)),
    checkboxInput("xlog", "Log x axis", FALSE),
    radioButtons("labels", "Label models on the plot",
                 choices = c("Outliers" = "outliers", "All" = "all", "None" = "none"), selected = "outliers", inline = TRUE),
    hr(),
    downloadButton("dl_outliers", "Outliers (CSV)", class = "btn-sm"),
    downloadButton("dl_points", "All points (CSV)", class = "btn-sm"),
    helpText("Use the camera icon on the plot to save a PNG.")
  ),
  uiOutput("notes"),
  card(full_screen = TRUE, card_header(textOutput("plot_title")),
       plotlyOutput("funnel", height = "720px")),
  layout_columns(
    col_widths = c(7, 5),
    card(card_header(textOutput("outlier_title")), tableOutput("outliers")),
    card(card_header("Cumulative excess revision by registry (outlier devices only)"), tableOutput("summary"),
         card_footer(helpText(
           tags$b("Only outlier devices are counted."), " An outlier device has a rate above the upper control limit ",
           "of its own registry at its own n. Devices inside the funnel add nothing, even if above the registry mean; ",
           "low outliers (green) are not counted either.", tags$br(),
           tags$b("Outlier devices:"), " number of devices above the upper limit.", tags$br(),
           tags$b("Procedures with an outlier device (n)"), " [procedures]: sum of n over outlier devices only, i.e. how many ",
           "patients received an outlier implant, revised or not. Exposure, not harm.", tags$br(),
           tags$b("Cumulative excess"), " [revisions]: sum over outlier devices only of n × (device rate − registry mean), the ",
           "approximate number of revisions above what the registry-mean rate would have given.", tags$br(),
           tags$i("Example: n = 10,000 at 5.0%, registry mean 4.0%, upper limit 4.6% → outlier: 10,000 procedures, ",
                  "100 excess revisions. The same n at 4.5% is inside the limit and adds 0. Screening numbers, not precise estimates."))))
  ),
  card(card_header("Registry means, 95% CI and funnel reference at this follow-up year"), tableOutput("means"),
       card_footer(helpText("CI: published (NJR, SIRIS) or pooled from the published CIs of the parts (LROI cemented + uncemented, ",
                            "AOANJRR male + female; approximate). EPRD has no published total, so it keeps its mean under every reference."))),
  card(card_header(textOutput("level_title")), tableOutput("level_counts"),
       card_footer(helpText("Points / rows kept (share of procedures kept) for the selected model. ",
                            "Rows without the detail a level needs are excluded, not regrouped.")))
))

ui <- page_navbar(
  title = "Registry revision funnel plots",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  fillable = FALSE,
  uka_ui,
  rare_ui()
)

# -------------------------------------------------------------------------------------------- server
server <- function(input, output, session) {

  rare_server(input, output, session)

  yr      <- reactive(as.numeric(input$year))
  levels_ <- reactive(if (input$level == "both") c(0.95, 0.998) else as.numeric(input$level))
  out_lvl <- reactive(max(levels_()))
  ref_sel <- reactive(if (is.null(input$reference)) "mean" else input$reference)
  # fixation filter only applies at detail levels 2 and 3 (where fixation is known)
  fix_sel <- reactive(if (is.null(input$fix_filter) || detail() < 2) "all" else input$fix_filter)
  label_sel <- reactive(if (is.null(input$labels)) "outliers" else if (isTRUE(input$labels)) "outliers"
                        else if (isFALSE(input$labels)) "none" else input$labels)
  labelled <- function(d) switch(label_sel(), all = d, none = d[0, ], d[d$outlier != "none", ])

  means <- reactive({
    registry_means(ext$device, ext$casemix, yr(), overrides) |>
      set_reference(ref_sel()) |>
      filter(registry %in% input$registries)
  })

  # all devices of the selected registries at this year (phi is estimated from all of them)
  detail <- reactive(as.integer(input$detail))

  # registry rows with their attributes (before pooling)
  rows <- reactive({
    device_points(ext$device, means(), yr(), families) |>
      filter(registry %in% input$registries, !is.na(p_ref)) |>
      derive_attributes(attr_rules, input$attr_src)
  })

  # points at the selected detail level (phi is estimated from all of them)
  all_points <- reactive({
    pts <- pool_to_level(rows(), detail(), input$attr3)
    left_join(pts, estimate_phi(pts), by = "registry")
  })

  points <- reactive({
    pts <- all_points()
    if (input$family != "__all__") pts <- filter(pts, family == input$family)
    if (fix_sel() != "all") pts <- filter(pts, fixation %in% fix_sel())
    pts <- classify_points(pts, levels_(), input$method)
    hi <- outlier_report(pts, out_lvl(), input$method)
    lo_lim <- purrr::pmap_dbl(list(pts$p_ref, pts$n_total, pts$phi),
                              function(p0, n, ph) funnel_limits(p0, n, out_lvl(), input$method, ph)$lower)
    up_lim <- purrr::pmap_dbl(list(pts$p_ref, pts$n_total, pts$phi),
                              function(p0, n, ph) funnel_limits(p0, n, out_lvl(), input$method, ph)$upper)
    pts |>
      mutate(upper_limit = up_lim, lower_limit = lo_lim,
             outlier = case_when(p > upper_limit ~ "high", p < lower_limit ~ "low", TRUE ~ "none"))
  })

  plot_regs <- reactive(intersect(means()$registry, unique(points()$registry)))

  curves <- reactive({
    req(nrow(points()) > 0)
    refs <- means() |>
      filter(registry %in% plot_regs()) |>
      left_join(distinct(all_points(), registry, phi), by = "registry") |>
      transmute(registry, p_ref, phi = coalesce(phi, 1), n_max = max(points()$n_total) * 1.1)
    build_limit_curves(refs, levels_(), input$method)
  })

  # funnels around the mean and both CI bounds (outlier level only), for the reference overlay
  ref_curves <- reactive({
    if (!isTRUE(input$ref_overlay %||% FALSE)) return(NULL)
    m <- filter(means(), registry %in% plot_regs())
    refs <- reference_refs(m, c("mean", "lcl", "ucl"), n_max = max(points()$n_total) * 1.1,
                           phi = distinct(all_points(), registry, phi))
    refs <- filter(refs, reference != ref_sel())   # the selected reference is already drawn
    if (!nrow(refs)) return(NULL)
    build_limit_curves(refs, out_lvl(), input$method)
  })
  REF_DASH <- c(mean = "longdash", lcl = "dot", ucl = "dashdot")

  cm_curves <- reactive({
    if (input$centre != "casemix" || is.null(ext$casemix)) return(NULL)
    st <- casemix_strata(ext$casemix, yr(), njr_level = "group") |>
      filter(registry %in% plot_regs())
    if (!nrow(st)) return(NULL)
    phis <- distinct(all_points(), registry, phi)
    refs <- st |>
      left_join(phis, by = "registry") |>
      transmute(registry, design, sex, age_group, stratum, p_ref = p, phi = coalesce(phi, 1),
                n_max = max(points()$n_total) * 1.1)
    build_limit_curves(refs, out_lvl(), input$method)
  })

  # LROI procedure-period funnels (K052B): one funnel per 2-year period around that period's own rate
  per_curves <- reactive({
    if (input$centre != "period" || !"LROI" %in% plot_regs()) return(NULL)
    per <- period_strata(ext$casemix, yr())
    if (is.null(per)) return(NULL)
    ph <- distinct(all_points(), registry, phi)
    per |>
      left_join(ph, by = "registry") |>
      transmute(registry, period, p_ref = p, phi = coalesce(phi, 1), n_max = max(points()$n_total) * 1.1) |>
      build_limit_curves(out_lvl(), input$method) |>
      mutate(period = factor(period, levels = levels(per$period)))
  })

  add_period_lines <- function(p, pc, ycols) {
    pcol <- period_colours(levels(pc$period))
    for (pr in levels(pc$period)) {
      c1 <- filter(pc, period == pr)
      if (!nrow(c1)) next
      for (side in c("upper", "lower")) {
        p <- add_lines(p, data = c1, x = ~n, y = c1[[ycols[[side]]]], opacity = 0.85,
                       line = list(color = pcol[[pr]], width = 1.3),
                       name = paste("LROI", pr), legendgroup = paste0("period_", pr), showlegend = side == "upper",
                       hovertemplate = paste0("LROI ", pr, " (", scales::percent(c1$p0[1], 0.01), ") ",
                                              level_label(out_lvl()), " ", side, "<br>n = %{x:,}<br>%{y:.2%}<extra></extra>"))
      }
    }
    p
  }

  # ---------------------------------------------------------------------------------------- notes
  output$notes <- renderUI({
    msgs <- character()
    missing <- setdiff(input$registries, plot_regs())
    if (length(missing)) msgs <- c(msgs, sprintf("No %s-year estimate for the selected model in: %s.",
                                                 yr(), paste(missing, collapse = ", ")))
    unver <- unique(points()$registry[points()$verification %in% "unverified"])
    if (length(unver)) msgs <- c(msgs, paste0("Values not yet verified: ", paste(unver, collapse = ", "), "."))
    fb <- means()$registry[grepl("FALLBACK", means()$mean_source)]
    if (length(fb)) msgs <- c(msgs, paste0("Mean is the n-weighted mean of listed devices (no registry total) for: ",
                                           paste(fb, collapse = ", "), "."))
    if (input$centre == "period") {
      pc <- per_curves()
      msgs <- c(msgs, if (is.null(pc)) paste0("No LROI procedure-period data at ", yr(), " years (or LROI not shown).")
        else paste0("LROI funnels by procedure period (Figure K052B): one ", level_label(out_lvl()),
          " funnel per 2-year period around that period's own ", yr(), "-year rate. ",
          if (input$layout == "overlay") "In the overlay each period funnel is centred on 0 (width only); use the panel layout to see them at their true rates. "
          else "", "K052B reports MAJOR revision (femur or tibia), the device tables any revision, so treat this as indicative. ",
          "Outliers (rings, tables) are still judged against the registry reference."))
    }
    if (input$centre == "casemix") {
      nocm <- setdiff(plot_regs(), CASEMIX_REGS)
      msgs <- c(msgs, paste0("Case-mix envelope available for ", paste(intersect(plot_regs(), CASEMIX_REGS), collapse = ", "),
                             if (length(nocm)) paste0(" only (not ", paste(nocm, collapse = ", "), ")") else "",
                             ". Outliers are still judged against the registry mean."))
    }
    if (detail() > 0) {
      ex <- attr(pool_to_level(rows(), detail(), input$attr3), "excluded")
      if (input$family != "__all__") ex <- filter(ex, family == input$family)
      if (nrow(ex)) msgs <- c(msgs, sprintf("%s: %d row(s) excluded for missing detail (%s procedures): %s.",
        level_title(detail(), input$attr3), nrow(ex), scales::comma(sum(ex$n_total)),
        { e <- sprintf("%s %s (%s)", ex$registry, ex$label, ex$excluded_because)
          paste0(paste(head(e, 8), collapse = "; "), if (length(e) > 8) sprintf("; ... and %d more", length(e) - 8) else "") }))
    }
    if (input$method == "overdispersed") {
      ph <- distinct(all_points(), registry, phi)
      msgs <- c(msgs, paste0("Overdispersion factor φ: ",
                             paste(sprintf("%s %.2f", ph$registry, ph$phi), collapse = ", "),
                             " (limits widened by √φ)."))
    }
    if (ref_sel() != "mean") {
      noci <- means()$registry[means()$reference_used == "mean"]
      msgs <- c(msgs, paste0("Funnels centred on the ", reference_label(ref_sel()), "; y and the funnels are standardised to ",
                             "this reference (0), and excess is measured from it.",
                             if (length(noci)) paste0(" No CI of the mean for ", paste(noci, collapse = ", "),
                                                      ": these use the mean.") else ""))
    }
    if (isTRUE(input$ref_overlay %||% FALSE)) msgs <- c(msgs, paste0("Thin extra funnels (", level_label(out_lvl()),
      "): long dash = around the registry mean, dotted = lower CI, dash-dot = upper CI of the mean. ",
      "Rings and tables use the selected reference only. ",
      if (input$layout == "overlay") "In the overlay every funnel is standardised to its own reference (all centred on 0), so the extra funnels show the change in width; points are placed relative to the selected reference. Use the panel layout (rate scale) to see the funnels at their true positions." else ""))
    if (fix_sel() != "all") msgs <- c(msgs, paste0("Showing ", fix_sel(), " models only (hybrid and SIRIS ",
      "'cemented/hybrid' rows are left out). Funnels are still centred on each registry's overall UKA mean, ",
      "not a ", fix_sel(), "-only mean."))
    if (!length(msgs)) return(NULL)
    div(class = "alert alert-secondary py-2 small", HTML(paste(msgs, collapse = "<br>")))
  })

  output$plot_title <- renderText({
    sprintf("%s: %s, %s, %d-year revision, %s limits around the %s (%s)", PROCEDURE,
            paste0(if (input$family == "__all__") "all models" else input$family,
                   if (fix_sel() != "all") paste0(" (", fix_sel(), " only)") else ""),
            level_title(detail(), input$attr3), yr(),
            paste(level_label(levels_()), collapse = " & "), reference_label(ref_sel()),
            names(METHODS)[METHODS == input$method])
  })

  # ----------------------------------------------------------------------------------------- plot
  output$funnel <- renderPlotly({
    req(nrow(points()) > 0)
    if (input$layout == "overlay") overlay_plot() else panel_plot()
  })

  hover_text <- function(d) {
    ci <- ifelse(is.na(d$lcl), if (detail() > 0) " (pooled, no CI)" else "",
                 sprintf(" (%s–%s)", scales::percent(d$lcl / 100, 0.01), scales::percent(d$ucl / 100, 0.01)))
    pooled <- if ("n_rows" %in% names(d)) ifelse(d$n_rows > 1, paste0("<br>Pooled rows: ", d$device_label), "") else ""
    sprintf(paste0("<b>%s</b> — %s<br>n = %s<br>%d-yr revision: %s%s%s<br>",
                   "Registry mean: %s<br>Reference (%s): %s<br>Δp vs reference: %s<br>Upper %s limit: %s<br>%s<br><i>%s %s, table %s p.%s</i>"),
            d$registry, d$label, scales::comma(d$n_total), yr(),
            scales::percent(d$p, 0.01), ci, pooled,
            scales::percent(d$p_mean, 0.01), d$reference_used, scales::percent(d$p_ref, 0.01),
            scales::percent(d$p - d$p_ref, 0.01),
            level_label(out_lvl()), scales::percent(d$upper_limit, 0.01),
            c(high = "<b style='color:#c0392b'>HIGH OUTLIER</b>", low = "<b style='color:#1e8449'>Low outlier</b>",
              none = "Within limits")[d$outlier],
            d$registry, d$report_year, d$table_id, d$pdf_page)
  }

  ring_col <- c(high = "#D62728", low = "#2CA02C")
  add_rings <- function(p, d, xcol, ycol) {
    for (kind in c("high", "low")) {
      r <- d[d$outlier == kind, ]
      if (!nrow(r)) next
      p <- add_trace(p, data = r, x = r[[xcol]], y = r[[ycol]], type = "scatter", mode = "markers",
                     marker = list(symbol = "circle-open", size = 18, color = ring_col[[kind]],
                                   line = list(width = 2.5, color = ring_col[[kind]])),
                     name = if (kind == "high") "Above upper limit" else "Below lower limit",
                     legendgroup = kind, hoverinfo = "skip", inherit = FALSE)
    }
    p
  }

  line_dash <- function(lv) if (lv >= 0.998) "solid" else "dash"

  overlay_plot <- function() {
    yc <- y_cols(input$yscale)
    cv <- curves(); pts <- points(); m <- means() |> filter(registry %in% plot_regs())
    labs <- registry_legend_labels(m, ci = TRUE)
    cols <- registry_colour(m$registry)
    p <- plot_ly()
    cm <- cm_curves()
    if (!is.null(cm)) {
      env <- cm |>
        group_by(registry, n) |>
        summarise(up_lo = min(.data[[yc$upper]]), up_hi = max(.data[[yc$upper]]),
                  lo_lo = min(.data[[yc$lower]]), lo_hi = max(.data[[yc$lower]]), .groups = "drop")
      for (reg in unique(env$registry)) {
        e <- filter(env, registry == reg)
        for (side in c("up", "lo")) {
          p <- add_ribbons(p, data = e, x = ~n, ymin = e[[paste0(side, "_lo")]], ymax = e[[paste0(side, "_hi")]],
                           fillcolor = scales::alpha(cols[[reg]], 0.2), line = list(width = 0),
                           name = paste(reg, "case-mix envelope"), legendgroup = paste0(reg, "_cm"),
                           showlegend = side == "up", hoverinfo = "skip")
        }
      }
    }
    pc <- per_curves()
    if (!is.null(pc)) p <- add_period_lines(p, pc, yc)
    for (reg in m$registry) {
      for (lv in levels_()) {
        c1 <- filter(cv, registry == reg, level == lv)
        for (side in c("upper", "lower")) {
          p <- add_lines(p, data = c1, x = ~n, y = c1[[yc[[side]]]],
                         line = list(color = cols[[reg]], width = 1.8, dash = line_dash(lv)),
                         name = labs[[reg]], legendgroup = reg,
                         showlegend = side == "upper" && lv == max(levels_()),
                         hovertemplate = paste0(reg, " ", level_label(lv), " ", side, " limit<br>n = %{x:,}<br>%{y:.2%}<extra></extra>"))
        }
      }
      rc <- ref_curves()
      if (!is.null(rc)) for (r in unique(rc$reference[rc$registry == reg])) {
        c1 <- filter(rc, registry == reg, reference == r)
        for (side in c("upper", "lower")) {
          p <- add_lines(p, data = c1, x = ~n, y = c1[[yc[[side]]]], opacity = 0.7,
                         line = list(color = cols[[reg]], width = 1.2, dash = REF_DASH[[r]]),
                         name = labs[[reg]], legendgroup = reg, showlegend = FALSE,
                         hovertemplate = paste0(reg, " ", level_label(out_lvl()), " ", side, " limit around the ",
                                                reference_label(r), "<br>n = %{x:,}<br>%{y:.2%}<extra></extra>"))
        }
      }
      d <- filter(pts, registry == reg)
      if (nrow(d)) {
        p <- add_markers(p, data = d, x = ~n_total, y = d[[yc$point]], text = hover_text(d), hoverinfo = "text",
                         marker = list(color = cols[[reg]], size = 9, line = list(color = "white", width = 1)),
                         name = labs[[reg]], legendgroup = reg, showlegend = FALSE)
      }
    }
    p <- add_rings(p, pts, "n_total", yc$point)
    {
      o <- labelled(pts)
      if (nrow(o)) p <- add_annotations(p, x = if (input$xlog) log10(o$n_total) else o$n_total, y = o[[yc$point]],
                                        text = o$label, showarrow = TRUE, arrowhead = 0, ax = 25, ay = -25,
                                        font = list(size = 10, color = unname(cols[o$registry])))
    }
    rng <- dynamic_ranges(bind_rows(cv, ref_curves(), per_curves()), pts, yc, "both", NULL, isTRUE(input$xlog))
    layout(p,
      xaxis = list(title = "Procedure volume (n)", type = if (input$xlog) "log" else "linear",
                   range = if (input$xlog) log10(rng$x) else rng$x, tickformat = ",", zeroline = FALSE),
      yaxis = list(title = yc$lab, range = rng$y, tickformat = if (input$yscale == "ratio") ".2f" else ".1%",
                   zeroline = TRUE, zerolinecolor = "#333"),
      legend = list(title = list(text = "<b>Registry (mean, n)</b>"), orientation = "v"),
      hoverlabel = list(align = "left")) |>
      config(toImageButtonOptions = list(format = "png", filename = plot_filename(), scale = 3), displaylogo = FALSE)
  }

  panel_plot <- function() {
    cv <- curves(); pts <- points(); m <- means() |> filter(registry %in% plot_regs())
    labs <- registry_legend_labels(m, ci = TRUE); cols <- registry_colour(m$registry)
    cm <- cm_curves()
    panels <- lapply(m$registry, function(reg) {
      c1 <- filter(cv, registry == reg); d <- filter(pts, registry == reg)
      p <- plot_ly()
      if (!is.null(cm) && reg %in% cm$registry) {
        e <- cm |> filter(registry == reg) |> group_by(n) |>
          summarise(up_lo = min(upper), up_hi = max(upper), lo_lo = min(lower), lo_hi = max(lower), .groups = "drop")
        for (side in c("up", "lo")) p <- add_ribbons(p, data = e, x = ~n, ymin = e[[paste0(side, "_lo")]],
                                                     ymax = e[[paste0(side, "_hi")]], line = list(width = 0),
                                                     fillcolor = scales::alpha(cols[[reg]], 0.2),
                                                     showlegend = FALSE, hoverinfo = "skip")
      }
      for (lv in levels_()) {
        c2 <- filter(c1, level == lv)
        for (side in c("upper", "lower")) p <- add_lines(p, data = c2, x = ~n, y = c2[[side]],
            line = list(color = cols[[reg]], width = 1.8, dash = line_dash(lv)), showlegend = FALSE,
            hovertemplate = paste0(level_label(lv), " ", side, "<br>n = %{x:,}<br>%{y:.2%}<extra></extra>"))
      }
      rc <- ref_curves()
      if (!is.null(rc)) for (r in unique(rc$reference[rc$registry == reg])) {
        c3 <- filter(rc, registry == reg, reference == r)
        for (side in c("upper", "lower")) p <- add_lines(p, data = c3, x = ~n, y = c3[[side]], opacity = 0.7,
            line = list(color = cols[[reg]], width = 1.2, dash = REF_DASH[[r]]), showlegend = FALSE,
            hovertemplate = paste0(level_label(out_lvl()), " ", side, " around the ", reference_label(r),
                                   "<br>n = %{x:,}<br>%{y:.2%}<extra></extra>"))
      }
      pc <- per_curves()
      if (!is.null(pc) && reg == "LROI") p <- add_period_lines(p, pc, list(upper = "upper", lower = "lower"))
      p <- add_lines(p, x = range(c1$n), y = rep(m$p_mean[m$registry == reg], 2),
                     line = list(color = cols[[reg]], width = 1), showlegend = FALSE, hoverinfo = "skip")
      if (m$reference_used[m$registry == reg] != "mean")
        p <- add_lines(p, x = range(c1$n), y = rep(m$p_ref[m$registry == reg], 2),
                       line = list(color = cols[[reg]], width = 1, dash = "dot"), showlegend = FALSE,
                       hovertemplate = paste0("Reference (", reference_label(ref_sel()), "): %{y:.2%}<extra></extra>"))
      if (nrow(d)) {
        p <- add_markers(p, data = d, x = ~n_total, y = ~p, text = hover_text(d), hoverinfo = "text",
                         marker = list(color = cols[[reg]], size = 9, line = list(color = "white", width = 1)),
                         showlegend = FALSE)
        p <- add_rings(p, d, "n_total", "p")
        o <- labelled(d)
        if (nrow(o)) p <- add_text(p, data = o, x = ~n_total, y = ~p, text = ~label, textposition = "top right",
                                   textfont = list(size = 10, color = cols[[reg]]), showlegend = FALSE,
                                   hoverinfo = "skip", inherit = FALSE)
      }
      p
    })
    # axis ranges are set after subplot(), per axis, so each panel keeps its own window
    axes <- list()
    for (i in seq_along(m$registry)) {
      reg <- m$registry[i]
      rng <- dynamic_ranges(bind_rows(filter(cv, registry == reg),
                                      if (reg == "LROI") per_curves()), filter(pts, registry == reg),
                            y_cols("absolute"), "both", NULL, isTRUE(input$xlog))
      sfx <- if (i == 1) "" else i
      axes[[paste0("xaxis", sfx)]] <- list(type = if (input$xlog) "log" else "linear", tickformat = ",",
                                           range = if (input$xlog) log10(rng$x) else rng$x)
      axes[[paste0("yaxis", sfx)]] <- list(tickformat = ".1%", range = c(max(0, rng$y[1]), rng$y[2]))
    }
    nr <- ceiling(length(panels) / 3)
    titles <- lapply(seq_along(m$registry), function(i) {
      list(text = paste0("<b>", labs[[m$registry[i]]], "</b>"), showarrow = FALSE, font = list(size = 12),
           xref = paste0("x", if (i == 1) "" else i, " domain"), yref = paste0("y", if (i == 1) "" else i, " domain"),
           x = 0.5, y = 1.08, xanchor = "center")
    })
    sp <- subplot(panels, nrows = nr, shareX = FALSE, shareY = FALSE, titleX = FALSE, titleY = FALSE,
                  margin = c(0.04, 0.04, 0.09, 0.06))
    do.call(layout, c(list(sp, showlegend = !is.null(per_curves()), hoverlabel = list(align = "left"), annotations = titles), axes)) |>
      config(toImageButtonOptions = list(format = "png", filename = plot_filename(), scale = 3), displaylogo = FALSE)
  }

  plot_filename <- function() {
    sprintf("%s_%dyr_L%s%s_%s_%s_%s_ref-%s", tolower(PROCEDURE), yr(), detail(),
            if (detail() == 0) "" else paste0("_", input$attr_src, if (detail() == 3) paste0("_", input$attr3) else ""),
            paste0(if (input$family == "__all__") "all" else gsub("[^a-z0-9]+", "_", tolower(input$family)),
                   if (fix_sel() != "all") paste0("_", fix_sel()) else ""),
            gsub("\\.", "", input$level), input$method, ref_sel())
  }

  # --------------------------------------------------------------------------------------- tables
  outliers_hi <- reactive({
    pts <- points()
    outlier_report(pts, out_lvl(), input$method)
  })

  output$outlier_title <- renderText(sprintf("High outliers (above their own registry's %s upper limit): %d",
                                             level_label(out_lvl()), nrow(outliers_hi())))

  output$outliers <- renderTable({
    o <- outliers_hi()
    if (!nrow(o)) return(data.frame(Result = "No devices above the upper limit."))
    o |> transmute(Registry = registry, 
                   Device = label, 
                   `n (procedures)` = scales::comma(n_total),
                   `Revision Rate (%)` = scales::percent(p, 0.01), 
                   `Reference Revision Rate (%)` = scales::percent(p_ref, 0.01),
                   `Upper limit (%)` = scales::percent(upper_limit, 0.01), 
                   `Delta p (% pts)` = scales::percent(delta_p, 0.01),
                   `Excess (revisions; n x delta p)` = sprintf("%.1f", excess))
  }, striped = TRUE, spacing = "s")

  output$summary <- renderTable({
    outlier_summary(outliers_hi(), points()) |>
      transmute(Registry = registry, `Devices plotted` = as.integer(n_devices), `Outlier devices` = as.integer(n_outliers),
                `Procedures with an outlier device (n)` = scales::comma(n_implants),
                `Cumulative excess (revisions)` = sprintf("%.1f", cumulative_excess))
  }, striped = TRUE, spacing = "s")

  output$means <- renderTable({
    m <- means() |> left_join(distinct(all_points(), registry, phi, n_devices_phi), by = "registry")
    m |> transmute(Registry = registry, Report = report_year, Metric = metric_type,
                   Mean = scales::percent(p_mean, 0.01),
                   `95% CI` = ifelse(is.na(p_lcl), "–", sprintf("%s–%s", scales::percent(p_lcl, 0.01), scales::percent(p_ucl, 0.01))),
                   `CI source` = ci_source,
                   Reference = ifelse(reference_used == "mean", "mean", reference_label(reference_used)),
                   `Reference rate` = scales::percent(p_ref, 0.01), n = scales::comma(n_total),
                   # `φ (overdispersion)` = sprintf("%.2f", phi), Source = mean_source)
                   Source = mean_source)
  }, striped = TRUE, spacing = "s")

  output$level_title <- renderText(sprintf("Detail levels: rows kept per registry (%s attributes%s)",
    input$attr_src, if (detail() == 3) paste0(", level 3 = ", input$attr3) else ""))
  output$level_counts <- renderTable({
    r <- rows()
    if (input$family != "__all__") r <- filter(r, family == input$family)
    if (!nrow(r)) return(NULL)
    level_counts(r, attr_rules, 1:3, input$attr3, input$attr_src) |>
      mutate(cell = sprintf("%d / %d (%s)", points, rows_kept, scales::percent(pct_procedures_kept, 1)),
             Level = level_title(level, input$attr3)) |>
      select(Level, registry, cell) |>
      tidyr::pivot_wider(names_from = registry, values_from = cell)
  }, striped = TRUE, spacing = "s")

  output$dl_outliers <- downloadHandler(
    filename = function() paste0(plot_filename(), "_outliers.csv"),
    content = function(file) readr::write_csv(
      mutate(outliers_hi(), detail_level = detail(), attribute_source = input$attr_src) |>
        select(detail_level, attribute_source, registry, family, point = label, device_label, n_total, p, p_mean, reference_used, p_ref,
               upper_limit, delta_p, excess, report_year, table_id, pdf_page), file))
  output$dl_points <- downloadHandler(
    filename = function() paste0(plot_filename(), "_points.csv"),
    content = function(file) readr::write_csv(
      mutate(points(), detail_level = detail(), attribute_source = input$attr_src) |>
        select(detail_level, attribute_source, registry, family, point = label, device_label, n_total, p, lcl, ucl, p_mean, reference_used, p_ref, phi,
             lower_limit, upper_limit, outlier, report_year, table_id, pdf_page, verification), file))
}

shinyApp(ui, server)
