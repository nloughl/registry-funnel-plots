# Rare joints tab (total ankle, total elbow) ------------------------------------------------------------
#
# Sourced by app.R. Same statistics as the UKA tab and the report (R/rare.R, R/limits.R): registry means
# from config/rare_joint_means.csv, model names from config/rare_device_families.csv. Input ids start "rj_".

RARE_CFG <- load_rare_config(ROOT)
RARE_EXT <- purrr::compact(stats::setNames(lapply(RARE_JOINTS, function(j) load_rare(DATA_DIR, j)), RARE_JOINTS))
RARE_AVAIL <- RARE_JOINTS[RARE_JOINTS %in% names(RARE_EXT)]

rare_ui <- function() {
  nav_panel(
    "Rare joints",
    if (!length(RARE_AVAIL)) {
      div(class = "alert alert-warning m-3",
          "No rare-joint outputs found next to the UKA data (", DATA_DIR, "). Run regextract with ",
          tags$code("--procedure ANKLE"), " and ", tags$code("--procedure ELBOW"), ".")
    } else layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 320,
        radioButtons("rj_joint", "Joint", choices = RARE_AVAIL, selected = RARE_AVAIL[1], inline = TRUE),
        uiOutput("rj_year_ui"),
        uiOutput("rj_family_ui"),
        uiOutput("rj_registries_ui"),
        radioButtons("rj_level", "Control limit",
                     choices = c("99.8%" = "0.998", "95%" = "0.95", "Both" = "both"), selected = "0.998", inline = TRUE),
        selectInput("rj_method", "Funnel distribution", choices = METHODS),
        radioButtons("rj_reference", "Funnel reference", choices = REFERENCES, selected = "mean"),
        conditionalPanel("input.rj_joint == 'ELBOW'",
          radioButtons("rj_dx", "AOANJRR elbow models (printed per diagnosis)",
                       choices = c("Pooled per stem (fracture + OA + RA)" = "pooled",
                                   "Per diagnosis, vs that diagnosis's ET9 rate" = "dx"),
                       selected = "pooled")),
        radioButtons("rj_yscale", "Y scale", choices = c("Difference from reference" = "difference", "Ratio to reference" = "ratio"),
                     selected = "difference", inline = TRUE),
        checkboxInput("rj_xlog", "Log x axis", TRUE),
        radioButtons("rj_labels", "Label models on the plot",
                     choices = c("Outliers" = "outliers", "All" = "all", "None" = "none"), selected = "all", inline = TRUE),
        hr(),
        downloadButton("rj_dl_outliers", "Outliers (CSV)", class = "btn-sm"),
        downloadButton("rj_dl_points", "All points (CSV)", class = "btn-sm")
      ),
      uiOutput("rj_notes"),
      card(full_screen = TRUE, card_header(textOutput("rj_title")), plotlyOutput("rj_funnel", height = "680px")),
      layout_columns(
        col_widths = c(7, 5),
        card(card_header(textOutput("rj_outlier_title")), tableOutput("rj_outliers")),
        card(card_header("Cumulative excess revision by registry (outlier models only)"), tableOutput("rj_summary"),
             card_footer(helpText("Same definitions as the knee tab: only models above their own registry's upper limit ",
                                  "count. Procedures = sum of n [procedures]; excess = sum of n × (rate − reference) [revisions].")))
      ),
      card(card_header("Registry means at this follow-up year"), tableOutput("rj_means"),
           card_footer(helpText(
             tags$b("Scope and means (decisions 2026-10):"), " elbow = total elbow replacement only. ",
             "Ankle means: AOANJRR all A11 diagnoses pooled (same population as the A15 models); NJR 3.A3 All cases; ",
             "LROI A014B all primary ankles. Elbow means: AOANJRR ET9 total elbow (fracture + OA + RA) pooled; ",
             "NJR 3.E6 total elbow (acute trauma + elective) pooled; LROI E017B total elbow (n approximate, from E001). ",
             "LROI reports no model-level rates (funnel only). LROI values are a manual transcription, not yet verified.")))
    )
  )
}

rare_server <- function(input, output, session) {
  if (!length(RARE_AVAIL)) return(invisible())

  rj_ext <- reactive(RARE_EXT[[req(input$rj_joint)]])
  rj_years <- reactive(rare_years(rj_ext()))
  output$rj_year_ui <- renderUI({
    yrs <- rj_years()
    selectInput("rj_year", "Follow-up year", choices = yrs, selected = if (3 %in% yrs) 3 else yrs[1])
  })
  rj_yr <- reactive(as.numeric(req(input$rj_year)))
  rj_levels <- reactive(if (input$rj_level == "both") c(0.95, 0.998) else as.numeric(input$rj_level))
  rj_out_lvl <- reactive(max(rj_levels()))

  rj_all_means <- reactive(rare_means(rj_ext(), RARE_CFG, rj_yr()))
  output$rj_registries_ui <- renderUI({
    regs <- sort(unique(rj_all_means()$registry))
    checkboxGroupInput("rj_registries", "Registries", choices = regs, selected = regs, inline = TRUE)
  })
  rj_means <- reactive({
    req(input$rj_registries)
    rj_all_means() |> set_reference(input$rj_reference) |> filter(registry %in% input$rj_registries)
  })
  rj_all_points <- reactive({
    pts <- rare_points(rj_ext(), rj_means(), RARE_CFG, rj_yr(), identical(input$rj_dx, "dx"))
    if (!nrow(pts)) return(pts)
    left_join(pts, estimate_phi(pts), by = "registry") |> mutate(phi = coalesce(phi, 1))
  })
  output$rj_family_ui <- renderUI({
    fams <- sort(unique(rj_all_points()$family))
    sel <- isolate(input$rj_family)
    selectInput("rj_family", "Model", choices = c("All models" = "__all__", fams),
                selected = if (!is.null(sel) && sel %in% fams) sel else "__all__")
  })
  rj_points <- reactive({
    pts <- rj_all_points()
    fam <- input$rj_family %||% "__all__"
    if (fam != "__all__") pts <- filter(pts, family == fam)
    if (!nrow(pts)) return(pts)
    pts <- classify_points(pts, rj_levels(), input$rj_method)
    lim <- purrr::pmap(list(pts$p_ref, pts$n_total, pts$phi),
                       function(p0, n, ph) funnel_limits(p0, n, rj_out_lvl(), input$rj_method, ph))
    pts |>
      mutate(upper_limit = purrr::map_dbl(lim, "upper"), lower_limit = purrr::map_dbl(lim, "lower"),
             outlier = case_when(p > upper_limit ~ "high", p < lower_limit ~ "low", TRUE ~ "none"))
  })
  rj_curves <- reactive({
    m <- rj_means()
    req(nrow(m) > 0)
    n_max <- max(c(rj_points()$n_total, 100), na.rm = TRUE) * 1.1
    ph <- if (nrow(rj_all_points())) distinct(rj_all_points(), registry, phi) else tibble::tibble(registry = character(), phi = numeric())
    m |>
      left_join(ph, by = "registry") |>
      transmute(registry, p_ref, phi = coalesce(phi, 1), n_max = n_max) |>
      build_limit_curves(rj_levels(), input$rj_method)
  })
  rj_outliers <- reactive({
    pts <- rj_points()
    if (!nrow(pts)) return(pts)
    outlier_report(pts, rj_out_lvl(), input$rj_method)
  })

  output$rj_title <- renderText({
    sprintf("%s: %s, %d-year revision, %s limits around the %s (%s)", names(RARE_JOINTS)[RARE_JOINTS == input$rj_joint],
            if ((input$rj_family %||% "__all__") == "__all__") "all models" else input$rj_family, rj_yr(),
            paste(level_label(rj_levels()), collapse = " & "), reference_label(input$rj_reference),
            names(METHODS)[METHODS == input$rj_method])
  })

  output$rj_notes <- renderUI({
    msgs <- character()
    m <- rj_means(); pts <- rj_points()
    only <- rare_mean_only(m, pts)
    if (length(only)) msgs <- c(msgs, paste0("Funnel only (no model-level rates", if ((input$rj_family %||% "__all__") != "__all__") " for this model" else "",
                                             "): ", paste(only, collapse = ", "), "."))
    if (identical(input$rj_joint, "ELBOW") && "AOANJRR" %in% pts$registry)
      msgs <- c(msgs, if (identical(input$rj_dx, "dx"))
        "AOANJRR stems shown per diagnosis; each is judged against its diagnosis's ET9 total-elbow rate (reference options apply to the pooled mean only)."
        else "AOANJRR stems pooled across fracture, OA and RA (n-weighted, no CI), marked [3 dx].")
    if ("LROI" %in% m$registry) msgs <- c(msgs, "LROI values are a manual transcription, not yet verified.")
    if (input$rj_method == "overdispersed") {
      ph <- distinct(rj_all_points(), registry, phi)
      msgs <- c(msgs, paste0("Overdispersion factor φ: ", paste(sprintf("%s %.2f", ph$registry, ph$phi), collapse = ", "),
                             " (few models per registry: φ is imprecise)."))
    }
    if (!length(msgs)) return(NULL)
    div(class = "alert alert-secondary py-2 small", HTML(paste(msgs, collapse = "<br>")))
  })

  rj_hover <- function(d) {
    ci <- ifelse(is.na(d$lcl), " (pooled, no CI)",
                 sprintf(" (%s–%s)", scales::percent(d$lcl / 100, 0.01), scales::percent(d$ucl / 100, 0.01)))
    dx <- ifelse(is.na(d$diagnosis), "", paste0("<br>Diagnosis: ", d$diagnosis))
    sprintf(paste0("<b>%s</b> — %s<br><i>%s</i>%s<br>n = %s<br>%d-yr revision: %s%s<br>",
                   "Reference (%s): %s<br>Δp vs reference: %s<br>Upper %s limit: %s<br>%s<br><i>%s %s, table %s p.%s</i>"),
            d$registry, d$label, d$device_label, dx, scales::comma(d$n_total), rj_yr(),
            scales::percent(d$p, 0.01), ci, d$reference_used, scales::percent(d$p_ref, 0.01),
            scales::percent(d$p - d$p_ref, 0.01), level_label(rj_out_lvl()), scales::percent(d$upper_limit, 0.01),
            c(high = "<b style='color:#c0392b'>HIGH OUTLIER</b>", low = "<b style='color:#1e8449'>Low outlier</b>",
              none = "Within limits")[d$outlier], d$registry, d$report_year, d$table_id, d$pdf_page)
  }

  output$rj_funnel <- renderPlotly({
    m <- rj_means(); cv <- rj_curves(); pts <- rj_points()
    yc <- y_cols(input$rj_yscale)
    labs <- registry_legend_labels(m, ci = TRUE); cols <- registry_colour(m$registry)
    p <- plot_ly()
    for (reg in m$registry) {
      for (lv in rj_levels()) {
        c1 <- filter(cv, registry == reg, level == lv)
        for (side in c("upper", "lower")) {
          p <- add_lines(p, data = c1, x = ~n, y = c1[[yc[[side]]]],
                         line = list(color = cols[[reg]], width = 1.8, dash = if (lv >= 0.998) "solid" else "dash"),
                         name = labs[[reg]], legendgroup = reg, showlegend = side == "upper" && lv == max(rj_levels()),
                         hovertemplate = paste0(reg, " ", level_label(lv), " ", side, " limit<br>n = %{x:,}<br>%{y:.2%}<extra></extra>"))
        }
      }
      d <- filter(pts, registry == reg)
      if (nrow(d)) p <- add_markers(p, data = d, x = ~n_total, y = d[[yc$point]], text = rj_hover(d), hoverinfo = "text",
                                    marker = list(color = cols[[reg]], size = 9, line = list(color = "white", width = 1)),
                                    name = labs[[reg]], legendgroup = reg, showlegend = FALSE)
    }
    for (kind in c("high", "low")) {
      r <- pts[pts$outlier %in% kind, ]
      if (nrow(r)) p <- add_trace(p, x = r$n_total, y = r[[yc$point]], type = "scatter", mode = "markers",
                                  marker = list(symbol = "circle-open", size = 18, color = c(high = "#D62728", low = "#2CA02C")[[kind]],
                                                line = list(width = 2.5, color = c(high = "#D62728", low = "#2CA02C")[[kind]])),
                                  name = if (kind == "high") "Above upper limit" else "Below lower limit",
                                  hoverinfo = "skip", inherit = FALSE)
    }
    lab_pts <- switch(input$rj_labels, all = pts, none = pts[0, ], pts[pts$outlier != "none", ])
    if (nrow(lab_pts)) p <- add_annotations(p, x = if (isTRUE(input$rj_xlog)) log10(lab_pts$n_total) else lab_pts$n_total,
                                            y = lab_pts[[yc$point]], text = lab_pts$label, showarrow = TRUE, arrowhead = 0,
                                            ax = 22, ay = -22, font = list(size = 10, color = unname(cols[lab_pts$registry])))
    rng <- dynamic_ranges(cv, if (nrow(pts)) pts else NULL, yc, "both", if (nrow(pts)) NULL else 50, isTRUE(input$rj_xlog))
    layout(p,
      xaxis = list(title = "Procedure volume (n)", type = if (isTRUE(input$rj_xlog)) "log" else "linear",
                   range = if (isTRUE(input$rj_xlog)) log10(rng$x) else rng$x, tickformat = ",", zeroline = FALSE),
      yaxis = list(title = sub("registry mean", "reference", yc$lab), range = rng$y,
                   tickformat = if (input$rj_yscale == "ratio") ".2f" else ".1%", zeroline = TRUE, zerolinecolor = "#333"),
      legend = list(title = list(text = "<b>Registry (mean, 95% CI, n)</b>")), hoverlabel = list(align = "left")) |>
      config(toImageButtonOptions = list(format = "png", scale = 3,
                                         filename = sprintf("%s_%dyr_%s", tolower(input$rj_joint), rj_yr(), input$rj_reference)),
             displaylogo = FALSE)
  })

  output$rj_outlier_title <- renderText(sprintf("High outliers (above their own registry's %s upper limit): %d",
                                                level_label(rj_out_lvl()), nrow(rj_outliers())))
  output$rj_outliers <- renderTable({
    o <- rj_outliers()
    if (!nrow(o)) return(data.frame(Result = "No models above the upper limit."))
    o |> transmute(Registry = registry, Model = label, `Printed name` = device_label, `n (procedures)` = scales::comma(n_total),
                   `Revision Rate (%)` = scales::percent(p, 0.01), `Reference Revision Rate (%)` = scales::percent(p_ref, 0.01),
                   `Upper limit (%)` = scales::percent(upper_limit, 0.01), `Delta p (% pts)` = scales::percent(delta_p, 0.01),
                   `Excess (revisions; n x delta p)` = sprintf("%.1f", excess))
  }, striped = TRUE, spacing = "s")
  output$rj_summary <- renderTable({
    pts <- rj_points()
    if (!nrow(pts)) return(NULL)
    outlier_summary(rj_outliers(), pts) |>
      transmute(Registry = registry, `Models plotted` = as.integer(n_devices), `Outlier models` = as.integer(n_outliers),
                `Procedures with an outlier model (n)` = scales::comma(n_implants),
                `Cumulative excess (revisions)` = sprintf("%.1f", cumulative_excess))
  }, striped = TRUE, spacing = "s")
  output$rj_means <- renderTable({
    m <- rj_means(); pts <- rj_all_points()
    m |> transmute(Registry = registry, Report = report_year, Metric = metric_type, Mean = scales::percent(p_mean, 0.01),
                   `95% CI` = ifelse(is.na(p_lcl), "–", sprintf("%s–%s", scales::percent(p_lcl, 0.01), scales::percent(p_ucl, 0.01))),
                   `Reference rate` = scales::percent(p_ref, 0.01), n = scales::comma(n_total),
                   Models = purrr::map_int(registry, ~ sum(pts$registry == .x)), Source = mean_source)
  }, striped = TRUE, spacing = "s")

  rj_file <- function(what) sprintf("%s_%dyr_%s_%s.csv", tolower(input$rj_joint), rj_yr(),
                                    gsub("[^a-z0-9]+", "_", tolower(input$rj_family %||% "all")), what)
  output$rj_dl_outliers <- downloadHandler(
    filename = function() rj_file("outliers"),
    content = function(file) readr::write_csv(
      select(rj_outliers(), registry, family, model = label, device_label, diagnosis, n_total, p, p_ref,
             upper_limit, delta_p, excess, report_year, table_id, pdf_page), file))
  output$rj_dl_points <- downloadHandler(
    filename = function() rj_file("points"),
    content = function(file) readr::write_csv(
      select(rj_points(), registry, family, model = label, device_label, diagnosis, n_total, p, lcl, ucl, p_ref,
             reference_used, lower_limit, upper_limit, outlier, report_year, table_id, pdf_page, verification), file))
}
