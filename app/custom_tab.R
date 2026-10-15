# Custom tab: funnel plot from the user's own data only ------------------------------------------------
#
# Sourced by app.R (after user_data.R, whose parser and plot helpers it reuses). Input ids start "cu_".
# Starts with one empty data set; "+ Add data set" adds more (up to CU_MAX_SETS). Each set has
#   - a name and a colour (fixed per set, see CU_COLOURS)
#   - model points, "name, revision rate (%), n" (typed, pasted or loaded from a CSV)
#   - a funnel reference: the n-weighted mean of its points (sum of p x n / sum of n), a rate the user
#     enters, or none (points only, judged against other sets' funnels)
# A set with an entered rate and no points draws a funnel only (e.g. a registry or national average).
# Points are judged against their own set's funnel, or against every funnel ("Judge each set's points
# against"). Limits and outliers use the same functions as the registry tabs (R/limits.R).
# Nothing is stored: the values live in the browser session only.

CU_MAX_SETS <- 8
# one colour per set, none of them a registry colour (R/palette.R)
CU_COLOURS <- c("#222222", "#7B3294", "#1B7837", "#B35806", "#2166AC", "#C51B7D", "#8C510A", "#5E3C99")

cu_set_panel <- function(i) {
  id <- function(x) paste0("cu_", x, "_", i)
  accordion_panel(
    title = paste("Set", i), value = paste0("set_", i),
    icon = tags$span(style = paste0("display:inline-block;width:12px;height:12px;border-radius:2px;background:",
                                    CU_COLOURS[(i - 1) %% length(CU_COLOURS) + 1])),
    textInput(id("name"), "Name", placeholder = "e.g. Hospital A, NJR average"),
    textAreaInput(id("points"), "Model points: name, revision rate (%), n", rows = 5, resize = "vertical",
                  placeholder = "e.g.  Oxford, 4.8, 1200"),
    fileInput(id("csv"), NULL, accept = c(".csv", ".tsv", ".txt"), buttonLabel = "Load CSV",
              placeholder = "or paste rows above"),
    uiOutput(id("msg")),
    radioButtons(id("ref_src"), "Funnel reference",
                 choices = c("n-weighted mean of the points" = "pooled", "A rate I enter" = "entered", "None (points only)" = "none"),
                 selected = "pooled"),
    conditionalPanel(sprintf("input.%s == 'entered'", id("ref_src")),
      layout_columns(col_widths = c(6, 6), gap = "0.5rem",
        numericInput(id("mean"), "Reference (%)", value = NA, min = 0, max = 100, step = 0.01),
        numericInput(id("mean_n"), "Procedures (n)", value = NA, min = 1, step = 1))),
    if (i > 1) actionLink(id("remove"), "Remove this set", class = "small text-danger")
  )
}

custom_ui <- function() {
  nav_panel(
    "Custom",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 360,
        tags$h6(class = "mb-1", "Data sets"),
        helpText("Each set gets its own colour and funnel. A set with a reference rate and no points draws a ",
                 "funnel only, e.g. a registry average."),
        accordion(id = "cu_sets", multiple = TRUE, open = "set_1", cu_set_panel(1)),
        div(class = "d-flex justify-content-between align-items-center my-2",
            actionButton("cu_add", tagList(icon("plus"), "Add data set"), class = "btn-sm btn-outline-primary"),
            downloadLink("cu_template", "CSV template", class = "small link-primary")),
        hr(),
        radioButtons("cu_compare", "Judge each set's points against",
                     choices = c("Its own funnel" = "own", "Every funnel" = "all"), selected = "own", inline = TRUE),
        radioButtons("cu_level", "Control limit",
                     choices = c("99.8%" = "0.998", "95%" = "0.95", "Both" = "both"), selected = "both", inline = TRUE),
        selectInput("cu_method", "Funnel distribution", choices = METHODS),
        radioButtons("cu_yscale", "Y scale",
                     choices = c("Revision rate (funnels at their own references)" = "absolute",
                                 "Standardised: deviation from own reference (all funnels at 0)" = "difference",
                                 "Standardised: ratio to own reference (all funnels at 1)" = "ratio"),
                     selected = "absolute"),
        conditionalPanel("input.cu_yscale != 'absolute' && input.cu_compare == 'all'",
          radioButtons("cu_std_points", "Points on the standardised plot",
                       choices = c("Once, vs their own set's reference" = "own",
                                   "Once per funnel (every comparison)" = "each"),
                       selected = "own")),
        checkboxInput("cu_xlog", "Log x axis", FALSE),
        radioButtons("cu_labels", "Label models on the plot",
                     choices = c("Outliers" = "outliers", "All" = "all", "None" = "none"), selected = "all", inline = TRUE),
        hr(),
        downloadButton("cu_dl", "Points and limits (CSV)", class = "btn-sm"),
        helpText("Use the camera icon on the plot to save a PNG.")
      ),
      uiOutput("cu_notes"),
      card(full_screen = TRUE, card_header(textOutput("cu_title")),
           plotlyOutput("cu_funnel", height = "640px")),
      card(card_header(textOutput("cu_table_title")), tableOutput("cu_table"),
           card_footer(helpText(uiOutput("cu_table_note", inline = TRUE)))),
      card(card_header("Summary by data set"), tableOutput("cu_summary"),
           card_footer(helpText(
             "Against each set's own funnel. An outlier is above the upper control limit at its own n. ",
             "Procedures = sum of n over high outliers; excess = sum of n x (rate - reference) over high outliers. ",
             "Screening numbers, not precise estimates.")))
    )
  )
}

custom_server <- function(input, output, session) {

  set_ids <- reactiveVal(1L)       # ids of the sets currently shown
  next_id <- reactiveVal(2L)
  cu_colour <- function(i) CU_COLOURS[(i - 1) %% length(CU_COLOURS) + 1]
  cu_in <- function(x, i) input[[paste0("cu_", x, "_", i)]]

  # per-set wiring: CSV load, parse message, panel title, remove link
  wire_set <- function(i) {
    observeEvent(cu_in("csv", i), {
      f <- cu_in("csv", i)$datapath
      if (is.null(f) || file.size(f) > 1e6) return()
      updateTextAreaInput(session, paste0("cu_points_", i),
                          value = paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
    })
    output[[paste0("cu_msg_", i)]] <- renderUI({
      pr <- parse_user_points(cu_in("points", i)); msgs <- pr$errors
      m <- cu_in("mean", i)
      if (identical(cu_in("ref_src", i), "entered") && !is.na(m %||% NA) && !(m > 0 && m < 100))
        msgs <- c("Reference must be between 0 and 100 (%).", msgs)
      ok <- if (nrow(pr$points)) sprintf("%d model point(s) read.", nrow(pr$points))
      if (!length(msgs) && is.null(ok)) return(NULL)
      div(class = paste("small mb-2", if (length(msgs)) "text-danger" else "text-body-secondary"),
          HTML(paste(c(ok, if (length(msgs) > 6) c(head(msgs, 6), sprintf("... and %d more.", length(msgs) - 6)) else msgs),
                     collapse = "<br>")))
    })
    observeEvent(cu_in("name", i), {
      nm <- trimws(cu_in("name", i) %||% "")
      accordion_panel_update("cu_sets", paste0("set_", i), title = if (nzchar(nm)) nm else paste("Set", i))
    }, ignoreInit = TRUE)
    if (i > 1) observeEvent(cu_in("remove", i), {
      accordion_panel_remove("cu_sets", paste0("set_", i))
      set_ids(setdiff(set_ids(), i))
    }, ignoreInit = TRUE, once = TRUE)
  }
  wire_set(1L)

  observeEvent(input$cu_add, {
    if (length(set_ids()) >= CU_MAX_SETS) {
      showNotification(sprintf("Up to %d data sets.", CU_MAX_SETS), type = "warning")
      return()
    }
    i <- next_id()
    accordion_panel_close("cu_sets", values = TRUE)   # keep the sidebar short: only the new set open
    accordion_panel_insert("cu_sets", cu_set_panel(i))
    accordion_panel_open("cu_sets", paste0("set_", i))
    set_ids(c(set_ids(), i)); next_id(i + 1L)
    wire_set(i)
  })

  output$cu_template <- downloadHandler(
    filename = "custom_models_template.csv",
    content = function(file) writeLines(c("model,revision_rate_pct,n", "Example model A,4.8,1200",
                                          "Example model B,3.1,450", "Example model C,7.9,300"), file))

  # where points go on the standardised scales; "each" only applies with Every funnel
  std_points <- reactive(if (input$cu_compare == "all") input$cu_std_points %||% "own" else "own")
  cu_levels  <- reactive(if (input$cu_level == "both") c(0.95, 0.998) else as.numeric(input$cu_level))
  cu_out_lvl <- reactive(max(cu_levels()))

  # one entry per shown set: name, colour, parsed points, reference (or NULL)
  sets <- reactive({
    lapply(set_ids(), function(i) {
      nm <- trimws(cu_in("name", i) %||% ""); nm <- if (nzchar(nm)) nm else paste("Set", i)
      pts <- parse_user_points(cu_in("points", i))$points
      src <- cu_in("ref_src", i) %||% "pooled"
      m <- cu_in("mean", i); mn <- cu_in("mean_n", i)
      ref <- if (src == "entered" && isTRUE(is.finite(m) && m > 0 && m < 100)) {
        list(p_ref = m / 100, n_total = if (isTRUE(is.finite(mn) && mn > 0)) round(mn) else NA_real_, source = "entered")
      } else if (src == "pooled" && nrow(pts) && sum(pts$p * pts$n_total) > 0) {
        list(p_ref = pool(pts$p, pts$n_total), n_total = sum(pts$n_total), source = "pooled")
      }
      list(id = i, name = nm, colour = cu_colour(i), src = src, points = pts, ref = ref)
    })
  })

  # set names must be unique: they key the funnels, legend and table columns
  set_names <- reactive({
    s <- sets(); nm <- make.unique(vapply(s, `[[`, "", "name"), sep = " #")
    stats::setNames(nm, vapply(s, function(x) as.character(x$id), ""))
  })

  refs <- reactive({
    s <- sets(); nm <- set_names()
    dplyr::bind_rows(lapply(s, function(x) if (!is.null(x$ref))
      tibble::tibble(set_id = x$id, registry = nm[[as.character(x$id)]], p_ref = x$ref$p_ref,
                     n_total = x$ref$n_total, source = x$ref$source, colour = x$colour, phi = 1)))
  })

  # one row per model x reference it is judged against
  cu_points <- reactive({
    s <- sets(); r <- refs(); nm <- set_names()
    if (!nrow(r)) return(NULL)
    rows <- dplyr::bind_rows(lapply(s, function(x) {
      if (!nrow(x$points)) return(NULL)
      comp <- if (input$cu_compare == "all") r else dplyr::filter(r, set_id == x$id)
      if (!nrow(comp)) return(NULL)
      tidyr::crossing(dplyr::mutate(x$points, set = nm[[as.character(x$id)]], set_id = x$id, set_colour = x$colour),
                      dplyr::transmute(comp, against_id = set_id, against = registry, p_ref, reference_used = source,
                                       colour, phi))
    }))
    if (!nrow(rows)) return(NULL)
    lim <- purrr::map2(rows$p_ref, rows$n_total, function(p0, n) funnel_limits(p0, n, cu_out_lvl(), input$cu_method))
    rows |>
      dplyr::mutate(own = set_id == against_id, registry = against,
                    upper_limit = purrr::map_dbl(lim, "upper"), lower_limit = purrr::map_dbl(lim, "lower"),
                    outlier = dplyr::case_when(p > upper_limit ~ "high", p < lower_limit ~ "low", TRUE ~ "none"),
                    dev = p - p_ref, ratio = p / p_ref, delta_p = p - p_ref,
                    excess = ifelse(outlier == "high", n_total * delta_p, NA_real_))
  })

  # points of sets without a funnel of their own, in "own" mode: shown on the rate scale, not judged
  unjudged <- reactive({
    s <- sets(); nm <- set_names()
    dplyr::bind_rows(lapply(s, function(x) {
      if (!nrow(x$points)) return(NULL)
      if (!is.null(x$ref) && input$cu_compare == "own") return(NULL)
      if (input$cu_compare == "all" && nrow(refs())) return(NULL)
      dplyr::mutate(x$points, set = nm[[as.character(x$id)]], set_colour = x$colour)
    }))
  })

  # right end of the funnels: the largest model n (funnel-only plots run to 10,000)
  cu_n_max <- reactive({
    n <- unlist(lapply(sets(), function(x) x$points$n_total))
    if (length(n)) max(n, 100) * 1.1 else 1e4
  })
  cu_curves <- reactive({
    r <- refs()
    if (!nrow(r)) return(NULL)
    build_limit_curves(dplyr::transmute(r, set_id, registry, p_ref, phi, n_max = cu_n_max()), cu_levels(), input$cu_method)
  })

  output$cu_notes <- renderUI({
    s <- sets(); r <- refs(); msgs <- character()
    n_pts <- sum(vapply(s, function(x) nrow(x$points), 0L))
    if (!nrow(r) && !n_pts) msgs <- "Enter model points and/or a reference rate in a data set on the left to draw a funnel."
    else if (!nrow(r)) msgs <- "No funnel yet: give at least one data set a reference (the mean of its points or a rate you enter)."
    for (x in s) {
      if (identical(x$ref$source, "pooled") && nrow(x$points) < 5)
        msgs <- c(msgs, paste0(htmltools::htmlEscape(x$name), ": funnel centred on the n-weighted mean of only ",
                               nrow(x$points), " model(s) (", scales::percent(x$ref$p_ref, 0.01),
                               "); it is largely set by the biggest ones. An external reference rate may be more meaningful."))
      if (x$src == "none" && nrow(x$points) && input$cu_compare == "own" && input$cu_yscale == "absolute")
        msgs <- c(msgs, paste0(htmltools::htmlEscape(x$name), ": no funnel of its own, so its points are only judged ",
                               "with 'Every funnel'. They are shown as open diamonds on the rate scale."))
      if (is.null(x$ref) && nrow(x$points) && input$cu_yscale != "absolute" && std_points() == "own")
        msgs <- c(msgs, paste0(htmltools::htmlEscape(x$name), ": no reference of its own, so its points cannot be ",
                               "standardised and are not plotted here (see the rate scale or the table)."))
    }
    if (input$cu_yscale != "absolute" && nrow(r)) {
      msgs <- c(msgs, if (std_points() == "own") paste0(
        "Standardised: every funnel is centred on ", if (input$cu_yscale == "ratio") "1" else "0",
        " (drawn relative to its own reference, so only the widths differ), and each model is plotted as its rate ",
        if (input$cu_yscale == "ratio") "divided by" else "minus", " its own set's reference.",
        if (input$cu_compare == "all") " The comparisons with the other funnels are in the table." else "")
      else "Every funnel: each model is drawn once per funnel, relative to that funnel's reference, in that funnel's colour.")
    }
    if (!length(msgs)) return(NULL)
    div(class = "alert alert-secondary py-2 small", HTML(paste(msgs, collapse = "<br>")))
  })

  output$cu_title <- renderText({
    r <- refs()
    sprintf("Custom: %d funnel(s), %s limits (%s)", nrow(r), paste(level_label(cu_levels()), collapse = " & "),
            names(METHODS)[METHODS == input$cu_method])
  })

  cu_hover <- function(d) {
    sprintf(paste0("<b>%s</b> (%s)<br>n = %s<br>Revision: %s<br>Judged against: %s (%s, %s)<br>",
                   "\u0394p vs reference: %s<br>%s limits: %s \u2013 %s<br>%s"),
            htmltools::htmlEscape(d$label), htmltools::htmlEscape(d$set), scales::comma(d$n_total),
            scales::percent(d$p, 0.01), htmltools::htmlEscape(d$against),
            ifelse(d$reference_used == "pooled", "n-weighted mean", "entered"), scales::percent(d$p_ref, 0.01),
            scales::percent(d$delta_p, 0.01), level_label(cu_out_lvl()),
            scales::percent(d$lower_limit, 0.01), scales::percent(d$upper_limit, 0.01),
            c(high = "<b style='color:#c0392b'>HIGH OUTLIER</b>", low = "<b style='color:#1e8449'>Low outlier</b>",
              none = "Within limits")[d$outlier])
  }

  output$cu_funnel <- renderPlotly({
    r <- refs(); cv <- cu_curves(); pts <- cu_points(); un <- unjudged()
    yc <- y_cols(input$cu_yscale); xlog <- isTRUE(input$cu_xlog); absolute <- input$cu_yscale == "absolute"
    cfg <- function(p) config(p, displaylogo = FALSE,
                              toImageButtonOptions = list(format = "png", scale = 3, filename = "custom_funnel_plot"))
    if (is.null(cv) && (!nrow(un) || !absolute)) {
      return(plotly_empty(type = "scatter", mode = "markers") |>
        layout(xaxis = list(title = "Procedure volume (n)", range = c(0, 1000), tickformat = ",", showgrid = TRUE),
               yaxis = list(title = sub("registry mean", "reference", yc$lab), tickformat = if (absolute) ".1%" else if (input$cu_yscale == "ratio") ".2f" else ".1%",
                            range = if (absolute) c(0, 0.1) else if (input$cu_yscale == "ratio") c(0, 2) else c(-0.05, 0.05)),
               annotations = list(list(text = "Your funnels will appear here", showarrow = FALSE, x = 0.5, y = 0.5,
                                       xref = "paper", yref = "paper", font = list(size = 14, color = "#888")))) |>
        cfg())
    }
    p <- plot_ly()
    dash <- function(lv) if (lv >= 0.998) "solid" else "dash"
    for (k in seq_len(nrow(r))) {
      rk <- r[k, ]
      c1 <- dplyr::filter(cv, set_id == rk$set_id)
      for (lv in cu_levels()) {
        c2 <- dplyr::filter(c1, level == lv)
        for (side in c("upper", "lower")) {
          p <- add_lines(p, x = c2$n, y = c2[[yc[[side]]]], line = list(color = rk$colour, width = 2, dash = dash(lv)),
                         name = htmltools::htmlEscape(user_legend_label(rk)), legendgroup = paste0("set", rk$set_id),
                         showlegend = side == "upper" && lv == max(cu_levels()),
                         hovertemplate = paste0(htmltools::htmlEscape(rk$registry), " ", level_label(lv), " ", side,
                                                " limit<br>n = %{x:,}<br>%{y:.2%}<extra></extra>"))
        }
      }
      if (absolute) p <- add_lines(p, x = range(c1$n), y = rep(rk$p_ref, 2), line = list(color = rk$colour, width = 1),
                                   legendgroup = paste0("set", rk$set_id), showlegend = FALSE,
                                   hovertemplate = paste0(htmltools::htmlEscape(rk$registry), " reference: ",
                                                          scales::percent(rk$p_ref, 0.01), "<extra></extra>"))
    }
    # which point rows to draw: on the rate scale a model has one position, so draw it once (against its
    # own funnel if it has one). Standardised: once, as the deviation from its own set's reference (the
    # default), or once per funnel ("each", Every funnel mode only)
    shown <- NULL
    each <- !absolute && std_points() == "each"
    if (!is.null(pts)) {
      shown <- if (absolute) {
        pts |> dplyr::group_by(set_id, label) |> dplyr::arrange(dplyr::desc(own), .by_group = TRUE) |>
          dplyr::slice(1) |> dplyr::ungroup()
      } else if (each) pts else dplyr::filter(pts, own)
      if (!nrow(shown)) shown <- NULL
    }
    if (!is.null(shown)) {
      fill <- if (each) shown$colour else shown$set_colour
      p <- add_markers(p, x = shown$n_total, y = shown[[yc$point]], text = cu_hover(shown), hoverinfo = "text",
                       name = if (absolute || each) "Models (\u25c6)" else "Models (\u25c6), vs own reference", legendgroup = "pts",
                       marker = list(symbol = "diamond", size = 12, color = fill,
                                     line = list(color = ifelse(shown$own, "white", shown$set_colour), width = 1.5)))
      for (kind in c("high", "low")) {
        rr <- shown[shown$outlier == kind, ]
        if (nrow(rr)) p <- add_trace(p, x = rr$n_total, y = rr[[yc$point]], type = "scatter", mode = "markers",
                                     marker = list(symbol = "circle-open", size = 20, color = c(high = "#D62728", low = "#2CA02C")[[kind]],
                                                   line = list(width = 2.5, color = c(high = "#D62728", low = "#2CA02C")[[kind]])),
                                     name = if (kind == "high") "Above upper limit" else "Below lower limit",
                                     hoverinfo = "skip", inherit = FALSE)
      }
      lab <- switch(input$cu_labels, all = shown, none = shown[0, ], shown[shown$outlier != "none", ])
      if (nrow(lab)) p <- add_annotations(p, x = if (xlog) log10(lab$n_total) else lab$n_total, y = lab[[yc$point]],
                                          text = htmltools::htmlEscape(lab$label), showarrow = TRUE, arrowhead = 0,
                                          ax = 22, ay = -22, font = list(size = 11, color = unname(lab$set_colour)))
    }
    # points of sets without a funnel to be judged against (rate scale only)
    if (nrow(un) && absolute) {
      p <- add_markers(p, x = un$n_total, y = un$p, name = "Not judged (no funnel)", legendgroup = "unjudged",
                       text = sprintf("<b>%s</b> (%s)<br>n = %s<br>Revision: %s<br>Not judged: no funnel",
                                      htmltools::htmlEscape(un$label), htmltools::htmlEscape(un$set),
                                      scales::comma(un$n_total), scales::percent(un$p, 0.01)), hoverinfo = "text",
                       marker = list(symbol = "diamond-open", size = 12, color = un$set_colour,
                                     line = list(color = un$set_colour, width = 2)))
      if (input$cu_labels == "all") p <- add_annotations(p, x = if (xlog) log10(un$n_total) else un$n_total, y = un$p,
                                                         text = htmltools::htmlEscape(un$label), showarrow = TRUE, arrowhead = 0,
                                                         ax = 22, ay = -22, font = list(size = 11, color = unname(un$set_colour)))
    }
    rng_pts <- dplyr::bind_rows(shown, if (absolute && nrow(un)) dplyr::mutate(un, dev = NA_real_, ratio = NA_real_))
    rng <- dynamic_ranges(if (is.null(cv)) tibble::tibble(n = c(20, cu_n_max()), lower = 0, upper = max(un$p)) else cv,
                          if (nrow(rng_pts)) rng_pts else NULL, yc, "both", if (nrow(rng_pts)) NULL else 100, xlog)
    if (absolute || input$cu_yscale == "ratio") rng$y[1] <- max(0, rng$y[1])   # rates and ratios are >= 0
    # the line all standardised funnels are centred on (0 for the difference, 1 for the ratio)
    if (!absolute) p <- add_lines(p, x = rng$x * c(if (xlog) 1 else 0, 1), y = rep(yc$ref, 2), inherit = FALSE,
                                  line = list(color = "#333", width = 1.5), showlegend = FALSE, hoverinfo = "skip")
    layout(p,
      xaxis = list(title = "Procedure volume (n)", type = if (xlog) "log" else "linear",
                   range = if (xlog) log10(rng$x) else rng$x, tickformat = ",", zeroline = FALSE),
      yaxis = list(title = switch(input$cu_yscale, absolute = "Cumulative revision",
                                  difference = if (each) "Deviation from each funnel's reference" else "Deviation from own reference (rate \u2212 reference)",
                                  ratio = if (each) "Ratio to each funnel's reference" else "Ratio to own reference"),
                   range = rng$y, tickformat = if (input$cu_yscale == "ratio") ".2f" else ".1%",
                   zeroline = FALSE),
      legend = list(title = list(text = if (absolute) "<b>Funnel (reference, n)</b>"
                                        else paste0("<b>Funnel (reference, n)</b><br><i>all centred on ",
                                                    if (input$cu_yscale == "ratio") "1" else "0", "</i>"))),
      hoverlabel = list(align = "left")) |>
      cfg()
  })

  # ---------------------------------------------------------------------------------------- tables
  output$cu_table_title <- renderText({
    pts <- cu_points()
    if (is.null(pts)) return("Models")
    sprintf("Models vs %s limits: %d above the upper limit, %d below the lower limit (of %d comparisons)",
            level_label(cu_out_lvl()), sum(pts$outlier == "high"), sum(pts$outlier == "low"), nrow(pts))
  })
  output$cu_table_note <- renderUI({
    if (input$cu_compare == "all") "One column per funnel: status and the upper limit at the model's n. The CSV has the lower limits, delta p and excess."
    else "Each model is judged against its own set's funnel. Switch to 'Every funnel' to compare it with the other sets as well."
  })
  output$cu_table <- renderTable({
    pts <- cu_points()
    if (is.null(pts)) return(data.frame(Result = "No models judged yet."))
    status <- c(high = "<b style='color:#c0392b'>HIGH</b>", low = "<b style='color:#1e8449'>low</b>", none = "within")
    if (input$cu_compare == "all") {
      return(pts |>
        dplyr::mutate(col = sprintf("%s<br><span class='fw-normal text-muted'>ref %s</span>", htmltools::htmlEscape(against),
                                    scales::percent(p_ref, 0.01)),
                      cell = sprintf("%s%s<br><span class='text-muted'>limit %s</span>", status[outlier],
                                     ifelse(own, " (own)", ""), scales::percent(upper_limit, 0.01))) |>
        dplyr::select(Set = set, Model = label, `n (procedures)` = n_total, `Revision Rate (%)` = p, col, cell) |>
        tidyr::pivot_wider(names_from = col, values_from = cell, values_fill = "") |>
        dplyr::mutate(Set = htmltools::htmlEscape(Set), Model = htmltools::htmlEscape(Model),
                      `n (procedures)` = scales::comma(`n (procedures)`),
                      `Revision Rate (%)` = scales::percent(`Revision Rate (%)`, 0.01)))
    }
    pts |>
      dplyr::arrange(set_id, dplyr::desc(outlier == "high"), dplyr::desc(delta_p)) |>
      dplyr::transmute(Set = htmltools::htmlEscape(set), Model = htmltools::htmlEscape(label),
                       `n (procedures)` = scales::comma(n_total), `Revision Rate (%)` = scales::percent(p, 0.01),
                       `Reference (%)` = paste0(scales::percent(p_ref, 0.01), ifelse(reference_used == "pooled", " (mean)", "")),
                       `Lower limit (%)` = scales::percent(lower_limit, 0.01),
                       `Upper limit (%)` = scales::percent(upper_limit, 0.01),
                       `Delta p (% pts)` = scales::percent(delta_p, 0.01),
                       `Excess (revisions)` = ifelse(is.na(excess), "-", sprintf("%.1f", excess)),
                       Status = c(high = "<b style='color:#c0392b'>HIGH OUTLIER</b>",
                                  low = "<b style='color:#1e8449'>Low outlier</b>", none = "Within limits")[outlier])
  }, striped = TRUE, spacing = "s", sanitize.text.function = identity)

  output$cu_summary <- renderTable({
    pts <- cu_points(); r <- refs(); s <- sets(); nm <- set_names()
    if (!length(s)) return(NULL)
    own <- if (is.null(pts)) tibble::tibble(set_id = integer(), outlier = character(), n_total = numeric(), excess = numeric())
           else dplyr::filter(pts, own)
    rows <- lapply(s, function(x) {
      d <- dplyr::filter(own, set_id == x$id)
      hi <- d[d$outlier %in% "high", ]
      data.frame(Set = htmltools::htmlEscape(nm[[as.character(x$id)]]),
                 Reference = if (is.null(x$ref)) "none" else paste0(scales::percent(x$ref$p_ref, 0.01),
                                                                    if (x$ref$source == "pooled") " (mean of points)" else " (entered)"),
                 `Models` = nrow(x$points), `High outliers` = if (is.null(x$ref)) "-" else as.character(nrow(hi)),
                 `Low outliers` = if (is.null(x$ref)) "-" else as.character(sum(d$outlier %in% "low")),
                 `Procedures with a high-outlier model (n)` = if (is.null(x$ref)) "-" else scales::comma(sum(hi$n_total)),
                 `Cumulative excess (revisions)` = if (is.null(x$ref)) "-" else sprintf("%.1f", sum(hi$excess)),
                 check.names = FALSE)
    })
    do.call(rbind, rows)
  }, striped = TRUE, spacing = "s", sanitize.text.function = identity)

  output$cu_dl <- downloadHandler(
    filename = function() "custom_funnel_points.csv",
    content = function(file) {
      pts <- cu_points()
      if (is.null(pts)) return(readr::write_csv(tibble::tibble(), file))
      readr::write_csv(dplyr::transmute(pts, set, model = label, n_total, p, judged_against = against, own_funnel = own,
                                        reference = reference_used, p_ref, method = input$cu_method, level = cu_out_lvl(),
                                        lower_limit, upper_limit, delta_p, excess, outlier), file)
    })
}
