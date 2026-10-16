# "Your data": user-supplied reference rate and model points, shown alongside the registry data ---------
#
# Sourced by app.R. Used by both tabs (ids "ud" and "rj_ud"). Nothing is stored: the values live in the
# browser session only.
#
# A user can give
#   - a reference: a mean revision rate (%) and, optionally, the number of procedures it is based on.
#     A funnel is drawn around it (black), with the same limits and distribution as the registries.
#   - model points: one per line, "name, revision rate (%), n" (typed, pasted from a spreadsheet, or
#     loaded from a CSV). Each point is judged against every comparator ticked in "Judge your models
#     against": the user's own reference and/or any of the selected registries' references (the registry
#     mean, or the CI bound chosen under "Funnel reference"), at the point's own n.
# User points never enter the registry outlier tables or the cumulative excess summary; they get their
# own table.

USER_COLOUR <- "#222222"   # not a registry colour (see R/palette.R)

# Control limits offered on the Rare joints and Custom tabs. Several can be drawn; outliers are judged
# at the widest one ticked. Plotly line style per level: 99.8% solid, 95% dashed, 90% dotted.
LIMIT_CHOICES <- c("90%" = "0.9", "95%" = "0.95", "99.8%" = "0.998")
limit_dash <- function(lv) if (lv >= 0.998) "solid" else if (lv >= 0.95) "dash" else "dot"
limit_levels <- function(x) { lv <- sort(as.numeric(x)); if (length(lv)) lv else 0.998 }   # none ticked: 99.8%
limit_input <- function(id, selected) {
  tagList(checkboxGroupInput(id, "Control limits", choices = LIMIT_CHOICES, selected = selected, inline = TRUE),
          helpText(class = "mt-n2", "99.8% solid, 95% dashed, 90% dotted. Outliers are judged at the widest limit ticked."))
}
USER_OWN    <- "__own__"   # comparator id for the user's own reference

# --------------------------------------------------------------------------------------------- UI
user_data_sidebar <- function(id) {
  ns <- NS(id)
  sidebar(
    id = ns("sidebar"), position = "right", open = FALSE, width = 340, title = "Your data",
    checkboxInput(ns("show"), "Show your data on the plot", TRUE),
    textInput(ns("label"), "Name", placeholder = "e.g. Hospital A"),
    tags$h6(class = "mt-2 mb-1", "Reference rate (optional)"),
    layout_columns(col_widths = c(6, 6), gap = "0.5rem",
      numericInput(ns("mean"), "Mean revision (%)", value = NA, min = 0, max = 100, step = 0.01),
      numericInput(ns("mean_n"), "Procedures (n)", value = NA, min = 1, step = 1)),
    helpText("Draws a black funnel around this rate, with the limits and distribution chosen on the left."),
    tags$h6(class = "mt-2 mb-1", "Model points (optional)"),
    textAreaInput(ns("points"), "One model per line: name, revision rate (%), n", rows = 6, resize = "vertical",
                  placeholder = "Oxford, 4.8, 1200\nZUK, 3.1, 450"),
    fileInput(ns("csv"), NULL, accept = c(".csv", ".tsv", ".txt"), buttonLabel = "Load CSV",
              placeholder = "or paste rows above"),
    div(class = "small mb-2", downloadLink(ns("template"), "Download a CSV template", class = "link-primary")),
    uiOutput(ns("parse_msg")),
    checkboxGroupInput(ns("against"), "Judge your models against", choices = character()),
    helpText("Each model is compared with each ticked reference at its own n. Registry references follow ",
             "the 'Funnel reference' setting. Your points are not added to the registry tables."),
    actionLink(ns("clear"), "Clear your data", class = "small text-danger")
  )
}

#' Header link that opens / closes the "Your data" sidebar.
user_data_toggle <- function(id) {
  actionLink(NS(id, "toggle"), tagList(icon("plus"), "Add your data"), class = "small text-nowrap link-primary")
}

# ----------------------------------------------------------------------------------------- parsing
#' Parse "name, rate (%), n" lines. Separators: comma, semicolon or tab (pasted from a spreadsheet).
#' With semicolon or tab separators a decimal comma is accepted. The last two fields are the rate and n,
#' so a model name may itself contain commas. A first line whose rate is not a number is read as a header.
parse_user_points <- function(txt) {
  empty <- tibble::tibble(label = character(), p = numeric(), n_total = numeric())
  lines <- trimws(strsplit(txt %||% "", "\r?\n")[[1]])
  keep <- nzchar(lines) & !startsWith(lines, "#")
  lines <- lines[keep]; line_no <- which(keep)
  if (!length(lines)) return(list(points = empty, errors = character()))
  num <- function(x, dec_comma) {
    x <- gsub("[%\\s]", "", x, perl = TRUE)
    if (dec_comma) x <- sub(",", ".", x, fixed = TRUE) else x <- gsub(",", "", x, fixed = TRUE)
    suppressWarnings(as.numeric(x))
  }
  rows <- lapply(seq_along(lines), function(i) {
    l <- lines[i]
    sep <- if (grepl("\t", l)) "\t" else if (grepl(";", l)) ";" else ","
    f <- trimws(strsplit(l, sep, fixed = TRUE)[[1]])
    if (length(f) < 3) return(list(err = sprintf("Line %d: needs name, rate (%%), n.", line_no[i])))
    k <- length(f)
    rate <- num(f[k - 1], sep != ",")
    n <- num(gsub("[.,](?=\\d{3}(\\D|$))", "", f[k], perl = TRUE), sep != ",")   # 1,234 / 1.234 = thousands
    name <- paste(f[seq_len(k - 2)], collapse = if (sep == ",") ", " else " ")
    if (is.na(rate)) return(list(err = if (i == 1) NA_character_   # header line
                                 else sprintf("Line %d: revision rate '%s' is not a number.", line_no[i], f[k - 1])))
    if (rate < 0 || rate >= 100) return(list(err = sprintf("Line %d: revision rate must be between 0 and 100 (%%).", line_no[i])))
    if (is.na(n) || n < 1) return(list(err = sprintf("Line %d: n '%s' must be a positive number.", line_no[i], f[k])))
    if (!nzchar(name)) name <- sprintf("Model %d", line_no[i])
    list(row = tibble::tibble(label = name, p = rate / 100, n_total = round(n)))
  })
  errs <- unlist(lapply(rows, `[[`, "err"))
  pts <- dplyr::bind_rows(empty, lapply(rows, `[[`, "row"))
  pts$label <- make.unique(pts$label, sep = " #")   # the same name twice stays two points
  list(points = pts, errors = errs[!is.na(errs)])
}

# ------------------------------------------------------------------------------------------ server
#' means   reactive: the tab's registry means (registry, p_ref, reference_used), already filtered
#' phis    reactive: registry, phi (overdispersion factors; missing = 1)
#' level   reactive: the outlier control level (e.g. 0.998);  method: reactive funnel distribution
#' Returns a reactive list: label, ref (one-row tibble or NULL), points (one row per model x comparator).
user_data_server <- function(id, means, phis, level, method) {
  moduleServer(id, function(input, output, session) {

    observeEvent(input$toggle, bslib::toggle_sidebar("sidebar"))

    observeEvent(input$csv, {
      f <- input$csv$datapath
      if (is.null(f) || file.size(f) > 1e6) return()
      updateTextAreaInput(session, "points", value = paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
    })

    output$template <- downloadHandler(
      filename = "your_models_template.csv",
      content = function(file) writeLines(c("model,revision_rate_pct,n", "Example model A,4.8,1200",
                                            "Example model B,3.1,450"), file))

    observeEvent(input$clear, {
      updateTextInput(session, "label", value = "")
      updateNumericInput(session, "mean", value = NA)
      updateNumericInput(session, "mean_n", value = NA)
      updateTextAreaInput(session, "points", value = "")
    })

    label <- reactive({ l <- trimws(input$label %||% ""); if (nzchar(l)) l else "Your data" })
    has_mean <- reactive(isTRUE(is.finite(input$mean) && input$mean > 0 && input$mean < 100))
    parsed <- reactive(parse_user_points(input$points))

    # comparator choices: own reference (when given) + the tab's selected registries
    choices <- reactive({
      regs <- means()$registry
      c(if (has_mean()) stats::setNames(USER_OWN, paste0(label(), " (your reference)")),
        stats::setNames(regs, regs))
    })
    observe({
      ch <- choices()
      sel <- isolate(input$against)
      sel <- if (is.null(sel)) unname(ch) else intersect(sel, ch)   # first time: everything ticked
      updateCheckboxGroupInput(session, "against", choices = ch, selected = sel, inline = TRUE)
    })
    # a reference typed in later gets ticked automatically
    observeEvent(has_mean(), if (has_mean())
      updateCheckboxGroupInput(session, "against", selected = union(input$against, USER_OWN)), ignoreInit = TRUE)

    output$parse_msg <- renderUI({
      pr <- parsed(); msgs <- pr$errors
      if (!is.na(input$mean %||% NA) && !has_mean()) msgs <- c("Mean revision must be between 0 and 100 (%).", msgs)
      ok <- if (nrow(pr$points)) sprintf("%d model point(s) read.", nrow(pr$points))
      if (!length(msgs) && is.null(ok)) return(NULL)
      div(class = paste("small mb-2", if (length(msgs)) "text-danger" else "text-body-secondary"),
          HTML(paste(c(ok, if (length(msgs) > 6) c(head(msgs, 6), sprintf("... and %d more.", length(msgs) - 6)) else msgs),
                     collapse = "<br>")))
    })

    ref <- reactive({
      if (!isTRUE(input$show) || !has_mean()) return(NULL)
      tibble::tibble(registry = label(), p_ref = input$mean / 100, p_mean = input$mean / 100,
                     n_total = if (isTRUE(is.finite(input$mean_n) && input$mean_n > 0)) round(input$mean_n) else NA_real_,
                     phi = 1)
    })

    pts <- reactive({
      up <- parsed()$points
      against <- intersect(input$against %||% character(), choices())
      if (!isTRUE(input$show) || !nrow(up) || !length(against)) return(NULL)
      m <- means()
      ph <- phis()
      comps <- dplyr::bind_rows(
        if (USER_OWN %in% against && has_mean())
          tibble::tibble(against = label(), comparator = USER_OWN, p_ref = input$mean / 100,
                         reference_used = "your mean", phi = 1, colour = USER_COLOUR),
        m |> dplyr::filter(registry %in% against) |>
          dplyr::left_join(if (is.null(ph) || !nrow(ph)) tibble::tibble(registry = character(), phi = numeric()) else ph,
                           by = "registry") |>
          dplyr::transmute(against = registry, comparator = registry, p_ref, reference_used,
                           phi = dplyr::coalesce(phi, 1), colour = unname(registry_colour(registry))))
      if (!nrow(comps)) return(NULL)
      out <- tidyr::crossing(up, comps) |> dplyr::filter(!is.na(p_ref))
      lim <- purrr::pmap(list(out$p_ref, out$n_total, out$phi),
                         function(p0, n, ph) funnel_limits(p0, n, level(), method(), ph))
      out |>
        dplyr::mutate(upper_limit = purrr::map_dbl(lim, "upper"), lower_limit = purrr::map_dbl(lim, "lower"),
                      outlier = dplyr::case_when(p > upper_limit ~ "high", p < lower_limit ~ "low", TRUE ~ "none"),
                      dev = p - p_ref, ratio = p / p_ref, delta_p = p - p_ref,
                      excess = ifelse(outlier == "high", n_total * delta_p, NA_real_),
                      registry = against) |>
        dplyr::arrange(label, against)
    })

    reactive(list(label = label(), ref = ref(), points = pts(), show = isTRUE(input$show)))
  })
}

# ------------------------------------------------------------------------------------ plot helpers
#' Registries used as comparators (their funnels must be drawn even without registry points).
user_comparator_regs <- function(ud) {
  if (is.null(ud$points)) character() else setdiff(unique(ud$points$comparator), USER_OWN)
}

user_legend_label <- function(ref) {
  sprintf("%s: %s%s", ref$registry, scales::percent(ref$p_ref, 0.01),
          if (is.na(ref$n_total)) "" else paste0(", n = ", scales::comma(ref$n_total)))
}

user_hover <- function(d, yr, level) {
  sprintf(paste0("<b>%s</b> (your data)<br>n = %s<br>%d-yr revision: %s<br>Judged against: %s (%s)<br>",
                 "Reference rate: %s<br>\u0394p vs reference: %s<br>Upper %s limit: %s<br>%s"),
          htmltools::htmlEscape(d$label), scales::comma(d$n_total), yr, scales::percent(d$p, 0.01),
          ifelse(d$comparator == USER_OWN, "your reference", htmltools::htmlEscape(d$against)), d$reference_used,
          scales::percent(d$p_ref, 0.01), scales::percent(d$delta_p, 0.01),
          level_label(level), scales::percent(d$upper_limit, 0.01),
          c(high = "<b style='color:#c0392b'>HIGH OUTLIER</b>", low = "<b style='color:#1e8449'>Low outlier</b>",
            none = "Within limits")[d$outlier])
}

#' Limit curves around the user's reference.
user_curves <- function(ud, levels, method, n_max) {
  if (is.null(ud$ref)) return(NULL)
  build_limit_curves(dplyr::transmute(ud$ref, registry, p_ref, phi, n_max = n_max), levels, method)
}

#' Funnel lines for the user's reference (black), in the style of the registry funnels.
add_user_funnel <- function(p, uc, ycols, levels, out_level, show_legend_label, line_dash) {
  if (is.null(uc)) return(p)
  for (lv in levels) {
    c1 <- dplyr::filter(uc, level == lv)
    for (side in c("upper", "lower")) {
      p <- add_lines(p, data = c1, x = ~n, y = c1[[ycols[[side]]]],
                     line = list(color = USER_COLOUR, width = 2, dash = line_dash(lv)),
                     name = htmltools::htmlEscape(show_legend_label), legendgroup = "user_ref",
                     showlegend = !is.na(show_legend_label) && side == "upper" && lv == max(levels),
                     hovertemplate = paste0(htmltools::htmlEscape(c1$registry[1]), " ", level_label(lv), " ", side,
                                            " limit<br>n = %{x:,}<br>%{y:.2%}<extra></extra>"))
    }
  }
  p
}

#' User model points: diamonds with a dark outline, filled in the comparator's colour.
add_user_points <- function(p, d, ycol, yr, level, show_legend = TRUE) {
  if (is.null(d) || !nrow(d)) return(p)
  add_markers(p, x = d$n_total, y = d[[ycol]], text = user_hover(d, yr, level), hoverinfo = "text",
              marker = list(symbol = "diamond", size = 12, color = d$colour,
                            line = list(color = USER_COLOUR, width = 1.5)),
              name = "Your models (\u25c6, colour = reference)", legendgroup = "user_pts",
              showlegend = show_legend, inherit = FALSE)
}

#' Compact table of user points: one row per model, one column per reference it is judged against.
#' Each cell gives the upper limit at the model's n and the status (HIGH in red). Full detail is in the CSV.
user_points_table <- function(d) {
  d |>
    dplyr::mutate(label_html = htmltools::htmlEscape(label),
                  ref_col = ifelse(comparator == USER_OWN, paste(htmltools::htmlEscape(against), "(yours)"), against),
                  ref_col = sprintf("%s<br><span class='fw-normal text-muted'>ref %s</span>", ref_col,
                                    scales::percent(p_ref, 0.01)),
                  cell = sprintf("%s<br><span class='text-muted'>limit %s</span>",
                                 c(high = "<b style='color:#c0392b'>HIGH</b>", low = "<b style='color:#1e8449'>low</b>",
                                   none = "within")[outlier], scales::percent(upper_limit, 0.01))) |>
    dplyr::select(Model = label_html, `n (procedures)` = n_total, `Revision Rate (%)` = p, ref_col, cell) |>
    tidyr::pivot_wider(names_from = ref_col, values_from = cell) |>
    dplyr::mutate(`n (procedures)` = scales::comma(`n (procedures)`),
                  `Revision Rate (%)` = scales::percent(`Revision Rate (%)`, 0.01))
}

user_points_csv <- function(d) {
  dplyr::transmute(d, model = label, n_total, p, judged_against = ifelse(comparator == USER_OWN, "your reference", against),
                   reference_used, p_ref, phi, lower_limit, upper_limit, delta_p, excess, outlier)
}

#' Card with the user's points, shown only when there are any.
user_points_card <- function(ns_table, ns_dl, ns_title) {
  card(card_header(class = "d-flex justify-content-between align-items-center",
                   textOutput(ns_title, inline = TRUE),
                   downloadButton(ns_dl, "Your models (CSV)", class = "btn-sm")),
       tableOutput(ns_table),
       card_footer(helpText("Your models are judged with the same limits and distribution as the registry devices, ",
                            "at their own n. Each cell: status and the upper limit at that n. The CSV adds the lower limit, delta p and excess. ",
                            "Your models are not counted in the registry outlier tables or the cumulative excess.")))
}
