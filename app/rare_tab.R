# Rare joints tab (total ankle, total elbow) ------------------------------------------------------------
#
# Sourced by app.R. Same statistics as the UKA tab and the report (R/rare.R, R/limits.R): registry means
# from config/rare_joint_means.csv, model names from config/rare_device_families.csv. Input ids start "rj_".

RARE_CFG <- load_rare_config(ROOT)
RARE_EXT <- purrr::compact(stats::setNames(lapply(RARE_JOINTS, function(j) load_rare(DATA_DIR, j)), RARE_JOINTS))
RARE_AVAIL <- RARE_JOINTS[RARE_JOINTS %in% names(RARE_EXT)]

# Leave-one-out (R/loo.R, same settings as manuscript_paper1.Rmd): random-effects meta-analysis of the model
# rates (REML, Knapp-Hartung) on the complementary log-log scale, registries with at least LOO_K_MIN models.
# Used two ways: "Lowest leave-one-out mean" as the funnel reference, and leaving out one chosen model per
# registry (click it on the funnel or in the Influence view). Runs are cached per joint and year.
LOO_K_MIN <- 3
LOO_SCALE <- "cloglog"
RJ_HAS_METAFOR <- requireNamespace("metafor", quietly = TRUE)
RJ_LOO_CACHE <- new.env()
rj_loo_run <- function(joint, yr) {
  key <- paste(joint, yr)
  if (!exists(key, envir = RJ_LOO_CACHE, inherits = FALSE))
    assign(key, tryCatch(suppressWarnings(loo_run(RARE_EXT[[joint]], RARE_CFG, yr, k_min = LOO_K_MIN, scale = LOO_SCALE)),
                         error = function(e) NULL), envir = RJ_LOO_CACHE)
  get(key, envir = RJ_LOO_CACHE, inherits = FALSE)
}
RJ_LOO_METHODS <- c("Random effects (metafor), as in the manuscript" = "re",
                    "Subtraction from the published total" = "sub")
RJ_GHOST <- "#9E9E9E"      # funnel around the published mean, drawn when the reference is a leave-one-out mean
RJ_LOO_KEY_SEP <- "||"     # leave-out keys are "REGISTRY||model"

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
        limit_input("rj_level", selected = "0.998"),
        selectInput("rj_method", "Funnel distribution", choices = METHODS),
        radioButtons("rj_reference", "Funnel reference",
                     choices = c(REFERENCES, "Lowest leave-one-out mean" = "loo"), selected = "mean"),
        radioButtons("rj_loo_method", "Leave-one-out mean", choices = RJ_LOO_METHODS, selected = "re"),
        helpText(class = "mt-n2", "Used by the lowest leave-one-out reference and when you leave a model out. ",
                 sprintf("Registries with at least %d models.", LOO_K_MIN)),
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
      div(class = "d-flex align-items-center gap-2 flex-wrap mb-2",
          tags$span(class = "small fw-semibold text-nowrap", "Leave out of the mean:"),
          div(style = "flex: 1; min-width: 300px;", class = "rj-leave",
              selectizeInput("rj_leave_out", NULL, choices = NULL, multiple = TRUE, width = "100%",
                             options = list(placeholder = "none (click a model on the plot or in Influence)",
                                            plugins = list("remove_button")))),
          actionLink("rj_leave_reset", "Reset", class = "small link-primary"),
          tags$style(".rj-leave .form-group { margin-bottom: 0; }
                      #rj_view .nav-link:not(.active) { color: var(--bs-secondary-color); }")),
      navset_card_underline(
        id = "rj_view", full_screen = TRUE,
        nav_panel("Funnel", value = "funnel",
          div(class = "small fw-semibold mb-1", textOutput("rj_title", inline = TRUE)),
          layout_sidebar(fillable = FALSE, border = FALSE, sidebar = user_data_sidebar("rj_ud"),
                         plotlyOutput("rj_funnel", height = "680px"))),
        nav_panel("Influence", value = "influence",
          uiOutput("rj_infl_intro"),
          uiOutput("rj_infl_plot_ui"),
          tableOutput("rj_infl_table")),
        nav_spacer(),
        nav_item(user_data_toggle("rj_ud"))),
      uiOutput("rj_ud_card"),
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
             "LROI reports no model-level rates (funnel only). LROI values are a manual transcription of image tables, checked against the report.")))
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
  rj_levels <- reactive(limit_levels(input$rj_level))
  rj_out_lvl <- reactive(max(rj_levels()))

  rj_all_means <- reactive(rare_means(rj_ext(), RARE_CFG, rj_yr()))

  # ------------------------------------------------------------------------------- leave-one-out
  rj_dx_mode <- reactive(identical(input$rj_joint, "ELBOW") && identical(input$rj_dx, "dx"))
  rj_loo_meth <- reactive(input$rj_loo_method %||% "re")
  rj_run <- reactive(if (RJ_HAS_METAFOR) rj_loo_run(input$rj_joint, rj_yr()) else NULL)
  # one row per model with a leave-one-out mean (chosen method), for the selected registries.
  # AOANJRR elbow models are judged against their diagnosis's rate in per-diagnosis mode, so not there.
  rj_loo_models <- reactive({
    run <- rj_run()
    if (is.null(run) || is.null(run$models) || !nrow(run$models)) return(NULL)
    mm <- filter(run$models, registry %in% (input$rj_registries %||% character()))
    if (rj_dx_mode()) mm <- filter(mm, registry != "AOANJRR")
    if (!nrow(mm)) return(NULL)
    re <- rj_loo_meth() == "re"
    mm |>
      left_join(select(run$summary, registry, k, p_pub, p_re), by = "registry") |>
      mutate(key = paste(registry, model, sep = RJ_LOO_KEY_SEP),
             loo = if (re) loo_re else loo_sub,
             loo_lb = if (re) loo_re_lb else pmax(loo_sub - 1.96 * loo_sub_se, 0),
             loo_ub = if (re) loo_re_ub else loo_sub + 1.96 * loo_sub_se)
  })
  # lowest leave-one-out mean per registry, and the model whose removal gives it
  rj_loo_summary <- reactive({
    mm <- rj_loo_models()
    if (is.null(mm)) return(NULL)
    mm |> group_by(registry) |> slice_min(loo, n = 1, with_ties = FALSE) |> ungroup() |>
      transmute(registry, loo_min = loo, drop_for_min = model, k, p_pub, p_re)
  })

  # leave-out selection: at most one model per registry (the latest one chosen wins)
  rj_key_reg <- function(key) sub("\\|\\|.*$", "", key)
  observe({
    mm <- rj_loo_models()
    ch <- if (is.null(mm)) list() else
      lapply(split(mm, mm$registry), function(d) {
        d <- arrange(d, model); stats::setNames(as.list(d$key), paste0(d$model, " (", d$registry, ")")) })
    sel <- intersect(isolate(input$rj_leave_out) %||% character(), unlist(ch))
    updateSelectizeInput(session, "rj_leave_out", choices = ch, selected = sel, server = FALSE)
  })
  observeEvent(input$rj_leave_out, {
    cur <- input$rj_leave_out
    keep <- cur[!duplicated(rj_key_reg(cur), fromLast = TRUE)]
    if (!identical(keep, cur)) updateSelectizeInput(session, "rj_leave_out", selected = keep)
  }, ignoreNULL = TRUE)
  observeEvent(input$rj_leave_reset, updateSelectizeInput(session, "rj_leave_out", selected = character()))
  rj_toggle_leave <- function(key) {
    if (is.null(key) || !length(key) || !nzchar(key[1])) return(invisible())
    key <- key[1]
    mm <- rj_loo_models()
    if (is.null(mm) || !key %in% mm$key) {
      showNotification(paste0("No leave-one-out mean for ", rj_key_reg(key), " here (it needs at least ", LOO_K_MIN,
                              " models with rates", if (rj_dx_mode()) ", and not in per-diagnosis mode" else "", ")."),
                       type = "message", duration = 5)
      return(invisible())
    }
    cur <- input$rj_leave_out %||% character()
    new <- if (key %in% cur) setdiff(cur, key) else c(cur[rj_key_reg(cur) != rj_key_reg(key)], key)
    updateSelectizeInput(session, "rj_leave_out", selected = new)
  }
  # (warnings suppressed: event_data() warns until the plot it listens to has been drawn)
  rj_click <- function(src) suppressWarnings(event_data("plotly_click", source = src))
  observeEvent(rj_click("rj_funnel"), rj_toggle_leave(rj_click("rj_funnel")$customdata))
  observeEvent(rj_click("rj_infl"), rj_toggle_leave(rj_click("rj_infl")$customdata))
  rj_left <- reactive({
    mm <- rj_loo_models(); sel <- input$rj_leave_out %||% character()
    if (is.null(mm) || !length(sel)) return(mm[0, ])
    filter(mm, key %in% sel)
  })
  output$rj_registries_ui <- renderUI({
    regs <- sort(unique(rj_all_means()$registry))
    checkboxGroupInput("rj_registries", "Registries", choices = regs, selected = regs, inline = TRUE)
  })
  rj_means <- reactive({
    req(input$rj_registries)
    ref <- input$rj_reference
    m <- rj_all_means() |> set_reference(if (ref == "loo") "mean" else ref) |> filter(registry %in% input$rj_registries) |>
      mutate(ref_kind = "standard", left_out = NA_character_)
    meth <- if (rj_loo_meth() == "re") "random effects" else "subtraction"
    if (ref == "loo") {
      s <- rj_loo_summary()
      if (!is.null(s)) {
        i <- match(m$registry, s$registry); has <- !is.na(i)
        m$p_ref[has] <- s$loo_min[i[has]]
        m$reference_used[has] <- paste0("lowest leave-one-out mean (", meth, "; without ", s$drop_for_min[i[has]], ")")
        m$ref_kind[has] <- "loo"
      }
    }
    lo <- rj_left()
    if (!is.null(lo) && nrow(lo)) {
      i <- match(m$registry, lo$registry); has <- !is.na(i)
      m$p_ref[has] <- lo$loo[i[has]]
      m$reference_used[has] <- paste0("mean without ", lo$model[i[has]], " (", meth, ")")
      m$ref_kind[has] <- "left_out"; m$left_out[has] <- lo$model[i[has]]
    }
    m
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
  rj_ud <- user_data_server("rj_ud", means = rj_means,
                            phis = reactive(if (nrow(rj_all_points())) distinct(rj_all_points(), registry, phi) else NULL),
                            level = rj_out_lvl, method = reactive(input$rj_method))
  rj_n_max <- reactive(max(c(rj_points()$n_total, rj_ud()$points$n_total, 100), na.rm = TRUE) * 1.1)
  rj_ud_curves <- reactive(user_curves(rj_ud(), rj_levels(), input$rj_method, rj_n_max()))

  rj_curves <- reactive({
    m <- rj_means()
    req(nrow(m) > 0)
    n_max <- rj_n_max()
    ph <- if (nrow(rj_all_points())) distinct(rj_all_points(), registry, phi) else tibble::tibble(registry = character(), phi = numeric())
    m |>
      left_join(ph, by = "registry") |>
      transmute(registry, p_ref, phi = coalesce(phi, 1), n_max = n_max) |>
      build_limit_curves(rj_levels(), input$rj_method)
  })
  # funnel around the published mean, for registries whose reference is a leave-one-out mean
  rj_ghost <- reactive({
    m <- filter(rj_means(), ref_kind != "standard")
    if (!nrow(m)) return(NULL)
    ph <- if (nrow(rj_all_points())) distinct(rj_all_points(), registry, phi) else tibble::tibble(registry = character(), phi = numeric())
    m |>
      left_join(ph, by = "registry") |>
      transmute(registry, p_centre = p_ref, p_ref = p_mean, phi = coalesce(phi, 1), n_max = rj_n_max()) |>
      build_limit_curves(rj_out_lvl(), input$rj_method)
  })

  rj_outliers <- reactive({
    pts <- rj_points()
    if (!nrow(pts)) return(pts)
    outlier_report(pts, rj_out_lvl(), input$rj_method)
  })

  output$rj_title <- renderText({
    sprintf("%s: %s, %d-year revision, %s limits around the %s (%s)", names(RARE_JOINTS)[RARE_JOINTS == input$rj_joint],
            if ((input$rj_family %||% "__all__") == "__all__") "all models" else input$rj_family, rj_yr(),
            paste(level_label(rj_levels()), collapse = " & "),
            paste0(if (input$rj_reference == "loo") "lowest leave-one-out mean" else reference_label(input$rj_reference),
                   if (nrow(rj_left() %||% data.frame())) ", with a model left out" else ""),
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
    ex <- rj_ext()
    unver <- unique(c(ex$device$registry[ex$device$verification %in% "unverified"],
                      ex$casemix$registry[ex$casemix$verification %in% "unverified"]))
    unver <- intersect(unver, m$registry)
    if (length(unver)) msgs <- c(msgs, paste0("Values not yet verified: ", paste(unver, collapse = ", "), "."))
    if (!RJ_HAS_METAFOR && (input$rj_reference == "loo"))
      msgs <- c(msgs, "Leave-one-out needs the metafor package: install.packages(\"metafor\"). Using the published mean.")
    s <- rj_loo_summary(); lo <- rj_left()
    if (input$rj_reference == "loo" && RJ_HAS_METAFOR) {
      lo_regs <- if (!is.null(lo)) lo$registry else character()
      ss <- if (is.null(s)) s else filter(s, !registry %in% lo_regs)
      if (!is.null(ss) && nrow(ss))
        msgs <- c(msgs, paste0("Funnels centred on the lowest leave-one-out mean (",
          if (rj_loo_meth() == "re") paste0("random effects, ", LOO_SCALE, " scale") else "subtraction from the published total", "): ",
          paste(sprintf("%s %s without %s (published %s)", ss$registry, scales::percent(ss$loo_min, 0.01),
                        htmltools::htmlEscape(ss$drop_for_min), scales::percent(ss$p_pub, 0.01)), collapse = "; "),
          ". It is the tightest reference, so it favours flagging."))
      none <- setdiff(m$registry, c(if (!is.null(s)) s$registry, lo_regs))
      if (length(none)) msgs <- c(msgs, paste0("No leave-one-out mean (fewer than ", LOO_K_MIN, " models with rates",
        if (rj_dx_mode()) ", or AOANJRR per diagnosis" else "", "): ", paste(none, collapse = ", "), " use the published mean."))
    }
    if (!is.null(lo) && nrow(lo))
      msgs <- c(msgs, paste0("Left out of the mean: ",
        paste(sprintf("%s <b>%s</b> (mean without it %s, published %s)", lo$registry, htmltools::htmlEscape(lo$model),
                      scales::percent(lo$loo, 0.01), scales::percent(lo$p_pub, 0.01)), collapse = "; "),
        ". Click the model again or use Reset to put it back."))
    if (any(rj_means()$ref_kind != "standard"))
      msgs <- c(msgs, paste0("Grey = the ", level_label(rj_out_lvl()), " funnel around the published mean, for comparison. ",
                             "Black square = model left out of the mean."))
    if (rj_dx_mode() && RJ_HAS_METAFOR && (input$rj_reference == "loo" || !is.null(lo)))
      msgs <- c(msgs, "AOANJRR per diagnosis: each stem keeps its diagnosis's ET9 rate; leave-one-out applies to the pooled view only.")
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
    p <- plot_ly(source = "rj_funnel")
    gh <- rj_ghost()
    if (!is.null(gh)) for (reg in unique(gh$registry)) {
      g1 <- filter(gh, registry == reg)
      for (side in c("upper", "lower")) {
        p <- add_lines(p, data = g1, x = ~n, y = g1[[yc[[side]]]], line = list(color = RJ_GHOST, width = 1.3),
                       name = "Funnel around the published mean", legendgroup = "ghost",
                       showlegend = side == "upper" && reg == unique(gh$registry)[1],
                       hovertemplate = paste0(reg, " ", level_label(rj_out_lvl()), " ", side,
                                              " limit around the published mean<br>n = %{x:,}<br>%{y:.2%}<extra></extra>"))
      }
    }
    for (reg in m$registry) {
      for (lv in rj_levels()) {
        c1 <- filter(cv, registry == reg, level == lv)
        for (side in c("upper", "lower")) {
          p <- add_lines(p, data = c1, x = ~n, y = c1[[yc[[side]]]],
                         line = list(color = cols[[reg]], width = 1.8, dash = limit_dash(lv)),
                         name = labs[[reg]], legendgroup = reg, showlegend = side == "upper" && lv == max(rj_levels()),
                         hovertemplate = paste0(reg, " ", level_label(lv), " ", side, " limit<br>n = %{x:,}<br>%{y:.2%}<extra></extra>"))
        }
      }
      d <- filter(pts, registry == reg)
      if (nrow(d)) p <- add_markers(p, data = d, x = ~n_total, y = d[[yc$point]], text = rj_hover(d), hoverinfo = "text",
                                    customdata = paste(d$registry, d$label, sep = RJ_LOO_KEY_SEP),
                                    marker = list(color = cols[[reg]], size = 9, line = list(color = "white", width = 1)),
                                    name = labs[[reg]], legendgroup = reg, showlegend = FALSE)
    }
    lo <- rj_left()
    if (!is.null(lo) && nrow(lo) && nrow(pts)) {
      lp <- semi_join(pts, lo, by = c("registry", "label" = "model"))
      if (nrow(lp)) p <- add_trace(p, x = lp$n_total, y = lp[[yc$point]], type = "scatter", mode = "markers",
                                   marker = list(symbol = "square-open", size = 22, color = "#000000",
                                                 line = list(width = 2.5, color = "#000000")),
                                   name = "Left out of the mean", hoverinfo = "skip", inherit = FALSE)
    }
    u <- rj_ud(); uc <- rj_ud_curves()
    if (!is.null(u$ref)) p <- add_user_funnel(p, uc, yc, rj_levels(), rj_out_lvl(), user_legend_label(u$ref),
                                              limit_dash)
    p <- add_user_points(p, u$points, yc$point, rj_yr(), rj_out_lvl())
    ring_pts <- bind_rows(if (nrow(pts)) select(pts, n_total, all_of(yc$point), outlier),
                          if (!is.null(u$points)) select(u$points, n_total, all_of(yc$point), outlier))
    for (kind in c("high", "low")) {
      r <- ring_pts[ring_pts$outlier %in% kind, ]
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
    uo <- u$points
    if (!is.null(uo) && input$rj_labels != "none") {
      if (input$rj_labels == "outliers") uo <- uo[uo$outlier != "none", ]
      if (nrow(uo)) p <- add_annotations(p, x = if (isTRUE(input$rj_xlog)) log10(uo$n_total) else uo$n_total,
                                         y = uo[[yc$point]], text = paste0("<b>", htmltools::htmlEscape(uo$label), "</b>"), showarrow = TRUE,
                                         arrowhead = 0, ax = -22, ay = -22, font = list(size = 10, color = USER_COLOUR))
    }
    all_pts <- bind_rows(if (nrow(pts)) pts, u$points)
    rng <- dynamic_ranges(bind_rows(cv, uc, gh), if (nrow(all_pts)) all_pts else NULL, yc, "both",
                          if (nrow(all_pts)) NULL else 50, isTRUE(input$rj_xlog))
    layout(p,
      xaxis = list(title = "Procedure volume (n)", type = if (isTRUE(input$rj_xlog)) "log" else "linear",
                   range = if (isTRUE(input$rj_xlog)) log10(rng$x) else rng$x, tickformat = ",", zeroline = FALSE),
      yaxis = list(title = sub("registry mean", "reference", yc$lab), range = rng$y,
                   tickformat = if (input$rj_yscale == "ratio") ".2f" else ".1%", zeroline = TRUE, zerolinecolor = "#333"),
      legend = list(title = list(text = "<b>Registry (mean, 95% CI, n)</b>")), hoverlabel = list(align = "left")) |>
      config(toImageButtonOptions = list(format = "png", scale = 3,
                                         filename = sprintf("%s_%dyr_%s", tolower(input$rj_joint), rj_yr(), input$rj_reference)),
             displaylogo = FALSE) |>
      event_register("plotly_click")
  })

  # ------------------------------------------------------------------------------- Influence view
  output$rj_infl_intro <- renderUI({
    run <- rj_run(); mm <- rj_loo_models()
    if (!RJ_HAS_METAFOR) return(div(class = "alert alert-warning small",
      "The Influence view needs the metafor package: install.packages(\"metafor\")."))
    skipped <- setdiff(input$rj_registries %||% character(), if (is.null(mm)) character() else unique(mm$registry))
    div(class = "small text-muted mb-2",
      tags$p(class = "mb-1", "Each dot is the registry mean recomputed without that model (",
             if (rj_loo_meth() == "re") paste0("random-effects meta-analysis, REML with Knapp-Hartung, ", LOO_SCALE,
                                                " scale, 95% CI") else "published total minus the model, binomial 95% CI",
             "). Black line = published mean", if (rj_loo_meth() == "re") "; blue dashed = random-effects mean of all models" else "",
             ". The pink ring marks the lowest leave-one-out mean. Dot size ~ the model's share of the registry's procedures."),
      tags$p(class = "mb-0", tags$b("Click a dot"), " to leave that model out; the funnel moves to that mean (see the Funnel tab).",
             if (length(skipped)) paste0(" No leave-one-out for ", paste(skipped, collapse = ", "), " (fewer than ", LOO_K_MIN,
                                         " models with rates", if (rj_dx_mode()) ", or AOANJRR per diagnosis" else "", ").") else ""))
  })
  output$rj_infl_plot_ui <- renderUI({
    mm <- rj_loo_models()
    if (is.null(mm)) return(NULL)
    h <- 70 * length(unique(mm$registry)) + 24 * nrow(mm) + 60
    plotlyOutput("rj_infl_plot", height = paste0(max(h, 280), "px"))
  })
  output$rj_infl_plot <- renderPlotly({
    mm <- rj_loo_models(); req(!is.null(mm))
    regs <- sort(unique(mm$registry)); cols <- registry_colour(regs)
    sel <- input$rj_leave_out %||% character()
    re <- rj_loo_meth() == "re"
    panels <- lapply(regs, function(reg) {
      d <- filter(mm, registry == reg) |> arrange(desc(loo))
      d$model_f <- factor(d$model, levels = d$model)
      low <- d[which.min(d$loo), ]
      share <- d$share_n
      p <- plot_ly(source = "rj_infl")
      yl <- c(-0.5, nrow(d) - 0.5)
      p <- add_lines(p, x = rep(d$p_pub[1], 2), y = yl, line = list(color = "#222222", width = 1.5), showlegend = FALSE,
                     hovertemplate = paste0(reg, " published mean: ", scales::percent(d$p_pub[1], 0.01), "<extra></extra>"))
      if (re) p <- add_lines(p, x = rep(d$p_re[1], 2), y = yl, line = list(color = "#2166AC", width = 1.5, dash = "dash"),
                             showlegend = FALSE,
                             hovertemplate = paste0(reg, " random-effects mean: ", scales::percent(d$p_re[1], 0.01), "<extra></extra>"))
      p <- add_markers(p, x = d$loo, y = as.numeric(d$model_f) - 1, customdata = d$key,
                       error_x = list(type = "data", symmetric = FALSE, array = d$loo_ub - d$loo, arrayminus = d$loo - d$loo_lb,
                                      color = cols[[reg]], thickness = 1.2, width = 3),
                       marker = list(size = 7 + 18 * sqrt(pmin(share, 1)), color = ifelse(d$key %in% sel, "#000000", cols[[reg]]),
                                     line = list(color = ifelse(d$key %in% sel, "#000000", "white"), width = 1)),
                       text = sprintf(paste0("<b>%s</b> (%s)<br>n = %s (%s of the registry)<br>Model rate: %s<br>",
                                             "Mean without it: %s (%s\u2013%s)<br>Change vs published mean: %s pp<br>",
                                             "Cook's distance: %.2f%s<br><i>Click to %s</i>"),
                                      htmltools::htmlEscape(d$model), reg, scales::comma(d$n_total), scales::percent(share, 1),
                                      scales::percent(d$p, 0.01), scales::percent(d$loo, 0.01), scales::percent(d$loo_lb, 0.01),
                                      scales::percent(d$loo_ub, 0.01), sprintf("%+.2f", 100 * (d$loo - d$p_pub)),
                                      d$cook_d, ifelse(d$influential, " (influential)", ""),
                                      ifelse(d$key %in% sel, "put it back", "leave it out")),
                       hoverinfo = "text", showlegend = FALSE)
      p <- add_markers(p, x = low$loo, y = as.numeric(low$model_f) - 1, hoverinfo = "skip", showlegend = FALSE,
                       marker = list(symbol = "circle-open", size = 7 + 18 * sqrt(pmin(low$share_n, 1)) + 10,
                                     color = "#E7298A", line = list(color = "#E7298A", width = 2.5)))
      layout(p, yaxis = list(tickmode = "array", tickvals = seq_len(nrow(d)) - 1, ticktext = htmltools::htmlEscape(d$model),
                             range = yl, title = list(text = paste0("<b>", reg, "</b> (k = ", nrow(d), ")"), font = list(size = 12)),
                             automargin = TRUE),
             xaxis = list(tickformat = ".1%", title = list(text = "Registry mean without the model", font = list(size = 11))))
    })
    ks <- unname(vapply(regs, function(r) sum(mm$registry == r), 0L))   # unnamed: subplot() fails on named heights
    subplot(panels, nrows = length(panels), heights = (ks + 2) / sum(ks + 2), shareX = FALSE, titleX = TRUE, titleY = TRUE,
            margin = c(0, 0, 0.07, 0.03)) |>
      layout(hoverlabel = list(align = "left")) |>
      config(displaylogo = FALSE, toImageButtonOptions = list(format = "png", scale = 3,
             filename = sprintf("%s_%dyr_leave_one_out_%s", tolower(input$rj_joint), rj_yr(), rj_loo_meth()))) |>
      event_register("plotly_click")
  })
  output$rj_infl_table <- renderTable({
    mm <- rj_loo_models(); req(!is.null(mm))
    sel <- input$rj_leave_out %||% character()
    low <- rj_loo_summary()
    mm |>
      arrange(registry, loo) |>
      mutate(is_low = paste(registry, model) %in% paste(low$registry, low$drop_for_min)) |>
      transmute(Registry = registry,
                Model = paste0(ifelse(key %in% sel, "<b>", ""), htmltools::htmlEscape(model), ifelse(key %in% sel, "</b> (left out)", ""),
                               ifelse(is_low, " <span style='color:#E7298A'>&#9679; lowest</span>", "")),
                `n (procedures)` = scales::comma(n_total), `Share of registry n` = scales::percent(share_n, 0.1),
                `Model rate (%)` = scales::percent(p, 0.01),
                `Mean without it (%)` = scales::percent(loo, 0.01),
                `95% CI` = sprintf("%s\u2013%s", scales::percent(loo_lb, 0.01), scales::percent(loo_ub, 0.01)),
                `Change vs published mean (pp)` = sprintf("%+.2f", 100 * (loo - p_pub)),
                `Cook's distance (random effects)` = sprintf("%.2f", cook_d),
                `Influential (metafor)` = ifelse(influential, "yes", ""))
  }, striped = TRUE, spacing = "s", sanitize.text.function = identity)

  output$rj_ud_card <- renderUI(if (!is.null(rj_ud()$points)) user_points_card("rj_ud_table", "rj_ud_dl", "rj_ud_title"))
  output$rj_ud_title <- renderText({
    d <- rj_ud()$points
    sprintf("Your models vs %s limits: %d above the upper limit (of %d comparisons)",
            level_label(rj_out_lvl()), sum(d$outlier == "high"), nrow(d))
  })
  output$rj_ud_table <- renderTable(user_points_table(req(rj_ud()$points)), striped = TRUE, spacing = "s",
                                    sanitize.text.function = identity)
  output$rj_ud_dl <- downloadHandler(
    filename = function() rj_file("your_models"),
    content = function(file) readr::write_csv(user_points_csv(rj_ud()$points), file))

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
    s <- rj_loo_summary()
    lo_txt <- if (is.null(s)) rep("\u2013", nrow(m)) else {
      i <- match(m$registry, s$registry)
      ifelse(is.na(i), "\u2013", sprintf("%s (without %s)", scales::percent(s$loo_min[i], 0.01), s$drop_for_min[i]))
    }
    m |> transmute(Registry = registry, Report = report_year, Metric = metric_type, Mean = scales::percent(p_mean, 0.01),
                   `95% CI` = ifelse(is.na(p_lcl), "–", sprintf("%s–%s", scales::percent(p_lcl, 0.01), scales::percent(p_ucl, 0.01))),
                   `Lowest leave-one-out mean` = lo_txt,
                   `Reference rate` = scales::percent(p_ref, 0.01), Reference = reference_used, n = scales::comma(n_total),
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
