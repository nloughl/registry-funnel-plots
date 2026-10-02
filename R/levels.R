# Implant detail levels ------------------------------------------------------------------------------
#
# Registries list devices at different levels of detail (NJR splits Oxford by peg and fixation, SIRIS
# only by fixation, LROI by cemented/uncemented table...). To compare like with like, device rows are
# pooled to a common level of detail within each registry (after the RSA approach of
# doi:10.2340/17453674.2026.45293, but rows that lack the needed detail are EXCLUDED rather than
# moved to an "other" category):
#
#   level 0  as listed      one point per registry row (no pooling)
#   level 1  brand          one point per brand per registry (config/device_families.csv)
#   level 2  + fixation     cemented / cementless / hybrid; rows with unknown fixation are excluded
#   level 3  + attribute    level 2 + bearing (fixed / mobile), femoral material or tibial type;
#                           rows missing either are excluded
#
# Attributes come from config/device_attributes.csv. With attribute_source = "stated" only details
# printed in the report count; "design" also fills gaps from general implant knowledge (source = design).
#
# Pooling: n = sum of the rows' n; p = n-weighted mean of the rows' KM estimates. This is an
# approximation (the rows' censoring patterns differ) and the pooled point has no CI.

DETAIL_LEVELS <- c("As listed" = 0, "1: Brand" = 1, "2: Brand + fixation" = 2, "3: Brand + fixation + attribute" = 3)
LEVEL3_ATTRS  <- c("Bearing (fixed / mobile)" = "bearing", "Femoral material" = "material",
                   "Tibial type (all-poly / metal-backed)" = "tibia")

load_attribute_rules <- function(path) {
  readr::read_csv(path, comment = "#", show_col_types = FALSE, col_types = readr::cols(.default = "c")) |>
    dplyr::mutate(field = dplyr::coalesce(field, ""))
}

# first matching rule's value for each element of `text` (NA if none)
first_match <- function(text, rules) {
  out <- rep(NA_character_, length(text))
  for (i in seq_len(nrow(rules))) {
    hit <- is.na(out) & !is.na(text) & grepl(rules$pattern[i], text, ignore.case = TRUE, perl = TRUE)
    out[hit] <- rules$value[i]
  }
  out
}

combine_fixation <- function(fem, tib) {
  dplyr::case_when(
    is.na(fem) | is.na(tib) ~ NA_character_,
    fem == tib ~ fem,
    TRUE ~ "hybrid")
}

#' Adds brand, fixation, bearing, material, tibia and a *_source column for each ("stated" / "design").
#' `d` needs device_label, femoral, tibial, section, table_key and `family` (from assign_family()).
derive_attributes <- function(d, rules, source = c("stated", "design")) {
  source <- match.arg(source)
  fem_txt <- dplyr::coalesce(d$femoral, d$device_label)
  tib_txt <- dplyr::coalesce(d$tibial, d$device_label)
  d$brand <- d$family
  for (att in unique(rules$attribute)) {
    r <- dplyr::filter(rules, attribute == att)
    st <- dplyr::filter(r, source == "stated"); dz <- dplyr::filter(r, source == "design")
    comp_scopes <- intersect(c("femoral", "tibial"), unique(r$scope))

    # 1. row-level statements (section heading, table, whole label)
    val <- rep(NA_character_, nrow(d))
    for (i in which(st$scope == "row")) {
      txt <- d[[st$field[i]]]
      hit <- is.na(val) & !is.na(txt) & grepl(st$pattern[i], txt, ignore.case = TRUE, perl = TRUE)
      val[hit] <- st$value[i]
    }
    src <- ifelse(is.na(val), NA_character_, "stated")

    # 2. component names
    comp <- function(rr, scope, txt) first_match(txt, dplyr::filter(rr, scope == !!scope))
    fem <- comp(st, "femoral", fem_txt); tib <- comp(st, "tibial", tib_txt)
    fem_s <- !is.na(fem); tib_s <- !is.na(tib)
    if (source == "design") {
      fem <- dplyr::coalesce(fem, comp(dz, "femoral", fem_txt))
      tib <- dplyr::coalesce(tib, comp(dz, "tibial", tib_txt))
      brand_val <- first_match(d$brand, dplyr::filter(dz, scope == "brand"))
      if (length(comp_scopes) == 2) {             # fixation: fill each missing component from the brand
        fem <- dplyr::coalesce(fem, brand_val); tib <- dplyr::coalesce(tib, brand_val)
      }
    }
    cv <- if (length(comp_scopes) == 2) combine_fixation(fem, tib) else
          if (identical(comp_scopes, "femoral")) fem else if (identical(comp_scopes, "tibial")) tib else
          rep(NA_character_, nrow(d))
    cs <- if (length(comp_scopes) == 2) fem_s & tib_s else
          if (identical(comp_scopes, "femoral")) fem_s else if (identical(comp_scopes, "tibial")) tib_s else FALSE
    fill <- is.na(val) & !is.na(cv)
    val[fill] <- cv[fill]
    src[fill] <- ifelse(cs[fill], "stated", "design")

    # 3. brand-level design knowledge for anything still missing
    if (source == "design") {
      bv <- first_match(d$brand, dplyr::filter(dz, scope == "brand"))
      fill <- is.na(val) & !is.na(bv)
      val[fill] <- bv[fill]; src[fill] <- "design"
    }
    d[[att]] <- val
    d[[paste0(att, "_source")]] <- src
  }
  d
}

#' The key each row is pooled on at a level; NA = row lacks the detail and is excluded.
level_key <- function(d, level, level3_attr = "bearing") {
  if (level == 0) return(d$label)
  k1 <- d$brand
  if (level == 1) return(k1)
  k2 <- ifelse(is.na(d$fixation), NA, paste(k1, d$fixation, sep = " | "))
  if (level == 2) return(k2)
  a3 <- d[[level3_attr]]
  ifelse(is.na(k2) | is.na(a3), NA, paste(k2, a3, sep = " | "))
}

why_excluded <- function(d, level, level3_attr) {
  dplyr::case_when(
    is.na(d$brand) ~ "brand not matched",
    level >= 2 & is.na(d$fixation) ~ "fixation not known",
    level >= 3 & is.na(d[[level3_attr]]) ~ paste(level3_attr, "not known"),
    TRUE ~ NA_character_)
}

#' Pool device points (output of device_points()) to a detail level.
#' Returns the pooled points (same columns the plots and outlier tables use) with attribute
#' `excluded` = the rows dropped for missing detail.
pool_to_level <- function(points, level = 1, level3_attr = "bearing", rules = NULL, source = "stated") {
  if (!"brand" %in% names(points)) points <- derive_attributes(points, rules, source)
  if (level == 0) {
    out <- dplyr::mutate(points, detail_level = 0L, n_rows = 1L, key = label)
    attr(out, "excluded") <- points[0, ]
    return(out)
  }
  points$key <- level_key(points, level, level3_attr)
  points$excluded_because <- why_excluded(points, level, level3_attr)
  keep <- dplyr::filter(points, !is.na(key))
  src_col <- c(NA, NA, "fixation_source", paste0(level3_attr, "_source"))[level + 1]
  out <- keep |>
    dplyr::group_by(registry, key) |>
    dplyr::summarise(
      family = dplyr::first(brand), brand = dplyr::first(brand),
      fixation = if (level >= 2) dplyr::first(fixation) else NA_character_,
      attr3 = if (level >= 3) dplyr::first(.data[[level3_attr]]) else NA_character_,
      detail_source = if (level >= 2) {
        s <- c(fixation_source, if (level >= 3) .data[[src_col]])
        if (any(s == "design", na.rm = TRUE)) "stated + design" else "stated"
      } else "stated",
      p = sum(p * n_total) / sum(n_total),
      lcl = if (dplyr::n() == 1) dplyr::first(lcl) else NA_real_,
      ucl = if (dplyr::n() == 1) dplyr::first(ucl) else NA_real_,
      n_total = sum(n_total),
      n_rows = dplyr::n(),
      device_label = paste(unique(device_label), collapse = " ; "),
      p_ref = dplyr::first(p_ref),
      p_mean = dplyr::first(p_mean),
      reference_used = dplyr::first(reference_used),
      report_year = dplyr::first(report_year),
      metric_type = dplyr::first(metric_type),
      table_id = paste(unique(table_id), collapse = "/"),
      pdf_page = paste(sort(unique(pdf_page)), collapse = "/"),
      verification = if (any(verification %in% "unverified")) "unverified" else dplyr::first(verification),
      low_at_risk = all(low_at_risk %in% TRUE),
      .groups = "drop") |>
    dplyr::mutate(
      estimate = p * 100,
      label = gsub(" | ", " · ", key, fixed = TRUE),
      label = ifelse(n_rows > 1, paste0(label, " [", n_rows, "]"), label),
      dev = p - p_ref, ratio = p / p_ref, detail_level = as.integer(level))
  attr(out, "excluded") <- dplyr::filter(points, is.na(key))
  out
}

#' Table-2-style counts: per registry, rows and procedures kept / excluded at each level.
level_counts <- function(points, rules, levels = 1:3, level3_attr = "bearing",
                         sources = c("stated", "design")) {
  purrr::map_dfr(sources, function(src) {
    att <- derive_attributes(points, rules, src)
    purrr::map_dfr(levels, function(lv) {
      key <- level_key(att, lv, level3_attr)
      att |>
        dplyr::mutate(kept = !is.na(key), key = key) |>
        dplyr::group_by(registry) |>
        dplyr::summarise(
          rows = dplyr::n(), rows_kept = sum(kept),
          points = dplyr::n_distinct(key[kept]),
          procedures = sum(n_total), procedures_kept = sum(n_total[kept]),
          .groups = "drop") |>
        dplyr::mutate(level = lv, attribute_source = src)
    })
  }) |>
    dplyr::mutate(pct_procedures_kept = procedures_kept / procedures) |>
    dplyr::relocate(attribute_source, level, registry)
}

level_title <- function(level, level3_attr = "bearing") {
  c("rows as listed", "level 1 (brand)", "level 2 (brand + fixation)",
    paste0("level 3 (brand + fixation + ", level3_attr, ")"))[level + 1]
}
