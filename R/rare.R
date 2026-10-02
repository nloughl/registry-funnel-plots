# Rare joints: total ankle and total elbow replacement ---------------------------------------------
#
# Data: the regextract outputs/<JOINT>/<JOINT>_long.csv (model rows) and _casemix_long.csv (overall,
# diagnosis, period, class rows). Decisions (2026-10, see README):
#   * elbow = TOTAL elbow only (radial head / hemi / unconfirmed rows are dropped here)
#   * registry means = config/rare_joint_means.csv (n-weighted pools of the listed rows, CI pooled)
#   * AOANJRR elbow models are printed per diagnosis (ET6-ET8): pooled to one point per humeral stem by
#     default; elbow_by_diagnosis = TRUE keeps them apart, each against its own diagnosis's ET9 rate
#   * model names harmonised with config/rare_device_families.csv
# Everything downstream (limits, outliers, plots) is the same code as the UKA report.

RARE_JOINTS <- c("Total ankle" = "ANKLE", "Total elbow" = "ELBOW")

#' Rare-joint extraction outputs. `data_dir` may be the UKA output folder (its parent holds ANKLE/, ELBOW/).
load_rare <- function(data_dir, joint) {
  root <- if (basename(normalizePath(data_dir, mustWork = FALSE)) %in% c("UKA", "TKA", "ANKLE", "ELBOW"))
    dirname(data_dir) else data_dir
  d <- file.path(root, joint)
  if (!file.exists(file.path(d, paste0(joint, "_long.csv")))) return(NULL)
  ext <- load_extraction(d, joint)
  ext$joint <- joint
  ext
}

load_rare_config <- function(root = ".") {
  list(
    means    = readr::read_csv(file.path(root, "config", "rare_joint_means.csv"), comment = "#",
                               show_col_types = FALSE, col_types = readr::cols(.default = "c")),
    families = readr::read_csv(file.path(root, "config", "rare_device_families.csv"), comment = "#",
                               show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  )
}

# rows in scope for the funnels (total elbow only)
rare_scope <- function(d, joint) {
  if (joint == "ELBOW" && "implant_class" %in% names(d))
    d <- dplyr::filter(d, is.na(implant_class) | implant_class == "Total elbow")
  d
}

#' Follow-up years with a model estimate in at least `min_regs` registries.
rare_years <- function(ext, min_regs = 2) {
  ext$device |>
    dplyr::filter(row_type == "device", has_estimate(dplyr::pick(dplyr::everything()))) |>
    dplyr::distinct(registry, time_yr) |>
    dplyr::count(time_yr) |>
    dplyr::filter(n >= min_regs) |>
    dplyr::pull(time_yr) |>
    sort()
}

#' Registry means for a joint at `time_yr`, from config/rare_joint_means.csv. Same columns as registry_means().
rare_means <- function(ext, cfg, time_yr) {
  spec <- dplyr::filter(cfg$means, joint == ext$joint)
  out <- purrr::pmap_dfr(spec, function(joint, registry, kind, table_key, row_types, label_regex, note) {
    src <- if (kind == "device") ext$device else ext$casemix
    lab <- if (kind == "device") src$device_label else src$stratum_label
    rows <- src[src$registry == registry & src$table_key == table_key & src$time_yr == time_yr &
                  src$row_type %in% strsplit(row_types, ";")[[1]] &
                  grepl(label_regex, dplyr::coalesce(lab, ""), ignore.case = TRUE, perl = TRUE), ]
    rows <- rows[!is.na(rows$estimate) & !is.na(rows$n_total), ]
    if (!nrow(rows)) return(NULL)
    mean_row(registry, pool(rows$estimate / 100, rows$n_total), sum(rows$n_total),
             paste0(note, " [", rows$table_id[1], "]"),
             pool_ci(rows$estimate, rows$lcl, rows$ucl, rows$n_total),
             if (nrow(rows) == 1) "published" else "pooled from published CIs (approx.)") |>
      dplyr::mutate(report_year = paste(unique(rows$report_year), collapse = "/"),
                    metric_type = paste(unique(rows$metric_type), collapse = "/"))
  })
  if (!nrow(out)) return(out)
  set_reference(out, "mean")
}

#' Diagnosis-specific means (AOANJRR ET9 total-elbow rows), for elbow_by_diagnosis.
rare_diagnosis_means <- function(ext, time_yr) {
  ext$casemix |>
    dplyr::filter(table_key == "elbow_class_diagnosis", implant_class == "Total Elbow", time_yr == !!time_yr,
                  !is.na(estimate)) |>
    dplyr::transmute(registry, diagnosis, p_dx = estimate / 100)
}

assign_rare_family <- function(d, cfg, joint) {
  fam <- dplyr::filter(cfg$families, .data$joint == !!joint)
  out <- rep(NA_character_, nrow(d))
  for (i in seq_len(nrow(fam))) {
    hit <- is.na(out) & grepl(fam$pattern[i], d$device_label, ignore.case = TRUE, perl = TRUE)
    out[hit] <- fam$family[i]
  }
  d$family <- dplyr::coalesce(out, d$device_label)
  d
}

#' Model points for a joint at `time_yr` (same columns as device_points()).
#' elbow_by_diagnosis: FALSE = AOANJRR ET6-ET8 pooled per stem (n-weighted, no CI);
#'                     TRUE  = one point per stem x diagnosis, judged against that diagnosis's ET9 rate.
rare_points <- function(ext, means, cfg, time_yr, elbow_by_diagnosis = FALSE) {
  d <- ext$device |>
    rare_scope(ext$joint) |>
    dplyr::filter(row_type == "device", time_yr == !!time_yr, !is.na(n_total)) |>
    dplyr::filter(has_estimate(dplyr::pick(dplyr::everything()))) |>
    dplyr::mutate(p = estimate / 100) |>
    assign_rare_family(cfg, ext$joint)
  if (!nrow(d)) return(d)
  by_dx <- !is.na(d$diagnosis) & d$registry == "AOANJRR" & ext$joint == "ELBOW"
  if (any(by_dx) && !elbow_by_diagnosis) {
    pooled <- d[by_dx, ] |>
      dplyr::group_by(registry, family) |>
      dplyr::summarise(
        device_label = paste(unique(device_label), collapse = " ; "),
        p = sum(p * n_total) / sum(n_total), lcl = if (dplyr::n() == 1) dplyr::first(lcl) else NA_real_,
        ucl = if (dplyr::n() == 1) dplyr::first(ucl) else NA_real_, n_total = sum(n_total),
        n_rows = dplyr::n(), diagnosis = paste(unique(diagnosis), collapse = " + "),
        table_id = paste(unique(table_id), collapse = "/"), pdf_page = paste(unique(pdf_page), collapse = "/"),
        report_year = dplyr::first(report_year), metric_type = dplyr::first(metric_type),
        verification = dplyr::first(verification), low_at_risk = FALSE, .groups = "drop") |>
      dplyr::mutate(estimate = p * 100)
    d <- dplyr::bind_rows(d[!by_dx, ] |> dplyr::mutate(n_rows = 1L, pdf_page = as.character(pdf_page)), pooled)
  } else {
    d <- dplyr::mutate(d, n_rows = 1L, pdf_page = as.character(pdf_page))
  }
  d <- d |>
    dplyr::left_join(dplyr::select(means, registry, p_ref, p_mean, reference_used), by = "registry") |>
    dplyr::filter(!is.na(p_ref))
  if (any(by_dx) && elbow_by_diagnosis) {
    dx <- rare_diagnosis_means(ext, time_yr)
    d <- d |>
      dplyr::left_join(dx, by = c("registry", "diagnosis")) |>
      dplyr::mutate(p_ref = dplyr::coalesce(p_dx, p_ref), p_mean = dplyr::coalesce(p_dx, p_mean),
                    reference_used = ifelse(!is.na(p_dx), "diagnosis (ET9)", reference_used)) |>
      dplyr::select(-p_dx)
  }
  d |>
    dplyr::group_by(registry, family, diagnosis_key = if (elbow_by_diagnosis) diagnosis else NA) |>
    dplyr::mutate(dup = dplyr::n() > 1) |>
    dplyr::ungroup() |>
    dplyr::select(-diagnosis_key) |>
    dplyr::mutate(
      family_label = ifelse(dup, device_label, family),
      label = ifelse(!is.na(diagnosis) & registry == "AOANJRR" & elbow_by_diagnosis & ext$joint == "ELBOW",
                     paste0(family_label, " (", sub("Fracture/Dislocation", "fracture", sub("Rheumatoid Arthritis", "RA",
                            sub("Osteoarthritis", "OA", diagnosis))), ")"),
                     ifelse(n_rows > 1, paste0(family_label, " [", n_rows, " dx]"), family_label)),
      dev = p - p_ref, ratio = p / p_ref) |>
    dplyr::select(-dup, -family_label)
}

#' Registries without model-level rates (e.g. LROI): listed so plots/tables can say so.
rare_mean_only <- function(means, points) setdiff(means$registry, unique(points$registry))
