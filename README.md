# registry-funnel-plots

Funnel plots of implant revision across national joint registries, built from the tables extracted by
[`regitry-outlier-detection`](https://github.com/nloughl/regitry-outlier-detection) (`regextract`).
This replaces `registry_funnel_plots_v1.Rmd` and `normal_binom_dist_comparison.R`, which read hand-built Excel sheets.

## Inputs
The extraction repo's CSVs are read directly from `params$data_dir`. The default is `../regitry-outlier-detection/outputs/UKA`, the sibling repo in `code/`.
- `UKA_long.csv`: one row per device × time point. Gives the device points and the registry summary rows.
- `UKA_casemix_long.csv`: AOANJRR KP22 (sex × age) and NJR 3.K6 (sex × age × fixation). Gives the age-sex classes.

## Run
Open `registry-funnel-plots.Rproj` in RStudio and knit `funnel_plots.Rmd`. Alternatively:
```r
rmarkdown::render("funnel_plots.Rmd")                                   # defaults: 3-yr, 95% + 99.8%
rmarkdown::render("funnel_plots.Rmd", params = list(time_yr = 5, limit_method = "exact"))
```
Packages: dplyr, tidyr, purrr, readr, ggplot2, ggrepel, scales, stringr, rmarkdown, knitr; metafor for `loo_influence.Rmd`.
Figures are written to `figures/<procedure>_<t>yr_<slug>.png`.

| Parameter | Default | Meaning |
|---|---|---|
| `time_yr` | 3 | Time point (every registry reports 1, 3, 5 and 10 years) |
| `levels` | 0.95, 0.998 | Control-limit levels available to the plots |
| `limit_method` | normal | `normal` (as v1) or `exact` (binomial with Spiegelhalter 2005 interpolation) |
| `y_scale` | difference | `difference`: registry mean set to 0. `ratio`: registry mean set to 1 |
| `families` | Oxford, ZUK, Journey, BalanSys | Device families that get their own plot |
| `detail_level` | 0 | Points in the main catalogue: 0 = rows as listed, 1 = brand, 2 = + fixation, 3 = + attribute |
| `detail_levels` | 1, 2, 3 | Levels compared in the "Detail levels" section |
| `level3_attr` | bearing | Level-3 attribute: `bearing`, `material` (femoral) or `tibia` (all-poly / metal-backed) |
| `fixation_filter` | all | At detail levels 2–3, show only `cemented` or `cementless` models (hybrid rows are left out) |
| `label_points` | outliers | Point labels: `outliers` (outside the plot's limit), `all`, or `none` |
| `reference` | mean | Funnel centre: `mean`, or `lcl` / `ucl` = lower / upper 95% CI of the registry mean |
| `attribute_source` | stated | `stated`: only details printed in the report. `design`: also fills gaps from implant knowledge |

## Reference: registry mean or its 95% CI
The funnels can be centred on the registry mean or on the lower or upper 95% CI of that mean. The limits are rebuilt around the chosen value.
- **Upper CI:** conservative, so fewer high outliers.
- **Lower CI:** sensitive, so more high outliers.

Everything is standardised to the chosen reference: y = rate − reference, so the funnel stays centred on 0 and the points shift by (bound − mean). Outlier status and excess revisions, n × (rate − reference), use the chosen reference.

Where each registry's CI comes from:
- **NJR, SIRIS, EPRD and AOANJRR:** published (EPRD 2025 Table 57, all UKA; AOANJRR 2025 Table KP18, all UKA for OA, modern prostheses).
- **LROI (cemented + uncemented):** pooled from the parts' published CIs. SEᵢ = (ucl − lcl)/3.92 and SE = √Σ(wᵢ·SEᵢ)², with n weights. This is approximate.
- Registries with no CI (none at present) keep their mean under every reference.
- **Overrides:** `config/registry_means_override.csv` accepts `lcl_pct` / `ucl_pct`.

The report's "Reference" section shows, for each family:
- a per-registry panel plot with the funnels around the mean and both CI bounds
- the outliers / procedures / excess under each reference
- the devices whose outlier status changes between references

The app has a reference selector and a checkbox to draw all three funnels.

## LROI funnels by procedure period
LROI Figure K052B gives the cumulative major revision of all primary UKAs by 2-year procedure period, from 2009-2010 to 2021-2022. It's extracted into `UKA_casemix_long.csv` with a `procedure_period` column. Each period gets its own funnel around its own rate.
- **Report:** the "LROI: funnels by procedure period" section has:
  - a rate-scale plot with the funnels at their true rates, the LROI registry-mean funnel dashed, and the LROI devices
  - a standardised plot, with each period funnel centred on 0
  - a table showing which period funnels each LROI device falls outside
- **App:** Funnel → "Registry mean + LROI funnels by procedure period".
- **Caveat:** K052B reports *major* revision (femur or tibia), while the device tables report any revision.

## Case-mix adjustment (`casemix_adjustment.Rmd`)
A separate report. The case-mix envelope in `funnel_plots.Rmd` only shows how funnel **width** changes with patient risk. Case mix also moves the **expected rate**. For example, with a registry mean of 3% and a highest-risk class at 5%, a model's deviation should be judged against 5%, i.e. 2 points lower. That shift is much larger than the change in width, which is only about 1.28× in this example. The report shows two ways to put the shift on the plot. Each registry has three references: its mean, its lowest-risk age-sex class and its highest-risk age-sex class.

- **Option A, adjusted funnels:** points stay at rate − registry mean (the observed data). Each reference gets its own funnel, centred on (p_ref − p_mean) with p_ref's width. A shaded band spans the shifted funnels of all age-sex classes.
- **Option B, point ranges:** funnels stay centred on 0. Each model gets one deviation per reference: p − p_mean, p − p_low (top) and p − p_high (bottom). They are drawn as three markers (**B1**) or as the main point with a shaded case-mix band (**B2**). Each is read against the funnel of its own reference.
- A and B are the same comparison drawn differently: p > upper(p_ref, n) ⇔ p − p_ref > upper(p_ref, n) − p_ref. Both are classified by `classify_casemix()`, and the report checks that they flag the same models.
- **Status per model** (point colour): high under every reference (case mix cannot explain it) / high vs mean but not vs highest-risk class / high only vs lowest-risk class / within / the low mirror images. The report also gives excess revisions against each reference; the excess against the highest-risk class is a lower bound.
- **Decisions:**
  - Registries: AOANJRR (KP22) and NJR (3.K6) only, as the other registries' extracted reports have no age-sex revision table.
  - NJR's cemented and uncemented/hybrid blocks are pooled per class, n-weighted, because fixation is a model attribute, not case mix.
  - Classes below `min_class_n` (500) are ignored.
- **Caveats:**
  - These are bounds, not an adjustment. The proper method is indirect standardisation with each model's own age-sex mix, which the reports don't publish.
  - Only age and sex are covered.
  - KP22 covers OA only.
- Code: `R/casemix_adjust.R`. Figures: `figures/casemix_adj_<proc>_<t>yr_*.png`. Parameters: `time_yr`, `level` (flagging level, 0.998), `levels`, `limit_method`, `detail_level`, `families`, `min_class_n`, `label_points`, `x_log`.

```r
rmarkdown::render("casemix_adjustment.Rmd")
```

## CRR over follow-up time across registries (`crr_meta_analysis.Rmd`)
A separate report that pools each model's published cumulative revision rates across registries at every follow-up year. It adapts the statistical methods of Combescure et al. 2025 (EFORT Open Rev, doi:10.1530/EOR-2024-0020) from hip cups and stems to knee brands.
- **Inputs:** `UKA_long.csv` (2025 reports). The brand comes from detail level 1 (`config/device_families.csv`), and each registry row is kept as a variant with its own CI.
- **Method:**
  - CRR is put on the cloglog scale, and each SE is backed out of the published 95% CI.
  - Per-year REML random-effects pooling, with variant nested in registry. CIs use a t-distribution with containment df.
  - Heterogeneity: Cochran's Q, multilevel I² and σ² per level.
  - Table 3/4 analogues: multivariable meta-regressions of the 5-year CRR on registry + fixation (Table 3) and on registry + design variant (Table 4), also at 1, 3 and 10 years.
  - A single model across all follow-up years (cloglog quadratic in log t, with random slopes).
  - Sensitivity analyses.
- **Outputs:**
  - Figure 1 analogue (pooled CRR over time, overall and by fixation).
  - Figure 2 analogue (registry variants over time, one line per registry row).
  - A one-line-per-registry version (rows collapsed per registry), with each registry's weight under the nested, independent and collapsed approaches.
  - Every trajectory figure in two versions: all follow-up years, and the core years 1/3/5/10 (files ending `_core`).
  - Table 2 analogue: `tables/crr_ma_<proc>_<models>_table2.csv`. Table 3/4 analogues: `..._table3_4.csv`.
  - Forest plots.
  - Figures are written to `figures/crr_ma_*`.
- **Parameters:** `models` (default Oxford; several brands give the full Figure 1 / Table 2 layout), `by_fixation`, `table34_year` (5), `core_years` (1, 3, 5, 10), `k_reg_min` (3), `test` (`t` or `z`), `attribute_source` (`design`).
- Code: `R/crr_meta.R`. Needs metafor and Matrix.

```r
rmarkdown::render("crr_meta_analysis.Rmd")
rmarkdown::render("crr_meta_analysis.Rmd", params = list(models = c("Oxford", "ZUK / Physica ZUK", "Journey")))
```

## Funnel trajectories (`funnel_trajectories.Rmd`)
This report runs the standard funnel at 1, 3, 5 and 10 years and follows each model through the funnels. Nothing in the funnel changes: each registry is centred on its overall registry mean at that year, and the limits come from `funnel_limits()`.
- **Side-by-side funnels at each year**, in two versions: the ratio scale with one y-axis for all panels, and the difference scale with the y-axis free per panel.
- **Funnel z-score**, z = (p − mean)/√(mean(1 − mean)/n), plotted over years. |z| > 1.96 means outside the 95% funnel and |z| > 3.09 outside the 99.8% funnel.
- **Trajectory labels:** persistent (beyond 95% at ≥ 2 consecutive years), early, late or isolated, each high or low.
- **Censoring sensitivity analysis:** n is replaced by the effective sample size n_eff = p(1 − p)/SE², with SE from the published KM 95% CI. The KM estimate handles censoring but binomial limits with total n do not, so they are too narrow at long follow-up.
- **Parameters:** `models` (Oxford), `detail_level` (1 = brand per registry, n-weighted; 0 = rows), `years`, `levels`, `limit_method`, `x_log`.
- **Code:** `R/funnel_traj.R`. Figures: `figures/ftraj_*`. Tables: `tables/ftraj_*`.

`funnel_over_time.Rmd` is kept as an **exploratory** report (rest-of-registry comparators, interval HRs, cross-registry pooling, joint model, false-alarm simulation).

```r
rmarkdown::render("funnel_trajectories.Rmd")
```

## Funnels over follow-up time, exploratory (`funnel_over_time.Rmd`)
This report keeps the funnels' within-registry comparison and adds follow-up time and cross-registry pooling from Combescure et al. 2025.
- **Comparison:** each model row is compared with the **rest of its own registry** (registry mean minus the model's rows). The comparison is a log cumulative HR on the cloglog scale with a z-score; |z| > 3.09 is outside the 99.8% funnel. SEs are taken from the published KM CIs.
- **Outputs:**
  - one standardised funnel per year (HR vs precision, all registries on one funnel);
  - z-score trajectories over 1/3/5/10 years (brand per registry and variants);
  - interval (conditional) HRs for 0–1, 1–3, 3–5 and 5–10 years, to separate early from late failure;
  - pooled HR per year and per interval (variants nested in registry, with a shared-comparator covariance);
  - a joint model across years, log HR = α + β log t, where β < 0 means the excess fades and β > 0 means it grows;
  - a trajectory classification (persistent / early / late / isolated; high / low);
  - simulated false-alarm rates of each flagging rule;
  - a comparison with the standard 3-year funnel.
- **Parameters:** `models` (Oxford), `core_years` (1, 3, 5, 10), `late_interval` (5–10), `k_reg_min` (3), `max_unsubtracted` (0.10; drops registry-years where the comparator would still contain much of the model), `n_sim`.
- **Code:** `R/funnel_time.R` (uses `R/crr_meta.R`). Figures: `figures/fot_*`. Tables: `tables/fot_*`.

```r
rmarkdown::render("funnel_over_time.Rmd")
```

## Random-effects leave-one-out report (`loo_metafor.Rmd`)
A focused version of `loo_influence.Rmd` containing only the metafor (random-effects) methods, for sharing.
- Contents: model inputs, the random-effects meta-analysis (`rma`, `predict`, `weights`), metafor's own `forest()` and `plot(influence())` figures, `leave1out()` means, influence diagnostics, standardised residuals, Wald vs exact limits, and 90% / 95% funnels around the published mean and the lowest leave-one-out mean.
- A table near the top lists which figures and tables come from metafor, and each figure has a "Source" line.
- **Summary** (at `summary_level` = 90%, `summary_method` = `normal`, i.e. Wald): the models driving the mean (metafor influential flag, and which way each pulls the mean), the recommended least-biased mean (the lowest random-effects leave-one-out mean, with CI and the model left out), and the outliers vs the published mean and vs that mean (counts, names, newly flagged / no longer flagged).

```r
rmarkdown::render("loo_metafor.Rmd")
```

## Leave-one-out and influence analysis (`loo_influence.Rmd`)
A separate, exploratory report for the rare joints (total ankle, total elbow; AOANJRR and NJR, the registries with model-level rates). With 4–13 models per registry and the largest holding 30–50% of procedures, one model can drive the registry mean. It can hide itself (self-masking) or hide other models (masking of others). The app is unchanged.

- **Subtraction leave-one-out:** p₍₋ᵢ₎ = (N·p_pub − nᵢ·pᵢ)/(N − nᵢ), taken from the published total. It is patient-weighted, like the published mean, and keeps the unlisted models. It is approximate, because KM rates are not exactly additive.
- **Meta-analytic leave-one-out (`metafor`):**
  - Random-effects `rma()`: logit scale, REML τ², Knapp-Hartung CI. Each model's SE comes from its **published KM 95% CI**. AOANJRR elbow stems are pooled across diagnoses, with SE from the delta method.
  - `leave1out()` gives the means without each model.
  - `influence()` gives the externally standardised residual (rstudent), DFFITS, Cook's distance, covariance ratio, hat value, τ² deleted, and metafor's "influential" flag (Viechtbauer & Cheung 2010).
  - The random-effects mean is model-weighted (a typical model), not patient-weighted.
- **Reference options compared per model:** published mean / own subtraction leave-one-out / random-effects mean / own random-effects leave-one-out / lowest leave-one-out mean (subtraction) / lowest leave-one-out mean (metafor). The lowest metafor leave-one-out mean is the least upward-biased mean, following the supervisor guidance. It is slightly low when there is no outlier, which suits screening.
- **Limits for small-volume models:**
  - Exact binomial limits (Spiegelhalter 2005) are the default in this report. A calibration table shows that normal limits flag models with the reference rate 2–7× too often below n ≈ 250.
  - The rstudent residual is shown as a third, meta-analytic test. It uses each model's own CI, cannot mask itself and allows for between-model variation (τ²).
- **Supervisor guidance (2026-10):**
  - The default flagging level is **90%** (`level`). A two-sided 90% limit is a one-sided 95% test, in line with non-inferiority testing.
  - The **lowest metafor leave-one-out mean** (`min_loo_re`) is a reference in its own right, as the least upward-biased mean when there are few models.
  - A **"Lower control limits for screening"** section compares 90%, 95% and 99.8% funnels around the published mean and the lowest leave-one-out mean, with power and false-alarm tables.
  - The supervisor's references are listed by DOI (full citations to be added).
- Every metric is explained in the report, which also has a decisions table and caveats (few models, two estimands, KM non-additivity, multiple comparisons).
- Code: `R/loo.R`. Figures: `figures/loo_*.png`. Parameters: `joints`, `time_yr`, `time_yrs`, `level` (0.90), `screen_levels` (0.90, 0.95, 0.998), `screen_refs`, `limit_method`, `ma_test`, `k_min`.

```r
rmarkdown::render("loo_influence.Rmd")     # needs the metafor package
```

## Rare joints: total ankle and total elbow
These are read from the regextract outputs `outputs/ANKLE/` and `outputs/ELBOW/`, the sibling folders of `data_dir`. Run `python -m regextract extract --procedure ANKLE` (and `ELBOW`) first.
- **Report:** the "Rare joints" section shows, for each joint:
  - the registry means
  - an all-models funnel plot
  - a plot for each model family reported by at least 2 registries
  - outlier and cumulative-excess tables under each plot

  All high outliers go to `tables/rare_joints_<t>yr_outliers_998.csv`.
- **Report parameters:**
  - `rare_joints`: default ANKLE, ELBOW
  - `rare_time_yr`: default 3. 1, 3 and 5 years are in every registry.
  - `elbow_by_diagnosis`
- **App:** a **Rare joints** tab next to the knee tab, with:
  - joint, follow-up year (the years reported by at least 2 registries), model, registries
  - control limits (any of 90%, 95% and 99.8%; outliers are judged at the widest one ticked), distribution and reference
  - AOANJRR elbow pooled vs per diagnosis
  - y scale, log x, labels

  Hover shows the model as printed, its diagnosis and its source table. There are CSV downloads.

  **Leave-one-out in the app.** The same leave-one-out analysis as `loo_influence.Rmd` and the manuscript (`R/loo.R`: random-effects meta-analysis, REML with Knapp-Hartung, cloglog scale, registries with at least 3 models), run per joint and year and cached:
  - **Funnel reference = Lowest leave-one-out mean:** each registry's funnel is centred on the lowest of its leave-one-out means (the model whose removal gives it is named in the note and the means table). It is the tightest reference, so it favours flagging. Registries without one (fewer than 3 models with rates, e.g. LROI) keep the published mean.
  - **Leave a model out:** click a model on the funnel (or a dot in the Influence view), or pick it in the "Leave out of the mean" box above the plot. That registry's funnel moves to the mean without that model, the model gets a black square, and the funnel around the published mean stays as a grey line for comparison. One model per registry at a time; click it again, remove its chip or use Reset to put it back.
  - **Leave-one-out mean:** random effects (default, as in the manuscript) or subtraction from the published total, (N p_pub - n_i p_i) / (N - n_i). It applies to both of the above.
  - **Influence** (sub-tab next to Funnel): per registry, the mean recomputed without each model with its 95% CI, against the published mean (black) and the random-effects mean (blue dashed); the lowest is ringed in pink and dot size follows the model's share of the registry's procedures. The table gives n, share, model rate, mean without it, the change from the published mean, Cook's distance and metafor's influence flag (as in Supplementary Table S7).
  - AOANJRR elbow per diagnosis: those stems are judged against their diagnosis's ET9 rate, so leave-one-out applies to the pooled view only.
  - Needs the `metafor` package; without it the app runs as before and says so.

### Decisions (2026-10)
- **Elbow = total elbow replacement only.** AOANJRR's supplement covers only total elbow. NJR radial head, distal humeral hemi and "unconfirmed" rows are dropped (`rare_scope()`).
- **Registry means** are set in `config/rare_joint_means.csv`. Each is an n-weighted pool of the listed rows, with the CI pooled from the published CIs.

  | Joint | AOANJRR | NJR | LROI |
  |---|---|---|---|
  | Ankle | A11, all diagnosis rows pooled (n = 5,379, same population as the A15 models) | 3.A3 All cases | A014B all primary ankles (n = 1,460 from the caption) |
  | Elbow | ET9 total elbow: fracture + OA + RA pooled | 3.E6 total elbow: acute trauma + elective pooled | E017B total elbow (n = 816, approximate, from Table E001) |

  Edit the CSV to change a definition. For example, to use AOANJRR A11 OA only, set `label_regex` to `^Osteoarthritis$`.
- **AOANJRR elbow models** come per diagnosis (ET6 fracture, ET7 OA, ET8 RA).
  - By default they're pooled into one point per humeral stem (n-weighted, no CI), labelled `[3 dx]`.
  - With `elbow_by_diagnosis` (report) or "Per diagnosis" (app), each stem × diagnosis is judged against that diagnosis's ET9 total-elbow rate.
- **Model names** are harmonised in `config/rare_device_families.csv`. Examples: BOX/Box, S.T.A.R/Star → STAR, Coonrad/Morrey and Coonrad Morrey → Coonrad-Morrey, the NJR and AOANJRR Latitude stems → Latitude.
  - Where a registry has several rows in one family, the printed name is kept on the plot, e.g. NJR "Latitude EV Stem" and "Latitude / Latitude EV Stem".
  - AOANJRR lists Salto and Salto Talaris separately. NJR's "Salto" is matched to Salto, but it may include Talaris.
- **LROI** gives no model-level rates for ankle or elbow, so its funnel is drawn without points. LROI values are a manual transcription, checked against the report (NL, 2026-10-05).
- Metrics differ between registries (AOANJRR CPR, NJR 1−KM, LROI KM percentage), so compare models against their own registry's funnel only.
- The AOANJRR ankle periods (A14) and the LROI total-ankle-for-OA periods (A015B) are extracted, but not plotted yet.

## Detail levels
Registry rows are pooled within each registry to a common level of detail. This follows the RSA approach in doi:10.2340/17453674.2026.45293, but rows that lack the detail a level needs are **excluded**, not regrouped:

| Level | One point per | Excluded when |
|---|---|---|
| 1 | brand (`config/device_families.csv`) | brand not matched |
| 2 | brand + fixation (cemented / cementless / hybrid) | fixation unknown |
| 3 | brand + fixation + bearing, femoral material or tibial type | either one unknown |

- Attributes come from the rules in `config/device_attributes.csv`. For each attribute, the first matching rule wins.
- `stated` rules read the report text: EPRD section headings, LROI's cemented/uncemented tables, and component names such as "(cless)", "Oxinium" or "All-Poly".
- `design` rules are general implant knowledge, e.g. ZUK is cemented and fixed bearing, and the Oxford tibia is cemented and mobile bearing. They are marked *verify*.
- A femoral and tibial component with different fixation is *hybrid*. SIRIS "Oxford cemented/hybrid" keeps its own category.
- A pooled point has n = sum of the rows' n and p = the n-weighted mean of the rows' KM estimates. This is approximate, and pooled points have no CI.
- The report's "Detail levels" section contains:
  - a table of points / rows kept per registry and level, for both attribute sources
  - the attributes of every row
  - an all-models plot with outlier tables for each level
  - a family × level outlier summary
- The report also writes `tables/<proc>_<t>yr_level_outliers_<source>_998.csv`.

## Interactive app
```r
install.packages(c("shiny", "bslib", "plotly"))   # once
shiny::runApp("app")                              # from the repo root; opens in your browser
```
It runs only on your computer and uses the same `R/` functions as the report. The controls are:
- **Follow-up year:** only years reported by at least 2 registries are listed. Registry means and device points both come from that year.
- **Model:** a family from `config/device_families.csv`, or all models.
- **Detail level:** as listed / 1 brand / 2 + fixation / 3 + attribute (pick bearing, femoral material or tibial type). You can also switch between stated attributes and stated + design knowledge. A note lists the rows excluded for missing detail, and a table at the bottom gives the rows kept per registry and level.
- **Registries:** tick boxes, all on by default.
- **Control limit:** 99.8%, 95%, or both.
- **Funnel distribution:**
  - Wald normal approximation.
  - Exact binomial (Spiegelhalter).
  - Overdispersion-adjusted. This widens each registry's limits by √φ, where φ is estimated from all of that registry's devices at that year using Spiegelhalter's multiplicative model with 10% winsorising.
- **Fixation shown** (detail levels 2–3): all, cemented only, or cementless only. The funnels stay centred on the registry's overall UKA mean.
- **Label models on the plot:** outliers, all, or none.
- **Funnel reference:** registry mean, lower CI or upper CI of the mean. An option draws all three funnels together (thin dotted/dash-dot lines). The rings and tables use the selected reference.
- **Funnel:** registry mean, or registry mean + a case-mix envelope. The envelope is the range of age-sex class limits, available for AOANJRR and NJR. Outliers are still judged against the registry mean.
- **Layout:** registries overlaid after standardising (difference or ratio), or one panel per registry on the rate scale.

Hover over a point to see:
- n
- the revision rate and its CI at that year
- the registry mean, Δp and the limit
- where the value came from (table and page)

Points get a red ring above the upper limit and a green ring below the lower limit.

Under the plot, the high-outlier table and the cumulative-excess summary update live. There are CSV downloads for both, and the camera icon on the plot saves a PNG.

To read the data from a different folder, set `REGISTRY_DATA_DIR` before starting the app.

### Your data
**+ Add your data** (top right of the plot, both tabs) opens a panel beside the plot where you can put your own numbers next to the registry data. Nothing is saved; the values live in the browser session only.
- **Reference rate (optional):** a mean revision rate (%) and, optionally, the number of procedures it is based on. A black funnel is drawn around it with the selected control limits and distribution. In the panel layout it gets its own panel.
- **Model points (optional):** one model per line as `name, revision rate (%), n`. Rows can be typed, pasted from a spreadsheet (tab- or semicolon-separated; decimal commas are accepted then) or loaded from a CSV (a template can be downloaded). Lines that cannot be read are listed under the box.
- **Judge your models against:** your own reference and/or any of the selected registries. Each model is compared with each ticked reference at its own n, using the same limits as the registry devices. Registry references follow the **Funnel reference** setting.

Your models are drawn as diamonds (◆) filled in the colour of the reference they are judged against (black = your reference), with the usual red/green rings. In the panel layout each one sits in its reference's panel. A table under the plot gives one row per model and one column per reference (status and upper limit at that n); its CSV adds the lower limit, Δp and excess. Your models are never counted in the registry outlier tables or the cumulative excess.

Small n gives wide funnels: with your points at a few hundred procedures and registry devices in the tens of thousands, the **Log x axis** option makes both readable.

### Custom tab
Funnel plots from your own numbers only, with no registry data. It starts with one empty data set; **+ Add data set** adds more (up to 8), and each set after the first has a **Remove this set** link. Each set has its own colour and funnel, and:
- **Name** (shown in the legend, tables and the set's header).
- **Model points:** one model per line as `name, revision rate (%), n`, typed, pasted from a spreadsheet or loaded from a CSV (template available), with the same parsing as the "Your data" panel.
- **Funnel reference:** the n-weighted mean of the set's points (sum of p x n / sum of n), a rate you enter (with an optional n for the legend), or none (points only). A set with an entered rate and no points draws a funnel only, e.g. a registry or national average. With fewer than 5 models the pooled mean is driven by the largest ones, and a note says so.

**Judge each set's points against:**
- **Its own funnel** (default): each model is compared with its own set's reference. Points of a set without a reference are drawn as open diamonds on the rate scale and not judged.
- **Every funnel:** each model is compared with every set's reference. The table becomes one row per model and one column per funnel (status and upper limit, "(own)" marks the model's own set).

**Y scale:**
- **Revision rate:** every funnel at its own reference rate, every model at its observed rate.
- **Standardised, deviation from own reference:** every funnel is drawn relative to its own reference, so all of them are centred on 0 and only their widths differ. Each model is plotted once, as its rate minus its own set's reference. Models of a set without a reference cannot be standardised and are left off (a note says so).
- **Standardised, ratio to own reference:** the same, as rate / reference, with every funnel centred on 1.
- With **Every funnel**, the standardised plot can instead show each model once per funnel, relative to that funnel's reference ("Points on the standardised plot"). By default it shows each model once against its own reference, and the other comparisons stay in the table.

Other settings: control limits (any of 90%, 95% and 99.8%, drawn dotted, dashed and solid; outliers are judged at the widest one ticked), funnel distribution, log x axis and labels. Points are diamonds in their set's colour, with red/green rings above/below the limits. The summary gives, per set and against its own funnel, the outlier counts, procedures with a high-outlier model and the cumulative excess, defined as in the registry tabs. The CSV download has one row per model and funnel, with the method and level used.

## Layout
| File | Contents |
|---|---|
| `R/palette.R` | **Registry colours** (one fixed colour per registry), legend labels (`NJR: 3.21%, n = 206,689`), line types per level, theme |
| `R/data.R` | Loads the CSVs, computes registry means, device points and case-mix strata |
| `R/limits.R` | `funnel_limits()` (normal / exact / overdispersed), `estimate_phi()`, `build_limit_curves()`, `classify_points()`, outlier reports |
| `app/app.R` | Interactive Shiny app |
| `R/plots.R` | `funnel_plot()` core function and the named variants |
| `R/rare.R` | Rare joints: `load_rare()`, `rare_means()`, `rare_points()` |
| `app/rare_tab.R` | Rare joints tab of the app |
| `app/user_data.R` | "Your data" panel (user reference rate and model points) used by both tabs |
| `app/custom_tab.R` | Custom tab: funnel plot from the user's own data only |
| `config/rare_joint_means.csv` | Which rows give each registry's rare-joint mean |
| `config/rare_device_families.csv` | Harmonised rare-joint model names |
| `R/casemix_adjust.R` | Case-mix adjustment report: class rates, shifted / centred reference funnels, `classify_casemix()`, option A and B plots |
| `casemix_adjustment.Rmd` | Case-mix adjustment report (adjusted funnels and model point ranges) |
| `R/crr_meta.R` | CRR meta-analysis across registries: cloglog inputs, per-year pooling, meta-regression, joint model across years |
| `crr_meta_analysis.Rmd` | Pooled CRR over follow-up time across registries (Combescure et al. 2025 approach) |
| `R/funnel_traj.R` | Funnel trajectories: points per year (rows or brand), funnel z-score, n_eff, trajectory labels |
| `funnel_trajectories.Rmd` | Standard funnel at several follow-up years, z-score trajectories, n_eff sensitivity analysis |
| `R/funnel_time.R` | Funnels over time: rest-of-registry comparators, log-HR contrasts, interval HRs, pooling, joint model, trajectory classification, false-alarm simulation |
| `funnel_over_time.Rmd` | Exploratory: outlier trajectories over follow-up, within registry and pooled across registries |
| `R/loo.R` | Leave-one-out / influence report: model inputs (SE from KM CI), `loo_fit()` (metafor), `loo_status()`, calibration, plots |
| `loo_metafor.Rmd` | Random-effects (metafor) leave-one-out report with summary of drivers, least-biased mean and outliers |
| `loo_influence.Rmd` | Leave-one-out and influence analysis of rare-joint registry means |
| `R/levels.R` | Detail levels: `derive_attributes()`, `pool_to_level()`, `level_counts()` |
| `config/device_families.csv` | Regex that groups each registry's device labels into families/brands (Oxford, ZUK, ...) |
| `config/device_attributes.csv` | Rules for fixation / bearing / material / tibia, each marked stated or design |
| `config/registry_means_override.csv` | Optional hand-set registry means, e.g. from a table not yet extracted |

### Several report years
`load_extraction(data_dir, procedure, report_years = "latest")` keeps the newest report per registry. That guards
against a regextract run with several editions (`--year 2024 2025`) mixing them. Use `report_years = "all"` or
specific years to keep them.

### Registry means
For each registry, the funnel centre at `time_yr` is taken from the first of these sources that is available:
1. A manual override in `config/registry_means_override.csv`.
2. The registry's own summary row: NJR "All unicompartmental", SIRIS "CH average", EPRD "Unicondylar knee arthroplasties" (2025 Table 57), AOANJRR "Unicompartmental" (2025 Table KP18: all primary UKA for OA, modern prostheses, n = 50,481), LROI cemented + uncemented pooled by n.
3. Pooled case-mix subtotals (e.g. AOANJRR KP22 male + female). Used only if a registry has no summary row. AOANJRR used this before KP18 was extracted; KP22 covers the same population (n = 50,481) and gave the same rates.

AOANJRR's model rows (Table KP17) cover modern prostheses for all diagnoses (TOTAL n = 50,910). KP18 is the OA subset, 99.2% of those procedures, so the mean and the models describe practically the same population.
4. **Fallback**: the n-weighted mean of the listed devices (no CI). No registry uses this at present; EPRD did with its 2024 report, which had no UKA total row.

`mean_source` in the report says which source each registry used.

### Axis ranges
Axis windows are worked out for each plot by `dynamic_ranges()` in `R/plots.R`:
- Every device point is inside the window.
- The y axis is symmetric about the reference (0 or 1), so both the upper and the lower side of each funnel show.
- Funnels are drawn down to the smallest device volume on the plot (`n_floor`).
- Plots with only funnels (the registry panels and the age-sex class panels) use a log x axis so the funnel shape is readable.
- To fix a window by hand, pass `xlim` / `ylim`.

### Case mix
Case-mix limits are drawn per **age-sex class**, e.g. `<55 Female` or `>=75 Male`. That gives 8 classes per registry or design block, in paired colours: one hue per age group, dark for female and light for male. The envelope here shows width only; see `casemix_adjustment.Rmd` for funnels shifted to each class's rate.

### Outlier reporting
Each device-family plot at 99.8% has two tables underneath:
1. **Outliers:** the devices above their own registry's upper 99.8% limit at their own n.
2. **Cumulative excess revision:** summed per registry and overall. It is Σ n × Δp, where Δp = device rate − registry mean. Only outlier devices count. These are screening numbers for generating hypotheses, not precise estimates.

The functions are `outlier_report()` and `outlier_summary()` in `R/limits.R`. All outlier rows are also written to `tables/<proc>_<t>yr_family_outliers_998.csv`.
To add tables under any other plot, give its row in the plot catalogue a `tables` function that returns `list(outliers, summary)`.

## Adding a plot
Every figure is one row in the `plots` table in the "Plot catalogue" chunk of `funnel_plots.Rmd`. Each row has a slug, a section, a size and a function that returns a ggplot.
Most variations are `funnel_plot()` called with different options:

```r
"njr_lroi_ratio_log", "Standardised overlay", 10, 6, function()
  funnel_plot(filter(reg_curves, registry %in% c("NJR", "LROI")),
              points = filter(points, registry %in% c("NJR", "LROI")),
              means = means, levels = 0.998, y = "ratio", x_log = TRUE, sides = "upper",
              title = "NJR vs LROI, 99.8% upper limits")
```

Options for `funnel_plot()`:

| Option | What it does |
|---|---|
| `levels` | Which control-limit levels to draw |
| `sides` | `"both"` or `"upper"` |
| `y` | `"difference"`, `"ratio"` or `"absolute"` |
| `facet` | Split the plot into panels |
| `x_log` | Log-scale x axis |
| `ribbon` | Shade the band between the limits |
| `label_filter` | Which points get labels |
| `xlim`, `ylim` | Axis windows |
| `colour_by` + `colour_scale` | Colour by something other than registry |

To add a registry, give it a colour in `REGISTRY_COLOURS` (`R/palette.R`). Plots stop with an error if a registry has no colour, so no registry is ever silently given a different colour.

## Caveats
- The limits treat a KM cumulative-revision estimate as a binomial proportion and ignore censoring. They are approximate, and more so when few patients remain at risk. NJR estimates with 250 or fewer at risk are drawn as open points.
- Metrics differ between registries: CPR, 1−KM or failure rate, with different populations and revision definitions. Standardising to each registry's own mean removes level differences between registries, not definitional ones.
- LROI values come from a manual transcription of image tables, checked against the report (NL, 2026-10-05). If a registry has rows marked `unverified` (e.g. a new report year), the captions and app notes name it.
- Fixation-segregated funnels (v1 "by fixation" section) are not included yet.
