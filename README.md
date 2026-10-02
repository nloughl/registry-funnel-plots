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
Packages: dplyr, tidyr, purrr, readr, ggplot2, ggrepel, scales, rmarkdown, knitr.
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
- **NJR and SIRIS:** published.
- **LROI (cemented + uncemented) and AOANJRR (male + female KP22):** pooled from the parts' published CIs. SEᵢ = (ucl − lcl)/3.92 and SE = √Σ(wᵢ·SEᵢ)², with n weights. This is approximate.
- **EPRD:** no CI (no published total), so it keeps its mean under every reference.
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
  - control limit, distribution and reference
  - AOANJRR elbow pooled vs per diagnosis
  - y scale, log x, labels

  Hover shows the model as printed, its diagnosis and its source table. There are CSV downloads.

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
- **LROI** gives no model-level rates for ankle or elbow, so its funnel is drawn without points. LROI values are a manual transcription that has not yet been verified.
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
| `config/rare_joint_means.csv` | Which rows give each registry's rare-joint mean |
| `config/rare_device_families.csv` | Harmonised rare-joint model names |
| `R/levels.R` | Detail levels: `derive_attributes()`, `pool_to_level()`, `level_counts()` |
| `config/device_families.csv` | Regex that groups each registry's device labels into families/brands (Oxford, ZUK, ...) |
| `config/device_attributes.csv` | Rules for fixation / bearing / material / tibia, each marked stated or design |
| `config/registry_means_override.csv` | Optional hand-set registry means, e.g. from a table not yet extracted |

### Registry means
For each registry, the funnel centre at `time_yr` is taken from the first of these sources that is available:
1. A manual override in `config/registry_means_override.csv`.
2. The registry's own summary row: NJR "All unicompartmental", SIRIS "CH average", LROI cemented + uncemented pooled by n.
3. Pooled case-mix subtotals: AOANJRR KP22 male + female.
4. **Fallback**: the n-weighted mean of the listed devices. EPRD currently uses this, because the 2024 table has no UKA total row.

`mean_source` in the report says which source each registry used.

### Axis ranges
Axis windows are worked out for each plot by `dynamic_ranges()` in `R/plots.R`:
- Every device point is inside the window.
- The y axis is symmetric about the reference (0 or 1), so both the upper and the lower side of each funnel show.
- Funnels are drawn down to the smallest device volume on the plot (`n_floor`).
- Plots with only funnels (the registry panels and the age-sex class panels) use a log x axis so the funnel shape is readable.
- To fix a window by hand, pass `xlim` / `ylim`.

### Case mix
Case-mix limits are drawn per **age-sex class**, e.g. `<55 Female` or `>=75 Male`. That gives 8 classes per registry or design block, in paired colours: one hue per age group, dark for female and light for male.

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
- LROI values come from a manual transcription that has not yet been verified. The figure captions say so until it is.
- Fixation-segregated funnels (v1 "by fixation" section) are not included yet.
