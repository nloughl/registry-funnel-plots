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
| `attribute_source` | stated | `stated`: only details printed in the report. `design`: also fills gaps from implant knowledge |

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
