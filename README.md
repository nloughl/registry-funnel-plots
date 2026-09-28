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

## Layout
| File | Contents |
|---|---|
| `R/palette.R` | **Registry colours** (one fixed colour per registry), legend labels (`NJR: 3.21%, n = 206,689`), line types per level, theme |
| `R/data.R` | Loads the CSVs, computes registry means, device points and case-mix strata |
| `R/limits.R` | `funnel_limits()` (normal / exact), `build_limit_curves()`, `classify_points()` |
| `R/plots.R` | `funnel_plot()` core function and the named variants |
| `config/device_families.csv` | Regex that groups each registry's device labels into families (Oxford, ZUK, ...) |
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
