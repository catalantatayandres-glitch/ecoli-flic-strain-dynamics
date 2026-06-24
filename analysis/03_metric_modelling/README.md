# 03 — Metric Modelling & Parameter Selection

Second-stage benchmarking of the DADA2 parameter grid. Computes quality
metrics on mock community samples, models the effect of each parameter,
and ranks all 96 parameter combinations to identify the optimal pipeline
configuration.

## Contents

| File | Description |
|------|-------------|
| `03_case_metric_modeling.Rmd` | Main analysis notebook — knit to reproduce all results |
| `03_case_metric_modeling_functions.R` | Companion functions sourced by the Rmd |

## What this analysis does

**Part 1 — Parameter effects on benchmark metrics**

1. Loads all 96 DADA2 case outputs and immediately rarefies each to 1 500 reads/sample
2. Computes 6 quality metrics on mockmix samples only (chimeric reads, total ASVs, mock-matching accuracy, read recovery)
3. Fits regression models per metric (quasibinomial GLM for proportions, Poisson GLM for counts, OLS for reads)
4. Computes observed marginal means per parameter level

**Part 2 — Composite score ranking**

5. Ranks all 96 parameter combinations using a composite score:
   `score = −z(noise) + z(recovery)`
   where noise = % ASVs not matching mock and recovery = % reads from mock-matching ASVs
6. Identifies the top 3 distinct score levels (tie-aware)
7. Produces scatter plots, bar charts, and a full ranking figure

## Inputs

```
data/
├── mock/
│   └── mock_trimmed_4.fasta     — 4 theoretical fliC reference sequences
└── cases/
    ├── case_001__F-Strict__Err-Binned__Pool-FALSE__O1e-200__B32.rds
    └── ...   (96 RDS files total)
```

Each RDS file must contain a list with at minimum:
- `$st` — seqtab before chimera removal (samples × ASVs matrix)
- `$st_nochim` — seqtab after chimera removal
- `$bim` — named logical vector: `TRUE` if ASV is chimeric
- `$track` — data.frame with read tracking columns

Mock samples are identified by `"Mock_mix"` in their rownames.

## Outputs

```
results/
└── 03_metric_modelling/
    ├── tables/
    │   ├── 01_case_metrics.csv
    │   ├── 02_metric_summary_stats.csv
    │   ├── 03_all_model_coefficients.csv
    │   ├── 04_marginal_means.csv
    │   ├── 05_all_cases_ranked.csv
    │   └── 06_top_parameter_combinations.csv
    ├── figures/
    │   ├── 01_marginal_means_<metric>.png   (4 figures)
    │   ├── 02_boxplot_<metric>_by_<param>.png  (20 figures)
    │   ├── 03_coef_plot_<metric>.png        (4 figures)
    │   ├── 04_heatmap_overview.png
    │   ├── 05_ranking_scatter.png
    │   ├── 06_top3_composite_scores.png
    │   ├── 07_top3_dual_metrics.png
    │   └── 08_full_ranking_all_cases.png
    └── model_summaries/
        ├── pct_chimeric_reads_model_summary.txt
        ├── mean_reads_per_sample_model_summary.txt
        ├── pct_asvs_matching_mock_model_summary.txt
        └── total_asvs_model_summary.txt
```

## How to run

1. Set your R working directory to the **project root**
2. Ensure the conda environment is active:
   ```bash
   source start_env.sh
   ```
3. Open `03_case_metric_modeling.Rmd` in RStudio and click **Knit**, or run:
   ```r
   rmarkdown::render("analysis/03_metric_modelling/03_case_metric_modeling.Rmd")
   ```

## Dependencies

Managed via `env_dada2_final.yml`. Key packages:

| Package | Source |
|---------|--------|
| `ggplot2` | CRAN |
| `dplyr` | CRAN |
| `tidyr` | CRAN |
| `vegan` | CRAN |

## Notes

- IUPAC ambiguity codes in reference sequences are handled correctly (e.g. the `V` in ESBL15 is matched as A, C, or G).
- Rarefaction uses `set.seed(42)` for full reproducibility.
- The composite score weighting (noise vs recovery = 1:1) can be adjusted via `weight_noise` and `weight_recovery` in `rank_parameter_combinations()`.
- This step builds on the case outputs from `pipeline/dada2/` and is independent of step 02.
