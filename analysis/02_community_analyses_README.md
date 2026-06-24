# 02 — Community Analyses

Community-level analysis of all 96 DADA2 pipeline parameter combinations,
assessing how parameter choice shapes the detected ASV community through
ASV sharing analysis, distance matrices, PERMANOVA, and PCoA ordination.

## Contents

| File | Description |
|------|-------------|
| `02_community_analyses.Rmd` | Main analysis notebook — knit to reproduce all results |
| `02_community_analyses_functions.R` | Companion functions sourced by the Rmd |

## What this analysis does

1. **Loads all 96 DADA2 output cases** (2 filter × 3 pool × 4 omega × 2 band × 2 error model)
2. **Removes pure mock controls** (`S_KPC`, `S_OXA`, `S_E15`, `S_E25`) — only MOCKMIX community standards are kept as positive controls
3. **Rarefies** all samples to a fixed depth of 1 500 reads
4. **ASV sharing analysis** — classifies ASVs as shared across or unique to each parameter level
5. **Distance matrices** — computes Jaccard (presence/absence) and Aitchison (CLR-Euclidean) distances with disk caching
6. **PERMANOVA** — tests each of the 6 parameters independently using focused subsamples (999 permutations)
7. **PCoA ordination** — 6-panel plots for both distance metrics
8. **Dendrogram heatmaps** — hierarchical clustering with annotation strips for all parameters

## Inputs

```
data/
└── cases/
    ├── case_001__F-Strict__Err-Binned__Pool-FALSE__O1e-200__B32.rds
    ├── case_002__F-Strict__Err-Binned__Pool-FALSE__O1e-200__B64.rds
    └── ...   (96 RDS files total)
```

Each RDS file contains a `phyloseq` object output by the DADA2 pipeline (see `pipeline/dada2/`).

## Outputs

```
results/
└── 02_community_analyses/
    ├── 01_asv_sharing/
    │   ├── ASV_sharing_all_params.png
    │   ├── ASV_sharing_Filter_Mode.csv
    │   ├── ASV_sharing_Pool.csv
    │   ├── ASV_sharing_Omega.csv
    │   ├── ASV_sharing_Band_Size.csv
    │   └── ASV_sharing_Error_Model.csv
    ├── 02_pcoa/
    │   ├── PCoA_Jaccard_6panels.png
    │   └── PCoA_Aitchison_6panels.png
    ├── 03_dendrogram/
    │   ├── Dendrogram_Jaccard.png
    │   └── Dendrogram_Aitchison.png
    ├── 04_permanova/
    │   ├── PERMANOVA_results.png
    │   └── PERMANOVA_all.csv
    ├── 05_combined/
    │   └── Community_summary.png
    └── distances/
        ├── dist_Jaccard.rds   (cached)
        └── dist_Aitchison.rds (cached)
```

## How to run

1. Set your R working directory to the **project root**
2. Ensure the conda environment is active:
   ```bash
   source start_env.sh
   ```
3. Open `02_community_analyses.Rmd` in RStudio and click **Knit**, or run from the terminal:
   ```r
   rmarkdown::render("analysis/02_community_analyses/02_community_analyses.Rmd")
   ```

## Dependencies

Managed via `env_dada2_final.yml`. Key packages:

| Package | Source |
|---------|--------|
| `phyloseq` | Bioconductor |
| `ComplexHeatmap` | Bioconductor |
| `vegan` | CRAN |
| `ape` | CRAN |
| `ggplot2` | CRAN |
| `patchwork` | CRAN |

## Notes

- Distance matrices are cached to `results/02_community_analyses/distances/` on first run. Set `FORCE_RECOMPUTE <- TRUE` in the `compute-distances` chunk to recompute from scratch.
- Windows users: parallel loading is automatically disabled (`N_CORES = 1`).
- This step must be run before `03_metric_modelling/`, which uses the same 96 case outputs.
