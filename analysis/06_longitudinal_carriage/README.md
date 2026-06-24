# 06 — Longitudinal E. coli ASV Carriage

Longitudinal profiling of *E. coli* ASV carriage across repeated faecal samples
from the same individual. Compares detection patterns and persistence between the
BEST and STANDARD DADA2 pipelines using only ASVs with perfect BLAST support.

## Contents

| File | Description |
|------|-------------|
| `06_Longitudinal_Carriage.Rmd` | Main analysis notebook — knit to reproduce all results |
| `06_Longitudinal_Carriage_functions.R` | Companion functions sourced by the Rmd |

## What this analysis does

1. **Loads case phyloseq objects** for the BEST and STANDARD pipelines
2. **Attaches longitudinal metadata** (Subject ID, days since baseline, sample type) from the reference phyloseq
3. **Filters to individual samples** — removes mock, negative control, and blank samples
4. **Parses BLAST results** from step 04 and retains only ASVs with perfect *E. coli* hits (100% identity, 100% query coverage)
5. **Plots carriage** for each Subject × Pipeline combination: filled circles sized by relative abundance, with dashed lines marking all sampled timepoints
6. **Saves** all figures as PNG and PDF, and exports a master detection dataframe and summary table as CSVs

## BLAST filter criteria

Only ASVs satisfying all three of the following are retained:
- `pident == 100` (100% sequence identity)
- `query_cover == 100` (100% query coverage)
- `sscinames` contains "Escherichia coli"

## Inputs

```
data/
├── cases/
│   ├── case_049__F-Strict__Err-Binned__Pool-FALSE__O1e-200__B32.rds   ← BEST
│   └── case_031__F-Loose__Err-PacBio__Pool-FALSE__O1e-40__B32.rds     ← STANDARD
└── wgs/
    └── phyloseq_NCBInt_process_all_sub_OMEGA_A-200_poolFALSE_minQ30EE05.RDS
        ← longitudinal metadata source (saved with save(), not saveRDS())

results/
└── 04_best_vs_standard/
    └── blast/
        └── blast_df_parsed.csv   ← BLAST results from step 04 (required)
```

> ⚠️ **Step 04 must be completed first** — this analysis depends on `blast_df_parsed.csv`.

## Outputs

```
results/
└── 06_longitudinal_carriage/
    ├── longitudinal_carriage_<subject>_BEST.png/.pdf     (one pair per subject)
    ├── longitudinal_carriage_<subject>_STANDARD.png/.pdf (one pair per subject)
    ├── 06_carriage_df.csv      — master long-format detection dataframe
    └── 06_carriage_summary.csv — summary stats per Subject × Pipeline
```

## How to run

1. Ensure step 04 has been completed and `blast_df_parsed.csv` exists
2. Set your R working directory to the **project root**
3. Ensure the conda environment is active:
   ```bash
   source start_env.sh
   ```
4. Open `06_Longitudinal_Carriage.Rmd` in RStudio and click **Knit**, or run:
   ```r
   rmarkdown::render("analysis/06_longitudinal_carriage/06_Longitudinal_Carriage.Rmd")
   ```

## Dependencies

Managed via `env_dada2_final.yml`. Key packages:

| Package | Source |
|---------|--------|
| `phyloseq` | Bioconductor |
| `tidyverse` | CRAN |
| `ggplot2` | CRAN |

## Notes

- The longitudinal metadata phyloseq was saved with `save()` rather than `saveRDS()`. `load_florentin_phyloseq()` handles this automatically by loading into a temporary environment and extracting the first phyloseq object found.
- All prevalence thresholds are set to 0 by default (keep all *E. coli* ASVs). Tighten `prev_high_cutoff`, `prev_med_cutoff`, and associated `min_ab` parameters in the `run_longitudinal_carriage()` call to apply adaptive filtering.
- Point size encodes relative abundance in five bins (< 0.5%, 0.5–1%, 1–2%, 2–5%, > 5%). Bin boundaries and sizes are fully configurable.
