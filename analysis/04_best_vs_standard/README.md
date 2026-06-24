# 04 — Best vs Standard Pipeline Comparison

Direct comparison of the BEST pipeline (identified in step 03) against the
STANDARD DADA2 pipeline. Evaluates sequence-level divergence via edit distance
analysis, validates ASVs biologically using BLAST and WGS data, and visualises
the results through an integrated dendrogram and biological group composition plots.

## Contents

| File | Description |
|------|-------------|
| `04_BEST_vs_STANDARD.Rmd` | Main analysis notebook — knit to reproduce all results |
| `04_functions.R` | Companion functions sourced by the Rmd |

## What this analysis does

1. **Loads case outputs** for the Best and Standard pipelines by filename pattern matching
2. **Extracts ASV sequences** from each pipeline, prefixing names to guarantee uniqueness
3. **Computes edit distance matrix** (DECIPHER alignment + pwalign Hamming distance) with disk caching
4. **Figure 1** — Edit distance distributions: within-Best, within-Standard, and cross-pipeline
5. **Writes FASTA files** split into parts for submission to NCBI BLAST (manual step)
6. **Parses BLAST JSON results** and computes genus-level LCA per ASV
7. **Assigns trust tiers** (5-tier system) based on BLAST identity/coverage and WGS confirmation
8. **Figure 2** — Integrated dendrogram: clustering + trust tiers + abundance + alignment
9. **Figure 3** — Bio Group composition (E. coli / Ambiguous / Non E. coli) by ASV origin

## Pipeline parameters

| Parameter | BEST | STANDARD |
|-----------|------|----------|
| Filtering | Strict (maxEE = 0.5) | Loose (maxEE = 2) |
| Error model | Binned | PacBio |
| Pooling | FALSE | FALSE |
| Omega | 1e-200 | 1e-40 |
| Band size | 32 | 32 |

## Inputs

```
data/
├── cases/
│   ├── case_NNN__F-Strict__Err-Binned__Pool-FALSE__O1e-200__B32.rds   ← Best
│   ├── case_NNN__F-Loose__Err-PacBio__Pool-FALSE__O1e-40__B32.rds     ← Standard
│   └── ...
├── wgs/
│   └── ps_wgs_EcfliC.RDS    — WGS phyloseq reference
└── blast/
    └── JSON/
        ├── blast_part_01.json   ← downloaded manually from NCBI BLAST
        ├── blast_part_02.json
        └── ...
```

> ⚠️ **BLAST is a manual step.** Run section 3b first to generate split FASTA
> files, submit each to [NCBI BLAST](https://blast.ncbi.nlm.nih.gov), download
> results as "Single-file JSON2", and save them into `data/blast/JSON/` before
> running sections 6–10.

## Outputs

```
results/
└── 04_best_vs_standard/
    ├── 04_all_asvs.fasta
    ├── cached_edit_distance_matrix.rds          ← reused on re-runs
    ├── Figure1_Edit_Distance_Distribution.pdf/.png
    ├── Figure2_Integrated_Dendrogram.pdf/.png
    └── Figure3_BioGroup_by_Origin.pdf/.png

data/blast/
    ├── blast_part_01.fasta … blast_part_07.fasta   ← BLAST input
    ├── JSON/blast_part_01.json … (manual input)
    └── blast_df_parsed.csv                          ← BLAST cache
```

## How to run

1. Set your R working directory to the **project root**
2. Ensure the conda environment is active:
   ```bash
   source start_env.sh
   ```
3. Open `04_BEST_vs_STANDARD.Rmd` in RStudio and knit sections 0–5 to generate FASTA files
4. Submit FASTA parts to NCBI BLAST and download JSON results into `data/blast/JSON/`
5. Knit sections 6–11 to complete the analysis

## Dependencies

Managed via `env_dada2_final.yml`. Key packages:

| Package | Source |
|---------|--------|
| `phyloseq` | Bioconductor |
| `DECIPHER` | Bioconductor |
| `Biostrings` | Bioconductor |
| `pwalign` | Bioconductor |
| `ggtree` / `ggtreeExtra` | Bioconductor |
| `ggplot2` / `patchwork` | CRAN |
| `jsonlite` | CRAN |
| `ape` / `phangorn` | CRAN |

## Notes

- The edit distance matrix is cached to `results/04_best_vs_standard/cached_edit_distance_matrix.rds`. Delete this file to force recomputation.
- The parsed BLAST table is cached to `data/blast/blast_df_parsed.csv`. Delete to re-parse JSON files.
- Trust tier assignment requires both BLAST results and WGS data. If BLAST is unavailable, all ASVs are assigned Tier 5 as a placeholder.
- Figure 2 is large (28 × 20 inches) and may take several minutes to render.
