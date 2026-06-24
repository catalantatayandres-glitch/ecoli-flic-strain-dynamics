# 05 — Mock Distributions (Check Mock QC)

Quality control of mock community positive controls. Assesses how well the
experimental mock samples match the theoretical composition for both the BEST
and STANDARD pipelines using Bray-Curtis similarity and stacked composition plots.

## Contents

| File | Description |
|------|-------------|
| `05_Mock_Distributions.Rmd` | Main analysis notebook — knit to reproduce all results |
| `05_Mock_Distributions_Functions.R` | Companion functions sourced by the Rmd |

## What this analysis does

1. **Loads the theoretical mock reference** (`ps_mock_trimmed_4.rds`, 4 known *fliC* sequences: ESBL15, ESBL25, KPC_1, Oxa48_1)
2. **Loads the BEST and STANDARD pipeline cases** (one RDS file each)
3. **Extracts all mock samples** — Mock Mix (`S_Mock_mix*`) and Singular Mocks (`S_KPC1`, `S_Oxa`, `S_E15`, `S_E25`) — combined into a single phyloseq
4. **Assigns taxonomy** to experimental ASVs by pairwise alignment against the theoretical reference sequences (not via DECIPHER::IdTaxa)
5. **Computes Bray-Curtis similarity** between each experimental mock sample and the theoretical composition
6. **Plots composition** as stacked bars alongside BC similarity values
7. **Saves figures** for both pipelines as PNG and PDF


## Inputs

```
data/
├── mock/
│   └── ps_mock_trimmed_4.rds          — theoretical mock (4 fliC sequences)
└── cases/
    ├── case_049__F-Strict__Err-Binned__Pool-FALSE__O1e-200__B32.rds  ← BEST
    └── case_031__F-Loose__Err-PacBio__Pool-FALSE__O1e-40__B32.rds    ← STANDARD
```

## Outputs

```
results/
└── 05_mock_qc/
    ├── 05_Mock_CheckZymo_case049_BEST.png
    ├── 05_Mock_CheckZymo_case049_BEST.pdf
    ├── 05_Mock_CheckZymo_case031_STANDARD.png
    └── 05_Mock_CheckZymo_case031_STANDARD.pdf
```

## How to run

1. Set your R working directory to the **project root**
2. Ensure the conda environment is active:
   ```bash
   source start_env.sh
   ```
3. Open `05_Mock_Distributions.Rmd` in RStudio and click **Knit**, or run:
   ```r
   rmarkdown::render("analysis/05_mock_qc/05_Mock_Distributions.Rmd")
   ```

## Dependencies

Managed via `env_dada2_final.yml`. Key packages:

| Package | Source |
|---------|--------|
| `phyloseq` | Bioconductor |
| `Biostrings` | Bioconductor |
| `pwalign` | Bioconductor |
| `microbiome` | Bioconductor |
| `chkMocks` | GitHub (microsud/chkMocks) |
| `ggplot2` / `patchwork` | CRAN |
| `RColorBrewer` | CRAN |

## Notes

- `pwalign` is required as `pairwiseAlignment` was moved out of `Biostrings` in version ≥ 2.77.1.
- ASVs that do not match any theoretical reference at ≥ `MIN_SIM`% identity receive a `nonMock::` prefix label and are retained as distinct taxa in all plots.
- The shared colour palette (`build_shared_palette()`) ensures the same taxon always has the same colour across BEST and STANDARD plots.
- This step does not depend on steps 02–04 and can be run independently.
