# High-Resolution *E. coli* Strain Profiling from Faecal DNA: Benchmarking Long-Read *fliC* Amplicon Sequencing

> **Benchmarking DADA2 parameters for long-read *fliC* amplicon sequencing and applying the validated pipeline to longitudinal *E. coli* carriage profiling in a childcare cohort.**

**Author:** Andrés Catalán Tatay &nbsp;·&nbsp; **Supervisor:** Florentin Constancias &nbsp;·&nbsp; **Director:** Sonja Lehtinen  
**Evolutionary Epidemiology Group · Computational Biology · University of Lausanne · 2025–2026**

---

## Why this project?

*E. coli* is a major gut coloniser, AMR reservoir, and early-life commensal — but strain-level profiling at scale is hard:

| Method | Problem |
|--------|---------|
| WGS + culture | Too resource-intensive for large cohorts |
| 16S rRNA amplicon | Cannot distinguish *E. coli* strains |
| Shotgun metagenomics | *E. coli* is < 1% of gut microbiome; insufficient depth |

**Solution:** culture-independent amplicon sequencing of the *fliC* (flagellin) gene with PacBio long reads. *fliC* has conserved primer sites flanking a hypervariable central region (~700–1700 bp) that accumulates strain-specific mutations. PacBio CCS is required because the amplicon is too long for short-read platforms.

Raw reads are processed by **DADA2**, whose tuneable parameters had never been systematically evaluated for long-read *fliC* data. This is the central aim of this work.

---

## What was done

1. **96 DADA2 configurations** benchmarked (2 filter × 2 error model × 3 pool × 4 omega × 2 band size)
2. **Mock community** of 4 known *E. coli fliC* sequences used as ground truth
3. **Optimal pipeline** identified via composite score (noise vs. recovery trade-off)
4. **Best vs. Standard** pipeline compared at sequence, taxonomic, and biological level
5. **Longitudinal carriage** of *E. coli* ASVs profiled in 3 individuals from a childcare cohort

---

## Main results

- **Pooling strategy** and **filter mode** drive pipeline performance; band size has no effect
- **Optimal config:** `Strict filter · Binned error model · FALSE pooling`
- Standard pipeline recovers marginally more true *E. coli* but generates substantially more noise → fundamental sensitivity–specificity trade-off
- Over 90% of ASVs unique to the standard pipeline lack reliable *E. coli* BLAST support
- Both pipelines recover the same dominant persistent strains longitudinally; differences arise in low-abundance, low-prevalence detections
- Subsequent experiments suggest the **Kinnex library prep** introduced contamination, contributing to the excess noise

---

## Repository structure

```
ecoli-flic-strain-dynamics/
│
├── dependencies/                  ← environment files
│   ├── env_dada2_final.yml        ← conda environment specification
│   ├── start_env.sh               ← environment activation script
│   └── .gitignore
│
├── pipeline/                      ← Part 1 · HPC cluster (UNIL)
│   ├── sample_selection/          ← subsample to 10k reads/sample; select cohort samples
│   ├── primer_trimming/           ← cutadapt primer removal
│   └── dada2/                     ← run DADA2 across all 96 parameter combinations
│       ├── dada2_script.R
│       └── run_dada2_script.sh
│
└── analysis/                      ← Part 2 · local
    ├── 02_community_analyses/     ← ASV sharing · PERMANOVA · PCoA · dendrograms
    ├── 03_metric_modelling/       ← quality metrics · GLMs · composite score ranking
    ├── 04_best_vs_standard/       ← edit distances · BLAST taxonomy · bio group plots
    ├── 05_mock_qc/                ← Bray-Curtis similarity vs theoretical mock
    └── 06_longitudinal_carriage/  ← E. coli carriage bubble plots per subject
```

Each `analysis/` step folder contains:
```
0X_step_name/
├── Scripts/        ← .Rmd notebook + companion .R functions file
└── Figures/        ← key output figures
```

---

## Data

> Raw FASTQ files are not in this repo. Populate `data/` locally before running anything.

```
data/
├── cases/
│   └── case_NNN__F-<filter>__Err-<err>__Pool-<pool>__O<omega>__B<band>.rds
│       96 phyloseq RDS files — one per DADA2 parameter combination
│       output of pipeline/dada2/
│
├── mock/
│   ├── mock_trimmed_4.fasta     4 theoretical E. coli fliC sequences (ground truth)
│   └── ps_mock_trimmed_4.rds    same sequences as a phyloseq object
│
└── wgs/
    ├── ps_wgs_EcfliC.RDS
    │   fliC sequences from WGS-assembled E. coli isolates from the same cohort
    │   used in step 04 to confirm ASVs against independent WGS evidence
    │
    └── phyloseq_NCBInt_process_all_sub_OMEGA_A-200_poolFALSE_minQ30EE05.RDS
        longitudinal metadata phyloseq (Subject ID, days since baseline, sample type)
        used in step 06 for carriage plots
```

### Mock community strains

| Strain ID | Resistance gene |
|-----------|----------------|
| Oxa48_1 | OXA-48 carbapenemase |
| ESBL25 | Extended-spectrum β-lactamase |
| ESBL15 | Extended-spectrum β-lactamase |
| KPC_1 | KPC carbapenemase |

### DADA2 parameter grid

| Parameter | Levels |
|-----------|--------|
| Filter mode | `Strict` (maxEE 0.5/100 bp) · `Loose` (maxEE 2/100 bp) |
| Error model | `Binned` · `PacBio` |
| Pooling | `FALSE` · `pseudo` · `TRUE` |
| Omega (ω) | `1e-40` · `1e-70` · `1e-100` · `1e-200` |
| Band size | `32` · `64` |

---

## Setup

```bash
# Create and activate the conda environment
conda env create -f dependencies/env_dada2_final.yml
source dependencies/start_env.sh
```

Key dependencies:

| Package | Source |
|---------|--------|
| `dada2` | Bioconductor |
| `phyloseq` | Bioconductor |
| `DECIPHER` + `pwalign` | Bioconductor |
| `ComplexHeatmap` | Bioconductor |
| `ggtree` + `ggtreeExtra` | Bioconductor |
| `vegan` | CRAN |
| `ggplot2` + `patchwork` | CRAN |
| `chkMocks` | GitHub — [microsud/chkMocks](https://github.com/microsud/chkMocks) |

---

## Reproducing the analysis

```bash
# Step 1 — run pipeline on the HPC cluster
# Output: data/cases/ (96 RDS files)

# Step 2 — run analysis steps locally (set working dir to project root)
```

```r
rmarkdown::render("analysis/02_community_analyses/Scripts/02_community_analyses.Rmd")
rmarkdown::render("analysis/03_metric_modelling/Scripts/03_case_metric_modeling.Rmd")
rmarkdown::render("analysis/04_best_vs_standard/Scripts/04_BEST_vs_STANDARD.Rmd")
rmarkdown::render("analysis/05_mock_qc/Scripts/05_Mock_Distributions.Rmd")
rmarkdown::render("analysis/06_longitudinal_carriage/Scripts/06_Longitudinal_Carriage.Rmd")
```

> ⚠️ **Step 04 has a manual BLAST step.** Generate FASTA files in section 3b, submit to [NCBI BLAST](https://blast.ncbi.nlm.nih.gov), download results as *Single-file JSON2*, and place in `data/blast/JSON/`. See `analysis/04_best_vs_standard/` for full instructions.

> ⚠️ **Step 06 requires step 04** to be completed first (`blast_df_parsed.csv` must exist).

---

## Contact

Open an issue for questions about the code or data.  
**Andrés Catalán Tatay** — University of Lausanne, Computational Biology Department
