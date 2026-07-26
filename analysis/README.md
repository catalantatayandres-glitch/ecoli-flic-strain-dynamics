# Analysis

Downstream benchmarking and biological analysis of the 96 DADA2 pipeline outputs.
Scripts in this folder were executed **locally**.

> Requires `data/cases/` to be populated with the 96 RDS files produced by `pipeline/dada2/`.

---

## Structure

```
analysis/
├── 02_community_analyses/     ← community-level comparison of all 96 pipelines
├── 03_metric_modelling/       ← quality metrics, modelling, and pipeline ranking
├── 04_best_vs_standard/       ← sequence and taxonomic comparison of top pipelines
├── 05_mock_qc/                ← mock community quality control
└── 06_longitudinal_carriage/  ← E. coli carriage in real cohort samples
```

Each step folder follows the same structure:
```
0X_step_name/
├── README.md     ← methods, results, and figures for that step
├── Scripts/      ← .Rmd notebook + companion .R functions file
└── Figures/      ← key output figures
```

Steps are numbered to enforce execution order. **Step 06 depends on step 04.**

---

## Steps

### 02 · Community analyses

**Goal:** Determine which DADA2 parameters have the strongest influence on the detected ASV community, before assessing whether those differences reflect true biology or noise.

- ASV sharing analysis — classifies ASVs as shared or unique per parameter level across all 96 cases
- Jaccard (presence/absence) and Aitchison (abundance-weighted) distance matrices
- PERMANOVA — formally tests the variance explained by each parameter (999 permutations), with sample type included as a reference variable so pipeline effects can be compared directly against biological signal
- PCoA ordination — 6-panel plots per distance metric to visualise the patterns driving PERMANOVA variance

**Key result:** Pooling strategy and filter mode are the dominant pipeline-driven sources of variation. The pooling effect is much stronger under Jaccard than Aitchison, raising the possibility that `TRUE` pooling predominantly introduces low-abundance sequences that may not represent genuine biological variants.

---

### 03 · Metric modelling & pipeline ranking

**Goal:** Quantify the effect of each parameter on interpretable mock community quality metrics, and formally rank all 96 configurations to identify the optimal pipeline.

- Quality metrics on mock samples only: % chimeric reads, total ASVs detected, % ASVs matching mock, % reads from mock-matching ASVs
- Quasibinomial GLMs (proportions) and Poisson GLMs (counts) — controls for all parameters simultaneously
- Marginal means heatmap — visual overview of all parameter effects across all metrics
- **Composite score** = −z(noise) + z(recovery) → ranks all 96 configurations on two dimensions simultaneously

**Key result:** Optimal pipeline identified as `Strict filtering · Binned error model · FALSE pooling`. Omega and band size can be set freely within the tested ranges without affecting performance.

---

### 04 · Best vs. Standard pipeline comparison

**Goal:** Directly compare the optimal pipeline against the standard DADA2 defaults at the sequence level and validate the biological origin of their respective ASV sets.

> ⚠️ Includes a **manual BLAST step** — see `04_best_vs_standard/README.md`.

- Pairwise edit distance matrix (DECIPHER + pwalign Hamming) — reveals sequence-level divergence between pipelines
- Edit distance distribution — within-Best, within-Standard, and cross-pipeline comparisons
- NCBI BLAST taxonomy via Single-file JSON parsing — all ASVs submitted for taxonomic validation
- 5-tier trust system (Tier 1 = all perfect *E. coli* hits → Tier 5 = no *E. coli* hits) collapsed into E. coli / Ambiguous / Non-E. coli biological groups
- Integrated dendrogram combining clustering, trust tiers, abundance, prevalence, and alignment
- Bio group composition by ASV origin (Shared / Unique-to-Best / Unique-to-Standard)

**Key result:** Over 90% of ASVs unique to the Standard pipeline lack reliable *E. coli* biological validation. The shared ASV pool has the highest proportion of confirmed *E. coli* sequences, confirming the core signal is robustly captured by both pipelines.

---

### 05 · Mock QC

**Goal:** Assess how faithfully the experimental mock community positive controls reflect the theoretical composition of 4 known *E. coli fliC* sequences, and evaluate whether unexpected sequences originate from pipeline noise or biological contamination.

- Pairwise alignment of experimental ASVs against the 4 theoretical reference sequences (replaces `checkZymoBiomics()`, which uses the Zymo 8-species training set incompatible with this reference)
- Bray-Curtis similarity between each mock sample and the theoretical composition
- Stacked composition plots for BEST and STANDARD pipelines side by side

**Key result:** Both pipelines detect unexpected sequences beyond the 4 expected strains. The BEST pipeline fails to recover Oxa48_1 (likely due to low concentration in the mock preparation), while the STANDARD pipeline recovers it marginally. Subsequent experiments confirmed the **Kinnex library preparation** as the likely contamination source.

---

### 06 · Longitudinal carriage

**Goal:** Apply both validated pipelines to real longitudinal cohort data to profile *E. coli* strain carriage over time and assess whether pipeline choice affects biological interpretation.

- Attaches longitudinal metadata (Subject ID, days since baseline, sample type) from the reference phyloseq
- Retains only ASVs with perfect BLAST support (100% identity, 100% query coverage, *Escherichia coli*)
- Bubble plots per subject × pipeline: circle size encodes relative abundance at each timepoint; dashed lines mark all sampled timepoints
- Outputs a master detection dataframe and summary table per subject

**Key result:** Both pipelines consistently recover the same dominant persistent strains. Differences arise only in low-abundance, low-prevalence detections whose biological origin — transient colonisation, detection limit, or contamination — cannot be resolved without paired WGS data.

---

## Running all steps

```r
# Set working directory to the project root, then:
rmarkdown::render("analysis/02_community_analyses/Scripts/02_community_analyses.Rmd")
rmarkdown::render("analysis/03_metric_modelling/Scripts/03_case_metric_modeling.Rmd")
rmarkdown::render("analysis/04_best_vs_standard/Scripts/04_BEST_vs_STANDARD.Rmd")
rmarkdown::render("analysis/05_mock_qc/Scripts/05_Mock_Distributions.Rmd")
rmarkdown::render("analysis/06_longitudinal_carriage/Scripts/06_Longitudinal_Carriage.Rmd")
```

Each subfolder contains its own `README.md` with detailed methods, figures, and conclusions.
