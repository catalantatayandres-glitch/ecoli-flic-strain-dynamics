# Analysis

Downstream benchmarking and biological analysis of the 96 DADA2 pipeline outputs.
Scripts in this folder were executed **locally**.

> Requires `data/cases/` to be populated with the 96 RDS files produced by `pipeline/dada2/`.

---

## Steps

```
analysis/
├── 02_community_analyses/     ← community-level comparison of all 96 pipelines
├── 03_metric_modelling/       ← quality metrics, modelling, and pipeline ranking
├── 04_best_vs_standard/       ← sequence and taxonomic comparison of top pipelines
├── 05_mock_qc/                ← mock community quality control
└── 06_longitudinal_carriage/  ← E. coli carriage in real cohort samples
```

Steps are numbered to enforce execution order. **Step 06 depends on step 04.**

---

### 02 · Community analyses

Assesses how each of the 5 DADA2 parameters shapes the detected ASV community across all 96 configurations.

- ASV sharing analysis — which ASVs are shared vs. unique per parameter level
- Jaccard and Aitchison distance matrices
- PERMANOVA — variance explained by each parameter (999 permutations)
- PCoA ordination — 6-panel plots per distance metric
- Hierarchical clustering dendrogram with parameter annotation strips

---

### 03 · Metric modelling & pipeline ranking

Quantifies parameter effects on mock community quality and identifies the optimal configuration.

- Quality metrics computed on mock samples only: % chimeric reads, total ASVs, % ASVs matching mock, % reads matching mock
- Quasibinomial GLMs and Poisson GLMs per metric
- Marginal means heatmap across all parameters
- **Composite score** = −z(noise) + z(recovery) → ranks all 96 configurations
- Identifies the **optimal pipeline**: `Strict · Binned · FALSE pooling`

---

### 04 · Best vs. Standard pipeline comparison

> ⚠️ Includes a **manual BLAST step** — see `04_best_vs_standard/README.md`.

Direct sequence-level and biological comparison of the BEST and STANDARD pipelines.

- Pairwise edit distance matrix (DECIPHER + pwalign Hamming)
- Three-panel edit distance distribution plot (within-Best, within-Standard, cross-pipeline)
- NCBI BLAST taxonomy via Single-file JSON parsing
- 5-tier trust system: Tier 1 (all perfect *E. coli*) → Tier 5 (no *E. coli* hits)
- Integrated dendrogram: clustering + trust tiers + abundance + alignment
- Bio group composition by ASV origin (Shared / Unique-to-Best / Unique-to-Standard)

---

### 05 · Mock QC

Assesses how well the experimental mock samples match the theoretical composition.

- Pairwise alignment of experimental ASVs against the 4 theoretical reference sequences
- Manual taxonomy assignment (replaces `checkZymoBiomics()`, which uses the Zymo 8-species set)
- Bray-Curtis similarity between each mock sample and the theoretical mock
- Stacked composition plots for BEST and STANDARD pipelines side by side

---

### 06 · Longitudinal carriage

Applies the validated pipelines to real cohort data to profile *E. coli* strain carriage over time.

- Attaches longitudinal metadata (Subject ID, days since baseline) from the reference phyloseq
- Retains only ASVs with perfect BLAST support (100% identity, 100% coverage, *Escherichia coli*)
- Bubble plots per subject × pipeline: circle size encodes relative abundance at each timepoint
- Outputs a master detection dataframe and summary table per subject

---

## Running all steps

```r
# Set working directory to the project root, then:
rmarkdown::render("analysis/02_community_analyses/02_community_analyses.Rmd")
rmarkdown::render("analysis/03_metric_modelling/03_case_metric_modeling.Rmd")
rmarkdown::render("analysis/04_best_vs_standard/04_BEST_vs_STANDARD.Rmd")
rmarkdown::render("analysis/05_mock_qc/05_Mock_Distributions.Rmd")
rmarkdown::render("analysis/06_longitudinal_carriage/06_Longitudinal_Carriage.Rmd")
```

Each subfolder contains its own `README.md` with detailed input/output descriptions.
