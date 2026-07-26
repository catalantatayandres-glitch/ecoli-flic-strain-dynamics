# 06 — Longitudinal E. coli ASV Carriage

Longitudinal profiling of *E. coli* ASV carriage across repeated faecal samples from
the same individual, comparing detection patterns between the BEST and STANDARD pipelines.
Only ASVs with **perfect BLAST support** (100% identity, 100% query coverage,
*Escherichia coli* annotation) are retained.

## Contents

| File | Description |
|------|-------------|
| `Scripts/06_Longitudinal_Carriage.Rmd` | Main analysis notebook |
| `Scripts/06_Longitudinal_Carriage_functions.R` | Companion functions |

---

## What this analysis does

1. Attaches longitudinal metadata (Subject ID, days since baseline, sample type) from the reference phyloseq
2. Filters to individual samples only — removes mock, negative control, and blank samples
3. Parses BLAST results from step 04; retains only ASVs with perfect *E. coli* hits
4. Plots carriage per Subject × Pipeline: filled circles sized by relative abundance, dashed lines marking all sampled timepoints

> ⚠️ **Requires step 04** — `blast_df_parsed.csv` must exist before running.

## How to run

```bash
source start_env.sh
```
```r
rmarkdown::render("analysis/06_longitudinal_carriage/Scripts/06_Longitudinal_Carriage.Rmd")
```

> The longitudinal metadata phyloseq was saved with `save()` rather than `saveRDS()` — `load_florentin_phyloseq()` handles this automatically.

---

## Results

Subject 813 was selected for presentation — it has the longest follow-up (~650 days)
and the richest carriage pattern in the cohort.

To profile *E. coli* strain carriage over time, bubble plots are produced for each
Subject × Pipeline combination. Each row represents one *E. coli* ASV (identified by
its global identifier); each bubble indicates a detection at a given timepoint with
**bubble size proportional to relative abundance** among all *E. coli fliC* reads at
that timepoint. Vertical dashed lines mark every sampled timepoint — an empty row at
a dashed line means that ASV was absent at that sampling, not that the sample is missing.
Comparing BEST and STANDARD side by side directly reveals whether pipeline choice
affects which strains are detected, how persistently, and at what abundance.

### Subject 813 — BEST Pipeline

![Longitudinal carriage 813 BEST](Figures/longitudinal_carriage_813_BEST.png)

### Subject 813 — STANDARD Pipeline

![Longitudinal carriage 813 STANDARD](Figures/longitudinal_carriage_813_STANDARD.png)

---

### Key findings

**Both pipelines recover the same dominant persistent strains**
- A core set of high-abundance ASVs — most notably **gASV034** and **gASV033** — is consistently detected across multiple timepoints spanning the full follow-up period in both pipelines
- This is the most important result: the primary biological signal — *E. coli* strains that truly colonise and persist — is robustly captured regardless of pipeline choice
- Several other ASVs (e.g. gASV002, gASV012, gASV028, gASV029) appear at comparable timepoints and abundances in both pipelines, further confirming reproducibility of the core carriage pattern

**The STANDARD pipeline detects more ASVs and extends temporal coverage**
- The STANDARD pipeline recovers a higher number of distinct *E. coli* ASVs overall, and captures some strains across a greater number of timepoints
- Additional detections unique to STANDARD are predominantly low-abundance, isolated detections (tiny bubbles at single timepoints), consistent with its higher sensitivity identified in the mock benchmarking

**Interpretation of additional detections in the STANDARD pipeline**
- Their biological interpretation is uncertain — several explanations are possible and cannot be distinguished without paired WGS data: genuine but transient colonisation episodes, very brief carriage events at the limit of detection, or *E. coli* contamination from wet-lab steps (which would pass the perfect BLAST filter as it constitutes a real *E. coli fliC* sequence)
- Even strains genuinely present may not be consistently recovered across all samples due to uneven spatial distribution within the gut or stochastic variation during DNA extraction

**Pipeline choice does not affect the core biological conclusion**
- Both pipelines support the same interpretation: a small number of dominant strains are persistently carried, with lower-abundance ASVs appearing transiently
- Whether those low-abundance, low-prevalence variants represent true biology or noise cannot be resolved from amplicon data alone — validation against paired WGS data represents the most important direction for future work

---

### Pipeline recommendation

The choice between pipelines should be guided by the biological question:

- **High sensitivity needed** (detect rare or transient strains, maximise *E. coli* diversity) → **STANDARD pipeline**. It recovers more ASVs and broader temporal coverage, but requires additional validation of low-confidence detections and increases downstream computational and curation burden.

- **High specificity needed** (transmission studies, focus on dominant persistent strains, strict confidence in biological origin) → **BEST pipeline**. It generates fewer ASVs, streamlines the analytical workflow, and reduces the risk of biological misinterpretation.

---

*Full methods in `Scripts/06_Longitudinal_Carriage.Rmd` and the project report.*
