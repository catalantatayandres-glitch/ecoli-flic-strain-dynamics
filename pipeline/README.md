# Pipeline

Upstream processing of raw PacBio *fliC* amplicon reads into DADA2 ASV outputs.
Scripts in this folder were executed on the **UNIL HPC cluster**.

---

## Steps

```
pipeline/
├── sample_selection/   ← 1. subsample reads; select cohort & mock samples
├── primer_trimming/    ← 2. remove PCR primers
└── dada2/              ← 3. run DADA2 across all 96 parameter combinations
```

### 1 · `sample_selection/`
- Selects the three cohort individuals (IDs: 431, 644, 813) used in the analysis
- Retains **mock mix** positive controls and **negative controls** for downstream benchmarking

### 2 · `primer_trimming/`
- Removes PCR primer sequences from both ends of each read using **cutadapt**
- Primer-trimmed reads are required before DADA2 to avoid interference with error modelling and sequence inference

### 3 · `dada2/`
- Runs DADA2 over all **96 parameter combinations** (2 filter × 2 error model × 3 pool × 4 omega × 2 band size)
- Each run produces one phyloseq RDS file saved as `data/cases/case_NNN__<params>.rds`
- See the [project README](../README.md) for the full parameter grid

---

## Input

Raw PacBio CCS FASTQ files (not tracked in this repo).

## Output

```
data/cases/
└── case_NNN__F-<filter>__Err-<err>__Pool-<pool>__O<omega>__B<band>.rds
    (96 files — one per parameter combination)
```

These RDS files are the input for everything in `analysis/`.
