# 02B_trim_primers_dada2.R
rm(list = ls()); gc()

library(ShortRead)
library(dada2)
library(readr)
library(dplyr)
library(tibble)
library(Biostrings)

# ----------------------------
# PATHS (cluster)
# ----------------------------
base_dir     <- "/users/acatala3/first_step/pilot_analysis/samples_pilot"
manifest_in  <- file.path(base_dir, "manifest_raw.tsv")

raw_dir      <- file.path(base_dir, "raw_samples")
trim_dir     <- file.path(base_dir, "trimmed_samples")
dir.create(trim_dir, recursive = TRUE, showWarnings = FALSE)

manifest_out <- file.path(base_dir, "manifest_trimmed.tsv")

# ----------------------------
# PRIMERS
# ----------------------------
primer_fwd   <- "GGTCAGGCGATTGCTAACCG"
primer_rev   <- "GACACTTCGGTCGCGTAGTC"
max_mismatch <- 2

# PI-style: use reverse-complement for right-side trimming
primer_rev_rc <- dada2:::rc(primer_rev)

# ----------------------------
# 1) Read manifest
# ----------------------------
man <- read_tsv(manifest_in, show_col_types = FALSE)

if (!("SampleID" %in% colnames(man))) {
  stop("manifest_raw.tsv must contain a column named 'SampleID'.\n",
       "Your columns are: ", paste(colnames(man), collapse = ", "))
}

fn_in <- file.path(raw_dir, paste0(man$SampleID, ".fastq.gz"))

exists_in <- file.exists(fn_in)
fn_in_ok  <- fn_in[exists_in]
ids_ok    <- man$SampleID[exists_in]

message("Manifest rows: ", nrow(man))
message("FASTQs found in raw_samples: ", length(fn_in_ok))

if (length(fn_in_ok) == 0) {
  stop("No input FASTQs found in: ", raw_dir,
       "\nCheck SampleID naming (should match filenames without .fastq.gz).")
}

fn_out_ok <- file.path(trim_dir, basename(fn_in_ok))

# ----------------------------
# 2) QC: reads before
# ----------------------------
count_reads_1 <- function(f) {
  x <- ShortRead::countFastq(f)
  if (length(x) == 1) return(as.integer(x))
  as.integer(sum(x))
}
reads_in <- vapply(fn_in_ok, count_reads_1, integer(1))

# ----------------------------
# 3) Remove primers (single-end) — PI-style trimLRPatterns on FASTQ
# ----------------------------
# Trims BOTH sequences + quality.
# Reads without matching patterns are left unchanged.
# Mismatch policy:
#   forward: 1 mismatch (PI strict)
#   reverse: up to 2 mismatches

for (i in seq_along(fn_in_ok)) {
  in_f  <- fn_in_ok[i]
  out_f <- fn_out_ok[i]

  fq <- ShortRead::readFastq(in_f)

  fq_trim <- ShortRead::trimLRPatterns(
    Lpattern = primer_fwd,
    Rpattern = primer_rev_rc,
    subject = fq,
    max.Lmismatch = max_mismatch,
    max.Rmismatch = max_mismatch
  )

  ShortRead::writeFastq(fq_trim, out_f, compress = TRUE)
}

# ----------------------------
# 4) QC: reads after
# ----------------------------
reads_out <- vapply(fn_out_ok, count_reads_1, integer(1))

# ----------------------------
# 5) Write updated manifest
# ----------------------------
proc <- tibble(
  SampleID       = ids_ok,
  fastq_raw      = fn_in_ok,
  fastq_trimmed  = fn_out_ok,
  reads_in       = as.integer(reads_in),
  reads_out      = as.integer(reads_out),
  pct_kept       = ifelse(reads_in > 0, 100 * reads_out / reads_in, NA_real_)
)

man_out <- man %>% left_join(proc, by = "SampleID")
write_tsv(man_out, manifest_out)

message("Done ✅")
message("Trimmed FASTQs in: ", trim_dir)
message("Updated manifest: ", manifest_out)
