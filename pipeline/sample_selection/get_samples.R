# 02A_select_and_copy_fastq.R
rm(list = ls()); gc()

library(phyloseq)
library(dplyr)
library(tibble)
library(readr)

# ------------------------------------------------
# CONFIG (CLUSTER PATHS)
# ------------------------------------------------
ps_path   <- "/users/acatala3/first_step/phyloseq_objects/ps_metaphlan.rds"

# Folder where ALL fastq files are currently
fastq_dir <- "/users/acatala3/first_step/sub_sampled"

# Destination folder (new folder with only selected samples)
dest_dir  <- "/users/acatala3/first_step/pilot_analysis/samples_pilot"

manifest_path <- file.path(dest_dir, "selected_samples.tsv")

# ------------------------------------------------
# 1) Load phyloseq and select individuals by metadata
# ------------------------------------------------
ps_metaphlan <- readRDS(ps_path)

ps_metaphlan %>%
  subset_samples(Subject %in% c(813, 644, 431)) %>%
  sample_data() %>%
  data.frame(stringsAsFactors = FALSE) %>%
  pull(Sample) -> indiv_ids

indiv_ids <- unique(as.character(indiv_ids))

# ------------------------------------------------
# 2) Select mocks + negatives by filename
# ------------------------------------------------
fn_all <- list.files(fastq_dir, pattern="\\.fastq\\.gz$", full.names=FALSE)
ids_all <- sub("\\.fastq\\.gz$", "", fn_all)

# A) Mock mix non-diluted
mockmix_ids <- ids_all[
  grepl("^S_Mock_mix_", ids_all, ignore.case = TRUE) &
  !grepl("dil", ids_all, ignore.case = TRUE)
]

# B) Single-sequence mocks (non-diluted only)
single_mock_ids <- c("S_E15", "S_E25", "S_KPC1", "S_Oxa")
single_mock_ids <- intersect(ids_all, single_mock_ids)

# C) 3 Tris negatives (choose which 3 you want)
tris_ids <- c("S_Tris_2-13", "S_Tris_2-16", "S_Tris_2-19")
tris_ids <- intersect(ids_all, tris_ids)

# ------------------------------------------------
# 3) Combine all targets
# ------------------------------------------------
targets <- unique(c(indiv_ids, mockmix_ids, single_mock_ids, tris_ids))

# Full paths
source_files <- file.path(fastq_dir, paste0(targets, ".fastq.gz"))

# Keep only existing files
exists <- file.exists(source_files)
source_files <- source_files[exists]
targets <- targets[exists]

# ------------------------------------------------
# 4) Copy files
# ------------------------------------------------
file.copy(from = source_files,
          to   = file.path(dest_dir, basename(source_files)),
          overwrite = FALSE)

# ------------------------------------------------
# 5) Write manifest
# ------------------------------------------------
manifest <- tibble(
  SampleID = targets,
  SourcePath = source_files,
  CopiedTo   = file.path(dest_dir, basename(source_files))
)

write_tsv(manifest, manifest_path)

message("Copied files: ", length(source_files))
message("Destination: ", dest_dir)
message("Manifest written to: ", manifest_path)
