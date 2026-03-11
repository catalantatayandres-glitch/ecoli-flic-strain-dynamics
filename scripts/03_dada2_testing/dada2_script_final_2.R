#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
  library(phyloseq)
  library(Biostrings)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(optparse)
})

msg <- function(...) {
  message(format(Sys.time(), "%F %T"), " | ", sprintf(...))
}

safe_dir_create <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

make_grid <- function() {
  expand_grid(
    filter_mode = c("Strict", "Loose"),
    pool        = c(FALSE, "pseudo", TRUE),
    OMEGA_A     = c(1e-40, 1e-70, 1e-100, 1e-200),
    BAND_SIZE   = c(32, 64),
    err_model   = c("Binned", "PacBio")
  ) %>%
    mutate(
      minQ  = if_else(filter_mode == "Strict", 30L, 6L),
      maxEE = if_else(filter_mode == "Strict", 0.5, 2.0)
    ) %>%
    arrange(filter_mode, err_model, pool, OMEGA_A, BAND_SIZE) %>%
    mutate(case_id = row_number())
}

get_denoised_count <- function(x) {
  sum(getUniques(x))
}

run_case <- function(
  fastq_files,
  row,
  out_dir,
  n_cores,
  case_index,
  n_cases,
  nbases = 1e8,
  minLen = 650,
  maxLen = 3000,
  maxN = 0,
  rm.phix = FALSE,
  DETECT_SINGLETONS = FALSE,
  bimera_method = "per-sample",
  minFoldParentOverAbundance = 3.5,
  minParentAbundance = 8
) {
  # ---------------------------
  # Case name
  # ---------------------------
  case_name <- sprintf(
    "case_%03d__F-%s__Err-%s__Pool-%s__O%s__B%s",
    row$case_id,
    row$filter_mode,
    row$err_model,
    as.character(row$pool),
    format(row$OMEGA_A, scientific = TRUE),
    row$BAND_SIZE
  )
  case_name <- gsub("[^A-Za-z0-9_.=-]", "", case_name)

  cases_dir <- file.path(out_dir, "cases")
  safe_dir_create(cases_dir)

  out_rds <- file.path(cases_dir, paste0(case_name, ".rds"))

  if (file.exists(out_rds)) {
    msg("[%d/%d] Skipping (already exists): %s", case_index, n_cases, out_rds)
    return(invisible(NULL))
  }

  msg("[%d/%d] Running %s", case_index, n_cases, case_name)

  # ---------------------------
  # Filtering cache
  # ---------------------------
  filt_key <- sprintf(
    "minQ%s__maxEE%s__minLen%s__maxLen%s__maxN%s__rmphix%s",
    row$minQ, row$maxEE, minLen, maxLen, maxN, rm.phix
  )
  filt_key <- gsub("[^A-Za-z0-9_.=-]", "", filt_key)

  filtered_dir <- file.path(out_dir, "filtered_cache", filt_key)
  safe_dir_create(filtered_dir)

  filts <- file.path(filtered_dir, basename(fastq_files))
  track_rds <- file.path(filtered_dir, "track_filter.rds")

  # ---------------------------
  # 1) Filter & trim (cached)
  # ---------------------------
  if (!all(file.exists(filts)) || !file.exists(track_rds)) {
    msg("[%d/%d] Filtering -> %s", case_index, n_cases, filtered_dir)

    track <- filterAndTrim(
      fwd         = fastq_files,
      filt        = filts,
      minQ        = row$minQ,
      minLen      = minLen,
      maxLen      = maxLen,
      maxN        = maxN,
      rm.phix     = rm.phix,
      maxEE       = row$maxEE,
      multithread = n_cores
    )

    track <- as.data.frame(track, stringsAsFactors = FALSE)

    if (nrow(track) != length(fastq_files)) {
      stop(sprintf(
        "filterAndTrim returned %d rows, but there are %d FASTQ files",
        nrow(track), length(fastq_files)
      ))
    }

    rownames(track) <- basename(fastq_files)
    saveRDS(track, track_rds)
    msg("[%d/%d] Saved filter track -> %s", case_index, n_cases, track_rds)

  } else {
    msg("[%d/%d] Reusing filtered FASTQs -> %s", case_index, n_cases, filtered_dir)
    track <- readRDS(track_rds)

    if (nrow(track) != length(fastq_files)) {
      stop(sprintf(
        "Cached filter track has %d rows, but there are %d FASTQ files",
        nrow(track), length(fastq_files)
      ))
    }

    rownames(track) <- basename(fastq_files)
  }

  # ---------------------------
  # 2) Error learning (cached)
  # ---------------------------
  err_fun <- switch(
    row$err_model,
    "PacBio" = dada2::PacBioErrfun,
    "Binned" = NULL,
    NULL
  )

  if (is.null(err_fun)) {
    err_fun <- makeBinnedQualErrfun(c(3, 10, 17, 22, 27, 35, 40))
  }

  err_key <- sprintf(
    "minQ%s__maxEE%s__minLen%s__maxLen%s__maxN%s__rmphix%s__errModel%s__nbases%s",
    row$minQ, row$maxEE, minLen, maxLen, maxN, rm.phix, row$err_model, nbases
  )
  err_key <- gsub("[^A-Za-z0-9_.=-]", "", err_key)

  err_dir <- file.path(out_dir, "error_cache")
  safe_dir_create(err_dir)
  err_rds <- file.path(err_dir, paste0(err_key, ".rds"))

  if (file.exists(err_rds)) {
    msg("[%d/%d] Reusing learned errors -> %s", case_index, n_cases, err_rds)
    err <- readRDS(err_rds)
  } else {
    msg("[%d/%d] Learning errors -> %s", case_index, n_cases, err_rds)

    err <- learnErrors(
      filts,
      errorEstimationFunction = err_fun,
      nbases = nbases,
      multithread = n_cores,
      randomize = TRUE,
      verbose = FALSE
    )

    saveRDS(err, err_rds)
    msg("[%d/%d] Saved learned errors -> %s", case_index, n_cases, err_rds)
  }

  plot_error <- plotErrors(err)

  # ---------------------------
  # 3) Dereplication
  # ---------------------------
  drp <- derepFastq(filts, verbose = FALSE)

  if (length(drp) != length(fastq_files)) {
    stop(sprintf(
      "Dereplication length mismatch: %d derep objects for %d FASTQ files",
      length(drp), length(fastq_files)
    ))
  }

  names(drp) <- basename(fastq_files)

  # ---------------------------
  # 4) DADA denoising
  # ---------------------------
  pool_val <- row$pool
  if (is.character(pool_val) && pool_val %in% c("TRUE", "FALSE")) {
    pool_val <- as.logical(pool_val)
  }

  setDadaOpt(
    OMEGA_A = row$OMEGA_A,
    DETECT_SINGLETONS = DETECT_SINGLETONS
  )

  dd <- dada(
    drp,
    err = err,
    BAND_SIZE = row$BAND_SIZE,
    pool = pool_val,
    multithread = n_cores,
    verbose = FALSE
  )

  if (length(dd) != length(fastq_files)) {
    stop(sprintf(
      "DADA output length mismatch: %d results for %d FASTQ files",
      length(dd), length(fastq_files)
    ))
  }

  # ---------------------------
  # 5) Sequence table
  # ---------------------------
  st <- makeSequenceTable(dd)

  if (nrow(st) != length(fastq_files)) {
    stop(sprintf(
      "Sequence table row mismatch: nrow(st)=%d, fastq_files=%d",
      nrow(st), length(fastq_files)
    ))
  }

  rownames(st) <- basename(fastq_files)

  # ---------------------------
  # 6) Chimera detection + removal
  # ---------------------------
  bim <- isBimeraDenovo(
    st,
    minFoldParentOverAbundance = minFoldParentOverAbundance,
    minParentAbundance = minParentAbundance,
    multithread = n_cores,
    verbose = FALSE
  )

  st_nochim <- removeBimeraDenovo(
    st,
    minFoldParentOverAbundance = minFoldParentOverAbundance,
    method = bimera_method,
    multithread = n_cores,
    verbose = FALSE
  )

  if (nrow(st_nochim) != length(fastq_files)) {
    stop(sprintf(
      "Nonchimera table row mismatch: nrow(st_nochim)=%d, fastq_files=%d",
      nrow(st_nochim), length(fastq_files)
    ))
  }

  rownames(st_nochim) <- basename(fastq_files)

  # ---------------------------
  # 7) Tracking
  # ---------------------------
  denoised_counts <- vapply(dd, get_denoised_count, numeric(1))
  nonchim_counts  <- rowSums(st_nochim)

  if (nrow(track) != length(denoised_counts) ||
      nrow(track) != length(nonchim_counts) ||
      nrow(track) != length(fastq_files)) {
    stop(sprintf(
      paste0(
        "Tracking mismatch: nrow(track)=%d, denoised=%d, ",
        "nonchim=%d, fastq_files=%d"
      ),
      nrow(track),
      length(denoised_counts),
      length(nonchim_counts),
      length(fastq_files)
    ))
  }

  tracking <- cbind(
    track,
    denoised = denoised_counts,
    nonchim  = nonchim_counts
  )

  rownames(tracking) <- basename(fastq_files)

  # ---------------------------
  # 8) Phyloseq
  # ---------------------------
  ps <- phyloseq(
    otu_table(t(st_nochim), taxa_are_rows = TRUE),
    sample_data(as.data.frame(tracking))
  )

  asv_seq <- DNAStringSet(taxa_names(ps))
  names(asv_seq) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, asv_seq)

  taxa_names(ps) <- paste0(
    "ASV",
    str_pad(seq_len(ntaxa(ps)), width = nchar(ntaxa(ps)), pad = "0")
  )
  sample_names(ps) <- str_remove(sample_names(ps), "\\.fastq\\.gz$")

  # ---------------------------
  # 9) Save
  # ---------------------------
  res <- list(
    physeq     = ps,
    track      = tracking,
    bim        = bim,
    st         = st,
    st_nochim  = st_nochim,
    plot_error = plot_error,
    params     = as.list(row)
  )

  saveRDS(res, out_rds)
  msg("[%d/%d] Saved %s", case_index, n_cases, out_rds)

  # help GC a bit on long runs
  rm(track, err, plot_error, drp, dd, st, bim, st_nochim, tracking, ps, asv_seq, res)
  invisible(gc())
}

# ----------------------------
# CLI
# ----------------------------
option_list <- list(
  make_option("--input_dir", type = "character", help = "Directory with primer-trimmed FASTQ.gz files"),
  make_option("--out_dir", type = "character", default = "dada2_sweep_out", help = "Output directory"),
  make_option("--pattern", type = "character", default = "fastq\\.gz$", help = "Regex pattern for FASTQs"),
  make_option("--targets", type = "character", default = NULL, help = "Optional comma-separated patterns to keep"),
  make_option("--n_cores", type = "integer", default = 6, help = "Threads"),
  make_option("--nbases", type = "double", default = 1e8, help = "Bases for learnErrors"),
  make_option("--minLen", type = "integer", default = 650, help = "Min length"),
  make_option("--maxLen", type = "integer", default = 3000, help = "Max length"),
  make_option("--case_id", type = "integer", default = NULL, help = "Optional single case_id to run")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_dir)) {
  stop("You must provide --input_dir")
}

safe_dir_create(opt$out_dir)

fastq_files <- sort(list.files(
  opt$input_dir,
  pattern = opt$pattern,
  full.names = TRUE
))

if (length(fastq_files) == 0) {
  stop("No FASTQ files found in input_dir.")
}

if (!is.null(opt$targets)) {
  tgs <- str_split(opt$targets, ",")[[1]] %>% str_trim()
  keep <- grepl(paste(tgs, collapse = "|"), basename(fastq_files))
  fastq_files <- fastq_files[keep]
}

if (length(fastq_files) == 0) {
  stop("No FASTQs left after targets filtering.")
}

msg("FASTQs: %d", length(fastq_files))

grid <- make_grid()
write.csv(grid, file.path(opt$out_dir, "grid_used.csv"), row.names = FALSE)
msg("Grid size: %d cases", nrow(grid))

if (!is.null(opt$case_id)) {
  if (!(opt$case_id %in% grid$case_id)) {
    stop(sprintf("Requested case_id %d not found in grid", opt$case_id))
  }
  grid_to_run <- grid %>% filter(case_id == opt$case_id)
  msg("Running single case_id: %d", opt$case_id)
} else {
  grid_to_run <- grid
  msg("Running all cases")
}

n_cases <- nrow(grid_to_run)

for (i in seq_len(n_cases)) {
  run_case(
    fastq_files = fastq_files,
    row         = grid_to_run[i, ],
    out_dir     = opt$out_dir,
    n_cores     = opt$n_cores,
    case_index  = i,
    n_cases     = n_cases,
    nbases      = opt$nbases,
    minLen      = opt$minLen,
    maxLen      = opt$maxLen
  )
}

msg("ALL DONE.")
