# ==============================================================================
# 02_community_analyses_functions.R
#
# Project : fliC Amplicon Sequencing for High-Resolution Profiling of E. coli
#           in the Infant Gut Microbiome
# Author  : Andrés Catalán Tatay
# Part    : Analysis — Step 02 (companion functions)
# Purpose : All helper functions for 02_community_analyses.Rmd.
#           Covers data loading, metadata enrichment, mock removal,
#           rarefaction, distance computation, PERMANOVA, PCoA, ASV sharing
#           analysis, and all associated plots.
#
# Parameters tested:
#   Filtering   : Loose | Strict
#   Pooling     : FALSE | pseudo | TRUE
#   Omega       : O1e-40 | O1e-70 | O1e-100 | O1e-200
#   Band_Size   : 32 | 64
#   Error_Model : Binned | PacBio
#   Sample type : Individual | Positive | Negative
#
# Inputs  : data/cases/   — 96 RDS files (one per DADA2 parameter combination)
# Outputs : Called by 02_community_analyses.Rmd; outputs written by that script
#
# Usage   : source("02_community_analyses_functions.R")
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(vegan)
  library(ape)
  library(patchwork)
  library(RColorBrewer)
  library(stringr)
  library(ggrepel)
  library(scales)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(parallel)
})

# ==============================================================================
# SECTION 1 — FILE PARSING & LOADING
# ==============================================================================

#' Parse parameter metadata from an RDS filename
#'
#' Expected format:
#'   case_NNN__F-<filter>__Err-<err>__Pool-<pool>__O<omega>__B<batch>.rds
#'
#' @param filepath  Full or relative path
#' @return Named list: case, Filter_Mode, Error_Model, Pool, Omega, Band_Size
parse_filename_params <- function(filepath) {
  fn <- basename(filepath)
  list(
    case        = str_match(fn, "case_(\\d+)")[1L, 2L],
    Filter_Mode = str_match(fn, "__F-([^_]+)")[1L, 2L],
    Error_Model = str_match(fn, "__Err-([^_]+)")[1L, 2L],
    Pool        = str_match(fn, "__Pool-([^_]+)")[1L, 2L],
    Omega       = str_match(fn, "__(O1e-\\d+)")[1L, 2L],
    Band_Size   = str_match(fn, "__B(\\d+)\\.rds$")[1L, 2L]
  )
}

#' Internal: extract the first phyloseq object from any object
#' @noRd
extract_phyloseq <- function(obj, filepath = "") {
  if (inherits(obj, "phyloseq")) return(obj)
  if (is.list(obj)) {
    hits <- which(vapply(obj, inherits, logical(1L), "phyloseq"))
    if (length(hits) > 0L) return(obj[[hits[1L]]])
  }
  stop("Cannot find a phyloseq object in: ", basename(filepath))
}

#' Load all 96 cases from a directory of RDS files
#'
#' Loads every .rds file found in rds_dir without any filename filtering.
#' Expected result: 2 filter × 3 pool × 4 omega × 2 band × 2 error = 96 files.
#'
#' @param rds_dir  Directory with all case RDS files
#' @param n_cores  Cores for parallel loading (1 on Windows)
#' @return Named list; each element = list(ps, params, file)
load_all_cases <- function(rds_dir, n_cores = 1L) {

  files <- list.files(rds_dir, pattern = "\\.rds$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0L)
    stop("No RDS files found in: ", rds_dir)

  message("Found ", length(files), " RDS files in directory (loading all).")

  load_one <- function(f) {
    obj    <- readRDS(f)
    ps     <- extract_phyloseq(obj, f)
    params <- parse_filename_params(f)
    list(ps = ps, params = params, file = basename(f))
  }

  total <- length(files)
  if (n_cores > 1L && .Platform$OS.type != "windows") {
    message("Loading ", total, " RDS files in parallel (", n_cores, " cores) ...")
    cases <- parallel::mclapply(files, load_one, mc.cores = n_cores)
  } else {
    message("Loading ", total, " RDS files sequentially ...")
    cases <- vector("list", total)
    for (i in seq_len(total)) {
      cases[[i]] <- load_one(files[i])
      if (i %% 10L == 0L || i == total)
        message("  Loaded ", i, " / ", total, " cases")
    }
  }

  names(cases) <- vapply(cases, `[[`, character(1L), "file")
  message("Done – loaded ", length(cases), " cases.")
  cases
}

# ==============================================================================
# SECTION 2 — PURE MOCK REMOVAL
# ==============================================================================

#' Names (patterns) that identify pure mock control samples
#'
#' These are the individual-strain controls, NOT the MOCKMIX community.
#' Pattern is case-insensitive to handle S_KPC1, S_OXA, S_oxa, etc.
#'
#' *** VERIFY THESE PATTERNS AGAINST YOUR ACTUAL SAMPLE NAMES ***
#' If your filenames use different spellings, add them here.
PURE_MOCK_PATTERNS <- c(
  "^S_KPC",    # covers S_KPC1, S_KPC2, S_KPC…
  "^S_OXA",    # covers S_OXA, S_oxa, S_OXA48…
  "^S_E15$",   # exact match for S_E15
  "^S_E25$"    # exact match for S_E25
)

#' Test whether a sample name matches any pure mock pattern
#'
#' @param sample_names Character vector
#' @return Logical vector, TRUE = pure mock
is_pure_mock <- function(sample_names) {
  # Strip the __caseN suffix if present before matching
  base_names <- sub("__case\\d+$", "", sample_names)
  pattern    <- paste(PURE_MOCK_PATTERNS, collapse = "|")
  grepl(pattern, base_names, ignore.case = TRUE)
}

#' Remove pure mock samples from a phyloseq object
#'
#' @param ps  phyloseq
#' @return    phyloseq with pure mock samples removed
remove_pure_mocks <- function(ps) {
  sn      <- sample_names(ps)
  is_mock <- is_pure_mock(sn)
  n_rem   <- sum(is_mock)
  if (n_rem > 0L) {
    message("  Removing ", n_rem, " pure mock sample(s): ",
            paste(sn[is_mock], collapse = ", "))
    ps <- prune_samples(!is_mock, ps)
  }
  ps
}

# ==============================================================================
# SECTION 3 — SAMPLE CLASSIFICATION
# ==============================================================================

#' Classify samples using actual naming conventions
#'
#' After removing pure mocks, three types remain:
#'   Negative   : name contains "TRIS"              → S_TRIS (buffer only)
#'   Positive   : name contains MOCKMIX or MOCK     → community mock mix
#'   Individual : S_ followed by digits (real patient/env samples)
#'
#' NOTE: KPC / OXA / E15 / E25 are intentionally absent here because
#' those pure mock samples are removed BEFORE classification.
#' If any accidentally remain they will fall through to "Individual"
#' and trigger a warning at the verify step in the Rmd.
#'
#' @param sample_names Character vector
#' @return Factor with levels Individual, Positive, Negative
classify_sample_type <- function(sample_names) {
  type <- dplyr::case_when(
    str_detect(sample_names,
               regex("TRIS", ignore_case = TRUE))               ~ "Negative",
    str_detect(sample_names,
               regex("MOCKMIX|MOCK", ignore_case = TRUE))       ~ "Positive",
    str_detect(sample_names,
               regex("^S_\\d+$", ignore_case = TRUE))           ~ "Individual",
    TRUE                                                          ~ "Individual"
  )
  factor(type, levels = c("Individual", "Positive", "Negative"))
}

# ==============================================================================
# SECTION 4 — METADATA ENRICHMENT
# ==============================================================================

#' Attach case-level parameters to sample_data and set ordered factors
#'
#' Includes Error_Model and Band_Size in addition to the original parameters.
#'
#' @param ps     phyloseq object
#' @param params Named list from parse_filename_params()
#' @return phyloseq with updated sample_data
enrich_sample_data <- function(ps, params) {
  sd_raw <- as(sample_data(ps), "data.frame")
  sd     <- data.frame(
    lapply(sd_raw, function(x)
      if (is.factor(x) || is.ordered(x)) as.character(x) else x),
    stringsAsFactors = FALSE,
    row.names        = rownames(sd_raw)
  )

  sd$Filter_Mode  <- params$Filter_Mode
  sd$Pool         <- params$Pool
  sd$Omega        <- params$Omega
  sd$Error_Model  <- params$Error_Model
  sd$Band_Size    <- params$Band_Size
  sd$case         <- params$case

  if (!"sample_type" %in% colnames(sd))
    sd$sample_type <- classify_sample_type(rownames(sd))

  # Ordered factors for consistent axis/legend ordering
  sd$Omega        <- factor(sd$Omega,
                             levels = c("O1e-40","O1e-70","O1e-100","O1e-200"))
  sd$Pool         <- factor(sd$Pool,         levels = c("FALSE","pseudo","TRUE"))
  sd$Filter_Mode  <- factor(sd$Filter_Mode,  levels = c("Loose","Strict"))
  sd$Error_Model  <- factor(sd$Error_Model,  levels = c("Binned","PacBio"))
  sd$Band_Size    <- factor(sd$Band_Size,    levels = c("32","64"))
  sd$sample_type  <- factor(sd$sample_type,
                             levels = c("Individual","Positive","Negative"))

  sample_data(ps) <- sample_data(sd)
  ps
}

# ==============================================================================
# SECTION 5 — RAREFACTION
# ==============================================================================

#' Compute the rarefaction depth from all 96 cases
#'
#' Depth = minimum number of reads found across ALL Individual and Positive
#' (mockmix) samples in all 96 phyloseq objects.
#' Negative controls are excluded from depth calculation.
#'
#' @param cases_enriched  List of phyloseq objects (already enriched & mocks removed)
#' @return Integer scalar: rarefaction depth
compute_rarefy_depth <- function(cases_enriched) {
  message("Calculating minimum sequencing depth across all 96 cases, ",
          "excluding negative controls and pure mocks...")

  all_depths <- unlist(lapply(cases_enriched, function(ps) {
    sd   <- as.data.frame(sample_data(ps))
    keep <- sd$sample_type %in% c("Individual", "Positive")
    if (!any(keep)) return(numeric(0L))
    sums <- sample_sums(ps)[rownames(sd)[keep]]
    sums[sums > 0]
  }))

  if (length(all_depths) == 0L)
    stop("No Individual/Positive samples found across any case. ",
         "Check sample classification.")

  min_depth <- min(all_depths)
  message("Minimum depth found: ", min_depth, " reads")
  message("This value will be used as the rarefaction depth.")
  as.integer(min_depth)
}

#' Apply rarefy_even_depth() to a list of phyloseq objects
#'
#' Logic:
#'   1. Rarefy Individual + Positive samples to min_depth.
#'   2. Negative controls with reads >= min_depth are rarefied too.
#'   3. Negative controls with reads < min_depth are kept UNRAREFIED
#'      (a warning is issued; they are flagged with a metadata column).
#'
#' @param cases_enriched List of enriched phyloseq objects (pure mocks removed)
#' @param min_depth      Rarefaction depth from compute_rarefy_depth()
#' @param seed           Random seed for reproducibility (default 42)
#' @return List of rarefied phyloseq objects (same structure as input)
rarefy_cases <- function(cases_enriched, min_depth, seed = 42L) {
  message("Applying rarefy_even_depth() to normalize sequencing depth ...")
  message("  Rarefaction depth: ", min_depth, " reads  |  seed: ", seed)

  total        <- length(cases_enriched)
  rarefied_lst <- vector("list", total)
  names(rarefied_lst) <- names(cases_enriched)

  for (i in seq_len(total)) {
    ps    <- cases_enriched[[i]]
    sd    <- as.data.frame(sample_data(ps))
    sums  <- sample_sums(ps)

    # Identify ALL samples with too few reads (any sample type)
    low_mask      <- sums < min_depth
    low_names     <- rownames(sd)[low_mask]

    if (any(low_mask)) {
      warning("Case ", i, ": ", sum(low_mask), " sample(s) have < ",
              min_depth, " reads and will be kept unrarefied: ",
              paste(low_names, collapse = ", "))
    }

    # Samples to rarefy: everything with reads >= min_depth
    to_rarefy   <- rownames(sd)[!low_mask]
    to_skip     <- low_names

    ps_sub_rarefy <- prune_samples(to_rarefy, ps)
    # Suppress the default "X samples removed" messages from rarefy_even_depth
    suppressMessages(
      ps_sub_r <- phyloseq::rarefy_even_depth(
        ps_sub_rarefy,
        sample.size  = min_depth,
        rngseed      = seed,
        replace      = FALSE,
        trimOTUs     = TRUE,
        verbose      = FALSE
      )
    )

    if (length(to_skip) > 0L) {
      ps_skip <- prune_samples(to_skip, ps)
      # Mark these samples
      sd_skip <- as.data.frame(sample_data(ps_skip))
      sd_skip$rarefied <- FALSE
      sample_data(ps_skip) <- sample_data(sd_skip)

      sd_r <- as.data.frame(sample_data(ps_sub_r))
      sd_r$rarefied <- TRUE
      sample_data(ps_sub_r) <- sample_data(sd_r)

      ps_out <- merge_phyloseq(ps_sub_r, ps_skip)
    } else {
      sd_r <- as.data.frame(sample_data(ps_sub_r))
      sd_r$rarefied <- TRUE
      sample_data(ps_sub_r) <- sample_data(sd_r)
      ps_out <- ps_sub_r
    }

    rarefied_lst[[i]] <- ps_out

    if (i %% 10L == 0L || i == total)
      message("  Rarefied ", i, " / ", total, " cases")
  }

  message("Rarefaction complete.")
  rarefied_lst
}

# ==============================================================================
# SECTION 6 — MERGING
# ==============================================================================

#' Enrich, remove pure mocks, and return a list of ready phyloseq objects
#'
#' This wrapper applies enrich_sample_data() + remove_pure_mocks() to each
#' case and makes sample names unique with a __caseN suffix.
#'
#' @param cases   List from load_all_cases()
#' @param n_cores Cores for parallel enrichment
#' @return List of enriched, mock-removed phyloseq objects (not yet merged)
enrich_and_clean_cases <- function(cases, n_cores = 1L) {
  message("Enriching metadata and removing pure mock samples ...")

  enrich_one <- function(i) {
    ps_i <- enrich_sample_data(cases[[i]]$ps, cases[[i]]$params)
    # Remove pure mocks BEFORE appending case suffix
    ps_i <- remove_pure_mocks(ps_i)
    sample_names(ps_i) <- paste0(sample_names(ps_i), "__case", i)
    ps_i
  }

  total <- length(cases)
  if (n_cores > 1L && .Platform$OS.type != "windows") {
    ps_list <- parallel::mclapply(seq_len(total), enrich_one, mc.cores = n_cores)
  } else {
    ps_list <- lapply(seq_len(total), enrich_one)
  }

  names(ps_list) <- names(cases)

  n_removed <- vapply(ps_list, function(ps) {
    nsamples(cases[[1L]]$ps) - nsamples(ps)
  }, integer(1L))
  message("  Removed pure mock samples: ",
          sum(n_removed), " sample(s) removed in total across all cases.")
  message("  (approximately ", round(mean(n_removed), 1), " per case)")

  ps_list
}

#' Merge all 96 cases into one phyloseq, making sample names unique
#'
#' Call AFTER enrich_and_clean_cases() and rarefy_cases().
#'
#' @param ps_list  List of enriched, rarefied phyloseq objects
#' @return Single merged phyloseq
merge_all_cases <- function(ps_list) {
  message("Merging ", length(ps_list), " phyloseq objects ...")
  ps_all <- do.call(merge_phyloseq, ps_list)
  message("  Result: ", nsamples(ps_all), " samples | ", ntaxa(ps_all), " ASVs")
  ps_all
}

# ==============================================================================
# SECTION 7 — ASV ORIGIN ANALYSIS
# ==============================================================================

#' Identify ASVs unique to or shared across ALL levels of one parameter
#'
#' Works for any of the 5 parameters: Filter_Mode, Pool, Omega,
#' Band_Size, Error_Model.
#'
#' @param ps_all  Merged phyloseq (all 96 cases combined)
#' @param param   Column in sample_data
#' @return data.frame: ASV, n_levels_present, origin,
#'                     plus one logical column per level (<param>_<level>)
summarise_asv_sharing <- function(ps_all, param) {

  sd_raw <- as(sample_data(ps_all), "data.frame")
  sd     <- data.frame(
    lapply(sd_raw, function(x)
      if (is.factor(x) || is.ordered(x)) as.character(x) else x),
    stringsAsFactors = FALSE,
    row.names        = rownames(sd_raw)
  )

  if (!param %in% colnames(sd))
    stop("Column '", param, "' not found in sample_data.")

  lvls <- levels(sd[[param]])
  if (is.null(lvls)) lvls <- sort(unique(as.character(sd[[param]])))
  lvls <- lvls[!is.na(lvls)]

  if (length(lvls) < 2L)
    stop("'", param, "' has fewer than 2 non-NA levels – cannot compare.")

  message("  ASV sharing for '", param,
          "'  levels: [", paste(lvls, collapse = " | "), "]",
          "  (using ALL cases, no restriction on other parameters)")

  asv_sets <- lapply(lvls, function(lv) {
    row_sel    <- !is.na(sd[[param]]) & as.character(sd[[param]]) == lv
    keep_samps <- rownames(sd)[row_sel]
    if (length(keep_samps) == 0L) {
      warning("No samples found for ", param, " == ", lv); return(character(0L))
    }
    ps_sub <- prune_samples(keep_samps, ps_all)
    ps_sub <- filter_taxa(ps_sub, function(x) sum(x) > 0, prune = TRUE)
    taxa_names(ps_sub)
  })
  names(asv_sets) <- lvls

  all_asvs <- unique(unlist(asv_sets))

  if (length(all_asvs) == 0L) {
    warning("No ASVs found for param '", param, "'")
    return(data.frame())
  }

  pres_mat           <- vapply(asv_sets,
                                function(s) all_asvs %in% s,
                                logical(length(all_asvs)))
  rownames(pres_mat) <- all_asvs

  n_present <- rowSums(pres_mat)
  n_levels  <- length(lvls)

  origin <- vapply(seq_along(all_asvs), function(i) {
    row <- pres_mat[i, ]
    n   <- sum(row)
    if (n == n_levels)  return("Shared_all")
    if (n == 1L)        return(paste0("Unique_", lvls[which(row)]))
    paste0("Shared_", paste(lvls[which(row)], collapse = "_"))
  }, character(1L))

  result            <- data.frame(
    ASV              = all_asvs,
    n_levels_present = n_present,
    origin           = origin,
    stringsAsFactors = FALSE
  )
  pres_df           <- as.data.frame(pres_mat)
  colnames(pres_df) <- paste0(param, "_", lvls)

  cbind(result, pres_df)
}

# ==============================================================================
# SECTION 8 — DISTANCE MATRICES
# ==============================================================================

#' Binary Jaccard distance
#' @param ps phyloseq
#' @return dist
compute_jaccard <- function(ps) {
  phyloseq::distance(ps, method = "jaccard", binary = TRUE)
}

#' Aitchison distance (CLR-Euclidean; pseudocount handles zeros)
#' @param ps          phyloseq
#' @param pseudocount Added before log transform (default 1)
#' @return dist
compute_aitchison <- function(ps, pseudocount = 1) {
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  otu <- otu + pseudocount
  clr <- t(apply(otu, 1L, function(x) log(x) - mean(log(x))))
  dist(clr, method = "euclidean")
}

#' Compute Jaccard and Aitchison distances with save/load caching
#'
#' @param ps             phyloseq
#' @param cache_dir      Directory to read/write cached dist objects
#' @param n_cores        Cores for parallel computation (ignored if cached)
#' @param force_recompute If TRUE, always recompute even if cache exists
#' @return Named list: Jaccard, Aitchison (both dist objects)
compute_distances_parallel <- function(ps,
                                        cache_dir       = NULL,
                                        n_cores         = 2L,
                                        force_recompute = FALSE) {

  cache_paths <- if (!is.null(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    list(
      Jaccard   = file.path(cache_dir, "dist_Jaccard.rds"),
      Aitchison = file.path(cache_dir, "dist_Aitchison.rds")
    )
  } else {
    list(Jaccard = NULL, Aitchison = NULL)
  }

  dist_fns <- list(
    Jaccard   = function() compute_jaccard(ps),
    Aitchison = function() compute_aitchison(ps)
  )

  res <- vector("list", 2L)
  names(res) <- c("Jaccard", "Aitchison")

  for (nm in c("Jaccard", "Aitchison")) {
    cache_file <- cache_paths[[nm]]
    if (!force_recompute && !is.null(cache_file) && file.exists(cache_file)) {
      message("Loading cached ", nm, " distance from: ", cache_file)
      res[[nm]] <- readRDS(cache_file)
    } else {
      message("Computing ", nm, " distance ...")
      res[[nm]] <- dist_fns[[nm]]()
      if (!is.null(cache_file)) {
        saveRDS(res[[nm]], cache_file)
        message("  Saved cache: ", cache_file)
      }
    }
    message("  ", nm, " matrix: ",
            attr(res[[nm]], "Size"), " x ", attr(res[[nm]], "Size"))
  }

  res
}

# ==============================================================================
# SECTION 9 — PCoA
# ==============================================================================

#' Run PCoA and return a tidy data frame for ggplot2
#'
#' @param dist_obj   dist object
#' @param ps         phyloseq (provides sample_data)
#' @param dist_label Label for plots
#' @return data.frame: PC1–PC4 + all metadata columns; attr "pct"
run_pcoa <- function(dist_obj, ps, dist_label = "Distance") {
  pcoa_res <- ape::pcoa(dist_obj)

  eig <- pcoa_res$values$Rel_corr_eig
  if (is.null(eig) || all(is.na(eig)))
    eig <- pcoa_res$values$Relative_eig
  eig[is.na(eig)] <- 0

  n_axes <- min(4L, sum(eig > 0))
  pct    <- round(eig[seq_len(n_axes)] * 100, 1)

  coords           <- as.data.frame(
    pcoa_res$vectors[, seq_len(n_axes), drop = FALSE])
  colnames(coords) <- paste0("PC", seq_len(n_axes))

  sd_raw <- as(sample_data(ps), "data.frame")
  meta   <- data.frame(
    lapply(sd_raw, function(x)
      if (is.factor(x) || is.ordered(x)) as.character(x) else x),
    stringsAsFactors = FALSE,
    row.names        = rownames(sd_raw)
  )
  common <- intersect(rownames(coords), rownames(meta))

  if (length(common) == 0L)
    stop("No matching sample IDs between PCoA coords and phyloseq.")

  df             <- cbind(coords[common, , drop = FALSE],
                          meta[common,   , drop = FALSE])
  df$dist_method <- dist_label
  attr(df, "pct") <- pct
  df
}

# ==============================================================================
# SECTION 10 — PERMANOVA
# ==============================================================================

#' Run a single PERMANOVA (vegan::adonis2) on a pre-computed distance
#'
#' @param dist_obj     dist object
#' @param ps           phyloseq (sample_data as metadata)
#' @param rhs          RHS formula string
#' @param permutations Number of permutations
#' @return adonis2 result object
run_permanova <- function(dist_obj, ps,
                          rhs          = "sample_type + Filter_Mode + Pool + Omega + Band_Size + Error_Model",
                          permutations = 999L) {

  sd_raw <- as(sample_data(ps), "data.frame")
  meta   <- data.frame(
    lapply(sd_raw, function(x)
      if (is.factor(x) || is.ordered(x)) as.character(x) else x),
    stringsAsFactors = FALSE,
    row.names        = rownames(sd_raw)
  )

  ids     <- labels(dist_obj)
  missing <- setdiff(ids, rownames(meta))
  if (length(missing) > 0L)
    stop(length(missing), " sample(s) in dist not found in sample_data.")

  meta <- meta[ids, , drop = FALSE]

  terms_vec <- trimws(unlist(strsplit(rhs, "\\+")))
  keep <- vapply(terms_vec, function(t) {
    if (!t %in% colnames(meta)) {
      warning("PERMANOVA: '", t, "' not in metadata - skipped."); return(FALSE)
    }
    if (length(unique(na.omit(meta[[t]]))) < 2L) {
      warning("PERMANOVA: '", t, "' is constant - skipped."); return(FALSE)
    }
    TRUE
  }, logical(1L))

  if (!any(keep)) stop("No valid PERMANOVA terms remain.")
  rhs_ok <- paste(terms_vec[keep], collapse = " + ")
  message("  [PERMANOVA] pos_axes ~ ", rhs_ok, "  (n=", length(ids), ")")

  pcoa_res    <- ape::pcoa(dist_obj)
  vec_rows    <- rownames(pcoa_res$vectors)
  missing_vec <- setdiff(ids, vec_rows)
  if (length(missing_vec) > 0L)
    stop("PCoA vectors missing rows for ", length(missing_vec), " sample(s).")

  pos_axes <- pcoa_res$vectors[ids, , drop = FALSE]

  fmla <- stats::as.formula(paste("pos_axes ~", rhs_ok))
  set.seed(42L)
  vegan::adonis2(fmla,
                 data         = meta,
                 method       = "euclidean",
                 permutations = as.integer(permutations),
                 by           = "margin")
}

#' Run PERMANOVA per-parameter (includes Band_Size and Error_Model)
#'
#' For each of the 6 parameters (Filter_Mode, Pool, Omega, Band_Size,
#' Error_Model, sample_type), all other parameters are fixed at reference
#' levels so each biological sample appears exactly ONCE per level of the
#' target parameter.
#'
#' Reference levels (most common / intermediate):
#'   Pool        = "pseudo"
#'   Filter_Mode = "Loose"
#'   Omega       = "O1e-200"
#'   Band_Size   = "32"
#'   Error_Model = "Binned"
#'
#' @param dist_list    Named list: Jaccard, Aitchison (full merged distances)
#' @param ps           Merged phyloseq
#' @param permutations Permutations per test
#' @param ref_pool     Reference Pool level
#' @param ref_filter   Reference Filter_Mode level
#' @param ref_omega    Reference Omega level
#' @param ref_band     Reference Band_Size level
#' @param ref_err      Reference Error_Model level
#' @return Named list of tidy PERMANOVA result data.frames (combined)
run_permanova_by_param <- function(dist_list,
                                    ps,
                                    permutations = 999L,
                                    ref_pool     = "pseudo",
                                    ref_filter   = "Loose",
                                    ref_omega    = "O1e-200",
                                    ref_band     = "32",
                                    ref_err      = "Binned") {

  sd_raw <- as(sample_data(ps), "data.frame")
  meta   <- data.frame(
    lapply(sd_raw, function(x)
      if (is.factor(x) || is.ordered(x)) as.character(x) else x),
    stringsAsFactors = FALSE,
    row.names        = rownames(sd_raw)
  )

  subset_ps_dist <- function(mask, dist_obj) {
    keep_samps <- rownames(meta)[mask]
    ps_sub     <- prune_samples(keep_samps, ps)
    dist_sub   <- as.dist(as.matrix(dist_obj)[keep_samps, keep_samps])
    list(ps = ps_sub, dist = dist_sub, n = length(keep_samps))
  }

  results <- list()

  for (dist_nm in names(dist_list)) {
    dist_obj <- dist_list[[dist_nm]]
    message("\n── PERMANOVA [", dist_nm, "] ──────────────────────────────")

    # 1. Test POOL: fix Filter_Mode, Omega, Band_Size, Error_Model
    mask_pool <- as.character(meta$Filter_Mode)  == ref_filter &
                 as.character(meta$Omega)         == ref_omega  &
                 as.character(meta$Band_Size)     == ref_band   &
                 as.character(meta$Error_Model)   == ref_err
    sub_pool  <- subset_ps_dist(mask_pool, dist_obj)
    message("  Testing Pool (n=", sub_pool$n, ")")
    res_pool  <- run_permanova(sub_pool$dist, sub_pool$ps,
                               rhs = "sample_type + Pool",
                               permutations = permutations)
    results[[paste0(dist_nm, "_Pool")]] <-
      tidy_permanova(res_pool, dist_label = dist_nm) |>
      dplyr::mutate(tested_param = "Pool")

    # 2. Test FILTER_MODE: fix Pool, Omega, Band_Size, Error_Model
    mask_filt <- as.character(meta$Pool)          == ref_pool  &
                 as.character(meta$Omega)         == ref_omega  &
                 as.character(meta$Band_Size)     == ref_band   &
                 as.character(meta$Error_Model)   == ref_err
    sub_filt  <- subset_ps_dist(mask_filt, dist_obj)
    message("  Testing Filter_Mode (n=", sub_filt$n, ")")
    res_filt  <- run_permanova(sub_filt$dist, sub_filt$ps,
                               rhs = "sample_type + Filter_Mode",
                               permutations = permutations)
    results[[paste0(dist_nm, "_Filter_Mode")]] <-
      tidy_permanova(res_filt, dist_label = dist_nm) |>
      dplyr::mutate(tested_param = "Filter_Mode")

    # 3. Test OMEGA: fix Pool, Filter_Mode, Band_Size, Error_Model
    mask_omg  <- as.character(meta$Pool)          == ref_pool   &
                 as.character(meta$Filter_Mode)   == ref_filter &
                 as.character(meta$Band_Size)     == ref_band   &
                 as.character(meta$Error_Model)   == ref_err
    sub_omg   <- subset_ps_dist(mask_omg, dist_obj)
    message("  Testing Omega (n=", sub_omg$n, ")")
    res_omg   <- run_permanova(sub_omg$dist, sub_omg$ps,
                               rhs = "sample_type + Omega",
                               permutations = permutations)
    results[[paste0(dist_nm, "_Omega")]] <-
      tidy_permanova(res_omg, dist_label = dist_nm) |>
      dplyr::mutate(tested_param = "Omega")

    # 4. Test BAND_SIZE: fix Pool, Filter_Mode, Omega, Error_Model
    mask_band <- as.character(meta$Pool)          == ref_pool   &
                 as.character(meta$Filter_Mode)   == ref_filter &
                 as.character(meta$Omega)         == ref_omega  &
                 as.character(meta$Error_Model)   == ref_err
    sub_band  <- subset_ps_dist(mask_band, dist_obj)
    message("  Testing Band_Size (n=", sub_band$n, ")")
    res_band  <- run_permanova(sub_band$dist, sub_band$ps,
                               rhs = "sample_type + Band_Size",
                               permutations = permutations)
    results[[paste0(dist_nm, "_Band_Size")]] <-
      tidy_permanova(res_band, dist_label = dist_nm) |>
      dplyr::mutate(tested_param = "Band_Size")

    # 5. Test ERROR_MODEL: fix Pool, Filter_Mode, Omega, Band_Size
    mask_err  <- as.character(meta$Pool)          == ref_pool   &
                 as.character(meta$Filter_Mode)   == ref_filter &
                 as.character(meta$Omega)         == ref_omega  &
                 as.character(meta$Band_Size)     == ref_band
    sub_err   <- subset_ps_dist(mask_err, dist_obj)
    message("  Testing Error_Model (n=", sub_err$n, ")")
    res_err   <- run_permanova(sub_err$dist, sub_err$ps,
                               rhs = "sample_type + Error_Model",
                               permutations = permutations)
    results[[paste0(dist_nm, "_Error_Model")]] <-
      tidy_permanova(res_err, dist_label = dist_nm) |>
      dplyr::mutate(tested_param = "Error_Model")

    # 6. Test SAMPLE_TYPE: single reference case (all 5 params fixed)
    mask_type <- as.character(meta$Pool)          == ref_pool   &
                 as.character(meta$Filter_Mode)   == ref_filter &
                 as.character(meta$Omega)         == ref_omega  &
                 as.character(meta$Band_Size)     == ref_band   &
                 as.character(meta$Error_Model)   == ref_err
    sub_type  <- subset_ps_dist(mask_type, dist_obj)
    message("  Testing sample_type (n=", sub_type$n, ", reference case)")
    res_type  <- run_permanova(sub_type$dist, sub_type$ps,
                               rhs = "sample_type",
                               permutations = permutations)
    results[[paste0(dist_nm, "_sample_type")]] <-
      tidy_permanova(res_type, dist_label = dist_nm) |>
      dplyr::mutate(tested_param = "sample_type")
  }

  dplyr::bind_rows(results)
}

#' Tidy an adonis2 result into a plottable data frame
#'
#' @param adonis_res adonis2 result
#' @param dist_label Distance label
#' @return data.frame: Term, Df, R2, F_stat, p_value, Significance, label
tidy_permanova <- function(adonis_res, dist_label = "Distance") {
  df      <- as.data.frame(adonis_res)
  df$Term <- rownames(df)

  df <- df[!is.na(df$F) &
             !df$Term %in% c("Residual", "Total"), , drop = FALSE]

  if (nrow(df) == 0L) {
    warning("tidy_permanova: no estimable terms found.")
    return(data.frame())
  }

  names(df)[names(df) == "Pr(>F)"] <- "p_value"
  names(df)[names(df) == "F"]      <- "F_stat"

  df$Significance <- dplyr::case_when(
    df$p_value < 0.001 ~ "***",
    df$p_value < 0.01  ~ "**",
    df$p_value < 0.05  ~ "*",
    df$p_value < 0.1   ~ ".",
    TRUE               ~ "ns"
  )
  df$label       <- paste0("R\u00b2=", round(df$R2, 3),
                            "  F=",    round(df$F_stat, 2),
                            "  p",     df$Significance)
  df$dist_method <- dist_label
  rownames(df)   <- NULL
  df
}

# ==============================================================================
# SECTION 11 — COLOUR PALETTES & BASE THEME
# ==============================================================================

#' Consistent colour and shape palettes across all plots
#'
#' Includes Band_Size and Error_Model palettes.
get_palettes <- function() {
  list(
    sample_type = c(
      Individual = "#2166AC",
      Positive   = "#1B7837",
      Negative   = "#D6604D"
    ),
    sample_type_shape = c(
      Individual = 16L,   # filled circle
      Positive   = 17L,   # filled triangle
      Negative   = 15L    # filled square
    ),
    Pool = c(
      "FALSE"  = "#C6DBEF",
      "pseudo" = "#4393C3",
      "TRUE"   = "#08306B"
    ),
    Filter_Mode = c(
      Loose  = "#FDAE6B",
      Strict = "#8C2D04"
    ),
    Omega = c(
      "O1e-40"  = "#D9F0A3",
      "O1e-70"  = "#78C679",
      "O1e-100" = "#238443",
      "O1e-200" = "#004529"
    ),
    Band_Size = c(
      "32" = "#FEE08B",
      "64" = "#D73027"
    ),
    Error_Model = c(
      "Binned" = "#762A83",
      "PacBio" = "#1A9850"
    )
  )
}

#' Base ggplot2 theme for all community analysis plots
theme_community <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold",
                                      size = base_size + 2, hjust = 0.5),
      plot.subtitle    = element_text(size = base_size - 1, hjust = 0.5,
                                      colour = "grey40"),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 2),
      axis.title       = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank()
    )
}

# ==============================================================================
# SECTION 12 — ASV SHARING PLOTS
# ==============================================================================

#' Assign intuitive colours to every unique origin label
#' @noRd
.origin_palette <- function(origin_levels) {
  base_cols <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
                 "#FF7F00","#A65628","#F781BF","#999999",
                 "#66C2A5","#FC8D62","#8DA0CB","#E78AC3")
  pal <- setNames(
    colorRampPalette(base_cols)(length(origin_levels)),
    origin_levels
  )
  for (i in seq_along(origin_levels)) {
    lv <- origin_levels[i]
    if (lv == "Shared_all")          pal[i] <- "#4DAF4A"
    else if (grepl("^Unique_", lv))  pal[i] <- "#E41A1C"
    else if (grepl("^Shared_", lv))  pal[i] <- "#984EA3"
  }
  pal
}

#' Bar chart of ASV sharing for ONE parameter
#'
#' @param sharing_df data.frame from summarise_asv_sharing()
#' @param param      Parameter name (for title/subtitle)
#' @return ggplot object
plot_asv_sharing_bar <- function(sharing_df, param) {

  if (is.null(sharing_df) || nrow(sharing_df) == 0L)
    return(ggplot() +
             labs(title = paste("No data for", param)) +
             theme_community())

  pres_cols   <- grep(paste0("^", param, "_"), colnames(sharing_df), value = TRUE)
  level_names <- gsub(paste0("^", param, "_"), "", pres_cols)

  sharing_filt <- sharing_df[
    sharing_df$origin == "Shared_all" |
      grepl("^Unique_", sharing_df$origin), ]

  if (nrow(sharing_filt) == 0L)
    return(ggplot() +
             labs(title = paste("No Shared_all / Unique entries for", param)) +
             theme_community())

  counts <- sharing_filt |>
    dplyr::count(origin, name = "n") |>
    dplyr::mutate(
      pct   = n / nrow(sharing_df) * 100,
      label = paste0(n, "  (", round(pct, 1), "%)")
    ) |>
    dplyr::arrange(dplyr::desc(n))

  unique_cats <- grep("^Unique_", counts$origin, value = TRUE)
  n_u         <- length(unique_cats)
  unique_cols <- if (n_u > 0L)
    setNames(RColorBrewer::brewer.pal(max(3L, n_u), "Set1")[seq_len(n_u)],
             unique_cats)
  else character(0L)

  pal <- c("Shared_all" = "#4DAF4A", unique_cols)[counts$origin]

  ggplot(counts,
         aes(x = reorder(origin, n), y = n, fill = origin)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.4) +
    geom_text(aes(label = label), hjust = -0.08,
              size = 3.3, fontface = "bold") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
    coord_flip() +
    labs(
      title    = paste0("ASV origin \u2013 ", param),
      subtitle = paste0("Total ASVs: ", nrow(sharing_df),
                        "  |  levels: ", paste(level_names, collapse = " / ")),
      x = NULL,
      y = "Number of ASVs"
    ) +
    theme_community() +
    theme(axis.text.y        = element_text(size = 10, face = "bold"),
          panel.grid.major.y = element_blank())
}

#' Combined ASV sharing figure for all 5 parameters
#'
#' Displayed as two rows: top row = Filter_Mode / Pool / Omega;
#' bottom row = Band_Size / Error_Model.
#'
#' @param sharing_list Named list from summarise_asv_sharing()
#' @return patchwork object
plot_asv_sharing_all <- function(sharing_list) {
  # Top row: original 3 parameters
  top_params  <- intersect(c("Filter_Mode","Pool","Omega"), names(sharing_list))
  # Bottom row: additional parameters
  bot_params  <- intersect(c("Band_Size","Error_Model"),    names(sharing_list))

  make_row <- function(params) {
    if (length(params) == 0L) return(NULL)
    plots <- lapply(params, function(p) plot_asv_sharing_bar(sharing_list[[p]], p))
    wrap_plots(plots, nrow = 1L)
  }

  top_row <- make_row(top_params)
  bot_row <- make_row(bot_params)

  combined <- if (!is.null(top_row) && !is.null(bot_row)) {
    top_row / bot_row
  } else if (!is.null(top_row)) {
    top_row
  } else {
    bot_row
  }

  combined +
    plot_annotation(
      title    = "ASV Sharing Across DADA2 Parameter Levels",
      subtitle = "Shared across all levels (green) vs. unique to each level (red tones)",
      theme = theme(
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 10)
      )
    )
}

# ==============================================================================
# SECTION 13 — PCoA PLOTS
# ==============================================================================

#' Single PCoA scatter: colour = parameter, shape = sample type
#'
#' @param pcoa_df   data.frame from run_pcoa() (must carry "pct" attribute)
#' @param color_by  Variable for colour
#' @param perm_df   Tidy PERMANOVA df or NULL
#' @param pal       Named colour vector (NULL = auto from get_palettes)
#' @param shapes    Named integer shape vector
#' @param title     Plot title (NULL = auto)
#' @param pc_x      PC index for x-axis (default 1)
#' @param pc_y      PC index for y-axis (default 2)
#' @return ggplot object
plot_pcoa <- function(pcoa_df,
                      color_by = "Pool",
                      perm_df  = NULL,
                      pal      = NULL,
                      shapes   = NULL,
                      title    = NULL,
                      pc_x = 1L, pc_y = 2L) {

  pct   <- attr(pcoa_df, "pct")
  if (is.null(pct)) pct <- rep(NA_real_, 4L)

  fmt <- function(i) if (!is.na(pct[i])) paste0(" [", pct[i], "%]") else ""
  x_col <- paste0("PC", pc_x);  y_col <- paste0("PC", pc_y)
  xlab  <- paste0("PC", pc_x, fmt(pc_x))
  ylab  <- paste0("PC", pc_y, fmt(pc_y))

  pals <- get_palettes()
  if (is.null(pal))    pal    <- pals[[color_by]]
  if (is.null(shapes)) shapes <- pals$sample_type_shape
  if (is.null(title))
    title <- paste0("PCoA \u2013 ", unique(pcoa_df$dist_method),
                    "  |  Colour: ", color_by)

  pcoa_df[[color_by]] <- factor(pcoa_df[[color_by]])
  pcoa_df$sample_type <- factor(pcoa_df$sample_type,
                                 levels = c("Individual","Positive","Negative"))

  grp_n        <- table(pcoa_df[[color_by]])
  grps_ellipse <- names(grp_n)[grp_n >= 3L]

  p <- ggplot(pcoa_df,
              aes(x      = .data[[x_col]],
                  y      = .data[[y_col]],
                  colour = .data[[color_by]],
                  shape  = sample_type)) +
    geom_point(size = 2.6, alpha = 0.82, stroke = 0.35) +
    scale_colour_manual(values = pal,    name = color_by,    drop = FALSE) +
    scale_shape_manual(values  = shapes, name = "Sample type", drop = FALSE) +
    labs(title = title, x = xlab, y = ylab) +
    theme_community() +
    guides(
      colour = guide_legend(override.aes = list(size = 4, alpha = 1)),
      shape  = guide_legend(override.aes = list(size = 4, alpha = 1))
    )

  if (length(grps_ellipse) > 0L) {
    ell_df <- pcoa_df[as.character(pcoa_df[[color_by]]) %in% grps_ellipse, ]
    p <- p + stat_ellipse(
      data        = ell_df,
      aes(group   = .data[[color_by]], colour = .data[[color_by]]),
      type        = "t", level = 0.8,
      linetype    = 2, linewidth = 0.5, alpha = 0.6,
      show.legend = FALSE
    )
  }

  if (!is.null(perm_df) && nrow(perm_df) > 0L) {
    if ("tested_param" %in% colnames(perm_df)) {
      row <- perm_df[perm_df$Term == color_by &
                       perm_df$tested_param == color_by, , drop = FALSE]
      if (nrow(row) == 0L)
        row <- perm_df[perm_df$Term == color_by, , drop = FALSE]
    } else {
      row <- perm_df[perm_df$Term == color_by, , drop = FALSE]
    }
    if (nrow(row) >= 1L) {
      row <- row[1L, ]
      ann <- paste0("PERMANOVA\nR\u00b2=", round(row$R2, 3),
                    "  F=",               round(row$F_stat, 2),
                    "\np",                row$Significance)
      p <- p + annotate("label",
                         x = Inf, y = Inf, label = ann,
                         hjust = 1.05, vjust = 1.1,
                         size = 2.7, fill = "grey96",
                         label.size = 0.25, colour = "grey25")
    }
  }
  p
}

#' 6-panel PCoA: sample_type / Pool / Filter_Mode / Omega / Band_Size / Error_Model
#'
#' Displayed as two rows of three panels each.
#'
#' @param pcoa_df data.frame from run_pcoa()
#' @param perm_df Tidy PERMANOVA df or NULL
#' @param params  Panels to produce (default all 6)
#' @return patchwork object
plot_pcoa_multi <- function(pcoa_df,
                             perm_df = NULL,
                             params  = c("sample_type","Pool",
                                         "Filter_Mode","Omega",
                                         "Band_Size","Error_Model")) {
  pals  <- get_palettes()
  plots <- lapply(params, function(p) {
    pal_p <- if (!is.null(pals[[p]])) pals[[p]] else NULL
    plot_pcoa(pcoa_df, color_by = p, perm_df = perm_df, pal = pal_p)
  })
  # Two rows of three
  wrap_plots(plots, nrow = 2L) +
    plot_annotation(
      title = paste0("PCoA \u2013 ", unique(pcoa_df$dist_method), " distance"),
      theme = theme(plot.title = element_text(face = "bold",
                                              hjust = 0.5, size = 14))
    )
}

# ==============================================================================
# SECTION 14 — DENDROGRAM + ANNOTATION HEATMAP
# ==============================================================================

#' Hierarchical clustering dendrogram with coloured annotation strips
#'
#' Includes Band_Size and Error_Model annotation columns.
#'
#' @param dist_obj       dist object
#' @param ps             phyloseq (sample_data)
#' @param dist_label     Title string
#' @param perm_df        Tidy PERMANOVA df or NULL
#' @param cluster_method hclust linkage method (default "ward.D2")
#' @return ComplexHeatmap Heatmap object – render with draw()
plot_dendrogram_heatmap <- function(dist_obj, ps,
                                     dist_label     = "Jaccard",
                                     perm_df        = NULL,
                                     cluster_method = "ward.D2") {
  pals <- get_palettes()

  sd_raw <- as(sample_data(ps), "data.frame")
  meta   <- data.frame(
    lapply(sd_raw, function(x)
      if (is.factor(x) || is.ordered(x)) as.character(x) else x),
    stringsAsFactors = FALSE,
    row.names        = rownames(sd_raw)
  )
  ids  <- labels(dist_obj)
  meta <- meta[ids, , drop = FALSE]

  hc <- hclust(dist_obj, method = cluster_method)

  ann_df   <- meta[, c("sample_type","Pool","Filter_Mode",
                        "Omega","Band_Size","Error_Model"), drop = FALSE]
  ann_df[] <- lapply(ann_df, as.character)

  ha <- HeatmapAnnotation(
    "Sample\ntype" = ann_df$sample_type,
    "Pool"         = ann_df$Pool,
    "Filter"       = ann_df$Filter_Mode,
    "Omega"        = ann_df$Omega,
    "Band"         = ann_df$Band_Size,
    "Err model"    = ann_df$Error_Model,
    col = list(
      "Sample\ntype" = pals$sample_type,
      "Pool"         = pals$Pool,
      "Filter"       = pals$Filter_Mode,
      "Omega"        = pals$Omega,
      "Band"         = pals$Band_Size,
      "Err model"    = pals$Error_Model
    ),
    which                = "row",
    show_annotation_name = FALSE,
    simple_anno_size     = unit(6, "mm"),
    gap                  = unit(1.5, "mm"),
    show_legend          = TRUE,
    annotation_legend_param = list(
      "Sample\ntype" = list(title = "Sample type",
                            title_gp  = gpar(fontsize = 9, fontface = "bold"),
                            labels_gp = gpar(fontsize = 8)),
      "Pool"         = list(title = "Pooling",
                            title_gp  = gpar(fontsize = 9, fontface = "bold"),
                            labels_gp = gpar(fontsize = 8)),
      "Filter"       = list(title = "Filter mode",
                            title_gp  = gpar(fontsize = 9, fontface = "bold"),
                            labels_gp = gpar(fontsize = 8)),
      "Omega"        = list(title = "Omega",
                            title_gp  = gpar(fontsize = 9, fontface = "bold"),
                            labels_gp = gpar(fontsize = 8)),
      "Band"         = list(title = "Band size",
                            title_gp  = gpar(fontsize = 9, fontface = "bold"),
                            labels_gp = gpar(fontsize = 8)),
      "Err model"    = list(title = "Error model",
                            title_gp  = gpar(fontsize = 9, fontface = "bold"),
                            labels_gp = gpar(fontsize = 8))
    )
  )

  # PERMANOVA caption
  perm_caption <- ""
  if (!is.null(perm_df) && nrow(perm_df) > 0L) {
    if ("tested_param" %in% colnames(perm_df)) {
      perm_show <- perm_df[perm_df$Term == perm_df$tested_param |
                             perm_df$Term == "sample_type", , drop = FALSE]
      perm_show <- perm_show[!duplicated(perm_show$Term), , drop = FALSE]
    } else {
      perm_show <- perm_df
    }
    perm_caption <- paste0(
      "PERMANOVA (", unique(perm_df$dist_method)[1L],
      ", margin, 999 perm):  ",
      paste(
        paste0(perm_show$Term,
               "  R\u00b2=", round(perm_show$R2, 3),
               "  F=",       round(perm_show$F_stat, 2),
               "  p",        perm_show$Significance),
        collapse = "   |   "
      )
    )
  }

  col_title <- paste0("Hierarchical Clustering \u2013 ",
                       dist_label, " distance (", cluster_method, ")\n",
                       "Annotation (left \u2192 right): ",
                       "Sample type | Pool | Filter | Omega | Band | Err model")

  mat_dummy <- matrix(nrow = length(ids), ncol = 0L,
                      dimnames = list(ids, character(0L)))

  ht <- Heatmap(
    matrix               = mat_dummy,
    cluster_rows         = hc,
    cluster_columns      = FALSE,
    show_row_names       = FALSE,
    show_column_names    = FALSE,
    right_annotation     = ha,
    column_title         = col_title,
    column_title_gp      = gpar(fontsize = 11, fontface = "bold"),
    row_dend_width       = unit(4, "cm"),
    row_dend_gp          = gpar(col = "grey20", lwd = 0.8),
    border               = FALSE,
    heatmap_legend_param = list(
      title_gp  = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    )
  )
  attr(ht, "perm_caption") <- perm_caption
  ht
}

# ==============================================================================
# SECTION 15 — STANDALONE PERMANOVA PLOT
# ==============================================================================

#' Lollipop plot of PERMANOVA R² for Jaccard and Aitchison
#'
#' Shows all 6 parameters.
#'
#' @param perm_tidy Named list of tidy PERMANOVA dfs (combined)
#' @return ggplot object
plot_permanova_results <- function(perm_tidy) {
  if (nrow(perm_tidy) == 0L) stop("plot_permanova_results: no data to plot.")

  # All 6 focal parameters (each tested in its own focused subsample)
  focal_params <- c("Pool","Filter_Mode","Omega","Band_Size","Error_Model","sample_type")

  plot_df <- perm_tidy |>
    dplyr::filter(
      (Term == "Pool"         & tested_param == "Pool")         |
      (Term == "Filter_Mode"  & tested_param == "Filter_Mode")  |
      (Term == "Omega"        & tested_param == "Omega")        |
      (Term == "Band_Size"    & tested_param == "Band_Size")    |
      (Term == "Error_Model"  & tested_param == "Error_Model")  |
      (Term == "sample_type"  & tested_param == "sample_type")
    ) |>
    dplyr::group_by(Term, dist_method) |>
    dplyr::slice_max(order_by = R2, n = 1L, with_ties = FALSE) |>
    dplyr::ungroup()

  if (nrow(plot_df) == 0L) {
    warning("plot_permanova_results: no matching rows after filtering. ",
            "Check that tested_param values match Term names.")
    return(ggplot() + labs(title = "No PERMANOVA data") + theme_community())
  }

  term_order <- plot_df |>
    dplyr::group_by(Term) |>
    dplyr::summarise(mean_R2 = mean(R2, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(mean_R2) |>
    dplyr::pull(Term)

  plot_df$Term <- factor(plot_df$Term, levels = term_order)

  plot_df$p_label <- paste0(
    formatC(plot_df$p_value, format = "g", digits = 2),
    " ", plot_df$Significance
  )

  dist_colours <- c(Jaccard = "#2166AC", Aitchison = "#B2182B")

  p_main <- ggplot(plot_df,
                   aes(x = Term, y = R2, colour = dist_method)) +
    geom_segment(aes(xend = Term, yend = 0),
                 position  = position_dodge(width = 0.6),
                 linewidth = 1.1, alpha = 0.6) +
    geom_point(aes(size = F_stat),
               position = position_dodge(width = 0.6),
               alpha    = 0.95) +
    geom_text(aes(label = paste0("R\u00b2=", round(R2, 3))),
              position  = position_dodge(width = 0.6),
              hjust = -0.18, size = 3.0, fontface = "bold",
              show.legend = FALSE) +
    scale_colour_manual(values = dist_colours, name = "Distance") +
    scale_size_continuous(
      name  = "F-statistic", range = c(3, 10),
      guide = guide_legend(override.aes = list(colour = "grey40"))
    ) +
    scale_y_continuous(
      limits = c(0, max(plot_df$R2, na.rm = TRUE) * 1.5),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_flip() +
    labs(
      title    = "PERMANOVA Results",
      subtitle = paste0("R\u00b2 = proportion of variance explained  |  ",
                        "by='margin'  |  999 permutations\n",
                        "Each parameter tested in its own focused subsample"),
      x = NULL,
      y = "R\u00b2"
    ) +
    theme_community() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.4),
      plot.margin        = margin(5, 5, 5, 5)
    )

  p_tab <- ggplot(plot_df,
                  aes(x = Term, y = dist_method, label = p_label,
                      colour = dist_method)) +
    geom_text(size = 3.0, fontface = "bold") +
    scale_colour_manual(values = dist_colours, guide = "none") +
    scale_x_discrete(limits = levels(plot_df$Term)) +
    coord_flip() +
    labs(x = NULL, y = NULL, title = "p-value") +
    theme_void() +
    theme(
      plot.title       = element_text(size = 10, face = "bold",
                                      hjust = 0.5, vjust = 1),
      axis.text.y      = element_blank(),
      panel.background = element_rect(fill = "grey97", colour = NA),
      plot.margin      = margin(5, 8, 5, 0)
    )

  p_main + p_tab + patchwork::plot_layout(widths = c(4, 1), guides = "collect")
}

# ==============================================================================
# SECTION 16 — SAVE HELPERS
# ==============================================================================

#' Save a ggplot / patchwork object
#'
#' @param p        ggplot or patchwork
#' @param filename Full path with extension (.png / .pdf / .svg)
#' @param width    Width in inches
#' @param height   Height in inches
#' @param dpi      Resolution (default 300)
save_plot <- function(p, filename, width = 14, height = 8, dpi = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot = p,
                  width = width, height = height, dpi = dpi, bg = "white")
  message("  Saved: ", filename)
  invisible(filename)
}

#' Save a ComplexHeatmap object to PNG or PDF
#'
#' @param ht       ComplexHeatmap / HeatmapList
#' @param filename Full path (.png or .pdf)
#' @param width    Width in inches
#' @param height   Height in inches
#' @param res      PNG resolution (default 300)
save_heatmap <- function(ht, filename, width = 18, height = 13, res = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  caption     <- attr(ht, "perm_caption")
  has_caption <- !is.null(caption) && nchar(trimws(caption)) > 0L
  out_height  <- if (has_caption) height + 0.7 else height

  ext <- tolower(tools::file_ext(filename))
  if (ext == "pdf") {
    grDevices::pdf(filename, width = width, height = out_height)
  } else {
    grDevices::png(filename, width = width, height = out_height,
                   units = "in", res = res, bg = "white")
  }
  on.exit(grDevices::dev.off())

  if (has_caption) {
    hm_frac <- height / out_height
    grid::pushViewport(
      grid::viewport(y = 1 - hm_frac / 2, height = hm_frac, just = "centre")
    )
    ComplexHeatmap::draw(ht,
                          merge_legend            = TRUE,
                          heatmap_legend_side     = "right",
                          annotation_legend_side  = "right",
                          newpage                 = FALSE)
    grid::popViewport()
    grid::grid.text(
      label = caption,
      x     = grid::unit(0.5, "npc"),
      y     = grid::unit(0.3, "inches"),
      just  = "centre",
      gp    = grid::gpar(fontsize = 7.5, col = "grey25")
    )
  } else {
    ComplexHeatmap::draw(ht,
                          merge_legend            = TRUE,
                          heatmap_legend_side     = "right",
                          annotation_legend_side  = "right")
  }
  message("  Saved: ", filename)
  invisible(filename)
}
