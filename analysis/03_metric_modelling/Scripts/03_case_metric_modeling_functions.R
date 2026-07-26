# ==============================================================================
# 03_case_metric_modeling_functions.R
#
# Project : fliC Amplicon Sequencing for High-Resolution Profiling of E. coli
#           in the Infant Gut Microbiome
# Author  : Andrés Catalán Tatay
# Part    : Analysis — Step 03 (companion functions)
# Purpose : All helper functions for 03_case_metric_modeling.Rmd.
#           Covers FASTA parsing, IUPAC-aware sequence matching, case loading
#           with rarefaction, metric computation, regression modelling,
#           marginal means, composite score ranking, and all associated plots.
#
# Data structure assumed:
#   Each case is a list with:
#     $st        - seqtab BEFORE chimera removal (samples x ASVs matrix)
#     $st_nochim - seqtab AFTER chimera removal  (samples x ASVs matrix)
#     $bim       - named logical vector: ASV sequence -> TRUE if chimeric
#     $track     - data.frame: input, filtered, denoised, nonchim
#     $physeq    - phyloseq object (not used in this script)
#     $plot_error- ggplot object (not used in this script)
#
# Mock samples are identified by "Mock_mix" in their rownames.
# Theoretical mock reference is provided as a FASTA file with IUPAC codes.
#
# Inputs  : data/mock/mock_trimmed_4.fasta — theoretical mock reference
#            data/cases/                   — 96 RDS files (DADA2 outputs)
# Outputs : Called by 03_case_metric_modeling.Rmd
#
# Usage   : source("03_case_metric_modeling_functions.R")
# ==============================================================================


# ==============================================================================
# SECTION 1 — FASTA PARSING
# ==============================================================================

#' Parse a FASTA file into a named character vector.
#'
#' @param fasta_path Path to the .fasta file.
#' @return Named character vector: names = sequence IDs, values = sequences
#'         (uppercase, no spaces).
parse_fasta <- function(fasta_path) {
  stopifnot(file.exists(fasta_path))
  lines <- readLines(fasta_path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[lines != ""]  # drop blank lines

  seqs        <- list()
  current_id  <- NULL
  current_seq <- character(0)

  for (line in lines) {
    if (startsWith(line, ">")) {
      # Save the previous entry
      if (!is.null(current_id)) {
        seqs[[current_id]] <- paste(current_seq, collapse = "")
      }
      current_id  <- sub("^>\\s*", "", line)
      # Keep only the first word as the ID (ignore description after first space)
      current_id  <- strsplit(current_id, "\\s+")[[1]][1]
      current_seq <- character(0)
    } else {
      current_seq <- c(current_seq, toupper(line))
    }
  }
  # Save the last entry
  if (!is.null(current_id)) {
    seqs[[current_id]] <- paste(current_seq, collapse = "")
  }

  unlist(seqs)
}


# ==============================================================================
# SECTION 2 — IUPAC-AWARE SEQUENCE MATCHING
# ==============================================================================

# IUPAC ambiguity code lookup table.
# Each ambiguity code maps to the set of nucleotides it represents.
IUPAC_BASES <- list(
  A = "A", C = "C", G = "G", T = "T",
  R = c("A","G"),        # puRine
  Y = c("C","T"),        # pYrimidine
  S = c("G","C"),        # Strong
  W = c("A","T"),        # Weak
  K = c("G","T"),        # Keto
  M = c("A","C"),        # aMino
  B = c("C","G","T"),    # not A
  D = c("A","G","T"),    # not C
  H = c("A","C","T"),    # not G
  V = c("A","C","G"),    # not T
  N = c("A","C","G","T") # aNy
)

#' Test whether a single observed base matches an IUPAC reference base.
#'
#' @param obs  Single character: observed nucleotide (must be A/C/G/T).
#' @param ref  Single character: reference nucleotide (may be IUPAC code).
#' @return Logical.
iupac_base_match <- function(obs, ref) {
  ref_bases <- IUPAC_BASES[[ref]]
  if (is.null(ref_bases)) {
    # Unknown code: treat as no-match
    return(FALSE)
  }
  obs %in% ref_bases
}

#' Test whether an observed ASV sequence matches a reference sequence,
#' allowing for IUPAC ambiguity codes in the reference.
#'
#' Both sequences must be the SAME LENGTH for an exact (full-length) match.
#'
#' @param observed   Character: observed ASV sequence (uppercase, A/C/G/T only).
#' @param reference  Character: reference sequence (uppercase, may contain IUPAC).
#' @return Logical: TRUE if every position matches under IUPAC rules.
iupac_seq_match <- function(observed, reference) {
  if (nchar(observed) != nchar(reference)) return(FALSE)

  obs_chars <- strsplit(observed,  "")[[1]]
  ref_chars <- strsplit(reference, "")[[1]]

  all(mapply(iupac_base_match, obs_chars, ref_chars))
}

#' For each ASV in a character vector, determine whether it matches ANY
#' sequence in the theoretical mock reference (IUPAC-aware).
#'
#' @param asv_seqs   Character vector of ASV sequences.
#' @param mock_seqs  Named character vector of reference sequences (from parse_fasta).
#' @return Named logical vector (same length as asv_seqs):
#'         TRUE  = ASV matches at least one mock reference sequence.
#'         FALSE = no match.
match_asvs_to_mock <- function(asv_seqs, mock_seqs) {
  # For speed: pre-group mock sequences by length; only compare equal-length pairs.
  mock_by_len <- split(mock_seqs, nchar(mock_seqs))

  result <- setNames(logical(length(asv_seqs)), asv_seqs)

  for (i in seq_along(asv_seqs)) {
    asv     <- asv_seqs[i]
    asv_len <- nchar(asv)
    refs    <- mock_by_len[[as.character(asv_len)]]

    if (is.null(refs)) {
      result[i] <- FALSE
      next
    }

    # Check each reference of the same length
    matched <- any(sapply(refs, function(ref) iupac_seq_match(asv, ref)))
    result[i] <- matched
  }

  result
}


# ==============================================================================
# SECTION 3 — CASE FILE LOADING & PARAMETER PARSING
# ==============================================================================

#' Parse DADA2 parameters from a case filename.
#'
#' Expected filename format:
#'   case_001__F-Loose__Err-Binned__Pool-FALSE__O1e-200__B32.rds
#'   case_NNN__F-<filter>__Err-<error>__Pool-<pool>__O<omega>__B<band>.rds
#'
#' @param filename  Character: basename of the file (with or without .rds).
#' @return Named list with fields: case_id, filter_mode, error_model,
#'         pool_mode, omega_raw, omega_num, band_size.
parse_case_filename <- function(filename) {
  # Strip directory and extension
  base <- basename(filename)
  base <- sub("\\.rds$", "", base, ignore.case = TRUE)

  # Split on double underscore
  parts <- strsplit(base, "__")[[1]]

  # Initialise with NA so we return something even if parsing fails
  result <- list(
    case_id    = NA_character_,
    filter_mode= NA_character_,
    error_model= NA_character_,
    pool_mode  = NA_character_,
    omega_raw  = NA_character_,
    omega_num  = NA_real_,
    band_size  = NA_integer_
  )

  for (p in parts) {
    # Case ID: "case_NNN"
    if (grepl("^case_\\d+$", p)) {
      result$case_id <- p
    }
    # Filter mode: "F-<value>"
    else if (grepl("^F-", p)) {
      result$filter_mode <- sub("^F-", "", p)
    }
    # Error model: "Err-<value>"
    else if (grepl("^Err-", p)) {
      result$error_model <- sub("^Err-", "", p)
    }
    # Pooling mode: "Pool-<value>"
    else if (grepl("^Pool-", p)) {
      result$pool_mode <- sub("^Pool-", "", p)
    }
    # Omega: "O<value>" e.g. "O1e-200"
    else if (grepl("^O[0-9]", p)) {
      raw <- sub("^O", "", p)
      result$omega_raw <- raw
      result$omega_num <- as.numeric(raw)   # converts "1e-200" -> 1e-200
    }
    # Band size: "B<integer>"
    else if (grepl("^B\\d+$", p)) {
      result$band_size <- as.integer(sub("^B", "", p))
    }
  }

  result
}

#' Safely load a single case RDS file.
#'
#' @param path  Full path to the .rds file.
#' @return The loaded list object, or NULL if loading fails.
load_case <- function(path) {
  tryCatch(
    readRDS(path),
    error = function(e) {
      message("  [WARNING] Could not load: ", path, "\n  Reason: ", conditionMessage(e))
      NULL
    }
  )
}

#' Load all case RDS files from a directory.
#'
#' Immediately after loading, each case is rarefied to 1500 reads per sample
#' using vegan::rrarefy(). Samples below the threshold are kept unrarefied.
#'
#' @param cases_dir  Path to the directory containing case .rds files.
#' @param pattern    Regex pattern to select files (default: "^case_.*\\.rds$").
#' @return Named list of case objects (names = filenames without extension).
#'         Files that failed to load are excluded.
load_all_cases <- function(cases_dir, pattern = "^case_.*\\.rds$") {
  stopifnot(dir.exists(cases_dir))

  rds_files <- list.files(cases_dir, pattern = pattern,
                           full.names = TRUE, ignore.case = TRUE)

  if (length(rds_files) == 0) {
    stop("No case files found in: ", cases_dir)
  }

  message("Found ", length(rds_files), " case files. Loading...")

  cases <- vector("list", length(rds_files))
  names(cases) <- tools::file_path_sans_ext(basename(rds_files))

  for (i in seq_along(rds_files)) {
    message("  Loading (", i, "/", length(rds_files), "): ", basename(rds_files[i]))
    case_obj <- load_case(rds_files[i])

    # ── Rarefaction ───────────────────────────────────────────────────────────
    # Immediately after loading, before any analysis, rarefy $st_nochim to
    # 1500 reads per sample for ALL cases (Filtering = Loose and Strict).
    #
    # Approach: vegan::rrarefy() subsamples each row (sample) of the
    # samples-x-ASVs matrix independently to exactly 1500 reads without
    # replacement. Samples with fewer than 1500 reads are kept as-is
    # (not rarefied and not discarded) so no data are lost.
    # ASV columns that become all-zero after rarefaction are removed.
    if (!is.null(case_obj)) {
      case_name_i <- tools::file_path_sans_ext(basename(rds_files[i]))
      if (!is.null(case_obj$st_nochim) &&
          is.matrix(case_obj$st_nochim) &&
          nrow(case_obj$st_nochim) > 0) {

        st_raw        <- case_obj$st_nochim
        sample_depths <- rowSums(st_raw)

        # Identify samples below and above the rarefaction threshold
        low_idx  <- sample_depths < 1500
        high_idx <- sample_depths >= 1500

        low_samples <- names(sample_depths)[low_idx]
        if (length(low_samples) > 0) {
          message("  [RAREFACTION WARNING] Case '", case_name_i,
                  "': the following sample(s) have fewer than 1500 reads",
                  " and will be kept as-is (not rarefied): ",
                  paste(low_samples, collapse = ", "))
        }

        if (!any(high_idx)) {
          message("  [RAREFACTION WARNING] Case '", case_name_i,
                  "': no samples have >= 1500 reads -- $st_nochim left entirely unchanged.")
        } else {
          message("  [RAREFACTION] Case '", case_name_i,
                  "' -- rarefying ", sum(high_idx),
                  " sample(s) to 1500 reads; ",
                  sum(low_idx), " sample(s) below threshold kept as-is.")

          # Rarefy only the samples that are deep enough
          set.seed(42)   # reproducible subsampling
          st_rarefied_high <- vegan::rrarefy(
            st_raw[high_idx, , drop = FALSE], sample = 1500
          )

          # Recombine: rarefied rows on top, low-depth rows unchanged,
          # then restore original row order
          st_combined <- rbind(st_rarefied_high,
                               st_raw[low_idx,  , drop = FALSE])
          st_combined <- st_combined[rownames(st_raw), , drop = FALSE]

          # Remove ASV columns that are zero across ALL samples
          asv_present <- colSums(st_combined) > 0
          st_combined <- st_combined[, asv_present, drop = FALSE]

          case_obj$st_nochim <- st_combined
        }

      } else {
        message("  [RAREFACTION] Case '", case_name_i,
                "': $st_nochim is absent or",
                " not a matrix -- skipping rarefaction.")
      }
    }
    # ─────────────────────────────────────────────────────────────────────────

    cases[[i]] <- case_obj
  }

  # Remove failed loads
  failed <- sapply(cases, is.null)
  if (any(failed)) {
    message("[WARNING] ", sum(failed), " case(s) could not be loaded and will be skipped.")
    cases <- cases[!failed]
  }

  message("Successfully loaded ", length(cases), " cases.")
  cases
}


# ==============================================================================
# SECTION 4 — METRIC COMPUTATION
# ==============================================================================

# Pattern used throughout to identify mock community samples
MOCK_PATTERN <- "Mock_mix"

#' Identify mock-community samples in a seqtab matrix.
#'
#' @param st  A seqtab matrix (rows = samples, columns = ASV sequences).
#' @return Character vector of matching rownames.
get_mock_samples <- function(st) {
  grep(MOCK_PATTERN, rownames(st), value = TRUE)
}

#' Compute benchmark metrics for a single case, restricted to mock samples.
#'
#' Metrics computed:
#'
#'   1. pct_chimeric_reads
#'      Percentage of reads (in mock samples) assigned to chimeric ASVs
#'      relative to ALL reads before chimera removal ($st).
#'
#'   2. mean_reads_per_sample
#'      Mean of per-sample total reads across mock samples ($st_nochim).
#'
#'   3. total_asvs
#'      Number of distinct ASVs with at least one read in any mock sample
#'      ($st_nochim).
#'
#'   4. pct_asvs_matching_mock
#'      Percentage of mock-present ASVs whose sequence matches at least one
#'      theoretical mock reference (IUPAC-aware, exact full-length match).
#'
#'   5. pct_asvs_not_matching_mock (= 100 - pct_asvs_matching_mock)
#'      Noise metric used in Part 2 ranking.
#'
#'   6. pct_reads_matching_mock
#'      Percentage of total mock reads ($st_nochim) from mock-matching ASVs.
#'      Recovery metric used in Part 2 ranking.
#'
#' @param case_obj   A loaded case list object.
#' @param mock_seqs  Named character vector of theoretical mock sequences.
#' @param mock_match_cache  Optional pre-computed match vector (speeds up batch runs).
#' @return Named numeric vector of metrics, or NAs if the case is malformed.
compute_case_metrics <- function(case_obj, mock_seqs, mock_match_cache = NULL) {

  # Validate object structure
  required_slots <- c("st", "st_nochim", "bim", "track")
  missing_slots  <- setdiff(required_slots, names(case_obj))
  if (length(missing_slots) > 0) {
    warning("Case object missing slots: ", paste(missing_slots, collapse = ", "))
    return(rep(NA_real_, 6) |> setNames(
      c("pct_chimeric_reads","mean_reads_per_sample","total_asvs",
        "pct_asvs_matching_mock","pct_asvs_not_matching_mock",
        "pct_reads_matching_mock")))
  }

  st        <- case_obj$st
  st_nochim <- case_obj$st_nochim
  bim       <- case_obj$bim   # named logical: TRUE = chimeric ASV

  # Identify mock samples
  mock_samples <- get_mock_samples(st_nochim)
  if (length(mock_samples) == 0) {
    warning("No mock samples found (pattern: '", MOCK_PATTERN, "').")
    return(rep(NA_real_, 6) |> setNames(
      c("pct_chimeric_reads","mean_reads_per_sample","total_asvs",
        "pct_asvs_matching_mock","pct_asvs_not_matching_mock",
        "pct_reads_matching_mock")))
  }

  # Subset to mock samples
  st_mock        <- st[mock_samples, , drop = FALSE]
  st_nochim_mock <- st_nochim[mock_samples, , drop = FALSE]

  # Metric 1: % chimeric reads (computed from $st, before chimera zeroing)
  chim_asvs       <- names(bim)[bim == TRUE]
  chim_asvs_in_st <- intersect(chim_asvs, colnames(st_mock))

  total_reads_mock <- sum(st_mock)
  chim_reads_mock  <- if (length(chim_asvs_in_st) > 0) sum(st_mock[, chim_asvs_in_st]) else 0

  pct_chimeric_reads <- if (total_reads_mock > 0) {
    chim_reads_mock / total_reads_mock * 100
  } else {
    NA_real_
  }

  # Metric 2: Mean reads per sample (from $st_nochim)
  reads_per_sample      <- rowSums(st_nochim_mock)
  mean_reads_per_sample <- mean(reads_per_sample)

  # Identify ASVs present in mock samples
  asv_totals_mock <- colSums(st_nochim_mock)
  present_asvs    <- names(asv_totals_mock)[asv_totals_mock > 0]

  # Metric 3: Total ASVs
  total_asvs <- length(present_asvs)

  # IUPAC-aware matching of present ASVs to theoretical mock
  if (!is.null(mock_match_cache)) {
    asv_is_mock <- mock_match_cache[present_asvs]
  } else {
    asv_is_mock <- match_asvs_to_mock(present_asvs, mock_seqs)
  }

  n_matching <- sum(asv_is_mock, na.rm = TRUE)

  # Metric 4: % ASVs matching mock
  pct_asvs_matching_mock <- if (total_asvs > 0) {
    n_matching / total_asvs * 100
  } else {
    NA_real_
  }

  # Metric 5: % ASVs NOT matching mock (noise metric)
  pct_asvs_not_matching_mock <- 100 - pct_asvs_matching_mock

  # Metric 6: % reads represented by matching ASVs (recovery metric)
  matching_asvs                <- present_asvs[asv_is_mock]
  reads_in_nochim_mock_total   <- sum(st_nochim_mock)
  reads_from_matching_asvs     <- if (length(matching_asvs) > 0) {
    sum(st_nochim_mock[, matching_asvs, drop = FALSE])
  } else {
    0
  }

  pct_reads_matching_mock <- if (reads_in_nochim_mock_total > 0) {
    reads_from_matching_asvs / reads_in_nochim_mock_total * 100
  } else {
    NA_real_
  }

  # Return all metrics
  c(
    pct_chimeric_reads         = pct_chimeric_reads,
    mean_reads_per_sample      = mean_reads_per_sample,
    total_asvs                 = total_asvs,
    pct_asvs_matching_mock     = pct_asvs_matching_mock,
    pct_asvs_not_matching_mock = pct_asvs_not_matching_mock,
    pct_reads_matching_mock    = pct_reads_matching_mock
  )
}

#' Run compute_case_metrics() on a list of cases and return a tidy data frame.
#'
#' @param cases_list  Named list of case objects (output of load_all_cases).
#' @param mock_seqs   Named character vector of theoretical mock sequences.
#' @return A data frame with one row per case, containing parsed parameters
#'         and all computed metrics.
build_metrics_table <- function(cases_list, mock_seqs) {
  rows <- vector("list", length(cases_list))

  for (i in seq_along(cases_list)) {
    case_name <- names(cases_list)[i]
    case_obj  <- cases_list[[i]]

    message("  Computing metrics for: ", case_name)

    # Parse parameters from filename
    params <- parse_case_filename(case_name)

    # Compute metrics (with error protection)
    metrics <- tryCatch(
      compute_case_metrics(case_obj, mock_seqs),
      error = function(e) {
        message("  [ERROR] Metric computation failed for ", case_name,
                ": ", conditionMessage(e))
        rep(NA_real_, 6) |> setNames(
          c("pct_chimeric_reads","mean_reads_per_sample","total_asvs",
            "pct_asvs_matching_mock","pct_asvs_not_matching_mock",
            "pct_reads_matching_mock"))
      }
    )

    rows[[i]] <- data.frame(
      case_name  = case_name,
      case_id    = params$case_id,
      filter_mode= params$filter_mode,
      error_model= params$error_model,
      pool_mode  = params$pool_mode,
      omega_raw  = params$omega_raw,
      omega_num  = params$omega_num,
      band_size  = params$band_size,
      as.data.frame(t(metrics)),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}


# ==============================================================================
# SECTION 5 — FACTOR ORDERING UTILITIES
# ==============================================================================

#' Apply sensible factor level ordering to the parameter columns of the
#' metrics table. Ensures consistent, readable axis ordering in plots.
#'
#' @param df  Data frame output of build_metrics_table().
#' @return The same data frame with factor columns added/replaced.
factorise_parameters <- function(df) {

  # Filter mode: alphabetical (Loose / Strict or whatever levels exist)
  if (!is.null(df$filter_mode)) {
    lvls <- sort(unique(df$filter_mode))
    df$filter_mode <- factor(df$filter_mode, levels = lvls)
  }

  # Error model: alphabetical
  if (!is.null(df$error_mode)) {
    lvls <- sort(unique(df$error_model))
    df$error_model <- factor(df$error_model, levels = lvls)
  }
  if (!is.null(df$error_model)) {
    lvls <- sort(unique(df$error_model))
    df$error_model <- factor(df$error_model, levels = lvls)
  }

  # Pool mode: logical-like order FALSE / pseudo / TRUE
  if (!is.null(df$pool_mode)) {
    pool_order <- c("FALSE", "pseudo", "TRUE")
    lvls <- pool_order[pool_order %in% unique(df$pool_mode)]
    df$pool_mode <- factor(df$pool_mode, levels = lvls)
  }

  # Omega: numeric order
  if (!is.null(df$omega_raw)) {
    omega_levels <- unique(df$omega_raw[order(df$omega_num)])
    df$omega_raw <- factor(df$omega_raw, levels = omega_levels)
  }

  # Band size: numeric order
  if (!is.null(df$band_size)) {
    band_levels <- sort(unique(df$band_size))
    df$band_size <- factor(df$band_size, levels = band_levels)
  }

  df
}


# ==============================================================================
# SECTION 6 — REGRESSION MODELLING
# ==============================================================================

#' Fit a regression model for a single metric against the five parameters.
#'
#' Model choice rationale:
#'   - pct_chimeric_reads / pct_asvs_*_mock / pct_reads_matching_mock:
#'     Quasibinomial GLM on the proportion (metric/100), logit link.
#'     Handles bounded [0,1] nature and overdispersion without extra packages.
#'   - mean_reads_per_sample: OLS on log-transformed reads.
#'     Log transform stabilises right-skewed variance.
#'   - total_asvs: Poisson GLM (log link). Count outcome.
#'
#' All models treat the five DADA2 parameters as fixed categorical effects.
#' Reference level for filter_mode is set to "Strict".
#'
#' @param df      Data frame from build_metrics_table() (factorised).
#' @param metric  Character: column name of the metric to model.
#' @return A list with elements: model, summary, tidy_coefs, family_used,
#'         metric, n_obs.
fit_metric_model <- function(df, metric) {

  # Drop rows with NA in the metric
  df_clean <- df[!is.na(df[[metric]]), ]

  if (nrow(df_clean) < 10) {
    warning("Too few complete cases (", nrow(df_clean), ") for metric: ", metric)
    return(NULL)
  }

  # Predictors: all five parameters as factors
  predictors <- c("filter_mode", "error_model", "pool_mode", "omega_raw", "band_size")
  for (p in predictors) {
    if (!is.factor(df_clean[[p]])) {
      df_clean[[p]] <- factor(df_clean[[p]])
    }
  }

  # Set filter_mode reference level to "Strict"
  if ("Strict" %in% levels(df_clean$filter_mode)) {
    df_clean$filter_mode <- relevel(factor(df_clean$filter_mode), ref = "Strict")
  }

  formula_str <- paste(metric, "~", paste(predictors, collapse = " + "))
  rhs_formula <- as.formula(paste("~", paste(predictors, collapse = " + ")))

  # Choose appropriate model
  if (metric %in% c("pct_chimeric_reads",
                     "pct_asvs_matching_mock",
                     "pct_asvs_not_matching_mock",
                     "pct_reads_matching_mock")) {
    # Quasibinomial GLM on proportion (metric / 100)
    prop_col <- paste0(metric, "_prop")
    df_clean[[prop_col]] <- df_clean[[metric]] / 100
    # Clamp to avoid 0/1 boundary issues
    df_clean[[prop_col]] <- pmax(pmin(df_clean[[prop_col]], 0.9999), 0.0001)

    formula_glm <- as.formula(paste(prop_col, "~", paste(predictors, collapse = " + ")))
    model <- glm(formula_glm, data = df_clean, family = quasibinomial(link = "logit"))
    family_used <- "quasibinomial (logit link) — response is metric/100"

  } else if (metric == "mean_reads_per_sample") {
    # OLS on log-transformed reads
    log_col <- paste0("log_", metric)
    df_clean[[log_col]] <- log(df_clean[[metric]] + 1)
    formula_lm <- as.formula(paste(log_col, "~", paste(predictors, collapse = " + ")))
    model <- lm(formula_lm, data = df_clean)
    family_used <- "Gaussian OLS — response is log(mean_reads_per_sample + 1)"

  } else if (metric == "total_asvs") {
    # Poisson GLM (counts)
    formula_glm <- as.formula(paste(metric, "~", paste(predictors, collapse = " + ")))
    model <- glm(formula_glm, data = df_clean, family = poisson(link = "log"))
    family_used <- "Poisson GLM (log link) — check for overdispersion"

  } else {
    # Fallback: OLS
    formula_lm <- as.formula(paste(metric, "~", paste(predictors, collapse = " + ")))
    model <- lm(formula_lm, data = df_clean)
    family_used <- "Gaussian OLS (fallback)"
  }

  # Tidy coefficient table (base R, no broom needed)
  coef_summary <- summary(model)$coefficients
  tidy_coefs   <- as.data.frame(coef_summary)

  # Standardise column names across lm / glm
  colnames(tidy_coefs) <- sub("^Estimate$",              "estimate",  colnames(tidy_coefs))
  colnames(tidy_coefs) <- sub("^Std\\. Error$",          "std_error", colnames(tidy_coefs))
  colnames(tidy_coefs) <- sub("^(t value|z value)$",     "statistic", colnames(tidy_coefs))
  colnames(tidy_coefs) <- sub("^Pr\\(>\\|[tz]\\|\\)$",  "p_value",   colnames(tidy_coefs))

  tidy_coefs$term <- rownames(tidy_coefs)
  rownames(tidy_coefs) <- NULL

  list(
    model       = model,
    summary     = summary(model),
    tidy_coefs  = tidy_coefs,
    family_used = family_used,
    metric      = metric,
    n_obs       = nrow(df_clean)
  )
}

#' Fit models for all target metrics and return a named list.
#'
#' @param df  Data frame from build_metrics_table() (factorised).
#' @return Named list of model objects (output of fit_metric_model).
fit_all_models <- function(df) {
  target_metrics <- c(
    "pct_chimeric_reads",
    "mean_reads_per_sample",
    "pct_asvs_matching_mock",
    "total_asvs"
  )

  models <- lapply(target_metrics, function(m) {
    message("  Fitting model for: ", m)
    fit_metric_model(df, m)
  })
  names(models) <- target_metrics
  models
}

#' Extract a combined coefficient table from all models.
#'
#' @param models_list  Output of fit_all_models().
#' @return Data frame with one row per coefficient per model.
extract_all_coefficients <- function(models_list) {
  rows <- lapply(names(models_list), function(metric) {
    m <- models_list[[metric]]
    if (is.null(m)) return(NULL)
    tidy <- m$tidy_coefs
    tidy$metric <- metric
    tidy$family <- m$family_used
    tidy$n_obs  <- m$n_obs
    tidy
  })
  dplyr::bind_rows(rows[!sapply(rows, is.null)])
}


# ==============================================================================
# SECTION 7 — MARGINAL MEANS
# ==============================================================================

#' Compute observed marginal means for a metric, grouped by a single parameter.
#'
#' Simple group-level mean of the raw metric — interpretable without
#' back-transformation, appropriate for communicating parameter effects.
#'
#' @param df       Data frame (factorised).
#' @param metric   Character: metric column name.
#' @param param    Character: parameter column name.
#' @return Data frame: param_value, mean, sd, n, se, ci_low, ci_high.
marginal_means_by_param <- function(df, metric, param) {
  df_clean <- df[!is.na(df[[metric]]) & !is.na(df[[param]]), ]

  result <- do.call(rbind, lapply(levels(df_clean[[param]]), function(lv) {
    sub <- df_clean[df_clean[[param]] == lv, metric]
    n   <- length(sub)
    mn  <- mean(sub)
    sd  <- sd(sub)
    se  <- sd / sqrt(n)
    ci  <- qt(0.975, df = n - 1) * se
    data.frame(
      param       = param,
      param_value = lv,
      mean        = mn,
      sd          = sd,
      n           = n,
      se          = se,
      ci_low      = mn - ci,
      ci_high     = mn + ci,
      stringsAsFactors = FALSE
    )
  }))
  result
}

#' Compute marginal means for ALL metrics × ALL parameters.
#'
#' @param df  Factorised metrics data frame.
#' @return Data frame with metric, param, param_value, mean, sd, n, se, ci_low, ci_high.
all_marginal_means <- function(df) {
  target_metrics <- c(
    "pct_chimeric_reads",
    "pct_asvs_matching_mock",
    "total_asvs"
  )
  params <- c("filter_mode", "error_model", "pool_mode", "omega_raw", "band_size")

  rows <- list()
  for (m in target_metrics) {
    for (p in params) {
      mm <- marginal_means_by_param(df, m, p)
      mm$metric <- m
      rows[[length(rows) + 1]] <- mm
    }
  }
  do.call(rbind, rows)
}


# ==============================================================================
# SECTION 8 — RANKING (Part 2)
# ==============================================================================

#' Rank parameter combinations using a transparent composite score.
#'
#' Composite score = -z_noise + z_recovery, where:
#'   z_noise    = standardised % ASVs not matching mock  (lower = better)
#'   z_recovery = standardised % reads matching mock     (higher = better)
#'
#' Higher score = less noise AND more recovery = BETTER.
#' Both metrics are weighted equally by default (weight = 1).
#'
#' @param df              Factorised metrics data frame.
#' @param weight_noise    Weight for noise metric (default 1).
#' @param weight_recovery Weight for recovery metric (default 1).
#' @return Input data frame with added columns:
#'           z_noise, z_recovery, composite_score, rank.
rank_parameter_combinations <- function(df,
                                         weight_noise    = 1,
                                         weight_recovery = 1) {

  df_rank <- df[!is.na(df$pct_asvs_not_matching_mock) &
                !is.na(df$pct_reads_matching_mock), ]

  if (nrow(df_rank) == 0) stop("No complete cases for ranking.")

  # Standardise
  mean_A <- mean(df_rank$pct_asvs_not_matching_mock)
  sd_A   <- sd(df_rank$pct_asvs_not_matching_mock)
  mean_B <- mean(df_rank$pct_reads_matching_mock)
  sd_B   <- sd(df_rank$pct_reads_matching_mock)

  df_rank$z_noise     <- (df_rank$pct_asvs_not_matching_mock - mean_A) / sd_A
  df_rank$z_recovery  <- (df_rank$pct_reads_matching_mock    - mean_B) / sd_B

  # Composite: penalise noise, reward recovery
  df_rank$composite_score <- -weight_noise * df_rank$z_noise +
                              weight_recovery * df_rank$z_recovery

  # Rank descending (best = rank 1)
  df_rank$rank <- rank(-df_rank$composite_score, ties.method = "min")

  # Sort by rank
  df_rank[order(df_rank$rank), ]
}

#' Extract the top N parameter combinations, including all ties at the cut-off.
#'
#' @param ranked_df  Output of rank_parameter_combinations().
#' @param n          Nominal cut-off (default 5).
#' @return Data frame of top cases (may contain more than N rows if ties exist).
get_top_n <- function(ranked_df, n = 5) {
  sorted <- ranked_df[order(-ranked_df$composite_score), ]

  if (nrow(sorted) < n) return(sorted)

  # Score of the case sitting at nominal position n — the inclusion threshold
  cutoff_score <- sorted$composite_score[n]

  # Return all cases at or above the threshold (captures all ties)
  ranked_df[ranked_df$composite_score >= cutoff_score, ]
}

#' Classify all cases into rank groups based on the top 3 distinct scores.
#'
#' Labels: "1st", "2nd", "3rd", "Other". All ties within a score level
#' receive the same group label.
#'
#' @param ranked_df  Output of rank_parameter_combinations().
#' @return Same data frame with a new factor column rank_group.
assign_rank_groups <- function(ranked_df) {
  top3_scores <- sort(unique(ranked_df$composite_score), decreasing = TRUE)[1:3]

  ranked_df$rank_group <- dplyr::case_when(
    ranked_df$composite_score == top3_scores[1] ~ "1st",
    ranked_df$composite_score == top3_scores[2] ~ "2nd",
    ranked_df$composite_score == top3_scores[3] ~ "3rd",
    TRUE ~ "Other"
  )

  ranked_df$rank_group <- factor(ranked_df$rank_group,
                                  levels = c("1st", "2nd", "3rd", "Other"))
  # Store the 3 score values so plots can display them in the legend
  attr(ranked_df, "top3_scores") <- top3_scores
  ranked_df
}

#' Build a human-readable label for a case's parameter combination.
#'
#' @param row  A single-row data frame (from the ranked table).
#' @return Character label like "F:Loose | Err:Binned | Pool:FALSE | Ω:1e-200 | B:32"
make_case_label <- function(row) {
  paste0(
    "F:", row$filter_mode, " | ",
    "Err:", row$error_model, " | ",
    "Pool:", row$pool_mode, " | ",
    "\u03a9:", row$omega_raw, " | ",
    "B:", row$band_size
  )
}


# ==============================================================================
# SECTION 9 — PLOT THEME
# ==============================================================================

#' Return a consistent ggplot2 theme for all benchmark figures.
#'
#' @return A ggplot2 theme object.
theme_benchmark <- function(base_size = 13) {
  ggplot2::theme_bw(base_size = base_size) +
  ggplot2::theme(
    plot.title      = ggplot2::element_text(face = "bold", size = base_size + 3,
                                            hjust = 0, margin = ggplot2::margin(b = 8)),
    plot.subtitle   = ggplot2::element_text(size = base_size + 1, color = "grey30",
                                            hjust = 0, margin = ggplot2::margin(b = 10)),
    axis.title      = ggplot2::element_text(face = "bold", size = base_size + 1),
    axis.text       = ggplot2::element_text(size = base_size - 1),
    axis.text.x     = ggplot2::element_text(angle = 30, hjust = 1),
    legend.title    = ggplot2::element_text(face = "bold"),
    legend.position = "right",
    strip.text      = ggplot2::element_text(face = "bold", size = base_size),
    strip.background= ggplot2::element_rect(fill = "#EEF2F7", colour = "grey80"),
    panel.grid.major= ggplot2::element_line(colour = "grey88"),
    panel.grid.minor= ggplot2::element_blank(),
    panel.border    = ggplot2::element_rect(colour = "grey70"),
    plot.margin     = ggplot2::margin(12, 15, 12, 12)
  )
}

# Colour palette for parameters (up to 8 levels per parameter)
PARAM_COLOURS <- c(
  "#2E86AB", "#A23B72", "#F18F01", "#C73E1D",
  "#3B1F2B", "#44BBA4", "#E94F37", "#393E41"
)


# ==============================================================================
# SECTION 10 — PLOT FUNCTIONS — PART 1
# ==============================================================================

#' Plot marginal means (with 95% CI) for all parameters × one metric.
#'
#' @param mm_df   Output of all_marginal_means() or marginal_means_by_param().
#' @param metric  Metric to plot.
#' @param y_label Y-axis label string.
#' @return A ggplot object.
plot_marginal_means <- function(mm_df, metric, y_label = NULL) {
  sub_df <- mm_df[mm_df$metric == metric, ]
  if (nrow(sub_df) == 0) stop("No data for metric: ", metric)

  if (is.null(y_label)) y_label <- gsub("_", " ", metric)

  param_labels <- c(
    filter_mode  = "Filter mode",
    error_model  = "Error model",
    pool_mode    = "Pooling mode",
    omega_raw    = "Omega (\u03c9)",
    band_size    = "Band size"
  )

  sub_df$param_label <- param_labels[sub_df$param]
  sub_df$param_label <- factor(sub_df$param_label,
                                levels = param_labels)

  ggplot2::ggplot(sub_df,
         ggplot2::aes(x = param_value, y = mean,
                      ymin = ci_low, ymax = ci_high,
                      colour = param_label, fill = param_label)) +
    ggplot2::geom_errorbar(width = 0.25, linewidth = 0.8, alpha = 0.7) +
    ggplot2::geom_point(size = 3.5, shape = 21, colour = "white", stroke = 1.5) +
    ggplot2::facet_wrap(~ param_label, scales = "free_x", nrow = 1) +
    ggplot2::scale_colour_manual(values = PARAM_COLOURS, guide = "none") +
    ggplot2::scale_fill_manual(values = PARAM_COLOURS, guide = "none") +
    ggplot2::labs(
      title    = paste("Effect of DADA2 parameters on:", y_label),
      subtitle = "Observed marginal means \u00b1 95% CI across all cases",
      x        = "Parameter level",
      y        = y_label
    ) +
    theme_benchmark()
}

#' Boxplot of a metric by parameter level, with individual case points.
#'
#' @param df      Factorised metrics data frame.
#' @param metric  Metric column name.
#' @param param   Parameter column name.
#' @param y_label Y-axis label.
#' @param x_label X-axis label.
#' @return A ggplot object.
plot_boxplot_by_param <- function(df, metric, param,
                                   y_label = NULL,
                                   x_label = NULL) {
  df_clean <- df[!is.na(df[[metric]]) & !is.na(df[[param]]), ]
  if (nrow(df_clean) == 0) stop("No data for ", metric, " by ", param)

  if (is.null(y_label)) y_label <- gsub("_", " ", metric)
  if (is.null(x_label)) x_label <- gsub("_", " ", param)

  ggplot2::ggplot(df_clean,
         ggplot2::aes_string(x = param, y = metric, fill = param)) +
    ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
    ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, colour = "grey20") +
    ggplot2::scale_fill_manual(values = PARAM_COLOURS, guide = "none") +
    ggplot2::labs(
      x = x_label,
      y = y_label
    ) +
    theme_benchmark()
}

#' Create a coefficient plot from a model fit (showing effect sizes with CIs).
#'
#' Intercept is excluded. For GLMs, coefficients are on the link scale.
#'
#' @param model_obj  Output of fit_metric_model() for one metric.
#' @return A ggplot object.
plot_coefficient <- function(model_obj) {
  if (is.null(model_obj)) return(NULL)

  tidy <- model_obj$tidy_coefs
  if (!all(c("estimate", "std_error") %in% colnames(tidy))) {
    warning("Cannot identify 'estimate' / 'std_error' columns in coefficient table.")
    return(NULL)
  }

  tidy$ci_low  <- tidy$estimate - 1.96 * tidy$std_error
  tidy$ci_high <- tidy$estimate + 1.96 * tidy$std_error

  # Remove intercept
  tidy <- tidy[tidy$term != "(Intercept)", ]
  if (nrow(tidy) == 0) return(NULL)

  # Colour by parameter
  tidy$parameter <- sub("(filter_mode|error_model|pool_mode|omega_raw|band_size).*", "\\1", tidy$term)

  param_labels <- c(
    filter_mode  = "Filter mode",
    error_model  = "Error model",
    pool_mode    = "Pooling mode",
    omega_raw    = "Omega (\u03c9)",
    band_size    = "Band size"
  )
  tidy$param_label <- param_labels[tidy$parameter]
  tidy$param_label[is.na(tidy$param_label)] <- "Other"

  # Order by effect size
  tidy <- tidy[order(tidy$estimate), ]
  tidy$term <- factor(tidy$term, levels = tidy$term)

  ggplot2::ggplot(tidy,
         ggplot2::aes(x = estimate, y = term, colour = param_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_low, xmax = ci_high),
                            height = 0.3, linewidth = 0.8, alpha = 0.7) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(values = PARAM_COLOURS, name = "Parameter") +
    ggplot2::labs(
      title    = paste("Coefficient plot:", gsub("_", " ", model_obj$metric)),
      subtitle = paste0("Model: ", model_obj$family_used,
                        "  |  n = ", model_obj$n_obs),
      x        = "Coefficient estimate (link scale)",
      y        = "Term"
    ) +
    theme_benchmark() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
}


# ==============================================================================
# SECTION 11 — PLOT FUNCTIONS — PART 2
# ==============================================================================

#' Scatter plot of the two ranking metrics, with top-selected cases highlighted.
#'
#' Text labels are omitted from the plot interior to avoid clutter.
#'
#' @param ranked_df  Output of rank_parameter_combinations().
#' @param n          Nominal cut-off passed to get_top_n() (default 5).
#' @return A ggplot object.
plot_ranking_scatter <- function(ranked_df, n = 5) {

  # Identify top cases (tie-aware: may be more than n rows)
  top_df       <- get_top_n(ranked_df, n)
  n_actual     <- nrow(top_df)
  top_names    <- top_df$case_name

  ranked_df$is_top <- ranked_df$case_name %in% top_names

  df_other <- ranked_df[!ranked_df$is_top, ]
  df_top   <- ranked_df[ ranked_df$is_top, ]

  legend_label <- if (n_actual == n) {
    paste0("Top ", n, " cases")
  } else {
    paste0("Top ", n, " + ties (n = ", n_actual, ")")
  }

  ggplot2::ggplot(mapping = ggplot2::aes(
      x = pct_asvs_not_matching_mock,
      y = pct_reads_matching_mock)) +
    ggplot2::geom_point(data = df_other,
                        colour = "grey72", size = 2.2,
                        alpha = 0.7, shape = 16) +
    ggplot2::geom_point(data = df_top,
                        colour = "#9B2226", fill = "#C73E1D",
                        size = 4, alpha = 0.92,
                        shape = 21, stroke = 1.2) +
    ggplot2::geom_point(data = data.frame(
                          pct_asvs_not_matching_mock = NA_real_,
                          pct_reads_matching_mock    = NA_real_,
                          group = factor(c("other", "top"),
                                         levels = c("top", "other"))),
                        ggplot2::aes(colour = group),
                        size = 3, na.rm = TRUE) +
    ggplot2::scale_colour_manual(
      values = c("top" = "#C73E1D", "other" = "grey72"),
      labels = c("top" = legend_label, "other" = "Other cases"),
      name   = NULL,
      guide  = ggplot2::guide_legend(
        override.aes = list(
          size   = c(4, 2.5),
          shape  = c(21, 16),
          fill   = c("#C73E1D", NA),
          stroke = c(1.2, 0),
          alpha  = c(0.92, 0.7)
        )
      )
    ) +
    ggplot2::labs(
      title    = "Parameter combination ranking",
      subtitle = paste0(
        "Each point = one of the 96 parameter combinations.  ",
        "Ideal = top-left (low noise, high recovery).\n",
        "Highlighted cases: score \u2265 score of rank-", n, " case."
      ),
      x = "% ASVs not matching mock  (noise \u2014 lower is better)",
      y = "% reads from mock-matching ASVs\n(recovery \u2014 higher is better)"
    ) +
    theme_benchmark() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = ggplot2::element_rect(fill = "white",
                                                colour = "grey80",
                                                linewidth = 0.4),
      legend.margin   = ggplot2::margin(5, 8, 5, 8)
    )
}

#' Bar plot of composite scores for the top-3-score cases, coloured by rank group.
#'
#' @param ranked_df  Output of rank_parameter_combinations().
#' @return A ggplot object.
plot_top3_scores <- function(ranked_df) {
  # Add rank_group if not already present
  if (!"rank_group" %in% names(ranked_df)) {
    ranked_df <- assign_rank_groups(ranked_df)
  }

  # Keep only top-3 score groups
  top_df <- ranked_df[ranked_df$rank_group %in% c("1st", "2nd", "3rd"), ]

  top_df <- top_df[order(top_df$rank_group, top_df$composite_score,
                          decreasing = c(TRUE, FALSE), method = "radix"), ]
  top_df <- top_df[order(top_df$rank_group, decreasing = TRUE), ]

  top_df$label <- paste0(top_df$rank_group, ": ",
                          sapply(seq_len(nrow(top_df)),
                                 function(i) make_case_label(top_df[i, ])))
  top_df$label <- factor(top_df$label, levels = rev(top_df$label))

  rank_colors <- c("1st" = "#C73E1D", "2nd" = "#F4831F", "3rd" = "#F7C948")

  s <- attr(ranked_df, "top3_scores")
  if (is.null(s)) s <- sort(unique(ranked_df$composite_score), decreasing = TRUE)[1:3]
  rank_labels <- c(
    "1st" = sprintf("1st (score = %.3f)", s[1]),
    "2nd" = sprintf("2nd (score = %.3f)", s[2]),
    "3rd" = sprintf("3rd (score = %.3f)", s[3])
  )

  n_groups <- table(top_df$rank_group)
  subtitle_str <- paste0(
    "1st: n=", n_groups["1st"], " | ",
    "2nd: n=", n_groups["2nd"], " | ",
    "3rd: n=", n_groups["3rd"],
    " \u2014 all tied cases within each score level shown"
  )

  ggplot2::ggplot(top_df,
         ggplot2::aes(x = composite_score, y = label,
                      fill = rank_group)) +
    ggplot2::geom_col(width = 0.6, alpha = 0.9) +
    ggplot2::scale_fill_manual(
      values = rank_colors,
      labels = rank_labels,
      name   = "Rank group"
    ) +
    ggplot2::labs(
      title    = "Top 3 score levels \u2014 composite score",
      subtitle = subtitle_str,
      x        = "Composite score (standardised units)",
      y        = "Parameter combination"
    ) +
    theme_benchmark() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
                   axis.text.y = ggplot2::element_text(size = 10),
                   legend.position = "right")
}

#' Dual-metric bar chart for top cases: shows both raw metrics side by side.
#'
#' @param ranked_df  Output of rank_parameter_combinations().
#' @param n          Nominal cut-off passed to get_top_n() (default 5).
#' @return A ggplot object.
plot_top_n_dual_metrics <- function(ranked_df, n = 5) {
  top_df   <- get_top_n(ranked_df, n)
  n_actual <- nrow(top_df)
  top_df   <- top_df[order(top_df$composite_score), ]
  top_df$label <- paste0("#", top_df$rank, ": ",
                          sapply(seq_len(nrow(top_df)),
                                 function(i) make_case_label(top_df[i,])))
  top_df$label <- factor(top_df$label, levels = top_df$label)

  title_str <- if (n_actual == n) {
    paste0("Top ", n, " combinations \u2014 raw metric values")
  } else {
    paste0("Top ", n, " + tied cases (n = ", n_actual, ") \u2014 raw metric values")
  }

  # Reshape to long format
  long_df <- rbind(
    data.frame(label  = top_df$label,
               metric = "% ASVs not matching mock\n(noise \u2014 lower is better)",
               value  = top_df$pct_asvs_not_matching_mock,
               good   = "lower",
               stringsAsFactors = FALSE),
    data.frame(label  = top_df$label,
               metric = "% reads from mock-matching ASVs\n(recovery \u2014 higher is better)",
               value  = top_df$pct_reads_matching_mock,
               good   = "higher",
               stringsAsFactors = FALSE)
  )

  fill_colors <- c(
    "% ASVs not matching mock\n(noise \u2014 lower is better)"           = "#2E86AB",
    "% reads from mock-matching ASVs\n(recovery \u2014 higher is better)" = "#44BBA4"
  )

  ggplot2::ggplot(long_df,
         ggplot2::aes(x = value, y = label, fill = metric)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                      width = 0.6, alpha = 0.88) +
    ggplot2::scale_fill_manual(values = fill_colors, name = "Metric") +
    ggplot2::labs(
      title    = title_str,
      subtitle = "Both decision metrics shown side by side",
      x        = "Value (%)",
      y        = "Parameter combination"
    ) +
    theme_benchmark() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.position = "bottom",
      legend.direction = "vertical"
    )
}


# ==============================================================================
# SECTION 12 — SAVE UTILITIES
# ==============================================================================

#' Save a ggplot to a file at high resolution.
#'
#' @param plot    A ggplot object.
#' @param path    Full file path (extension determines format: .png, .pdf, .svg).
#' @param width   Width in inches (default 10).
#' @param height  Height in inches (default 6).
#' @param dpi     Resolution (default 300).
save_plot <- function(plot, path, width = 10, height = 6, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height,
                   dpi = dpi, bg = "white")
  message("  Saved: ", path)
}

#' Save a data frame as a CSV file.
#'
#' @param df    Data frame to save.
#' @param path  Full file path.
save_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE)
  message("  Saved: ", path)
}

#' Save a model summary as a plain-text file.
#'
#' @param model_obj  Output of fit_metric_model().
#' @param path       Full file path (.txt).
save_model_summary <- function(model_obj, path) {
  if (is.null(model_obj)) return(invisible(NULL))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  sink(path)
  cat("==========================================================\n")
  cat("Metric: ", model_obj$metric, "\n")
  cat("Family: ", model_obj$family_used, "\n")
  cat("N observations: ", model_obj$n_obs, "\n")
  cat("==========================================================\n\n")
  print(model_obj$summary)
  sink()
  message("  Saved: ", path)
}
