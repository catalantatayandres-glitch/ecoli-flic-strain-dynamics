# ==============================================================================
# 05_Mock_Distributions_Functions.R
#
# Project : fliC Amplicon Sequencing for High-Resolution Profiling of E. coli
#           in the Infant Gut Microbiome
# Author  : Andrés Catalán Tatay
# Part    : Analysis — Step 05 (companion functions)
# Purpose : All helper functions for 05_Mock_Distributions.Rmd.
#           Implements a custom Check Mock workflow to carry a project-specific reference-based
#           approach using pairwise alignment and Bray-Curtis similarity.
#

#
# Inputs  : data/cases/  — DADA2 output RDS files
#            data/mock/ps_mock_trimmed_4.rds  — theoretical mock reference
# Outputs : Called by 05_Mock_Distributions.Rmd
#
# Usage   : source("05_Mock_Distributions_Functions.R")
# ==============================================================================


# ==============================================================================
# SECTION 1 — PACKAGE INSTALLATION & LOADING
# ==============================================================================

.install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("[INSTALL] Missing: ", pkg, " -- installing...")
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

.install_bioc_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("[INSTALL] Missing Bioc package: ", pkg, " -- installing...")
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
}

.install_chkMocks <- function() {
  if (!requireNamespace("chkMocks", quietly = TRUE)) {
    message("[INSTALL] chkMocks not found -- installing from GitHub...")
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools", repos = "https://cloud.r-project.org")
    }
    devtools::install_github("microsud/chkMocks")
  }
}

.install_if_missing(c("devtools", "dplyr", "ggplot2", "patchwork", "RColorBrewer",
                       "tidyr", "tibble"))
.install_bioc_if_missing(c("phyloseq", "Biostrings", "microbiome", "DECIPHER", "pwalign"))
.install_chkMocks()

library(phyloseq)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(tidyr)
library(tibble)
library(microbiome)
library(pwalign)   # pairwiseAlignment moved here from Biostrings >= 2.77.1
library(chkMocks)

message("[OK] All packages loaded successfully.")


# ==============================================================================
# SECTION 2 — DATA LOADING
# ==============================================================================

#' Load a case RDS file and extract the phyloseq object.
#'
#' @param path  Full path to the .rds file.
#' @return phyloseq object.
load_case_phyloseq <- function(path) {
  message("[LOAD] Case file: ", path)
  if (!file.exists(path)) stop("[ERROR] File not found: ", path)
  obj <- readRDS(path)
  if (inherits(obj, "phyloseq")) {
    message("[LOAD] File is directly a phyloseq object.")
    return(obj)
  }
  candidates <- c("ps", "phyloseq", "physeq", "ps_object", "phyloseq_object",
                  "result", "output", "pseq", "ps_all", "ps_final",
                  "ps_dada2", "data")
  if (is.list(obj)) {
    for (nm in candidates) {
      if (!is.null(obj[[nm]]) && inherits(obj[[nm]], "phyloseq")) {
        message("[LOAD] Extracted phyloseq from element: '", nm, "'")
        return(obj[[nm]])
      }
    }
    for (nm in names(obj)) {
      if (inherits(obj[[nm]], "phyloseq")) {
        message("[LOAD] Extracted phyloseq from element: '", nm, "'")
        return(obj[[nm]])
      }
    }
    stop("[ERROR] No phyloseq found.\n  Elements: ",
         paste(names(obj), collapse = ", "))
  }
  stop("[ERROR] Not a phyloseq or list. Class: ", class(obj))
}

#' Load the theoretical mock reference phyloseq.
#'
#' @param path  Full path to the .rds file.
#' @return phyloseq object.
load_theoretical_mock <- function(path) {
  message("[LOAD] Theoretical mock: ", path)
  if (!file.exists(path)) stop("[ERROR] File not found: ", path)
  obj <- readRDS(path)
  if (inherits(obj, "phyloseq")) {
    message("[LOAD] Theoretical mock is directly a phyloseq object.")
    return(obj)
  }
  if (is.list(obj)) {
    for (nm in names(obj)) {
      if (inherits(obj[[nm]], "phyloseq")) {
        message("[LOAD] Extracted from element: '", nm, "'")
        return(obj[[nm]])
      }
    }
  }
  stop("[ERROR] No phyloseq found. Elements: ", paste(names(obj), collapse=", "))
}


# ==============================================================================
# SECTION 3 — DIAGNOSTICS
# ==============================================================================

#' Print a diagnostic summary of a phyloseq object.
#'
#' @param ps     phyloseq object.
#' @param label  Label string for the output header.
#' @return Invisibly returns ps.
diagnose_ps <- function(ps, label = "phyloseq diagnostics") {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("DIAGNOSTICS: ", label, "\n", sep = "")
  cat(strrep("-", 60), "\n")
  cat("  ntaxa            :", ntaxa(ps), "\n")
  cat("  nsamples         :", nsamples(ps), "\n")
  has_otu  <- !is.null(otu_table(ps, errorIfNULL = FALSE))
  has_tax  <- !is.null(access(ps, "tax_table", errorIfNULL = FALSE))
  has_ref  <- !is.null(refseq(ps, errorIfNULL = FALSE))
  has_tree <- !is.null(phy_tree(ps, errorIfNULL = FALSE))
  cat("  otu_table        :", if (has_otu)  "PRESENT" else "ABSENT", "\n")
  cat("  tax_table        :", if (has_tax)  "PRESENT" else "ABSENT", "\n")
  cat("  refseq()         :", if (has_ref)  "PRESENT" else "ABSENT", "\n")
  cat("  phy_tree         :", if (has_tree) "PRESENT" else "ABSENT", "\n")
  tn <- taxa_names(ps)
  cat("  taxa_names[1:3]  :\n    ",
      paste(substr(tn[seq_len(min(3, length(tn)))], 1, 70),
            collapse = "\n    "), "\n")
  cat("  taxa_names = DNA :", all(grepl("^[ACGTNacgtn]+$", tn)), "\n")
  if (has_ref) {
    rs <- refseq(ps)
    cat("  refseq n         :", length(rs), "\n")
    cat("  refseq in sync   :", setequal(names(rs), tn), "\n")
  }
  if (has_otu) cat("  zero-count taxa  :", sum(taxa_sums(ps) == 0), "\n")
  cat(strrep("=", 60), "\n\n")
  invisible(ps)
}


# ==============================================================================
# SECTION 4 — MOCK SAMPLE EXTRACTION
# ==============================================================================

#' Extract S_Mock_mix samples and remove zero-count taxa.
#'
#' @param ps  phyloseq object.
#' @return phyloseq subset containing only S_Mock_mix* samples.
extract_mock_samples <- function(ps) {
  message("[MOCK] Searching for samples matching '^S_Mock_mix'...")
  mock_samples <- grep("^S_Mock_mix", sample_names(ps), value = TRUE)
  if (length(mock_samples) == 0) {
    stop("[ERROR] No samples starting with 'S_Mock_mix' found.\n",
         "  All sample names: ", paste(sample_names(ps), collapse = ", "))
  }
  message("[MOCK] Found ", length(mock_samples), " mock sample(s): ",
          paste(mock_samples, collapse = ", "))
  ps_mocks <- prune_samples(mock_samples, ps)
  n_before <- ntaxa(ps_mocks)
  ps_mocks <- prune_taxa(taxa_sums(ps_mocks) > 0, ps_mocks)
  n_pruned <- n_before - ntaxa(ps_mocks)
  if (n_pruned > 0) message("[MOCK] Removed ", n_pruned, " all-zero taxa after subsetting.")
  message("[MOCK] Result: ", ntaxa(ps_mocks), " taxa x ", nsamples(ps_mocks), " samples.")
  ps_mocks
}

#' Extract Singular Mock samples (S_KPC1, S_Oxa, S_E15, S_E25) and remove
#' zero-count taxa.
#'
#' @param ps  phyloseq object.
#' @return phyloseq subset containing only the four singular mock samples.
extract_singular_mock_samples <- function(ps) {
  singular_mock_names <- c("S_KPC1", "S_Oxa", "S_E15", "S_E25")
  message("[SINGULAR MOCK] Searching for samples: ",
          paste(singular_mock_names, collapse = ", "), "...")
  found   <- intersect(singular_mock_names, sample_names(ps))
  missing <- setdiff(singular_mock_names, sample_names(ps))
  if (length(missing) > 0)
    warning("[SINGULAR MOCK] Not found: ", paste(missing, collapse = ", "))
  if (length(found) == 0)
    stop("[ERROR] None of the Singular Mock samples were found.\n",
         "  All sample names: ", paste(sample_names(ps), collapse = ", "))
  message("[SINGULAR MOCK] Found ", length(found), " sample(s): ",
          paste(found, collapse = ", "))
  ps_singular <- prune_samples(found, ps)
  n_before <- ntaxa(ps_singular)
  ps_singular <- prune_taxa(taxa_sums(ps_singular) > 0, ps_singular)
  n_pruned <- n_before - ntaxa(ps_singular)
  if (n_pruned > 0)
    message("[SINGULAR MOCK] Removed ", n_pruned, " all-zero taxa after subsetting.")
  message("[SINGULAR MOCK] Result: ", ntaxa(ps_singular), " taxa x ",
          nsamples(ps_singular), " samples.")
  ps_singular
}

#' Extract both Mock Mix and Singular Mock samples into one combined phyloseq.
#'
#' Handles the NA-filling bug in merge_phyloseq() by using rowSums(na.rm=TRUE)
#' for zero-taxon pruning, and re-attaches the refseq() slot which merge_phyloseq()
#' silently drops.
#'
#' @param ps  Full case phyloseq object.
#' @return Combined phyloseq with all mock samples.
extract_all_mock_samples <- function(ps) {
  ps_mix      <- extract_mock_samples(ps)
  ps_singular <- extract_singular_mock_samples(ps)
  message("[ALL MOCK] Merging Mock Mix and Singular Mock samples...")
  ps_all <- merge_phyloseq(ps_mix, ps_singular)
  n_before <- ntaxa(ps_all)

  # Use rowSums(na.rm = TRUE) to avoid merge_phyloseq() NA-filling bug:
  # merge_phyloseq() fills absent taxa with NA rather than 0. taxa_sums()
  # propagates NAs, so prune_taxa(taxa_sums > 0) would silently drop taxa
  # present in only one of the two input objects (e.g. taxa only in S_E15).
  otu_raw  <- as(otu_table(ps_all), "matrix")
  if (!taxa_are_rows(ps_all)) otu_raw <- t(otu_raw)
  keep     <- rowSums(otu_raw, na.rm = TRUE) > 0
  ps_all   <- prune_taxa(keep, ps_all)
  n_pruned <- n_before - ntaxa(ps_all)
  if (n_pruned > 0)
    message("[ALL MOCK] Removed ", n_pruned, " all-zero taxa after merging.")
  message("[ALL MOCK] Combined result: ", ntaxa(ps_all), " taxa x ",
          nsamples(ps_all), " samples.")

  # Re-attach refseq() slot: merge_phyloseq() silently drops it when merging
  # two objects that each carry their own refseq. assign_taxonomy_from_theoretical()
  # requires refseq() to be present for every surviving taxon.
  rs_full <- refseq(ps, errorIfNULL = FALSE)
  if (!is.null(rs_full)) {
    surviving <- taxa_names(ps_all)
    rs_keep   <- rs_full[names(rs_full) %in% surviving]
    if (length(rs_keep) > 0) {
      ps_all <- merge_phyloseq(ps_all, rs_keep)
      message("[ALL MOCK] refseq() re-attached for ", length(rs_keep), " taxa.")
    }
  } else {
    warning("[ALL MOCK] Source phyloseq has no refseq() slot -- ",
            "assign_taxonomy_from_theoretical() will fail if refseq is needed.")
  }

  ps_all
}


# ==============================================================================
# SECTION 5 — CUSTOM CHECK MOCK (PROJECT-SPECIFIC REFERENCE)
# ==============================================================================

#' Assign taxonomy to experimental ASVs by pairwise alignment against
#' the theoretical mock reference sequences.
#'
#' Similarity denominator: nmatch + nmismatch (base-vs-base positions only),
#' not any string length. This ensures a perfect-sequence ASV always scores
#' exactly 100% regardless of length differences between query and subject.
#'
#' Non-matching ASVs receive a unique "nonMock::<ASV_ID>" label so they
#' survive aggregate_taxa("species") as distinct taxa and can be individually
#' coloured in plots.
#'
#' @param ps_mocks       Experimental mock phyloseq (refseq slot required).
#' @param ps_theoretical Theoretical mock phyloseq (refseq slot required).
#' @param min_sim        Minimum % identity to assign a taxon (default 100).
#' @return ps_mocks with a tax_table added (column "species").
assign_taxonomy_from_theoretical <- function(ps_mocks,
                                             ps_theoretical,
                                             min_sim = 100) {

  message("[TAX] Assigning taxonomy by sequence matching against theoretical mock...")

  exp_rs <- refseq(ps_mocks, errorIfNULL = FALSE)
  if (is.null(exp_rs))
    stop("[ERROR] ps_mocks has no refseq() slot. Cannot assign taxonomy.")

  theo_rs <- refseq(ps_theoretical, errorIfNULL = FALSE)
  if (is.null(theo_rs))
    stop("[ERROR] ps_theoretical has no refseq() slot.")

  exp_seqs  <- as.character(exp_rs)
  exp_ids   <- taxa_names(ps_mocks)
  theo_seqs <- as.character(theo_rs)

  # Use taxa_names as strain labels (e.g. ESBL15, ESBL25, KPC_1).
  # Using the Species tax_table column would collapse all E. coli strains
  # into one "Escherichia coli" label, losing strain-level resolution.
  theo_species <- taxa_names(ps_theoretical)
  names(theo_species) <- taxa_names(ps_theoretical)
  message("[TAX] Strain labels: ", paste(theo_species, collapse = ", "))

  exp_dna    <- DNAStringSet(exp_seqs)
  theo_dna   <- DNAStringSet(theo_seqs)
  theo_names <- names(theo_seqs)

  message("[TAX] Aligning ", length(exp_dna), " experimental ASVs against ",
          length(theo_dna), " theoretical sequences...")

  assigned_species <- character(length(exp_ids))
  assigned_sim     <- numeric(length(exp_ids))

  for (i in seq_along(exp_ids)) {
    query  <- exp_dna[[i]]
    scores <- numeric(length(theo_dna))
    sims   <- numeric(length(theo_dna))

    for (j in seq_along(theo_dna)) {
      aln       <- pwalign::pairwiseAlignment(query, theo_dna[[j]],
                                              type = "global-local",
                                              scoreOnly = FALSE)
      scores[j] <- score(aln)
      # Correct similarity: nmatch / (nmatch + nmismatch)
      # avoids inflation by gap characters or length differences
      n_compared <- nmatch(aln) + nmismatch(aln)
      sims[j]    <- if (n_compared > 0) (nmatch(aln) / n_compared) * 100 else 0
    }

    best_j        <- which.max(scores)
    best_sim      <- sims[best_j]
    best_taxon_id <- theo_names[best_j]

    if (best_sim >= min_sim) {
      assigned_species[i] <- theo_species[best_taxon_id]
      message("  ", exp_ids[i], " -> ", assigned_species[i],
              " (", round(best_sim, 1), "% identity)")
    } else {
      # Unique nonMock:: label preserves each unmatched ASV as a distinct taxon
      assigned_species[i] <- paste0("nonMock::", exp_ids[i])
      message("  ", exp_ids[i], " -> non-mock (best match: ",
              round(best_sim, 1), "% < ", min_sim, "%) -- labelled nonMock::",
              exp_ids[i])
    }
    assigned_sim[i] <- best_sim
  }

  tax_mat <- matrix(
    NA_character_,
    nrow     = length(exp_ids),
    ncol     = 7,
    dimnames = list(exp_ids,
                    c("domain", "phylum", "class", "order",
                      "family", "genus", "species"))
  )
  tax_mat[, "species"] <- assigned_species
  tax_table(ps_mocks)  <- tax_table(tax_mat)

  n_matched <- sum(!startsWith(assigned_species, "nonMock::"))
  message("[TAX] Done. ", n_matched, " / ", length(exp_ids),
          " ASVs matched; ", length(exp_ids) - n_matched, " labelled nonMock.")
  ps_mocks
}


#' Build the theoretical reference as a single-sample compositional phyloseq.
#'
#' @param ps_theoretical  Theoretical mock phyloseq.
#' @return A single-sample phyloseq named "MockTheoretical" with % abundances.
build_theoretical_ps <- function(ps_theoretical) {

  message("[THEO] Building theoretical reference phyloseq...")

  species_labels <- taxa_names(ps_theoretical)
  names(species_labels) <- taxa_names(ps_theoretical)
  message("[THEO] Strain labels: ", paste(species_labels, collapse = ", "))

  otu_theo <- as(otu_table(ps_theoretical), "matrix")
  if (!taxa_are_rows(ps_theoretical)) otu_theo <- t(otu_theo)

  theo_df <- data.frame(species   = species_labels,
                        abundance = rowSums(otu_theo),
                        stringsAsFactors = FALSE)

  theo_agg <- theo_df %>%
    group_by(species) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    mutate(abundance = abundance / sum(abundance) * 100)

  otu_mat <- matrix(theo_agg$abundance, nrow = nrow(theo_agg), ncol = 1,
                    dimnames = list(theo_agg$species, "MockTheoretical"))

  tax_mat <- matrix(theo_agg$species, nrow = nrow(theo_agg), ncol = 7,
                    dimnames = list(theo_agg$species,
                                    c("domain","phylum","class","order",
                                      "family","genus","species")))

  sd_theo <- data.frame(ZymoType = "Theoretical", row.names = "MockTheoretical")

  ps_theo_out <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    tax_table(tax_mat),
    sample_data(sd_theo)
  )
  message("[THEO] Theoretical phyloseq: ", ntaxa(ps_theo_out),
          " species, 1 sample (MockTheoretical).")
  ps_theo_out
}


#' Run the full custom Check Mock analysis for one pipeline.
#'
#' Identical logic is applied to both BEST and STANDARD pipelines.
#' Steps:
#'   1. Assign taxonomy by sequence matching against ps_theoretical.
#'   2. Agglomerate to species level (manual, NA-safe); convert to %.
#'   3. Build a theoretical reference sample.
#'   4. Merge experimental + theoretical with merge_phyloseq.
#'   5. Compute Bray-Curtis similarity (NA cells treated as 0).
#'
#' Why manual aggregation instead of microbiome::aggregate_taxa():
#'   aggregate_taxa() relies on tax_glom() which silently discards taxa
#'   whose higher-rank columns are all NA — exactly what happens for
#'   nonMock:: ASVs and for strict-filter pipelines. Manual aggregation
#'   via the tax_table "species" column avoids that data loss.
#'
#' Why Bray-Curtis and not Spearman/Pearson:
#'   The theoretical mock has 4 strains at 25% each — all ranks are tied,
#'   so Spearman = 0. Bray-Curtis dissimilarity is appropriate for
#'   compositional data and handles equal abundances correctly.
#'
#' @param ps_mocks       Combined mock phyloseq (Mock Mix + Singular Mocks).
#' @param ps_theoretical Theoretical mock phyloseq.
#' @param min_sim        Min % identity for taxonomy assignment (default 100).
#' @return List: ps_asv, ps_species, corrTable.
run_custom_check_mock <- function(ps_mocks,
                                  ps_theoretical,
                                  min_sim = 100) {

  message("[RUN] Starting custom Check Mock analysis...")

  # Step 1: Assign taxonomy
  ps_asv <- assign_taxonomy_from_theoretical(ps_mocks, ps_theoretical, min_sim)
  diagnose_ps(ps_asv, "After taxonomy assignment (ASV level)")

  # Step 2: Agglomerate to species level (manual, NA-safe); convert to %
  message("[RUN] Agglomerating to species level (manual, NA-safe)...")

  otu_asv <- as(otu_table(ps_asv), "matrix")
  if (!taxa_are_rows(ps_asv)) otu_asv <- t(otu_asv)   # -> taxa x samples
  # Treat NA as 0: merge_phyloseq() fills absent taxa with NA rather than 0.
  # Without this replacement, colSums() and normalisation return NaN, causing
  # singular mock samples to appear empty in all plots.
  otu_asv[is.na(otu_asv)] <- 0

  sp_labels <- as.character(tax_table(ps_asv)[, "species"])
  names(sp_labels) <- taxa_names(ps_asv)

  unique_sp <- unique(sp_labels)
  agg_mat   <- matrix(0.0, nrow = length(unique_sp), ncol = ncol(otu_asv),
                      dimnames = list(unique_sp, colnames(otu_asv)))
  for (sp in unique_sp) {
    rows <- names(sp_labels)[sp_labels == sp]
    agg_mat[sp, ] <- if (length(rows) == 1)
                       otu_asv[rows, ]
                     else
                       colSums(otu_asv[rows, , drop = FALSE])
  }

  # Convert to compositional (%) per sample
  col_sums <- colSums(agg_mat)
  col_sums[col_sums == 0] <- 1   # guard against empty samples
  comp_mat <- sweep(agg_mat, 2, col_sums, "/") * 100

  tax_agg <- matrix(NA_character_, nrow = length(unique_sp), ncol = 7,
                    dimnames = list(unique_sp,
                                    c("domain","phylum","class","order",
                                      "family","genus","species")))
  tax_agg[, "species"] <- unique_sp

  ps_species_exp <- phyloseq(
    otu_table(comp_mat, taxa_are_rows = TRUE),
    tax_table(tax_agg),
    sample_data(ps_asv)
  )
  sample_data(ps_species_exp)$ZymoType <- "UserSamples"
  message("[RUN] Experimental: ", ntaxa(ps_species_exp), " species x ",
          nsamples(ps_species_exp), " samples.")

  # Step 3: Build theoretical reference sample
  ps_theo_ref      <- build_theoretical_ps(ps_theoretical)
  theo_sample_name <- sample_names(ps_theo_ref)
  message("[RUN] Theoretical sample name: '", theo_sample_name, "'")

  message("[RUN] Experimental taxa: ", paste(taxa_names(ps_species_exp), collapse=", "))
  message("[RUN] Theoretical taxa : ", paste(taxa_names(ps_theo_ref),  collapse=", "))

  # Step 4: Merge experimental + theoretical
  message("[RUN] Merging experimental and theoretical phyloseq objects...")
  ps_all <- merge_phyloseq(ps_species_exp, ps_theo_ref)
  message("[RUN] Merged: ", ntaxa(ps_all), " species x ",
          nsamples(ps_all), " samples.")
  message("[RUN] All sample names: ", paste(sample_names(ps_all), collapse=", "))

  # Step 5: Bray-Curtis similarity
  # merge_phyloseq fills cells with NA where a taxon is absent in one object.
  # Treating NA as 0 means nonMock:: taxa correctly contribute to dissimilarity.
  message("[RUN] Computing Bray-Curtis similarity with '", theo_sample_name, "'...")

  otu_mat <- as(otu_table(ps_all), "matrix")
  if (!taxa_are_rows(ps_all)) otu_mat <- t(otu_mat)

  message("[RUN] OTU matrix: ", nrow(otu_mat), " taxa x ", ncol(otu_mat), " samples.")
  message("[RUN] Taxa   : ", paste(rownames(otu_mat), collapse=", "))
  message("[RUN] Samples: ", paste(colnames(otu_mat), collapse=", "))

  if (!theo_sample_name %in% colnames(otu_mat))
    stop("[ERROR] Theoretical sample '", theo_sample_name, "' not found.\n",
         "  Columns present: ", paste(colnames(otu_mat), collapse=", "))

  theo_vec         <- otu_mat[, theo_sample_name]
  exp_sample_names <- setdiff(colnames(otu_mat), theo_sample_name)

  bc_similarity <- function(p, q) {
    p[is.na(p)] <- 0
    q[is.na(q)] <- 0
    denom <- sum(p + q)
    if (denom == 0) return(NA_real_)
    1 - sum(abs(p - q)) / denom
  }

  corr_vals <- sapply(exp_sample_names, function(s) {
    bc_similarity(otu_mat[, s], theo_vec)
  })

  for (s in exp_sample_names)
    message("  ", s, ": BC similarity = ", round(corr_vals[s], 3))

  sorted_exp   <- exp_sample_names[order(corr_vals, na.last = FALSE)]
  stable_order <- c(sorted_exp, theo_sample_name)

  tex_cor <- tibble::tibble(
    sample.chkmks   = c(exp_sample_names, theo_sample_name),
    ZymoTheoretical = c(corr_vals, 1)
  ) %>%
    dplyr::mutate(row.sample = factor(sample.chkmks, levels = stable_order))

  message("[RUN] Custom Check Mock analysis complete.")

  list(ps_asv     = ps_asv,
       ps_species = ps_all,
       corrTable  = tex_cor)
}


# ==============================================================================
# SECTION 6 — PREPARE PHYLOSEQ FOR chkMocks (reference / future use)
# ==============================================================================

#' Swap ASV ID taxa_names to DNA sequences.
#'
#' Only needed if using checkZymoBiomics() directly (not used in this workflow).
#'
#' @param ps  phyloseq object with refseq() slot.
#' @return Invisibly returns ps with taxa_names as DNA sequences.
prepare_ps_for_chkMocks <- function(ps) {
  message("[PREP] Checking taxa_names format...")
  tn <- taxa_names(ps)
  if (all(grepl("^[ACGTNacgtn]+$", tn))) {
    message("[PREP] taxa_names already contain DNA sequences -- no swap needed.")
    return(invisible(ps))
  }
  rs <- tryCatch(refseq(ps), error = function(e) NULL)
  if (is.null(rs)) stop("[ERROR] refseq() slot is empty.")
  if (!setequal(names(rs), tn)) stop("[ERROR] names(refseq(ps)) don't match taxa_names.")
  seqs_ordered <- as.character(rs[tn])
  ps_slots <- list(otu_table(ps), sample_data(ps))
  if (!is.null(access(ps, "tax_table", errorIfNULL = FALSE)))
    ps_slots <- c(ps_slots, list(tax_table(ps)))
  ps_no_ref <- do.call(merge_phyloseq, ps_slots)
  taxa_names(ps_no_ref) <- seqs_ordered
  new_rs <- DNAStringSet(seqs_ordered)
  names(new_rs) <- seqs_ordered
  ps <- merge_phyloseq(ps_no_ref, new_rs)
  message("[PREP] taxa_names swapped to DNA sequences. Object ready.")
  invisible(ps)
}


# ==============================================================================
# SECTION 7 — POST-PROCESSING & OUTPUT
# ==============================================================================

#' Rename corrTable columns for display.
#'
#' @param output  List from run_custom_check_mock().
#' @return Renamed tibble.
get_correlation_table <- function(output) {
  output$corrTable %>%
    dplyr::rename(
      MockSampleID         = sample.chkmks,
      BrayCurtisSimilarity = ZymoTheoretical,
      MockSampleID_ordered = row.sample
    )
}


#' Build a shared colour palette from one or more run_custom_check_mock outputs.
#'
#' Call once before plotting so the same taxon always gets the same colour
#' across BEST and STANDARD plots.
#'
#' @param ...  One or more lists returned by run_custom_check_mock().
#' @return Named character vector of hex colours.
build_shared_palette <- function(...) {

  outputs <- list(...)

  all_raw <- character(0)
  for (out in outputs) {
    sp_col  <- as.character(tax_table(out$ps_asv)[, "species"])
    all_raw <- union(all_raw, sp_col)
  }

  is_nonmock     <- startsWith(all_raw, "nonMock::")
  mock_raw       <- sort(unique(all_raw[!is_nonmock]))
  nonmock_raw    <- sort(unique(all_raw[is_nonmock]))
  mock_labels    <- mock_raw
  nonmock_labels <- sub("^nonMock::", "", nonmock_raw)
  all_labels     <- c(mock_labels, nonmock_labels)
  n_total        <- length(all_labels)

  make_colour_pool <- function(n_needed) {
    qpal_names <- c("Set1", "Dark2", "Set2", "Paired", "Set3", "Accent")
    pool <- unique(unlist(lapply(qpal_names, function(p) {
      max_n <- RColorBrewer::brewer.pal.info[p, "maxcolors"]
      RColorBrewer::brewer.pal(max_n, p)
    })))
    pool <- pool[!pool %in% c("#BEBEBE","#A6A6A6","#969696","#737373",
                               "#525252","#FFFFFF","#F7F7F7","#B3B3B3")]
    if (length(pool) < n_needed) {
      extra <- grDevices::hcl(
        h = seq(10, 350, length.out = n_needed - length(pool) + 2)[
              -c(1, n_needed - length(pool) + 2)],
        c = 75, l = 52)
      pool <- unique(c(pool, extra))
    }
    pool[seq_len(min(n_needed, length(pool)))]
  }

  colour_pool  <- make_colour_pool(n_total + 10)
  n_mock       <- length(mock_labels)
  n_non        <- length(nonmock_labels)
  mock_cols    <- colour_pool[seq_len(n_mock)]
  nonmock_cols <- colour_pool[n_mock + seq_len(n_non)]
  pal          <- setNames(c(mock_cols, nonmock_cols), all_labels)

  message("[PAL] Shared palette: ", n_mock, " mock taxa + ",
          n_non, " non-mock ASVs = ", n_total, " colours total.")
  pal
}


#' Plot stacked composition bars and Bray-Curtis similarity panel.
#'
#' Left panel: stacked barplot of ASV/taxon composition per sample (%).
#' Right panel: horizontal Bray-Curtis similarity bars vs. theoretical mock.
#' nonMock:: prefix is stripped from display labels.
#'
#' @param output      List from run_custom_check_mock().
#' @param shared_pal  Named colour vector from build_shared_palette().
#' @param title       Optional plot title string.
#' @return A patchwork ggplot object.
plot_check_mock <- function(output, shared_pal = NULL, title = NULL) {

  ps_sp   <- output$ps_species
  cor_dat <- output$corrTable

  # Long-format abundance table
  otu_mat <- as(otu_table(ps_sp), "matrix")
  if (!taxa_are_rows(ps_sp)) otu_mat <- t(otu_mat)
  otu_mat[is.na(otu_mat)] <- 0

  ldf <- as.data.frame(otu_mat) %>%
    tibble::rownames_to_column("species") %>%
    tidyr::pivot_longer(cols = -species, names_to = "sample", values_to = "abundance")

  # Sample display order from corrTable
  sample_order <- levels(cor_dat$row.sample)
  ldf$sample   <- factor(ldf$sample, levels = sample_order)

  # Correlation panel data
  cor_plot_dat <- cor_dat %>%
    dplyr::mutate(sample   = factor(sample.chkmks, levels = sample_order),
                  corr_val = ZymoTheoretical) %>%
    dplyr::filter(!is.na(sample))

  # Strip nonMock:: prefix for display labels
  ldf$species <- sub("^nonMock::", "", ldf$species)

  # Colour palette
  if (!is.null(shared_pal)) {
    pal <- shared_pal
    missing_labels <- setdiff(unique(ldf$species), names(pal))
    if (length(missing_labels) > 0) {
      warning("[plot_check_mock] Labels not in shared_pal (grey): ",
              paste(missing_labels, collapse = ", "))
      pal <- c(pal, setNames(rep("#CCCCCC", length(missing_labels)), missing_labels))
    }
  } else {
    all_sp <- unique(ldf$species)
    n_sp   <- length(all_sp)
    make_pool_local <- function(n) {
      qpal_names <- c("Set1", "Dark2", "Set2", "Paired", "Set3", "Accent")
      pool <- unique(unlist(lapply(qpal_names, function(p) {
        RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[p,"maxcolors"], p)
      })))
      pool <- pool[!pool %in% c("#BEBEBE","#A6A6A6","#969696","#737373",
                                 "#525252","#FFFFFF","#F7F7F7","#B3B3B3")]
      if (length(pool) < n) {
        extra <- grDevices::hcl(h = seq(10,350, length.out=n-length(pool)+2)[
                                      -c(1,n-length(pool)+2)], c=75, l=52)
        pool <- unique(c(pool, extra))
      }
      pool[seq_len(min(n, length(pool)))]
    }
    tax_sp_col <- as.character(tax_table(output$ps_asv)[, "species"])
    is_nm   <- startsWith(tax_sp_col, "nonMock::")
    mock_sp <- sort(unique(tax_sp_col[!is_nm]))
    nm_sp   <- sort(unique(sub("^nonMock::", "", tax_sp_col[is_nm])))
    mock_sp <- intersect(mock_sp, all_sp)
    nm_sp   <- intersect(nm_sp,   all_sp)
    pool    <- make_pool_local(n_sp + 10)
    pal     <- setNames(c(pool[seq_len(length(mock_sp))],
                          pool[length(mock_sp) + seq_len(length(nm_sp))]),
                        c(mock_sp, nm_sp))
    extra_sp <- setdiff(all_sp, names(pal))
    if (length(extra_sp) > 0)
      pal <- c(pal, setNames(pool[length(mock_sp)+length(nm_sp)+seq_len(length(extra_sp))],
                             extra_sp))
  }

  # Left panel: stacked composition bars
  p_bar <- ggplot(ldf, aes(x = abundance, y = sample, fill = species)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = pal) +
    labs(x = "Abundance (%)", y = "Samples", fill = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right", axis.text.y = element_text(size = 9))

  # Right panel: BC similarity bars
  cor_plot_dat <- cor_plot_dat %>%
    dplyr::mutate(
      corr_display = dplyr::coalesce(corr_val, 0),
      corr_label   = dplyr::if_else(is.na(corr_val), "NA", sprintf("%.2f", corr_val))
    )

  p_cor <- ggplot(cor_plot_dat, aes(x = corr_display, y = sample)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = corr_label), hjust = -0.1, size = 3, colour = "grey20") +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    labs(x = "Bray-Curtis Similarity\nwith Theoretical", y = NULL) +
    xlim(0, 1.15) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  # Combine panels
  p_combined <- p_bar + p_cor + plot_layout(widths = c(3, 1), guides = "collect")

  if (!is.null(title)) {
    p_combined <- p_combined +
      plot_annotation(
        title = title,
        theme = theme(plot.title = element_text(size = 15, face = "bold",
                                                hjust = 0.5, margin = margin(b = 8)))
      )
  }

  p_combined
}


#' Save a figure to PNG (300 dpi) and PDF.
#'
#' @param plot_obj  ggplot or patchwork object.
#' @param base_name Base filename (no extension).
#' @param out_dir   Output directory (created if it does not exist).
#' @param width     Width in inches (default 12).
#' @param height    Height in inches (default 7).
#' @return Invisibly returns a list with png and pdf paths.
save_figure <- function(plot_obj, base_name, out_dir = "results/05_mock_qc",
                        width = 12, height = 7) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    message("[SAVE] Created directory: ", out_dir)
  }
  png_path <- file.path(out_dir, paste0(base_name, ".png"))
  pdf_path <- file.path(out_dir, paste0(base_name, ".pdf"))
  ggsave(png_path, plot = plot_obj, width = width, height = height,
         dpi = 300, bg = "white")
  ggsave(pdf_path, plot = plot_obj, width = width, height = height)
  message("[SAVE] PNG -> ", png_path)
  message("[SAVE] PDF -> ", pdf_path)
  invisible(list(png = png_path, pdf = pdf_path))
}
