# ==============================================================================
# 04_functions.R
#
# Project : fliC Amplicon Sequencing for High-Resolution Profiling of E. coli
#           in the Infant Gut Microbiome
# Author  : Andrés Catalán Tatay
# Part    : Analysis — Step 04 (companion functions)
# Purpose : All helper functions for 04_BEST_vs_STANDARD.Rmd.
#           Covers phyloseq loading, sequence extraction, edit distance
#           computation, FASTA writing, BLAST JSON parsing, LCA computation,
#           trust tier assignment, and all three figures.
#
# Inputs  : data/cases/              — 96 RDS files (DADA2 outputs)
#            data/wgs/               — WGS phyloseq reference
#            data/blast/JSON/        — BLAST Single-file JSON results
# Outputs : Called by 04_BEST_vs_STANDARD.Rmd
#
# Usage   : source("04_functions.R"); load_libraries()
# ==============================================================================


# ==============================================================================
# SECTION 0 — LIBRARY LOADER
# ==============================================================================

load_libraries <- function() {
  suppressPackageStartupMessages({
    library(phyloseq)
    library(DECIPHER)
    library(Biostrings)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(ggtree)
    library(ggtreeExtra)
    library(ggExtra)
    library(reshape2)
    library(stringr)
    library(ape)
    library(phangorn)
    library(RColorBrewer)
    library(tibble)
    library(purrr)
    library(jsonlite)   # JSON parsing for BLAST Single-file JSON output
  })
}


# ==============================================================================
# SECTION 1 — LOAD PHYLOSEQ FROM CASE LIST
# ==============================================================================

#' Load a phyloseq object stored inside a named list (.rds file).
#'
#' @param path  Full path to the .rds file.
#' @return phyloseq object.
load_case_ps <- function(path) {
  obj <- readRDS(path)
  if (inherits(obj, "phyloseq")) return(obj)
  if (is.list(obj)) {
    ps_idx <- which(sapply(obj, inherits, what = "phyloseq"))
    if (length(ps_idx) == 0) stop("No phyloseq object found inside list: ", path)
    return(obj[[ps_idx[1]]])
  }
  stop("Unexpected object class: ", class(obj))
}


# ==============================================================================
# SECTION 1A — LOAD WGS PHYLOSEQ
# ==============================================================================

#' Load a named phyloseq object from an .RDS file.
#'
#' @param path         Full path to the .RDS file.
#' @param object_name  Name of the object to extract (default "ps_wgs").
#' @return phyloseq object.
load_wgs_ps <- function(path, object_name = "ps_wgs") {
  if (!file.exists(path)) stop("File does not exist: ", path)
  env <- new.env(parent = emptyenv())
  loaded_objects <- load(path, envir = env)
  message("Objects loaded: ", paste(loaded_objects, collapse = ", "))
  if (!object_name %in% loaded_objects)
    stop("Object '", object_name, "' not found. Available objects: ",
         paste(loaded_objects, collapse = ", "))
  ps <- get(object_name, envir = env)
  if (!inherits(ps, "phyloseq"))
    stop("Object '", object_name, "' is not a phyloseq object.")
  return(ps)
}


# ==============================================================================
# SECTION 2 — EXTRACT ASV SEQUENCES
# ==============================================================================

#' Extract ASV sequences from a phyloseq object.
#'
#' Uses refseq() if available, otherwise interprets taxa_names as sequences.
#'
#' @param ps  phyloseq object.
#' @return DNAStringSet of ASV sequences.
get_asv_sequences <- function(ps) {
  if (!is.null(refseq(ps))) {
    seqs <- refseq(ps)
    if (is.character(seqs)) seqs <- DNAStringSet(seqs)
    return(seqs)
  }
  taxa <- taxa_names(ps)
  if (all(grepl("^ASV", taxa, ignore.case = TRUE)))
    stop("taxa_names are ASV labels, not sequences, and refseq() is NULL.")
  DNAStringSet(setNames(taxa, taxa))
}


# ==============================================================================
# SECTION 3 — MERGE PHYLOSEQ OBJECTS
# ==============================================================================

#' Merge two phyloseq objects.
#'
#' @param ps_best  Best pipeline phyloseq.
#' @param ps_std   Standard pipeline phyloseq.
#' @return Merged phyloseq object.
merge_two_ps <- function(ps_best, ps_std) {
  merge_phyloseq(ps_best, ps_std)
}


# ==============================================================================
# SECTION 4 — COMPUTE PAIRWISE EDIT DISTANCE MATRIX
# ==============================================================================

#' Align sequences with DECIPHER::AlignSeqs then compute pairwise Hamming
#' distances using pwalign::stringDist(method = "hamming").
#'
#' @param seqs  DNAStringSet to align and compare.
#' @param ...   Additional arguments passed to AlignSeqs.
#' @return dist object of raw integer bp edit distances.
compute_edit_dist <- function(seqs, ...) {
  if (!requireNamespace("pwalign", quietly = TRUE))
    stop("pwalign required: BiocManager::install('pwalign')")
  seqs_aligned <- AlignSeqs(seqs, verbose = FALSE, ...)
  as.dist(as.matrix(pwalign::stringDist(seqs_aligned, method = "hamming")))
}


# ==============================================================================
# SECTION 5 — EDIT DISTANCES TO CLOSEST MATCH
# ==============================================================================

#' For each ASV in query_names, find the minimum edit distance to any ASV
#' in ref_names, explicitly excluding self-comparisons.
#'
#' @param dm_bp         dist object OR matrix (bp-level edit distances).
#' @param query_names   Character vector — ASVs to query.
#' @param ref_names     Character vector — ASVs to compare against.
#' @param pipeline_label  String label for the 'pipeline' output column.
#' @return data.frame: asv, min_edit_dist, pipeline.
closest_match_dist <- function(dm_bp, query_names, ref_names, pipeline_label) {

  empty_df <- data.frame(asv = character(), min_edit_dist = numeric(),
                         pipeline = character(), stringsAsFactors = FALSE)

  dm_mat <- tryCatch(as.matrix(dm_bp), error = function(e) NULL)
  if (is.null(dm_mat) || nrow(dm_mat) == 0) return(empty_df)

  q_present   <- intersect(query_names, rownames(dm_mat))
  ref_present <- intersect(ref_names,   colnames(dm_mat))

  if (length(q_present) == 0 || length(ref_present) == 0) return(empty_df)

  sub_mat <- dm_mat[q_present, ref_present, drop = FALSE]

  # Self-exclusion: set diagonal cells to NA
  common <- intersect(q_present, ref_present)
  if (length(common) > 0)
    for (nm in common) sub_mat[nm, nm] <- NA_real_

  results <- lapply(q_present, function(q) {
    row_vals   <- sub_mat[q, , drop = TRUE]
    valid_vals <- row_vals[!is.na(row_vals)]
    if (length(valid_vals) == 0)
      return(data.frame(asv = q, min_edit_dist = NA_real_,
                        pipeline = pipeline_label, stringsAsFactors = FALSE))
    data.frame(asv = q, min_edit_dist = min(valid_vals),
               pipeline = pipeline_label, stringsAsFactors = FALSE)
  })

  out <- bind_rows(results)
  if (nrow(out) == 0) return(empty_df)
  out
}


# ==============================================================================
# SECTION 6 — FIGURE 1: EDIT DISTANCE DISTRIBUTIONS
# ==============================================================================

#' Three-panel plot of edit distances to the closest match:
#'   Panel 1: Within Best pipeline
#'   Panel 2: Within Standard pipeline
#'   Panel 3: Cross-Pipeline (Best ASVs vs Standard ASVs)
#'
#' Axes use a positional encoding so that Zone A (0–5 bp) occupies ~35% of
#' the x-axis width and larger distances are progressively compressed.
#'
#' @param dm_bp     dist object of pairwise edit distances.
#' @param best_asvs Character vector of Best pipeline ASV names.
#' @param std_asvs  Character vector of Standard pipeline ASV names.
#' @param colors    Named list with elements best, std, cross (hex colours).
#' @return ggplot object.
figure1_edit_dist <- function(dm_bp,
                               best_asvs, std_asvs,
                               colors = list(best  = "#4472C4",
                                             std   = "#70AD47",
                                             cross = "#ED7D31")) {

  df_best  <- closest_match_dist(dm_bp, best_asvs, best_asvs, "1. Within Best")
  df_best$fill_color <- colors$best

  df_std   <- closest_match_dist(dm_bp, std_asvs,  std_asvs,  "2. Within Standard")
  df_std$fill_color  <- colors$std

  df_cross_b <- closest_match_dist(dm_bp, best_asvs, std_asvs, "3. Cross-Pipeline")
  df_cross_s <- closest_match_dist(dm_bp, std_asvs,  best_asvs,"3. Cross-Pipeline")
  df_cross   <- bind_rows(df_cross_b, df_cross_s)
  df_cross$fill_color <- colors$cross

  df_all <- bind_rows(df_best, df_std, df_cross) %>%
    filter(!is.na(min_edit_dist)) %>%
    mutate(
      pipeline = factor(pipeline, levels = c("1. Within Best",
                                             "2. Within Standard",
                                             "3. Cross-Pipeline")),
      min_edit_dist = as.integer(min_edit_dist)
    )

  max_dist <- max(df_all$min_edit_dist, na.rm = TRUE)

  # ── Compressed x-axis with four zones ─────────────────────────────────────
  # Zone A (0–5 bp):    one position per bp — full resolution
  # Zone B (6–20 bp):   one position per 2 bp — moderate compression
  # Zone C (21–100 bp): one position per 10 bp — strong compression
  # Zone D (>100 bp):   one position per 100 bp — maximum compression

  zone_a <- data.frame(
    bin_lo = 0:5, bin_hi = 0:5, centre = 0:5, pos = 0:5,
    label  = as.character(0:5), bar_w = 0.85, stringsAsFactors = FALSE
  )

  zone_b_starts <- seq(6, min(20, max_dist), by = 2)
  zone_b_ends   <- pmin(zone_b_starts + 1, max_dist)
  zone_b <- data.frame(
    bin_lo = zone_b_starts, bin_hi = zone_b_ends,
    centre = round((zone_b_starts + zone_b_ends) / 2),
    pos    = max(zone_a$pos) + seq_along(zone_b_starts),
    label  = as.character(zone_b_starts), bar_w = 0.85, stringsAsFactors = FALSE
  )

  if (max_dist > 20) {
    zone_c_starts <- seq(21, min(100, max_dist), by = 10)
    zone_c_ends   <- pmin(zone_c_starts + 9, max_dist)
    zone_c <- data.frame(
      bin_lo = zone_c_starts, bin_hi = zone_c_ends,
      centre = round((zone_c_starts + zone_c_ends) / 2),
      pos    = max(zone_b$pos) + seq_along(zone_c_starts),
      label  = as.character(zone_c_starts), bar_w = 0.85, stringsAsFactors = FALSE
    )
  } else {
    zone_c <- data.frame(bin_lo = integer(), bin_hi = integer(), centre = integer(),
                         pos = integer(), label = character(), bar_w = numeric(),
                         stringsAsFactors = FALSE)
  }

  if (max_dist > 100) {
    zone_d_starts <- seq(101, max_dist, by = 100)
    zone_d_ends   <- pmin(zone_d_starts + 49, max_dist)
    prev_max_pos  <- if (nrow(zone_c) > 0) max(zone_c$pos) else max(zone_b$pos)
    zone_d <- data.frame(
      bin_lo = zone_d_starts, bin_hi = zone_d_ends,
      centre = round((zone_d_starts + zone_d_ends) / 2),
      pos    = prev_max_pos + seq_along(zone_d_starts),
      label  = as.character(zone_d_starts), bar_w = 0.85, stringsAsFactors = FALSE
    )
  } else {
    zone_d <- data.frame(bin_lo = integer(), bin_hi = integer(), centre = integer(),
                         pos = integer(), label = character(), bar_w = numeric(),
                         stringsAsFactors = FALSE)
  }

  axis_map <- bind_rows(zone_a, zone_b, zone_c, zone_d)

  assign_pos <- function(d) {
    idx <- findInterval(d, axis_map$bin_lo, rightmost.closed = TRUE)
    idx <- pmax(1L, pmin(idx, nrow(axis_map)))
    axis_map$pos[idx]
  }

  df_all <- df_all %>% mutate(pos = assign_pos(min_edit_dist))

  df_binned <- df_all %>%
    group_by(pipeline, pos, fill_color) %>%
    summarise(count = n(), .groups = "drop") %>%
    left_join(axis_map %>% select(pos, bar_w, label), by = "pos")

  y_max         <- max(df_binned$count, na.rm = TRUE)
  break_positions <- axis_map$pos[axis_map$label != ""]
  break_labels    <- axis_map$label[axis_map$label != ""]
  zone_boundary   <- max(zone_a$pos) + 0.5   # separates 0–5 bp from compressed region

  p <- ggplot(df_binned,
              aes(x = pos, y = count, fill = fill_color, width = bar_w)) +
    geom_col(color = NA) +
    geom_vline(xintercept = zone_boundary,
               linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    facet_wrap(~ pipeline, ncol = 1, scales = "fixed", strip.position = "top") +
    scale_fill_identity() +
    scale_x_continuous(breaks = break_positions, labels = break_labels,
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(limits = c(0, y_max * 1.05),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title    = "Distribution of Edit Distances to Closest Match",
         subtitle = "Comparing Within-Pipeline diversity vs Cross-Pipeline artifacts",
         x        = "Edit Distance (bp)",
         y        = "Count of ASVs") +
    theme_bw(base_size = 12) +
    theme(strip.background   = element_rect(fill = "black"),
          strip.text         = element_text(color = "white", face = "bold", size = 11),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title         = element_text(face = "bold"),
          legend.position    = "none",
          panel.spacing      = unit(0.3, "lines"),
          axis.text.x        = element_text(size = 8, angle = 45, hjust = 1))
  p
}


# ==============================================================================
# SECTION 7 — WRITE FASTA
# ==============================================================================

#' Write sequences to FASTA, optionally splitting into n_parts balanced files.
#'
#' @param seqs      DNAStringSet to write.
#' @param out_path  Output file path.
#' @param chunk_size  Legacy splitting by chunk size (NULL = disabled).
#' @param n_parts     Split into this many balanced files (NULL = single file).
#' @return Invisibly: path(s) to the written FASTA file(s).
write_fasta <- function(seqs, out_path, chunk_size = NULL, n_parts = NULL) {

  if (!requireNamespace("Biostrings", quietly = TRUE))
    stop("Package 'Biostrings' is required for write_fasta().")

  out_dir <- dirname(out_path)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Split into n_parts balanced FASTA files
  if (!is.null(n_parts) && n_parts > 1) {
    n_parts  <- as.integer(n_parts)
    n_seqs   <- length(seqs)
    if (n_seqs == 0) stop("No sequences found in 'seqs'.")
    if (n_parts > n_seqs)
      warning("n_parts > number of sequences. Some files may be empty.")
    part_idx <- split(seq_len(n_seqs),
                      cut(seq_len(n_seqs), breaks = n_parts, labels = FALSE))
    part_paths <- file.path(out_dir, sprintf("blast_part_%02d.fasta", seq_len(n_parts)))
    message("Splitting ", n_seqs, " sequences into ", n_parts, " FASTA files:")
    for (i in seq_len(n_parts)) {
      idx       <- part_idx[[as.character(i)]]
      part_seqs <- if (is.null(idx)) seqs[0] else seqs[idx]
      Biostrings::writeXStringSet(part_seqs, filepath = part_paths[i])
      message("  Wrote ", basename(part_paths[i]), " | ", length(part_seqs), " sequences")
    }
    missing_files <- part_paths[!file.exists(part_paths)]
    if (length(missing_files) > 0)
      stop("Some FASTA parts not created:\n", paste(missing_files, collapse = "\n"))
    message("Done. ", length(part_paths), "/", n_parts, " FASTA files created.")
    return(invisible(part_paths))
  }

  # Legacy chunk_size splitting
  if (!is.null(chunk_size) && length(seqs) > chunk_size) {
    n_chunks    <- ceiling(length(seqs) / chunk_size)
    base        <- tools::file_path_sans_ext(out_path)
    ext         <- tools::file_ext(out_path)
    chunk_paths <- character(n_chunks)
    for (i in seq_len(n_chunks)) {
      idx <- ((i - 1) * chunk_size + 1):min(i * chunk_size, length(seqs))
      chunk_paths[i] <- sprintf("%s_chunk%02d.%s", base, i, ext)
      Biostrings::writeXStringSet(seqs[idx], filepath = chunk_paths[i])
    }
    message("Wrote ", n_chunks, " FASTA chunks to ", out_dir)
    return(invisible(chunk_paths))
  }

  # Single FASTA file
  Biostrings::writeXStringSet(seqs, filepath = out_path)
  message("Wrote single FASTA: ", out_path, " | ", length(seqs), " sequences")
  return(invisible(out_path))
}


# ==============================================================================
# SECTION 8 — LOAD BLAST RESULTS FROM SINGLE-FILE JSON
# ==============================================================================

# BLAST results are downloaded from NCBI as "Single-file JSON" (one JSON per
# FASTA chunk). All JSON files in file.path(BLAST_DIR, "JSON") are parsed and
# combined into a single data.frame.
#
# JSON hierarchy:
#   BlastOutput2[i]
#     .report.results.search
#       .query_title   — ASV name
#       .query_len     — query sequence length
#       .hits[j]
#         .description[0].accession / taxid / sciname / title
#         .hsps[k]
#           .bit_score / .evalue / .identity / .align_len / .gaps
#           .query_from / .query_to / .hit_from / .hit_to
#
# query_coverage = (query_to - query_from + 1) / query_len * 100
#
# Output columns match the standard BLAST tabular format:
#   qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend,
#   sstart, send, evalue, bitscore, query_cover, staxids, sscinames, stitle

#' Load and combine all BLAST JSON files from a directory.
#'
#' @param blast_dir         Root directory; JSON files in blast_dir/JSON/.
#' @param max_hits_per_asv  Max hits to retain per query (default 50).
#' @return data.frame, one row per hit.
load_blast_results <- function(blast_dir, max_hits_per_asv = 50) {

  JSON_DIR <- file.path(blast_dir, "JSON")

  if (!dir.exists(JSON_DIR))
    stop("JSON subdirectory not found: ", JSON_DIR,
         "\nExpected: file.path(BLAST_DIR, 'JSON')  containing blast_part_NN.json files.")

  json_files <- sort(list.files(JSON_DIR, pattern = "\\.json$",
                                full.names = TRUE, recursive = FALSE))

  if (length(json_files) == 0)
    stop("No .json files found in: ", JSON_DIR)

  message(">>> JSON BLAST directory: ", JSON_DIR)
  message(">>> JSON files found: ", length(json_files))

  all_parts <- lapply(json_files, function(f) {
    message("  Parsing: ", basename(f))
    part_df <- .parse_single_blast_json(f, max_hits_per_asv = max_hits_per_asv)
    if (is.null(part_df)) {
      message("    SKIPPED (parse error or 0 rows)")
      return(NULL)
    }
    message("    Rows: ", nrow(part_df),
            " | Unique ASVs: ", dplyr::n_distinct(part_df$qseqid))
    part_df
  })

  successful <- Filter(Negate(is.null), all_parts)
  if (length(successful) == 0)
    stop("All JSON files failed to parse.")

  blast_df <- dplyr::bind_rows(successful)

  # Warn about duplicated ASV IDs across files
  qseqid_per_file <- lapply(successful, function(d) unique(d$qseqid))
  seen <- dups <- character(0)
  for (ids in qseqid_per_file) {
    dups <- c(dups, intersect(ids, seen))
    seen <- c(seen, ids)
  }
  if (length(dups) > 0)
    warning("ASV IDs found in >1 JSON file (possible duplicates):\n  ",
            paste(unique(dups), collapse = ", "))

  blast_df <- .ensure_blast_columns(blast_df)
  .blast_load_summary(blast_df)
  blast_df
}


# ── Internal: parse a single BLAST JSON file ──────────────────────────────────
.parse_single_blast_json <- function(f, max_hits_per_asv = 50) {

  .n <- function(x) suppressWarnings(as.numeric(x))
  .i <- function(x) suppressWarnings(as.integer(x))

  raw <- tryCatch(jsonlite::fromJSON(f, simplifyVector = FALSE),
                  error = function(e) {
                    warning("Cannot parse JSON: ", basename(f), " — ", e$message)
                    NULL
                  })
  if (is.null(raw)) return(NULL)

  reports <- raw[["BlastOutput2"]]
  if (is.null(reports)) {
    warning("No 'BlastOutput2' key in: ", basename(f))
    return(NULL)
  }

  no_hit_queries <- character(0)

  rows <- lapply(reports, function(report_obj) {

    search <- tryCatch(
      report_obj[["report"]][["results"]][["search"]],
      error = function(e) NULL
    )
    if (is.null(search)) return(NULL)

    query_title <- trimws(as.character(search[["query_title"]] %||% "UNKNOWN_QUERY"))
    if (!nzchar(query_title)) query_title <- "UNKNOWN_QUERY"
    query_len   <- .i(search[["query_len"]])

    hits <- search[["hits"]]
    if (is.null(hits) || length(hits) == 0) {
      no_hit_queries <<- c(no_hit_queries, query_title)
      return(NULL)
    }

    hit_rows <- lapply(hits, function(hit) {

      descr_list <- hit[["description"]]
      if (is.null(descr_list) || length(descr_list) == 0) return(NULL)
      descr     <- descr_list[[1]]

      accession <- as.character(descr[["accession"]] %||% NA_character_)
      staxids   <- as.character(descr[["taxid"]]     %||% NA_character_)
      sscinames <- as.character(descr[["sciname"]]   %||% NA_character_)
      stitle    <- as.character(descr[["title"]]     %||% NA_character_)

      # Fallback accession from 'id' field
      if (is.na(accession) || !nzchar(accession)) {
        id_raw    <- as.character(descr[["id"]] %||% "")
        accession <- gsub(".*\\|([^|]+)\\|?$", "\\1", id_raw)
        accession <- sub("\\.\\d+$", "", accession)   # strip version
      }

      hsps <- hit[["hsps"]]
      if (is.null(hsps) || length(hsps) == 0) return(NULL)

      hsp_df <- tryCatch({
        do.call(rbind, lapply(hsps, function(hsp) {
          al    <- .i(hsp[["align_len"]])
          iden  <- .n(hsp[["identity"]])
          gaps  <- .i(hsp[["gaps"]])
          pi    <- if (!is.na(iden) && !is.na(al) && al > 0)
                     round(iden / al * 100, 3) else NA_real_
          mm    <- if (!is.na(al) && !is.na(iden) && !is.na(gaps))
                     as.integer(al - iden - gaps) else NA_integer_
          qfrom <- .i(hsp[["query_from"]])
          qto   <- .i(hsp[["query_to"]])
          qcov  <- if (!is.na(query_len) && query_len > 0 &&
                       !is.na(qfrom) && !is.na(qto))
                     round((qto - qfrom + 1L) / query_len * 100, 2) else NA_real_
          data.frame(
            pident   = pi,   length   = al,   mismatch = mm,
            gapopen  = gaps, qstart   = qfrom, qend     = qto,
            sstart   = .i(hsp[["hit_from"]]), send = .i(hsp[["hit_to"]]),
            evalue   = .n(hsp[["evalue"]]),  bitscore = .n(hsp[["bit_score"]]),
            query_cover = qcov, stringsAsFactors = FALSE
          )
        }))
      }, error = function(e) NULL)

      if (is.null(hsp_df) || nrow(hsp_df) == 0) return(NULL)

      # Best HSP only (highest bitscore)
      best_hsp <- hsp_df[which.max(hsp_df$bitscore), , drop = FALSE]

      data.frame(
        qseqid    = query_title,
        sseqid    = accession,
        pident    = best_hsp$pident,
        length    = best_hsp$length,
        mismatch  = best_hsp$mismatch,
        gapopen   = best_hsp$gapopen,
        qstart    = best_hsp$qstart,
        qend      = best_hsp$qend,
        sstart    = best_hsp$sstart,
        send      = best_hsp$send,
        evalue    = best_hsp$evalue,
        bitscore  = best_hsp$bitscore,
        query_cover = best_hsp$query_cover,
        staxids   = if (is.na(staxids)   || !nzchar(staxids))   NA_character_ else staxids,
        sscinames = if (is.na(sscinames) || !nzchar(sscinames)) NA_character_ else sscinames,
        stitle    = if (is.na(stitle)    || !nzchar(stitle))    NA_character_ else stitle,
        stringsAsFactors = FALSE
      )
    })

    hit_rows_clean <- Filter(Negate(is.null), hit_rows)
    if (length(hit_rows_clean) == 0) return(NULL)

    query_hits <- dplyr::bind_rows(hit_rows_clean)

    # Limit hits per ASV: keep top max_hits_per_asv by bitscore
    if (nrow(query_hits) > max_hits_per_asv) {
      query_hits <- query_hits %>%
        dplyr::arrange(dplyr::desc(bitscore)) %>%
        dplyr::slice_head(n = max_hits_per_asv)
    }

    query_hits
  })

  if (length(no_hit_queries) > 0)
    message("    No-hit queries: ", length(no_hit_queries),
            " (e.g., ", paste(head(no_hit_queries, 3), collapse = ", "), ")")

  result <- dplyr::bind_rows(Filter(Negate(is.null), rows))
  if (nrow(result) == 0) {
    warning("0 rows parsed from: ", basename(f))
    return(NULL)
  }
  result
}

# Null-coalescing helper
`%||%` <- function(a, b) if (!is.null(a)) a else b


# ── Internal: guarantee all required BLAST columns exist ──────────────────────
.ensure_blast_columns <- function(df) {
  required <- list(
    qseqid      = NA_character_, sseqid    = NA_character_,
    pident      = NA_real_,      length    = NA_integer_,
    mismatch    = NA_integer_,   gapopen   = NA_integer_,
    qstart      = NA_integer_,   qend      = NA_integer_,
    sstart      = NA_integer_,   send      = NA_integer_,
    evalue      = NA_real_,      bitscore  = NA_real_,
    query_cover = NA_real_,
    staxids     = NA_character_, sscinames = NA_character_,
    stitle      = NA_character_
  )
  for (col in names(required)) {
    if (!col %in% colnames(df)) {
      df[[col]] <- required[[col]]
      message("  Added missing column '", col, "' as NA.")
    }
  }
  df$pident      <- as.numeric(df$pident)
  df$query_cover <- as.numeric(df$query_cover)
  df$bitscore    <- as.numeric(df$bitscore)
  df$evalue      <- as.numeric(df$evalue)
  df
}


# ── Internal: print a standard BLAST load summary ─────────────────────────────
.blast_load_summary <- function(df) {
  qcov_present <- !all(is.na(df$query_cover))
  tax_present  <- !all(is.na(df$sscinames))
  cat("\n── BLAST Load Summary ──────────────────────────────────────────\n")
  cat("Total rows              :", nrow(df), "\n")
  cat("Unique qseqid (ASVs)    :", dplyr::n_distinct(df$qseqid), "\n")
  cat("query_cover available   :", if (qcov_present) "YES" else "NO (all NA)", "\n")
  cat("sscinames available     :", if (tax_present)  "YES" else "NO (all NA)", "\n")
  cat("pident range            :", round(min(df$pident, na.rm = TRUE), 1), "–",
      round(max(df$pident, na.rm = TRUE), 1), "\n")
  if (qcov_present)
    cat("query_cover range       :", round(min(df$query_cover, na.rm = TRUE), 1), "–",
        round(max(df$query_cover, na.rm = TRUE), 1), "\n")
  cat("Columns                 :", paste(colnames(df), collapse = ", "), "\n")
  cat("────────────────────────────────────────────────────────────────\n\n")
  if (all(is.na(df$pident)))
    stop("BLAST error: pident is all NA.")
  if (all(is.na(df$sscinames)))
    stop("BLAST error: sscinames all NA. Ensure JSON files contain sciname fields.")
  if (!"qseqid" %in% colnames(df))
    stop("BLAST error: qseqid column missing.")
}


# ==============================================================================
# SECTION 9 — COMPUTE LCA FROM BLAST
# ==============================================================================

# Strategy:
#   1. Filter hits to top_n by bitscore AND min_pident threshold
#   2. For each query, check genus-level consensus among qualifying hits
#   3. LCA is the shared genus name if all hits agree; "ambiguous" otherwise
#   4. All queries from blast_df are preserved (lca_name = NA, n_hits = 0
#      if their hits fail the pident filter)
#
# NOTE: top_pident is NOT returned here — it is computed from blast_df in the
# Rmd and stored in asv_base to avoid column collision in assign_trust_tiers().

#' Compute genus-level LCA from BLAST results.
#'
#' @param blast_df    data.frame from load_blast_results().
#' @param top_n       Max hits per query to consider (default 10).
#' @param min_pident  Minimum %identity threshold (default 97).
#' @return data.frame(qseqid, lca_name, n_hits).
compute_lca <- function(blast_df,
                        top_n      = 10,
                        min_pident = 97) {

  all_queries <- dplyr::distinct(blast_df, qseqid)

  lca_result <- blast_df %>%
    dplyr::mutate(
      pident    = as.numeric(pident),
      bitscore  = as.numeric(bitscore),
      sscinames = as.character(sscinames)
    ) %>%
    dplyr::filter(!is.na(pident), pident >= min_pident) %>%
    dplyr::group_by(qseqid) %>%
    dplyr::arrange(dplyr::desc(bitscore), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::summarise(
      lca_name = {
        names_vec       <- sscinames
        names_vec       <- names_vec[!is.na(names_vec) & names_vec != ""]
        if (length(names_vec) == 0) {
          NA_character_
        } else {
          # Genus-level LCA: if all top hits share the same genus → use it
          genera        <- sapply(strsplit(names_vec, " "), `[`, 1)
          unique_genera <- unique(na.omit(genera))
          if (length(unique_genera) == 1) unique_genera[1] else "ambiguous"
        }
      },
      n_hits = dplyr::n(),
      .groups = "drop"
    )

  # Restore ASVs whose hits were all below min_pident
  lca_result <- all_queries %>%
    dplyr::left_join(lca_result, by = "qseqid") %>%
    dplyr::mutate(
      lca_name = dplyr::if_else(is.na(lca_name), NA_character_, lca_name),
      n_hits   = dplyr::if_else(is.na(n_hits),   0L, as.integer(n_hits))
    )

  lca_result
}


# ==============================================================================
# SECTION 10 — ASSIGN TRUST TIERS
# ==============================================================================

# 5-tier system determined by three independent evidence sources:
#   1. BLAST identity & coverage (NCBI-derived)
#   2. LCA taxonomy
#   3. WGS confirmation (edit distance = 0)
#
# Tier definitions:
#   Tier 1: All hits are 100% id + 100% cov E. coli  OR  WGS exact match
#   Tier 2: At least one hit is 100% id + 100% cov E. coli (but not all)
#   Tier 3: No perfect hits, all hits are E. coli
#   Tier 4: No perfect hits, mixed E. coli and non-E. coli hits
#   Tier 5: No E. coli hits at all
#
# Dendrogram grouping:
#   Tiers 1–2 → "E. coli"
#   Tiers 3–4 → "Ambiguous"
#   Tier  5   → "Non E. coli"

#' Assign BLAST-based trust tiers to all ASVs.
#'
#' @param asv_df       data.frame: asv, top_pident, top_qcov, wgs_confirmed.
#' @param wgs_ps       WGS phyloseq (passed for signature compatibility).
#' @param lca_df       data.frame: qseqid, lca_name, n_hits (from compute_lca).
#' @param blast_df     Full BLAST result data.frame.
#' @param target_genus Genus name to test against (default "Escherichia").
#' @return Input data.frame with added columns: trust_tier, trust_label, trust_group.
assign_trust_tiers <- function(asv_df,
                                wgs_ps,
                                lca_df,
                                blast_df         = NULL,
                                target_genus     = "Escherichia") {

  # Guards
  for (col in c("asv", "top_pident", "wgs_confirmed")) {
    if (!col %in% colnames(asv_df))
      stop("assign_trust_tiers(): '", col, "' missing from asv_df.")
  }
  if (!"top_qcov" %in% colnames(asv_df)) {
    message("assign_trust_tiers(): 'top_qcov' not in asv_df — setting to 0.")
    asv_df$top_qcov <- 0
  }

  # Remove any top_pident/top_qcov that might exist in lca_df
  lca_df <- lca_df %>%
    dplyr::select(-dplyr::any_of(c("top_pident", "top_qcov")))

  for (req_col in c("qseqid", "lca_name", "n_hits")) {
    if (!req_col %in% colnames(lca_df))
      stop("assign_trust_tiers(): lca_df missing '", req_col, "'.")
  }

  # Pre-compute per-ASV hit-level summaries
  if (!is.null(blast_df) && nrow(blast_df) > 0 &&
      all(c("qseqid", "pident", "sscinames") %in% colnames(blast_df))) {

    qcov_col <- intersect(c("query_cover", "qcovs", "qcov"), colnames(blast_df))
    if (length(qcov_col) == 0) {
      blast_df$query_cover <- NA_real_
      qcov_col <- "query_cover"
    } else {
      qcov_col <- qcov_col[1]
    }

    hit_summary <- blast_df %>%
      dplyr::mutate(
        pident      = as.numeric(pident),
        qcov        = as.numeric(.data[[qcov_col]]),
        sscinames   = as.character(sscinames),
        is_ecoli    = grepl(target_genus, sscinames, ignore.case = TRUE) & !is.na(sscinames),
        is_perfect  = !is.na(pident) & pident >= 100 & !is.na(qcov) & qcov >= 100
      ) %>%
      dplyr::group_by(qseqid) %>%
      dplyr::summarise(
        ncbi_total_perfect_hits     = sum(is_perfect,             na.rm = TRUE),
        ncbi_ecoli_hits             = sum(is_ecoli & is_perfect,  na.rm = TRUE),
        ncbi_seen_in_ecoli          = any(is_ecoli & is_perfect,  na.rm = TRUE),
        ncbi_top_hits_contain_ecoli = any(is_ecoli,               na.rm = TRUE),
        all_ecoli_hits              = all(is_ecoli,               na.rm = TRUE),
        any_perfect_ecoli           = any(is_ecoli & is_perfect,  na.rm = TRUE),
        all_perfect_ecoli           = all(is_ecoli & is_perfect,  na.rm = TRUE),
        any_perfect_nonecoli        = any(!is_ecoli & is_perfect, na.rm = TRUE),
        n_hits_total                = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        ncbi_ecoli_pct = round(ncbi_ecoli_hits /
                               pmax(ncbi_total_perfect_hits, 1) * 100, 1)
      )

  } else {
    message("assign_trust_tiers(): blast_df not available — all ASVs → Tier 5.")
    hit_summary <- data.frame(
      qseqid = character(), ncbi_total_perfect_hits = integer(),
      ncbi_ecoli_hits = integer(), ncbi_seen_in_ecoli = logical(),
      ncbi_top_hits_contain_ecoli = logical(), all_ecoli_hits = logical(),
      any_perfect_ecoli = logical(), all_perfect_ecoli = logical(),
      any_perfect_nonecoli = logical(), n_hits_total = integer(),
      ncbi_ecoli_pct = numeric(), stringsAsFactors = FALSE
    )
  }

  # Join and classify
  result <- asv_df %>%
    dplyr::left_join(lca_df,      by = c("asv" = "qseqid")) %>%
    dplyr::left_join(hit_summary, by = c("asv" = "qseqid")) %>%
    dplyr::mutate(
      lca_name                    = dplyr::coalesce(as.character(lca_name), NA_character_),
      top_pident                  = as.numeric(top_pident),
      top_qcov                    = as.numeric(top_qcov),
      wgs_confirmed               = as.logical(wgs_confirmed),
      ncbi_total_perfect_hits     = dplyr::coalesce(as.integer(ncbi_total_perfect_hits), 0L),
      ncbi_ecoli_hits             = dplyr::coalesce(as.integer(ncbi_ecoli_hits), 0L),
      ncbi_seen_in_ecoli          = dplyr::coalesce(ncbi_seen_in_ecoli,          FALSE),
      ncbi_top_hits_contain_ecoli = dplyr::coalesce(ncbi_top_hits_contain_ecoli, FALSE),
      all_ecoli_hits              = dplyr::coalesce(all_ecoli_hits,              FALSE),
      any_perfect_ecoli           = dplyr::coalesce(any_perfect_ecoli,           FALSE),
      all_perfect_ecoli           = dplyr::coalesce(all_perfect_ecoli,           FALSE),
      any_perfect_nonecoli        = dplyr::coalesce(any_perfect_nonecoli,        FALSE),
      n_hits                      = dplyr::coalesce(as.integer(n_hits),          0L),

      trust_tier = dplyr::case_when(
        wgs_confirmed | all_perfect_ecoli                                    ~ 1L,
        any_perfect_ecoli & !all_perfect_ecoli                               ~ 2L,
        !any_perfect_ecoli & all_ecoli_hits & ncbi_top_hits_contain_ecoli    ~ 3L,
        !any_perfect_ecoli & ncbi_top_hits_contain_ecoli & !all_ecoli_hits   ~ 4L,
        TRUE                                                                 ~ 5L
      ),

      trust_label = dplyr::case_when(
        trust_tier == 1L ~ "Tier 1: All perfect E. coli (or WGS match)",
        trust_tier == 2L ~ "Tier 2: At least one perfect E. coli hit",
        trust_tier == 3L ~ "Tier 3: No perfect hits, all hits E. coli",
        trust_tier == 4L ~ "Tier 4: No perfect hits, mixed taxonomy",
        TRUE             ~ "Tier 5: No E. coli hits"
      ),

      trust_group = dplyr::case_when(
        trust_tier %in% 1:2 ~ "E. coli",
        trust_tier %in% 3:4 ~ "Ambiguous",
        TRUE                ~ "Non E. coli"
      )
    )

  result
}


# ==============================================================================
# SECTION 11 — FIGURE 2: INTEGRATED DENDROGRAM
# ==============================================================================

#' Integrated dendrogram combining clustering, trust tiers, biological groups,
#' abundance, prevalence, origin strip, and sequence alignment panel.
#'
#' @param dm_bp         dist object of pairwise edit distances.
#' @param asv_meta      data.frame from build_asv_meta().
#' @param seqs_aligned  Aligned DNAStringSet.
#' @param sample_abund  OTU abundance matrix (ASVs × samples).
#' @param merge_dist    Merge distance threshold (for title only, default 0).
#' @param hclust_method Linkage method for hclust (default "average").
#' @param out_width     Output width hint in inches (default 36).
#' @param out_height    Output height hint in inches (default 20).
#' @return patchwork ggplot object.
figure2_dendrogram <- function(dm_bp,
                                asv_meta,
                                seqs_aligned,
                                sample_abund,
                                merge_dist    = 0,
                                hclust_method = "average",
                                out_width     = 36,
                                out_height    = 20) {

  # Clustering
  hc        <- hclust(dm_bp, method = hclust_method)
  tree      <- as.phylo(hc)
  tip_order <- tree$tip.label

  meta <- asv_meta %>%
    filter(asv %in% tip_order) %>%
    mutate(asv = factor(asv, levels = tip_order))

  # Colour palettes
  trust_tier_colors <- c(
    "Tier 1: All perfect E. coli (or WGS match)"  = "#1B7837",
    "Tier 2: At least one perfect E. coli hit"    = "#74C476",
    "Tier 3: No perfect hits, all hits E. coli"   = "#E65100",
    "Tier 4: No perfect hits, mixed taxonomy"     = "#FFB74D",
    "Tier 5: No E. coli hits"                     = "#C62828"
  )

  trust_group_colors <- c(
    "E. coli"     = "#1B7837",
    "Ambiguous"   = "#E65100",
    "Non E. coli" = "#C62828"
  )

  origin_colors <- c(
    "Best_only"     = "#4B0082",
    "Shared"        = "#B0C4DE",
    "Standard_only" = "#FFD700"
  )

  # Tree
  p_tree <- ggtree(tree, layout = "rectangular", ladderize = TRUE) +
    geom_tiplab(size = 1.5, align = TRUE, linesize = 0.2) +
    theme_tree2() +
    ggtitle(sprintf("Integrated Dendrogram: Merge Distance %d", merge_dist))

  strip_theme <- function() {
    theme_void() +
      theme(legend.position = "right",
            legend.key.size = unit(0.35, "cm"),
            legend.text     = element_text(size = 7),
            legend.title    = element_text(size = 8, face = "bold"),
            axis.text.y     = element_blank())
  }

  # Origin strip
  p_origin <- ggplot(meta %>% select(asv, origin),
                     aes(x = 0.5, y = asv, fill = origin)) +
    geom_tile(width = 0.85) +
    scale_fill_manual(values = origin_colors, name = "Origin") +
    scale_x_continuous(expand = c(0, 0)) +
    strip_theme()

  # Trust tier strips
  trust_strip_df <- meta %>%
    select(asv, trust_label, trust_group) %>%
    mutate(
      trust_label = factor(trust_label, levels = names(trust_tier_colors)),
      trust_group = factor(trust_group, levels = c("E. coli", "Ambiguous", "Non E. coli"))
    )
  p_trust_group <- ggplot(trust_strip_df, aes(x = 0.5, y = asv, fill = trust_group)) +
    geom_tile(width = 0.85) +
    scale_fill_manual(values = trust_group_colors, name = "Bio Group", drop = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    strip_theme()

  # Abundance heatmap (sqrt-scaled)
  abund_df <- meta %>%
    select(asv, abund_pct) %>%
    mutate(abund_sqrt = sqrt(pmax(abund_pct, 0)))

  abund_max_sqrt <- max(abund_df$abund_sqrt, na.rm = TRUE)

  p_abund <- ggplot(abund_df, aes(x = 0.5, y = asv, fill = abund_sqrt)) +
    geom_tile(width = 0.85) +
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#9ECAE1", "#2171B5", "#08306B"),
      limits  = c(0, abund_max_sqrt),
      name    = "Abund %\n(\u221a-scaled)",
      labels  = function(b) round(b^2, 2)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    strip_theme()

  # Prevalence heatmap
  p_prev <- ggplot(meta %>% select(asv, prev_pct),
                   aes(x = 0.5, y = asv, fill = prev_pct)) +
    geom_tile(width = 0.85) +
    scale_fill_gradientn(
      colours = c("#FFFFFF", "#BDD7E7", "#6BAED6", "#2171B5", "#08306B"),
      limits  = c(0, 100), name = "Prev %"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    strip_theme()

  # Sequence alignment panel
  nt_colors <- c(a = "#33A02C", c = "#1F78B4",
                 g = "#FF7F00", t = "#E31A1C", `-` = "#D9D9D9")

  seq_mat <- as.matrix(seqs_aligned)
  seq_df  <- as.data.frame(seq_mat) %>%
    rownames_to_column("asv") %>%
    filter(asv %in% tip_order) %>%
    pivot_longer(-asv, names_to = "pos", values_to = "nt") %>%
    mutate(
      pos = as.integer(sub("V", "", pos)),
      nt  = tolower(nt),
      asv = factor(asv, levels = tip_order)
    )

  p_align <- ggplot(seq_df, aes(x = pos, y = asv, fill = nt)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_manual(values = nt_colors, name = "Nucleotide", na.value = "grey95") +
    scale_x_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "right",
          legend.key.size = unit(0.35, "cm"),
          legend.text     = element_text(size = 7),
          legend.title    = element_text(size = 8, face = "bold"),
          axis.text.y     = element_blank())

  # Assemble with patchwork
  p_final <- p_tree + p_origin + p_trust_group + p_abund + p_prev + p_align +
    plot_layout(
      nrow   = 1,
      widths = c(4, 0.6, 0.6, 0.6, 0.6, 12),
      guides = "collect"
    ) &
    theme(legend.position = "right")

  p_final
}


# ==============================================================================
# SECTION 12 — BUILD ASV METADATA FOR FIGURE 2
# ==============================================================================

#' Build a per-ASV metadata data.frame for figure2_dendrogram().
#'
#' @param ps_best        phyloseq (Best pipeline).
#' @param ps_std         phyloseq (Standard pipeline).
#' @param trust_df       data.frame with asv (prefixed), trust_label, trust_tier, trust_group.
#' @param best_asvs_can  Character — "Best_*" prefixed names.
#' @param std_asvs_can   Character — "Standard_*" prefixed names.
#' @param seqs_best_char Character vector of Best sequences (content).
#' @param seqs_std_char  Character vector of Standard sequences (content).
#' @return data.frame ready for figure2_dendrogram().
build_asv_meta <- function(ps_best, ps_std, trust_df,
                            best_asvs_can  = NULL,
                            std_asvs_can   = NULL,
                            seqs_best_char = NULL,
                            seqs_std_char  = NULL) {

  best_asvs <- if (!is.null(best_asvs_can)) best_asvs_can else
    paste0("Best_", taxa_names(ps_best))
  std_asvs  <- if (!is.null(std_asvs_can))  std_asvs_can  else
    paste0("Standard_", taxa_names(ps_std))
  all_asvs  <- c(best_asvs, std_asvs)

  # Detect shared sequences (same nucleotide content in both pipelines)
  if (!is.null(seqs_best_char) && !is.null(seqs_std_char)) {
    best_seq_map   <- setNames(seqs_best_char, best_asvs)
    std_seq_map    <- setNames(seqs_std_char,  std_asvs)
    is_shared_best <- best_seq_map[best_asvs] %in% seqs_std_char
    is_shared_std  <- std_seq_map[std_asvs]   %in% seqs_best_char
    names(is_shared_best) <- best_asvs
    names(is_shared_std)  <- std_asvs
  } else {
    is_shared_best <- setNames(rep(FALSE, length(best_asvs)), best_asvs)
    is_shared_std  <- setNames(rep(FALSE, length(std_asvs)),  std_asvs)
  }

  origin_best <- ifelse(is_shared_best, "Shared", "Best_only")
  origin_std  <- ifelse(is_shared_std,  "Shared", "Standard_only")

  # Abundance / prevalence from merged phyloseq
  ps_merged <- merge_phyloseq(ps_best, ps_std)
  ps_rel    <- transform_sample_counts(ps_merged, function(x) x / sum(x))
  otu_mat   <- as(otu_table(ps_rel), "matrix")
  if (!taxa_are_rows(ps_rel)) otu_mat <- t(otu_mat)

  abund_pct <- rowMeans(otu_mat, na.rm = TRUE) * 100
  prev_pct  <- rowSums(otu_mat > 0, na.rm = TRUE) / ncol(otu_mat) * 100

  strip_prefix <- function(x) sub("^(Best_|Standard_)", "", x)

  abund_best <- abund_pct[strip_prefix(best_asvs)]
  abund_std  <- abund_pct[strip_prefix(std_asvs)]
  prev_best  <- prev_pct[strip_prefix(best_asvs)]
  prev_std   <- prev_pct[strip_prefix(std_asvs)]

  meta <- data.frame(
    asv       = all_asvs,
    abund_pct = c(abund_best, abund_std),
    prev_pct  = c(prev_best,  prev_std),
    origin    = c(origin_best, origin_std),
    stringsAsFactors = FALSE
  )

  meta <- meta %>%
    left_join(trust_df %>%
                select(asv, trust_label, trust_tier,
                       dplyr::any_of("trust_group")),
              by = "asv") %>%
    mutate(
      trust_tier  = tidyr::replace_na(as.integer(trust_tier), 5L),
      trust_label = dplyr::coalesce(trust_label, paste0("Tier ", trust_tier, ": No E. coli hits")),
      trust_group = dplyr::coalesce(
        if ("trust_group" %in% names(.)) trust_group else NA_character_,
        dplyr::case_when(
          trust_tier %in% 1:2 ~ "E. coli",
          trust_tier %in% 3:4 ~ "Ambiguous",
          TRUE                ~ "Non E. coli"
        )
      )
    )

  meta
}


# ==============================================================================
# SECTION 13 — FIGURE 3: BIO GROUP DISTRIBUTION BY ASV ORIGIN
# ==============================================================================

#' Proportional + count stacked barplot comparing Bio Group composition across
#' Shared / Unique-to-Best / Unique-to-Standard ASV categories.
#'
#' @param asv_meta  data.frame from build_asv_meta() — must contain columns:
#'                  asv, origin, trust_group.
#' @return patchwork of two ggplots (proportions top, counts bottom).
figure3_biogroup_by_origin <- function(asv_meta) {

  trust_group_colors <- c(
    "E. coli"     = "#1B7837",
    "Ambiguous"   = "#E65100",
    "Non E. coli" = "#C62828"
  )

  df <- asv_meta %>%
    filter(!is.na(origin), !is.na(trust_group)) %>%
    mutate(
      category = dplyr::recode(origin,
                               "Shared"        = "Shared",
                               "Best_only"     = "Unique to Best",
                               "Standard_only" = "Unique to Standard"
      ),
      category    = factor(category,
                           levels = c("Unique to Standard", "Shared", "Unique to Best")),
      trust_group = factor(trust_group,
                           levels = c("E. coli", "Ambiguous", "Non E. coli"))
    )

  summary_df <- df %>%
    group_by(category, trust_group) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(category) %>%
    mutate(
      total      = sum(n),
      proportion = n / total * 100
    ) %>%
    ungroup()

  # Panel A: proportional stacked barplot
  p_prop <- ggplot(summary_df,
                   aes(x = category, y = proportion, fill = trust_group)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(proportion >= 5,
                                 paste0(round(proportion, 1), "%"), "")),
              position = position_stack(vjust = 0.5),
              size = 3.2, colour = "white", fontface = "bold") +
    scale_fill_manual(values = trust_group_colors, name = "Bio Group") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title    = "Bio Group composition by ASV origin",
         subtitle = "Proportion of E. coli / Ambiguous / Non-E. coli per category",
         x = NULL, y = "Proportion (%)") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.title         = element_text(face = "bold"),
          legend.position    = "right")

  # Panel B: absolute count stacked barplot
  p_count <- ggplot(summary_df,
                    aes(x = category, y = n, fill = trust_group)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(n >= 1, n, "")),
              position = position_stack(vjust = 0.5),
              size = 3.2, colour = "white", fontface = "bold") +
    scale_fill_manual(values = trust_group_colors, name = "Bio Group") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(subtitle = "Absolute counts",
         x = NULL, y = "Number of ASVs") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          legend.position    = "right")

  # Assemble
  p_prop / p_count +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
}
