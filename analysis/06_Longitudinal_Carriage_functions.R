# ==============================================================================
# 06_Longitudinal_Carriage_functions.R
#
# Project : fliC Amplicon Sequencing for High-Resolution Profiling of E. coli
#           in the Infant Gut Microbiome
# Author  : Andrés Catalán Tatay
# Part    : Analysis — Step 06 (companion functions)
# Purpose : All helper functions for 06_Longitudinal_Carriage.Rmd.
#           Covers data loading, metadata attachment, sample filtering, BLAST
#           parsing, and longitudinal carriage visualisation per subject and
#           pipeline. Adapted from supervisor's plot_carriage_blocks() function.
#
# Inputs  : data/cases/              — DADA2 output RDS files
#            data/wgs/               — longitudinal metadata phyloseq
#            results/04_best_vs_standard/blast/blast_df_parsed.csv
# Outputs : Called by 06_Longitudinal_Carriage.Rmd
#
# Usage   : source("06_Longitudinal_Carriage_functions.R")
# ==============================================================================

library(tidyverse)
library(phyloseq)


# ==============================================================================
# SECTION 1 — DATA LOADING
# ==============================================================================

#' Load a case RDS file and extract the phyloseq object.
#'
#' Handles both direct phyloseq objects and lists containing $physeq.
#'
#' @param path  Full path to the .rds file.
#' @return phyloseq object.
load_case_phyloseq <- function(path) {
  obj <- readRDS(path)
  if (inherits(obj, "phyloseq")) {
    message("  -> Loaded as direct phyloseq object.")
    return(obj)
  } else if (is.list(obj) && !is.null(obj$physeq) && inherits(obj$physeq, "phyloseq")) {
    message("  -> Loaded as list; extracting $physeq.")
    return(obj$physeq)
  } else {
    stop("RDS file does not contain a phyloseq object (neither directly nor under $physeq).")
  }
}


#' Load a phyloseq object saved with save() rather than saveRDS().
#'
#' Searches the loaded environment for any phyloseq object and returns the first.
#'
#' @param path  Full path to the .RDS file saved with save().
#' @return phyloseq object.
load_florentin_phyloseq <- function(path) {
  env          <- new.env()
  loaded_names <- load(path, envir = env)
  ps_candidates <- loaded_names[sapply(loaded_names, function(n) inherits(env[[n]], "phyloseq"))]
  if (length(ps_candidates) == 0) {
    stop("No phyloseq object found in the file loaded from: ", path)
  }
  if (length(ps_candidates) > 1) {
    warning("Multiple phyloseq objects found; using the first: ", ps_candidates[1])
  }
  message("  -> Loaded longitudinal metadata phyloseq object: '", ps_candidates[1], "'")
  return(env[[ps_candidates[1]]])
}


# ==============================================================================
# SECTION 2 — METADATA ATTACHMENT
# ==============================================================================

#' Attach longitudinal metadata from a reference phyloseq to a case phyloseq.
#'
#' Matches sample names between the case and reference, prunes unmatched samples,
#' and replaces the case sample_data with the specified metadata columns.
#'
#' @param ps_case       Case phyloseq object.
#' @param ps_florentin  Reference phyloseq carrying longitudinal metadata.
#' @param meta_cols     Character vector of metadata column names to transfer.
#' @return Case phyloseq with updated sample_data.
attach_longitudinal_metadata <- function(ps_case, ps_florentin,
                                         meta_cols = c("Subject", "days_since_baseline",
                                                        "weeks_since_baseline", "Type",
                                                        "Type2", "n_sample", "Time")) {
  flor_meta   <- as(sample_data(ps_florentin), "data.frame")
  case_snames <- sample_names(ps_case)
  flor_snames <- sample_names(ps_florentin)

  overlap <- intersect(case_snames, flor_snames)
  n_case  <- length(case_snames)
  n_over  <- length(overlap)

  message("  -> Case samples: ", n_case,
          " | Reference samples: ", length(flor_snames),
          " | Matched: ", n_over)

  if (n_over == 0) {
    stop("No sample names match between case phyloseq and reference phyloseq.\n",
         "  Case examples:      ", paste(head(case_snames, 5), collapse = ", "), "\n",
         "  Reference examples: ", paste(head(flor_snames, 5), collapse = ", "))
  }
  if (n_over < n_case) {
    warning(n_case - n_over, " case sample(s) not found in reference phyloseq and will be dropped.")
  }

  ps_case        <- prune_samples(overlap, ps_case)
  available_cols <- intersect(meta_cols, colnames(flor_meta))
  missing_cols   <- setdiff(meta_cols, colnames(flor_meta))
  if (length(missing_cols) > 0) {
    warning("These metadata columns are absent in the reference phyloseq: ",
            paste(missing_cols, collapse = ", "))
  }
  new_meta <- flor_meta[overlap, available_cols, drop = FALSE]
  sample_data(ps_case) <- sample_data(new_meta)
  return(ps_case)
}


# ==============================================================================
# SECTION 3 — SAMPLE FILTERING
# ==============================================================================

#' Filter a phyloseq object to retain only individual (patient) samples.
#'
#' Removes samples whose Type or Type2 metadata column matches any of the
#' following patterns: mock, ntc, neg, blank, control, singlestrain, dilut, standard.
#'
#' @param ps  phyloseq object with Type and/or Type2 in sample_data.
#' @return Filtered phyloseq object.
filter_individual_samples <- function(ps) {
  meta <- as(sample_data(ps), "data.frame")

  if ("Type"  %in% colnames(meta)) message("  -> Unique Type values:  ", paste(unique(meta$Type),  collapse = " | "))
  if ("Type2" %in% colnames(meta)) message("  -> Unique Type2 values: ", paste(unique(meta$Type2), collapse = " | "))

  exclude_patterns <- c("mock", "ntc", "neg", "blank", "control",
                        "singlestrain", "dilut", "standard")

  keep <- rep(TRUE, nsamples(ps))

  if ("Type" %in% colnames(meta)) {
    type_lower <- tolower(as.character(meta$Type))
    bad_type   <- Reduce(`|`, lapply(exclude_patterns, function(p) grepl(p, type_lower)))
    keep       <- keep & !bad_type
  }
  if ("Type2" %in% colnames(meta)) {
    type2_lower <- tolower(as.character(meta$Type2))
    bad_type2   <- Reduce(`|`, lapply(exclude_patterns, function(p) grepl(p, type2_lower)))
    keep        <- keep & !bad_type2
  }

  ps_filtered <- prune_samples(keep, ps)
  message("  -> Samples after individual filter: ", nsamples(ps_filtered),
          " (removed ", sum(!keep), ")")
  return(ps_filtered)
}


# ==============================================================================
# SECTION 4 — BLAST PARSING
# ==============================================================================

#' Parse the BLAST results CSV and identify E. coli ASVs per pipeline.
#'
#' Retains only ASVs with perfect BLAST hits:
#'   - pident == 100 (100% sequence identity)
#'   - query_cover == 100 (100% query coverage)
#'   - sscinames contains "Escherichia coli"
#'
#' Pipeline is inferred from the qseqid prefix ("Best_" or "Standard_").
#'
#' @param blast_file  Path to blast_df_parsed.csv (output of step 04).
#' @return data.frame with columns: Pipeline, ASV_ID.
parse_ecoli_asvs <- function(blast_file) {
  blast_df <- read.csv(blast_file, stringsAsFactors = FALSE)

  blast_df <- blast_df %>%
    mutate(
      Pipeline = case_when(
        stringr::str_detect(qseqid, "^Best_")     ~ "BEST",
        stringr::str_detect(qseqid, "^Standard_") ~ "STANDARD",
        TRUE                                       ~ NA_character_
      ),
      ASV_ID = stringr::str_remove(qseqid, "^(Best_|Standard_)")
    )

  unrecognised <- blast_df %>% filter(is.na(Pipeline))
  if (nrow(unrecognised) > 0) {
    warning(nrow(unrecognised),
            " BLAST row(s) have unrecognised qseqid prefix. They will be ignored.")
  }

  ecoli_asvs <- blast_df %>%
    filter(
      pident      == 100,
      query_cover == 100,
      stringr::str_detect(sscinames, "Escherichia coli")
    ) %>%
    distinct(Pipeline, ASV_ID)

  n_best     <- sum(ecoli_asvs$Pipeline == "BEST",     na.rm = TRUE)
  n_standard <- sum(ecoli_asvs$Pipeline == "STANDARD", na.rm = TRUE)
  message("  -> E. coli ASVs retained — BEST: ", n_best, " | STANDARD: ", n_standard)

  if ((n_best + n_standard) == 0) {
    stop("No E. coli ASVs passed the BLAST filter (pident=100, coverage=100, Escherichia coli).")
  }
  return(ecoli_asvs)
}


# ==============================================================================
# SECTION 5 — PHYLOSEQ FILTERING TO E. COLI ASVs
# ==============================================================================

#' Filter a phyloseq object to retain only E. coli ASVs for a given pipeline.
#'
#' @param ps             phyloseq object.
#' @param ecoli_asvs     data.frame from parse_ecoli_asvs() with Pipeline and ASV_ID.
#' @param pipeline_label Character: "BEST" or "STANDARD".
#' @return phyloseq subset containing only the E. coli ASVs detected in that pipeline.
filter_ecoli_ps <- function(ps, ecoli_asvs, pipeline_label) {
  target_asvs <- ecoli_asvs %>%
    filter(Pipeline == pipeline_label) %>%
    pull(ASV_ID)

  if (length(target_asvs) == 0) stop("No E. coli ASVs found for pipeline: ", pipeline_label)

  available <- intersect(target_asvs, taxa_names(ps))
  missing   <- setdiff(target_asvs, taxa_names(ps))

  if (length(missing) > 0) {
    warning("Pipeline ", pipeline_label, ": ", length(missing),
            " E. coli ASV(s) from BLAST not found in phyloseq taxa names:\n  ",
            paste(head(missing, 10), collapse = ", "))
  }
  if (length(available) == 0) {
    stop("None of the E. coli ASVs for pipeline '", pipeline_label,
         "' are present in the phyloseq object.")
  }

  ps_ecoli <- prune_taxa(available, ps)
  ps_ecoli <- prune_samples(sample_sums(ps_ecoli) > 0, ps_ecoli)
  message("  -> [", pipeline_label, "] E. coli ASVs in phyloseq: ", ntaxa(ps_ecoli),
          " | Samples with detections: ", nsamples(ps_ecoli))
  return(ps_ecoli)
}


# ==============================================================================
# SECTION 6 — SINGLE SUBJECT CARRIAGE PLOT
# ==============================================================================

#' Plot longitudinal E. coli carriage for one Subject × Pipeline combination.
#'
#' Visualisation: filled circles (shape 21) sized by relative abundance.
#' Only real detections (Abundance > 0) are plotted; absent timepoints show as
#' gaps. Vertical dashed lines mark all sampled timepoints.
#'
#' @param ps_ecoli           phyloseq filtered to E. coli ASVs (relative abundance).
#' @param subject            Character: Subject ID to plot.
#' @param pipeline_label     Character: pipeline label for plot title.
#' @param group_by           Metadata column to use as x-axis (default "days_since_baseline").
#' @param point_size_scale   Scalar multiplier for all point sizes (default 1.0).
#' @param abund_breaks       Numeric vector of abundance bin breakpoints.
#' @param abund_labels       Character vector of abundance bin labels.
#' @param abund_sizes        Named numeric vector of point sizes per abundance bin.
#' @param prev_high_cutoff   Global prevalence threshold for high-prevalence category.
#' @param prev_high_min_ab   Minimum abundance filter for high-prevalence ASVs.
#' @param prev_med_cutoff    Global prevalence threshold for medium-prevalence category.
#' @param prev_med_min_ab    Minimum abundance filter for medium-prevalence ASVs.
#' @param prev_low_min_ab    Minimum abundance filter for low-prevalence ASVs.
#' @param min_subject_samples Minimum number of detections required per ASV per subject.
#' @param peak_trust_abund   Minimum peak abundance to retain an ASV regardless of count.
#' @return Named list: plot (ggplot), df (detection data.frame), or NULL if no data.
plot_ecoli_carriage_single <- function(
    ps_ecoli,
    subject,
    pipeline_label,
    group_by            = "days_since_baseline",
    point_size_scale    = 1.0,
    abund_breaks        = c(-Inf, 0.5, 1, 2, 5, Inf),
    abund_labels        = c("< 0.5%", "0.5-1%", "1-2%", "2-5%", "> 5%"),
    abund_sizes         = c("< 0.5%" = 1.0, "0.5-1%" = 2.5,
                            "1-2%"   = 4.0, "2-5%"   = 6.0, "> 5%" = 8.0),
    prev_high_cutoff    = 0,
    prev_high_min_ab    = 0,
    prev_med_cutoff     = 0,
    prev_med_min_ab     = 0,
    prev_low_min_ab     = 0,
    min_subject_samples = 2,
    peak_trust_abund    = 0
) {

  # Convert to relative abundance
  ps_base <- prune_samples(sample_sums(ps_ecoli) > 0, ps_ecoli)
  ps_perc <- transform_sample_counts(ps_base, function(x) 100 * x / sum(x))

  # Global prevalence (dataset-wide) for adaptive filtering
  mat            <- as(otu_table(ps_base), "matrix")
  if (!taxa_are_rows(ps_base)) mat <- t(mat)
  asv_prevalence <- rowSums(mat > 0) / ncol(mat)

  prev_thresh_df <- tibble(
    OTU               = names(asv_prevalence),
    Global_Prevalence = asv_prevalence
  ) %>%
    mutate(
      Dynamic_Min_Abund = case_when(
        Global_Prevalence >= prev_high_cutoff ~ prev_high_min_ab,
        Global_Prevalence >= prev_med_cutoff  ~ prev_med_min_ab,
        TRUE                                  ~ prev_low_min_ab
      )
    )

  # Subset to subject
  sam_df    <- as(sample_data(ps_perc), "data.frame")
  keep_sams <- sam_df$Subject == subject
  keep_sams[is.na(keep_sams)] <- FALSE
  ps_sub    <- prune_samples(keep_sams, ps_perc)

  if (nsamples(ps_sub) == 0) {
    warning("No samples for Subject ", subject, " in pipeline ", pipeline_label, ". Skipping.")
    return(NULL)
  }
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  if (ntaxa(ps_sub) == 0) {
    warning("No taxa detected for Subject ", subject, " in pipeline ", pipeline_label, ". Skipping.")
    return(NULL)
  }

  meta             <- as(sample_data(ps_sub), "data.frame")
  meta[[group_by]] <- as.numeric(as.character(meta[[group_by]]))
  all_times        <- sort(unique(na.omit(meta[[group_by]])))

  # Melt and remove zeros — only real detections plotted
  df             <- psmelt(ps_sub)
  df[[group_by]] <- as.numeric(as.character(df[[group_by]]))
  df             <- df %>% filter(Abundance > 0)

  if (nrow(df) == 0) {
    warning("No detections (Abundance > 0) for Subject ", subject,
            " Pipeline ", pipeline_label, ". Skipping.")
    return(NULL)
  }

  # Adaptive prevalence filtering
  df_det <- df %>%
    left_join(prev_thresh_df, by = "OTU") %>%
    filter(Abundance >= Dynamic_Min_Abund) %>%
    group_by(OTU) %>%
    filter(n() >= min_subject_samples | max(Abundance) >= peak_trust_abund) %>%
    ungroup()

  if (nrow(df_det) == 0) {
    warning("No ASVs met filtering criteria for Subject ", subject,
            " Pipeline ", pipeline_label, ".")
    return(NULL)
  }

  # Labels and fill (one colour per ASV)
  df_det$LegendLabel <- df_det$OTU
  df_det$FillValue   <- df_det$OTU

  # Abundance binning for point size
  df_det <- df_det %>%
    mutate(
      Abund_Bin = cut(Abundance, breaks = abund_breaks, labels = abund_labels, right = FALSE)
    )

  final_sizes <- abund_sizes * point_size_scale

  # Factor ordering: alphabetical, reversed so ASV001 is at the top
  asv_levels         <- sort(unique(df_det$LegendLabel), decreasing = TRUE)
  df_det$LegendLabel <- factor(df_det$LegendLabel, levels = asv_levels)

  # Build plot: filled circles only
  p <- ggplot() +
    geom_vline(xintercept = all_times, linetype = "dashed", color = "grey60", alpha = 0.8) +
    geom_point(
      data  = df_det,
      aes(x    = .data[[group_by]],
          y    = LegendLabel,
          fill = FillValue,
          size = Abund_Bin),
      shape = 21,
      color = "black"
    ) +
    scale_size_manual(values = final_sizes, drop = FALSE) +
    guides(
      fill = "none",
      size = guide_legend(override.aes = list(fill = "grey50", shape = 21))
    ) +
    theme_linedraw() +
    labs(
      title    = paste0("Longitudinal E. coli Carriage — Subject: ", subject,
                        " | Pipeline: ", pipeline_label),
      subtitle = "Only ASVs with perfect BLAST support (100% identity, 100% coverage, Escherichia coli)",
      x        = "Days since baseline",
      y        = "E. coli ASV",
      size     = "Relative\nAbundance"
    ) +
    theme(
      panel.grid.major = element_line(color = "grey92"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 11, face = "bold"),
      plot.subtitle    = element_text(size = 8,  color = "grey40"),
      axis.text.y      = element_text(size = 8)
    )

  return(list(plot = p, df = df_det))
}


# ==============================================================================
# SECTION 7 — FULL ANALYSIS ACROSS ALL SUBJECTS AND PIPELINES
# ==============================================================================

#' Run the full longitudinal carriage analysis across all subjects and pipelines.
#'
#' For each Subject × Pipeline combination: generates and saves a carriage plot,
#' accumulates detection data into a master long-format dataframe, and computes
#' summary statistics. Mirrors the supervisor's carriage_df structure.
#'
#' @param ps_best            phyloseq — BEST pipeline (with longitudinal metadata).
#' @param ps_standard        phyloseq — STANDARD pipeline (with longitudinal metadata).
#' @param ecoli_asvs         data.frame from parse_ecoli_asvs().
#' @param subjects           Character vector of Subject IDs to process.
#' @param output_dir         Directory for saving figures and CSVs.
#' @param group_by           Metadata column for x-axis (default "days_since_baseline").
#' @param point_size_scale   Scalar multiplier for point sizes (default 1.0).
#' @param abund_breaks       Numeric vector of abundance bin breakpoints.
#' @param abund_labels       Character vector of abundance bin labels.
#' @param abund_sizes        Named numeric vector of point sizes per abundance bin.
#' @param prev_high_cutoff   Global prevalence threshold — high-prevalence tier.
#' @param prev_high_min_ab   Minimum abundance for high-prevalence ASVs.
#' @param prev_med_cutoff    Global prevalence threshold — medium-prevalence tier.
#' @param prev_med_min_ab    Minimum abundance for medium-prevalence ASVs.
#' @param prev_low_min_ab    Minimum abundance for low-prevalence ASVs.
#' @param min_subject_samples Minimum detections required per ASV per subject.
#' @param peak_trust_abund   Minimum peak abundance to retain an ASV unconditionally.
#' @param fig_width          Figure width in inches (default 10).
#' @param fig_height         Figure height in inches (default 6).
#' @return Named list:
#'   $individual_plots — named list of ggplot objects (Subject_Pipeline)
#'   $carriage_df      — master long-format detection data.frame
#'   $summary_df       — summary stats per Subject × Pipeline
run_longitudinal_carriage <- function(
    ps_best,
    ps_standard,
    ecoli_asvs,
    subjects,
    output_dir,
    group_by            = "days_since_baseline",
    point_size_scale    = 1.0,
    abund_breaks        = c(-Inf, 0.5, 1, 2, 5, Inf),
    abund_labels        = c("< 0.5%", "0.5-1%", "1-2%", "2-5%", "> 5%"),
    abund_sizes         = c("< 0.5%" = 1.0, "0.5-1%" = 2.5,
                            "1-2%"   = 4.0, "2-5%"   = 6.0, "> 5%" = 8.0),
    prev_high_cutoff    = 0,
    prev_high_min_ab    = 0,
    prev_med_cutoff     = 0,
    prev_med_min_ab     = 0,
    prev_low_min_ab     = 0,
    min_subject_samples = 2,
    peak_trust_abund    = 0,
    fig_width           = 10,
    fig_height          = 6
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pipelines <- list(BEST = ps_best, STANDARD = ps_standard)

  summary_rows  <- list()
  plot_list     <- list()
  carriage_list <- list()

  for (pipeline_label in names(pipelines)) {
    ps_full  <- pipelines[[pipeline_label]]
    ps_ecoli <- filter_ecoli_ps(ps_full, ecoli_asvs, pipeline_label)

    for (subj in subjects) {
      message("\nProcessing Subject: ", subj, " | Pipeline: ", pipeline_label)

      result <- plot_ecoli_carriage_single(
        ps_ecoli            = ps_ecoli,
        subject             = subj,
        pipeline_label      = pipeline_label,
        group_by            = group_by,
        point_size_scale    = point_size_scale,
        abund_breaks        = abund_breaks,
        abund_labels        = abund_labels,
        abund_sizes         = abund_sizes,
        prev_high_cutoff    = prev_high_cutoff,
        prev_high_min_ab    = prev_high_min_ab,
        prev_med_cutoff     = prev_med_cutoff,
        prev_med_min_ab     = prev_med_min_ab,
        prev_low_min_ab     = prev_low_min_ab,
        min_subject_samples = min_subject_samples,
        peak_trust_abund    = peak_trust_abund
      )

      plot_key <- paste0(subj, "_", pipeline_label)

      if (!is.null(result)) {
        p      <- result$plot
        df_det <- result$df %>% mutate(Subject = subj, Pipeline = pipeline_label)

        # Save PNG
        png_path <- file.path(output_dir,
                              paste0("longitudinal_carriage_", subj, "_", pipeline_label, ".png"))
        ggsave(png_path, plot = p, width = fig_width, height = fig_height, dpi = 300)
        message("  -> Saved PNG: ", png_path)

        # Save PDF
        pdf_path <- file.path(output_dir,
                              paste0("longitudinal_carriage_", subj, "_", pipeline_label, ".pdf"))
        ggsave(pdf_path, plot = p, width = fig_width, height = fig_height)
        message("  -> Saved PDF: ", pdf_path)

        plot_list[[plot_key]]     <- p
        carriage_list[[plot_key]] <- df_det

        # Summary stats
        n_asvs       <- length(unique(df_det$OTU))
        detected     <- paste(sort(unique(df_det$OTU)), collapse = ";")
        first_detect <- min(df_det[[group_by]], na.rm = TRUE)
        last_detect  <- max(df_det[[group_by]], na.rm = TRUE)
        n_timepoints <- length(unique(df_det[[group_by]]))

      } else {
        n_asvs <- 0; detected <- NA
        first_detect <- NA; last_detect <- NA; n_timepoints <- 0
      }

      summary_rows[[plot_key]] <- data.frame(
        Subject                = subj,
        Pipeline               = pipeline_label,
        n_ecoli_asvs           = n_asvs,
        asvs_detected          = detected,
        first_detection_day    = first_detect,
        last_detection_day     = last_detect,
        n_detection_timepoints = n_timepoints,
        stringsAsFactors       = FALSE
      )
    }
  }

  # Master carriage dataframe (mirrors supervisor's carriage_df structure)
  carriage_df <- bind_rows(carriage_list) %>%
    select(Subject, Pipeline, OTU, LegendLabel,
           all_of(group_by), Abundance, Abund_Bin, everything())

  carriage_csv <- file.path(output_dir, "06_carriage_df.csv")
  write.csv(carriage_df, carriage_csv, row.names = FALSE)
  message("\n-> Carriage dataframe saved: ", carriage_csv)

  # Summary CSV
  summary_df  <- bind_rows(summary_rows)
  summary_csv <- file.path(output_dir, "06_carriage_summary.csv")
  write.csv(summary_df, summary_csv, row.names = FALSE)
  message("-> Summary CSV saved: ", summary_csv)

  return(list(
    individual_plots = plot_list,
    carriage_df      = carriage_df,
    summary_df       = summary_df
  ))
}
