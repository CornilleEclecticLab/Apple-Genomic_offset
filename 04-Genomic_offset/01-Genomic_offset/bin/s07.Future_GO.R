#!/usr/bin/env Rscript

# ============================================================
# Script: s07.Future_GO.R
# Purpose: Run genomic offset (genetic.gap) for future climate
#          conditions
# ============================================================

# ---- Clean environment ----
rm(list = ls())
source("./bin/fonction.R")

# Load libraries silently to keep logs clean
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(broom.mixed)
  library(lme4)
  library(ggcorrplot)
  library(sf)
  library(raster)
  library(gstat)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(patchwork)
})

# ---- Input arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript bin/s07.Future_GO.R <project_name> <env_folder> <reference_info> <K_value>")
}

PROJECT    <- args[1]
env_folder <- args[2]
ref_info   <- args[3]
K          <- as.numeric(args[4])

# PROJECT <- "dom_admixed_rm"
# env_folder <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/dom_admixed_rm/"
# ref_info <- "./in/data_tree_dom_admixed_rm_V2.txt"
# K <- 5

message("🔧 Running Genomic Offset for project: ", PROJECT)

# ---- Define folders ----
go_main  <- file.path("./out", PROJECT, "s07.GO_future_V2")
file_dir <- file.path(go_main, "file")
fig_dir  <- file.path(go_main, "figure")
env_dir  <- file.path(go_main, "env")
lfmm2_main <- file.path("./out", PROJECT, "s02.lfmm2_V3")

invisible(lapply(c(go_main, file_dir, fig_dir, env_dir), create_directory))

# ---- Reference & Climate Prep ----
info_ref_all_SYL <- read.table(ref_info, header = TRUE, sep = "\t")
sorted_info <- info_ref_all_SYL[order(info_ref_all_SYL$Order_vcf), ]

# Helper to process .env files
make_env_file <- function(pattern, src_folder, out_folder, ref_df) {
  files <- list.files(src_folder, pattern = pattern, full.names = TRUE)
  for (f in files) {
    env_name <- basename(sub("\\.csv$", "", f))
    env_data <- read.csv(f, header = TRUE, dec = ".")
    # Match order to genotype matrix
    env_data <- env_data[match(ref_df$Seq_ID, env_data$Seq_ID), ]
    write.env(env_data[, -1], file.path(out_folder, paste0(env_name, ".env")))
  }
}

message("📂 Preparing environment files...")
make_env_file("climate_extract_current_*", env_folder, env_dir, sorted_info)
make_env_file("climate_extract_future_*", env_folder, env_dir, sorted_info)

X <- read.env(file.path(env_dir, paste0("climate_extract_current_", PROJECT, ".env")))
full_geno <- read.geno(file.path(lfmm2_main, "file", paste0("Full_chr_pvalue_", PROJECT, ".geno")))
full_geno_filt <- read.geno(file.path(lfmm2_main, "file", paste0("Full_chr_pvalue_filt_", PROJECT, ".geno")))

# ---- Core Function with Generalized Legend ----
Go_future_pred <- function(Y, X, K, future_files, sorted_info, candidate.loci = NULL, type = "Full", suffix = "") {
  
  all_results <- list()
  
  # 1. Data Processing Loop
  for (env_path in future_files) {
    fname <- basename(env_path)
    timeline <- str_extract(fname, "\\d{4}-\\d{4}")
    ssp      <- str_extract(fname, "ssp\\d+")
    
    message(paste("   > Calculating offset for:", timeline, ssp))
    
    X.pred <- read.env(env_path)
    
    # Genetic Gap Calculation
    g.gap <- genetic.gap(input = Y, env = X, pred.env = X.pred, 
                         scale = TRUE, K = K, candidate.loci = candidate.loci)

    Delta <- X - X.pred
    data_offset <- cbind(
      sorted_info[c("Seq_ID", "Country", "Site", "Latitude", "Longitude")],
      squared.env.dist = rowSums(Delta^2), 
      offset = g.gap$offset, 
      offset_Relstad = g.gap$distance
    )
    
    out_name <- sprintf("All_chr_offset_%s%s_%s_%s.csv", type, suffix, timeline, ssp)
    write.csv(data_offset, file.path(file_dir, out_name), row.names = FALSE)
    
    all_results[[paste(timeline, ssp, sep = "_")]] <- data_offset
  }

  # --- Visualization (Generalized Legend Logic) ---
  combined_df <- bind_rows(all_results, .id = "ID")
  
  # Calculate global limits to ensure color scales match exactly
  global_limits <- range(combined_df$offset, na.rm = TRUE)
  
  timelines <- sort(unique(str_extract(names(all_results), "\\d{4}-\\d{4}")))
  ssps      <- sort(unique(str_extract(names(all_results), "ssp\\d+")))
  
  plot_list <- list()
  for (t in timelines) {
    for (s in ssps) {
      pid <- paste(t, s, sep = "_")
      if (!is.null(all_results[[pid]])) {
        
        # Calculate stats for this specific map
        current_data <- all_results[[pid]]
        mean_val <- mean(current_data$offset, na.rm = TRUE)
        n_points <- nrow(current_data)
        
        # Create plot
        p <- Map_Go_extrapolated(
          df = current_data, 
          timeline = t, ssp = s, type = type,
          limits = global_limits
        )
        
        # FORCE identical scales and add stats as subtitle
        # Use '&' later to apply theme broadly, but scale must be set here
        plot_list[[pid]] <- p + 
          labs(
            title = paste(t, "|", s),
            subtitle = sprintf("N = %d | Mean GO = %.3f", n_points, mean_val)
          ) +
          # CRITICAL: The 'name' must be identical across all plots for patchwork to collect
          scale_color_viridis_c(limits = global_limits, name = "Genomic Offset") + 
          theme(
            axis.title = element_blank(), 
            axis.text = element_blank(),
            plot.title = element_text(size = 10, face = "bold"),
            plot.subtitle = element_text(size = 8)
          )
      }
    }
  }

  # Combine and Collect
  # The '&' operator applies the theme to all subplots simultaneously
  figure_3x3 <- wrap_plots(plot_list, ncol = 3) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "right")

  # Add global annotation
  figure_3x3 <- figure_3x3 + 
    plot_annotation(
      title = paste("Future Genetic Offset:", type),
      subtitle = paste("Correction:", suffix, "| Global Range:", 
                       round(global_limits[1],3), "-", round(global_limits[2],3)),
      theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
    )

  fig_name <- sprintf("Figure_3x3_future_offset_%s%s.png", type, suffix)
  ggsave(file.path(fig_dir, fig_name), figure_3x3, width = 15, height = 15, bg = "white")
}

# ---- Execution Loop ----
GEA_files <- list.files(path = file.path(lfmm2_main, "file"), pattern = "All_SNPs_.*\\.csv$", full.names = TRUE)
future_envs <- list.files(path = env_dir, pattern = ".*ssp.*\\.env$", full.names = TRUE)

for(file in GEA_files) {
  is_filtered <- grepl("Filtered", file)
  Y_geno      <- if(is_filtered) full_geno_filt else full_geno
  suffix      <- if(is_filtered) "_filtered" else ""
  
  gea_data <- read.csv(file)
  
  # 1. Bonferroni
  cand_bonf <- gea_data$LFMM_pos[gea_data$Bonf_signif == TRUE]
  if(length(cand_bonf) > 0) {
    Go_future_pred(Y_geno, X, K, future_envs, sorted_info, cand_bonf, "Bonf", suffix)
  }

  # 2. FDR
  cand_fdr <- gea_data$LFMM_pos[gea_data$FDR_signif == TRUE]
  if(length(cand_fdr) > 0) {
    Go_future_pred(Y_geno, X, K, future_envs, sorted_info, cand_fdr, "FDR", suffix)
  }

  # 3. Full Model (All Loci)
  Go_future_pred(Y_geno, X, K, future_envs, sorted_info, NULL, "Full", suffix)
  
  gc() # Garbage collection to free up RAM
}

message("✅ Script completed successfully.")