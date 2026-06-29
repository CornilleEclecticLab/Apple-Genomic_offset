#!/usr/bin/env Rscript

# ============================================================
# Script: s08.Past_climate_GO.R
# Purpose: Run genomic offset (genetic.gap) for future climate
#          conditions
# ============================================================

# ---- Clean environment ----
rm(list = ls())
source("./bin/fonction.R")

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

# PROJECT <- "dom_admixed_rm"
# env_folder <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/dom_admixed_rm/"
# ref_info <- "./in/data_tree_dom_admixed_rm_V2.txt"
# K <- 5

# ---- Input arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript bin/s07.Future_GO.R <project_name> <env_folder> <reference_info> <K_value>")
}

PROJECT    <- args[1]
env_folder <- args[2]
ref_info   <- args[3]
K          <- as.numeric(args[4])

message("🔧 Running Past Genomic Offset for project: ", PROJECT)

# ---- Define folders ----
go_main  <- file.path("./out", PROJECT, "s08.Past_climate_GO")
file_dir <- file.path(go_main, "file")
fig_dir  <- file.path(go_main, "figure")
env_dir  <- file.path(go_main, "env")
lfmm2_main <- file.path("./out", PROJECT, "s02.lfmm2_V3")

invisible(lapply(c(go_main, file_dir, fig_dir, env_dir), create_directory))

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

Go_past_pred <- function(Y, X, K, past_env_path, sorted_info, candidate.loci = NULL, type = "Full", suffix = "") {
  
  # 1. Extract Metadata from filename
  # This extracts 'LGM' or 'MidHolocene' from the filename
  fname <- basename(past_env_path)
  period <- str_extract(fname, "(?<=extract_)[^._]+") 
  if(is.na(period)) period <- "Past"
  
  message(paste("   > Calculating offset for period:", period))
  
  # 2. Data Processing
  X.pred <- read.env(past_env_path)
  
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
  
  # Save the CSV
  out_name <- sprintf("All_chr_offset_%s%s_%s.csv", type, suffix, period)
  write.csv(data_offset, file.path(file_dir, out_name), row.names = FALSE)
  
  # 3. Visualization (Single Map)
  # Calculate stats for the subtitle
  mean_val <- mean(data_offset$offset, na.rm = TRUE)
  n_points <- nrow(data_offset)
  
  # Generate map
  # Note: limits are set based on the current data since it's a single map
  p <- Map_Go_extrapolated(
    df = data_offset, 
    timeline = period, ssp = "", type = type,
    limits = range(data_offset$offset, na.rm = TRUE)
  )
  
  # Final styling for a single plot
  final_plot <- p + 
    labs(
      title = paste("Past Genetic Offset:", period),
      subtitle = sprintf("N = %d | Mean GO = %.3f | Type: %s%s", n_points, mean_val, type, suffix)
    ) +
    #scale_color_viridis_c(name = "Genomic Offset") + 
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "right"
    )

  # Save Plot
  fig_name <- sprintf("Figure_Past_offset_%s%s_%s.png", type, suffix, period)
  ggsave(file.path(fig_dir, fig_name), final_plot, width = 8, height = 7, bg = "white")
  
  return(data_offset) # Return data in case you want to use it further
}

message("📂 Preparing environment files...")
make_env_file("climate_extract_current_*", env_folder, env_dir, sorted_info)
make_env_file("climate_extract_LGM_*", env_folder, env_dir, sorted_info)

X <- read.env(file.path(env_dir, paste0("climate_extract_current_", PROJECT, ".env")))
full_geno <- read.geno(file.path(lfmm2_main, "file", paste0("Full_chr_pvalue_", PROJECT, ".geno")))
full_geno_filt <- read.geno(file.path(lfmm2_main, "file", paste0("Full_chr_pvalue_filt_", PROJECT, ".geno")))

GEA_files <- list.files(path = file.path(lfmm2_main, "file"), pattern = "All_SNPs_.*\\.csv$", full.names = TRUE)
past_envs <- list.files(path = env_dir, pattern = ".*LGM.*\\.env$", full.names = TRUE)[1] # there is only one file

for(file in GEA_files) {
  is_filtered <- grepl("Filtered", file)
  Y_geno      <- if(is_filtered) full_geno_filt else full_geno
  suffix      <- if(is_filtered) "_filtered" else ""
  
  gea_data <- read.csv(file)
  
  # 1. Bonferroni
  cand_bonf <- gea_data$LFMM_pos[gea_data$Bonf_signif == TRUE]
  if(length(cand_bonf) > 0) {
    Go_past_pred(Y_geno, X, K, past_envs, sorted_info, cand_bonf, "Bonf", suffix)
  }

  # 2. FDR
  cand_fdr <- gea_data$LFMM_pos[gea_data$FDR_signif == TRUE]
  if(length(cand_fdr) > 0) {
    Go_past_pred(Y_geno, X, K, past_envs, sorted_info, cand_fdr, "FDR", suffix)
  }

  # 3. Full Model (All Loci)
  Go_past_pred(Y_geno, X, K, past_envs, sorted_info, NULL, "Full", suffix)
  
  gc() # Garbage collection to free up RAM
}