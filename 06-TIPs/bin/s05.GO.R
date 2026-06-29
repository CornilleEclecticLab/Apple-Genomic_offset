#!/usr/bin/env Rscript
rm(list = ls())
source("./bin/fonction.R")
source("./bin/config.R")
library(ggplot2)
library(ggrepel)

# --- 1. Setup Directories ---
message("Folders creation...")
create_directory(go_folder)
work.directory <- paste0(go_folder, "file/")
fig_output     <- paste0(go_folder, "figure/")
env_output     <- paste0(go_folder, "env/")

lapply(c(work.directory, fig_output, env_output), create_directory)
message("...Done")

# --- 2. Load Genomic & Environmental Data ---
file_name <- basename(sub("\\.lfmm$", "", lfmm))
info_ref_all_SYL <- read.table(ref, header = TRUE, sep = '\t')
sorted_info      <- info_ref_all_SYL[order(info_ref_all_SYL$Order_vcf), ]

# Load Genotype (Y)
geno_path <- ifelse(
  file.exists(paste0(computational_folder, "file/", file_name, "_filt.lfmm_imputed.geno")),
  paste0(computational_folder, "file/", file_name, "_filt.lfmm_imputed.geno"),
  paste0(computational_folder, "file/", file_name, "_filt.geno")
)
Y <- read.geno(geno_path)

# Load Environmental (X and X.pred)
env_name_1 <- basename(sub("\\.csv$", "", env.1))
X          <- read.env(paste0(lfmm_folder, "env/", env_name_1, ".env"))

env_name_2 <- basename(sub("\\.csv$", "", env.2))
out_file   <- file.path(env_output, paste0(env_name_2, ".env"))
next_env   <- read.delim(env.2, header = FALSE, sep = '\t', dec = '.')
env.format.2 <- scale(next_env, center = FALSE, scale = FALSE)
write.env(env.format.2, out_file)
X.pred     <- read.env(out_file)

# --- 3. Define SNP Sets to Process ---
# Load candidate lists (adjust paths if necessary)
gea_bonf <- read.csv(paste0(lfmm_folder, "file/", file_name, "_TIPS_GEA_Bonf.csv"))$x
gea_fdr  <- read.csv(paste0(lfmm_folder, "file/", file_name, "_TIPS_GEA_FDR.csv"))$x

# Create a named list: NULL means "All SNPs"
tips_lists <- list(
  "all_SNPs" = NULL,
  "GEA_Bonferroni" = gea_bonf,
  "GEA_FDR" = gea_fdr
)

# --- 4. Calculation and Plotting Loop ---
# Pre-calculate environmental distance once
Delta <- X[, 1:19] - X.pred[, 1:19]
squared.env.dist <- rowSums(Delta^2)

for (set_name in names(tips_lists)) {
  message(paste("Processing set:", set_name))
  
  # Run genetic gap
  # candidate.loci = NULL defaults to all SNPs in the 'genetic.gap' function
  g_gap <- genetic.gap(
    input = Y, 
    env = X, 
    pred.env = X.pred, 
    scale = TRUE, 
    K = K, 
    candidate.loci = tips_lists[[set_name]]
  )
  
  # Combine results
  res_df <- data.frame(
    squared.env.dist = squared.env.dist,
    offset = g_gap$offset,
    offset_Relstad = g_gap$distance,
    Country = sorted_info$Country.code..ISO.3166.,
    Site = sorted_info$Site
  )
  print(cor(res_df$offset, res_df$squared.env.dist))
  # Save TSV
  write.csv(res_df[,1:3], 
              paste0(work.directory, file_name, "_", set_name, "_offset_prediction.csv"), 
              row.names = FALSE)
  
  # Generate Plots for both metrics
  for (metric in c("offset", "offset_Relstad")) {
    p <- ggplot(res_df, aes(
      x = log(squared.env.dist),
      y = .data[[metric]],
      color = Country,
      label = Site
    )) +
      geom_point(size = 5) +
      geom_text_repel(size = 3, max.overlaps = 40) +
      geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
      labs(title = paste(set_name, "-", metric),
           x = "Log(Environmental distance)", 
           y = metric) +
      theme_bw()
    
    # Save files
    base_path <- paste0(fig_output, file_name, "_", set_name, "_", metric)
    ggsave(paste0(base_path, ".png"), plot = p, width = 10, height = 10, dpi = 300)
    ggsave(paste0(base_path, ".svg"), plot = p, width = 10, height = 10)
  }
}

message("Processing complete.")