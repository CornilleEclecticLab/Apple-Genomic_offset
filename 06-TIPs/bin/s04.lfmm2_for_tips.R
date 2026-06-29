#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 1. SETUP & LIBRARIES
# ------------------------------------------------------------------------------
rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(qvalue)
  library(data.table)
  library(LEA) # For read.geno, read.env, lfmm2
})

# Load custom functions (ensure create_directory, plot_density, etc., are here)
source("./bin/fonction.R")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4){
  stop("Usage: Rscript bin/s02.postprocess_LFMM.R <project_name> <env_csv_file> <position_data> <K>")
}

PROJECT       <- args[1]
env_csv       <- args[2]
position_data <- args[3]
K             <- as.numeric(args[4])
# geno_LDpruned <- args[5] ## on va pas en avoir besoin ici

## for test
PROJECT       <- "tips_dom_admixed_rm_V3"
env_csv       <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/dom_admixed_rm/climate_extract_current_dom_admixed_rm.csv"
position_data <- "./in/data_tree_dom_admixed_rm_V2.txt"
K             <- 5

# ------------------------------------------------------------------------------
# 2. DIRECTORIES & OUTPUT SETUP
# ------------------------------------------------------------------------------
lfmm2_main <- file.path("./out", PROJECT, "s04.lfmm2_for_tips")
file_dir   <- file.path(lfmm2_main, "file")
fig_dir    <- file.path(lfmm2_main, "figure")
env_dir    <- file.path(lfmm2_main, "env")

invisible(lapply(c(file_dir, fig_dir, env_dir), create_directory))

# ------------------------------------------------------------------------------
# 3. ENVIRONMENTAL DATA PREPARATION
# ------------------------------------------------------------------------------
env_data  <- read.csv(env_csv, header = TRUE, sep = ",", dec = ".")
accession <- read.delim(position_data, header = TRUE, sep = "\t")

# Ensure environmental data matches the genomic individual order
env_ordered <- env_data[match(accession$Seq_ID, env_data$Seq_ID), ]
env_matrix  <- env_ordered[, -1, drop = FALSE] 
name_plots  <- colnames(env_matrix)

out_env_file <- file.path(env_dir, paste0(basename(tools::file_path_sans_ext(env_csv)), ".env"))
write.env(env_matrix, out_env_file)
X <- read.env(out_env_file)

# ------------------------------------------------------------------------------
# 4. RUN LATENT FACTOR MODEL (LFMM2)
# ------------------------------------------------------------------------------
comp_main <- file.path("./out", PROJECT, "s03.computation_for_tips")
tips_folder <- list.dirs(comp_main, recursive = FALSE) # only one

for(fold in tips_folder){
    if (file.exists(paste0(fold, "/file/GraffiTE.merged.genotypes_",PROJECT,"_filt.lfmm_imputed.geno"))){
        Y <- read.geno(paste0(fold, "/file/GraffiTE.merged.genotypes_",PROJECT,"_filt.lfmm_imputed.geno"))
    } else {
        Y <- read.geno(paste0(fold, "/file/GraffiTE.merged.genotypes_",PROJECT,"_filt.geno"))
    }
}

message(">>> Running LFMM2 (K=", K, ")...")
mod_lfmm2 <- lfmm2(input = Y, env = X, K = K)

# Plot Latent Factors
# png(file.path(fig_dir, "Latent_Factors_U.png"), width = 8, height = 6, units = "in", res = 300)
# plot(mod_lfmm2@U, col = "darkgrey", pch = 19, xlab = "Factor 1", ylab = "Factor 2", main = "LFMM2 Latent Factors")
# grid()
# dev.off()

df_lfmm <- data.frame(
  Factor_1 = mod_lfmm2@U[, 1],
  Factor_2 = mod_lfmm2@U[, 2]
)

# --- 2. Build the plot ---
p <- ggplot(df_lfmm, aes(x = Factor_1, y = Factor_2)) +
  # Adds points; shape 19 is a solid circle, color matches "darkgrey"
  geom_point(color = "darkgrey", shape = 19, size = 2) + 
  
  # Clean, minimal theme with grid lines that closely mimics base grid()
  theme_bw() + 
  
  # Titles and axis labels
  labs(
    title = "LFMM2 Latent Factors",
    x = "Factor 1",
    y = "Factor 2"
  ) +
  
  # Center align the title to match base R's style
  theme(plot.title = element_text(hjust = 0.5))

# --- 3. Save the plot ---
# ggsave automatically handles opening, writing, and closing the device
ggsave(
  filename = "Latent_Factors_U.png",
  plot = p,
  path = fig_dir,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------
# 5. GEA (ASSOCIATION TESTING)
# ------------------------------------------------------------------------------

# Using lists for speed (avoids O(n^2) rbind overhead)
full_res_list <- list()
filt_res_list <- list()
bio_sep_list  <- list()

for (fold in tips_folder) {
  tips_file_name <- basename(fold)
  message(">>> Processing Chromosome: ", tips_file_name)
  #chr_id <- sub(".*(chr[0-9]+).*", "\\1", chr_name)
  
  # Load SNP names and Imputed Genotypes
  tips_csv <- list.files(file.path(fold, "file"), "_TIPs_names_filtered.csv$", full.names = TRUE)[1]
  tips_names <- read.csv(tips_csv)[,1]
  
#   geno_f <- list.files(file.path(chr_folder, "file"), "_filt.lfmm_imputed.geno$", full.names = TRUE)[1]
#   Y_chr  <- read.geno(geno_f)
#   full_geno_list[[chr_name]] <- Y_chr

  # A. Full Model Test (All env variables combined)
  pv_full <- lfmm2.test(object = mod_lfmm2, input = Y, env = X, full = TRUE)
  df_file  <- data.frame(TIPs = tips_names, pvalues = as.numeric(pv_full$pvalues), file = tips_file_name)
  full_res_list <- df_file
  
  # B. P-value Filtering (Noise reduction for p >= 0.9)
  keep_idx <- which(df_file$pvalues < 0.9)
  Y_filt   <- Y[, keep_idx, drop = FALSE]
  #full_geno_list_filt[[tips_file_name]] <- Y_filt # Save for combined geno file
  
  pv_filt <- lfmm2.test(object = mod_lfmm2, input = Y_filt, env = X, full = TRUE)
  filt_res_list <- data.frame(TIPs = tips_names[keep_idx], pvalues = as.numeric(pv_filt$pvalues), file = tips_file_name)
  
  # C. Separated Bio Variable Test
  pv_sep  <- lfmm2.test(object = mod_lfmm2, input = Y_filt, env = X, full = FALSE)
  df_bio  <- data.frame(TIPs = tips_names[keep_idx], file = tips_file_name)
  
  for(i in 1:ncol(X)) {
    df_bio[[paste0("p_", name_plots[i])]] <- as.numeric(pv_sep$pvalues[i,])
  }
  bio_sep_list <- df_bio
}

# ------------------------------------------------------------------------------
# 6. GLOBAL MERGING & SIGNIFICANCE SCORING
# ------------------------------------------------------------------------------
message(">>> Merging results and calculating global significance...")
all_res    <- full_res_list
all_res_rm <- filt_res_list
all_bio    <- bio_sep_list

# Function to add significance columns
add_sig <- function(df, p_col, is_filtered = FALSE) {
  df$LFMM_pos <- rownames(df)
  p_vals <- df[[p_col]]
  # Use specific lambda for filtered data to account for truncation
   if(is_filtered){ lambda_seq <- seq(0.05, 0.89, 0.05) 
    q_obj <- qvalue::qvalue(p_vals, fdr.level = 0.05, lambda = lambda_seq)
   } else {
    q_obj <- qvalue::qvalue(p_vals, fdr.level = 0.05)
   }
  df$qvalue <- q_obj$qvalues
  df$Bonf_signif <- p_vals < (0.05 / nrow(df))
  df$FDR_signif  <- q_obj$significant
  return(df)
}

all_res    <- add_sig(all_res, "pvalues"); print(head(all_res))
all_res_rm <- add_sig(all_res_rm, "pvalues", is_filtered = TRUE); print(head(all_res_rm))

write.csv(all_res, file.path(file_dir, "All_TIPs_Full_Results.csv"), row.names = FALSE)
write.csv(all_res_rm, file.path(file_dir, "All_TIPs_Filtered_Results.csv"), row.names = FALSE)


# Combine and save the filtered genotype matrix
# final_geno <- do.call(cbind, full_geno_list)
# write.geno(final_geno, file.path(file_dir, paste0("Full_chr_pvalue_", PROJECT, ".geno")))

# final_geno_filt <- do.call(cbind, full_geno_list_filt)
# write.geno(final_geno_filt, file.path(file_dir, paste0("Full_chr_pvalue_filt_", PROJECT, ".geno")))

# ------------------------------------------------------------------------------
# 7. VISUALIZATION FUNCTIONS
# ------------------------------------------------------------------------------
# plot_geas_manhattan <- function(df, p_col, title, out_file) {
#   # Calculate cumulative position for Manhattan plot
#   df_plot <- df %>%
#     group_by(chromosome) %>%
#     summarise(chr_len = n()) %>%
#     mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
#     select(-chr_len) %>%
#     left_join(df, ., by = "chromosome") %>%
#     arrange(chromosome, SNPs) %>%
#     mutate(index = row_number() + (tot*0.05)) # get 5% in space between chr
  
#   axis_df <- df_plot %>% group_by(chromosome) %>% summarize(center = mean(index))
#   thresh <- -log10(0.05 / nrow(df))
  
#   p <- ggplot(df_plot, aes_string(x = "index", y = paste0("-log10(", p_col, ")"))) +
#     # Background points (alternating chromosome colors)
#     geom_point(aes(color = chromosome), alpha = 0.5, size = 0.7) +
#     # Highlight Significant SNPs
#     geom_point(data = subset(df_plot, FDR_signif & !Bonf_signif), color = "#A0B5E5", size = 1) +
#     geom_point(data = subset(df_plot, Bonf_signif), color = "#F5B0B0", size = 1.2) +
    
#     scale_color_manual(values = rep(c("grey60", "grey60"), 50)) +
#     scale_x_continuous(label = axis_df$chromosome, breaks = axis_df$center) +
#     geom_hline(yintercept = thresh, color = "black", linetype = "dashed") +
#     labs(title = title, subtitle = "Red: Bonferroni | Blue: FDR 5%", x = "Chromosome", y = "-log10(p)") +
#     theme_minimal() +
#     theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

#   ggsave(out_file, plot = p, width = 13, height = 10, dpi = 300)
# }

# ------------------------------------------------------------------------------
# 8. GENERATE ALL PLOTS
# ------------------------------------------------------------------------------
# message(">>> Generating plots...")

# # Global Manhattan
# plot_geas_manhattan(all_res, "pvalues", "Global GEA: Full Dataset", file.path(fig_dir, "Manhattan_Global_Full.png"))
# plot_geas_manhattan(all_res_rm, "pvalues", "Global GEA: P-Filtered Dataset", file.path(fig_dir, "Manhattan_Global_Filtered.png"))

# Per-Variable Manhattan
for(bio in name_plots) {
  p_col <- paste0("p_", bio)
  # Temp significance for the specific bio variable
  tmp_bio <- all_bio[, c("TIPs", p_col)]
  tmp_bio <- add_sig(tmp_bio, p_col)
  write.csv(tmp_bio, file.path(file_dir, paste0(bio, "_TIPs_Filtered_Results.csv")), row.names = FALSE)
  
#   plot_geas_manhattan(tmp_bio, p_col, paste("GEA Association:", bio), file.path(fig_dir, paste0("Manhattan_", bio, ".png")))
}

# Global Distributions
plot_histo(all_res, "Global P-value Distribution", file.path(fig_dir, "Pval_Histo_Global.png"))

message("✅ LFMM2 Post-processing complete for project: ", PROJECT)