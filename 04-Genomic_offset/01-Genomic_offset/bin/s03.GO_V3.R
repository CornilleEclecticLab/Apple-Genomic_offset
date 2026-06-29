#!/usr/bin/env Rscript

# ---- Setup & Arguments ----
rm(list = ls())
source("./bin/fonction.R")
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: Rscript s03.R <proj> <env_curr> <new_env_dir> <ref_info> <K>")

PROJECT    <- args[1]
ENV_CURR   <- args[2]
NEW_ENV_D  <- args[3]
REF_INFO   <- args[4]
K          <- as.numeric(args[5])

message(
  "PROJECT: ",PROJECT,"\n",     
  "ENV_CURR: ",ENV_CURR,"\n",
  "NEW_ENV_D: ",NEW_ENV_D,"\n",   
  "REF_INFO: ",REF_INFO,"\n",    
  "K: ",K,"\n"          
)

# PROJECT    <- "dom_admixed_rm"
# ENV_CURR   <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/dom_admixed_rm/climate_extract_current_dom_admixed_rm.csv"
# NEW_ENV_D  <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s03.orchard_generalization/dom_admixed_rm" 
# REF_INFO   <- "./in/data_tree_dom_admixed_rm_V2.txt"
# K          <- 5

# ---- Directories ----
go_main  <- file.path("./out", PROJECT, "s03.GO_V3")
sub_dirs <- setNames(nm = c("file", "figure", "env"))
lapply(file.path(go_main, sub_dirs), create_directory)

# ---- Data Loading ----
info_ref <- read.table(REF_INFO, header = TRUE, sep = "\t")
info_ref <- info_ref[order(info_ref$Order_vcf), ]

# Process Current Env
env_current <- read.delim(ENV_CURR, header = TRUE, sep = ",", dec = ".")
env_current_name <- basename(sub("\\.csv$", "", ENV_CURR))
env_current <- env_current[match(info_ref$Seq_ID, env_current$Seq_ID), ]
write.env(env_current[-1], file.path(go_main, "env", paste0(env_current_name, ".env")))
X <- read.env(file.path(go_main, "env", paste0(env_current_name, ".env"))) # Remove ID col

# Process Future Env Files
new_env_files <- list.files(NEW_ENV_D, pattern = "env_.*", full.names = TRUE)
print(new_env_files)
# for (f in new_env_files) {
#   print(f)
#   tmp <- read.delim(f, header = F, sep = "\t", dec = ".")
#   write.env(tmp, file.path(go_main, "env", paste0(basename(sub("\\.csv$", "", f)), ".env")))
# }

# ---- Plotting function ----
plot_go_env <- function(data.gap, data.info, type=NULL, fig.dir, data.type=c("GEA_SNPs", "Full_SNPs"), OUT_FILE ){

  for (metric in c("offset", "offset_Relstad")) {
  
    plot1 <- ggplot(data.gap, aes(x = log(squared.env.dist), y = data.gap[, metric],
                          color = data.info$Country.code..ISO.3166., label = data.info$Site)) +
      geom_point(size = 5) +
      geom_text_repel(size = 3, max.overlaps = 40) +
      geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
      labs(x = "Log(Environmental distance)", y = metric, color = "Country") +
      theme_bw()

    ggsave(file.path(fig.dir, paste0(OUT_FILE, metric, "_", data.type, "_", type, ".png")),
    plot = plot1,
    width = 6, height = 6, units = "in", dpi = 300)

    plot2 <- ggplot(data.gap, aes(
        x = log(squared.env.dist),
        y = .data[[metric]],
        color = data.info$Country.code..ISO.3166.,
        label = data.info$Seq_ID
      )) +
    geom_point(size = 4) +
    geom_text_repel(size = 2.5, max.overlaps = 150) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    facet_wrap(~ data.info$Country.code..ISO.3166., scales = "free") +
    labs(x = "Log(Environmental distance)",
        y = metric,
        color = "Country") +
    theme_bw()

    ggsave(
      file.path(fig.dir, paste0(OUT_FILE, metric, "_", data.type, "_", type, "_FACET.png")),
      plot = plot2,
      width = 12, height = 8, units = "in", dpi = 300
    )
  }
}

# ---- Helper: GO Calculation & Export ----
run_and_save_go <- function(geno, env_curr, env_pred, loci = NULL, label, suffix, info) {
  # Calculate Gap
  g <- genetic.gap(input = geno, env = env_curr, pred.env = env_pred, 
                   scale = TRUE, K = K, candidate.loci = loci)
  
  # Prepare Data
  df <- info
  df$squared.env.dist     <- rowSums((env_curr - env_pred)^2)
  df$offset          <- g$offset
  df$offset_Relstad  <- g$distance
  
  # Generalized File Writing
  write.csv(df, file.path(go_main, "file", paste0("GO_", label, "_", suffix, ".csv")), row.names = FALSE)
  
  # Plotting
  plot_go_env(df, info, type = suffix, file.path(go_main, "figure"), label, paste0("GO_", label, "_"))
}

# ---- Load Genotypes & Loci ----
# Define a list of "Genotype Scenarios" to iterate over
lfmm_path <- file.path("./out", PROJECT, "s02.lfmm2_V3/file")

scenarios <- list(
  Full = list(
    geno = read.geno(file.path(lfmm_path, paste0("Full_chr_pvalue_", PROJECT, ".geno"))),
    meta = read.csv(file.path(lfmm_path, "All_SNPs_Full_Results.csv"))
  ),
  Filt = list(
    geno = read.geno(file.path(lfmm_path, paste0("Full_chr_pvalue_filt_", PROJECT, ".geno"))),
    meta = read.csv(file.path(lfmm_path, "All_SNPs_Filtered_Results.csv"))
  )
)

# ---- Main Execution Loop ----
#env_targets <- list.files(file.path(go_main, "env"), pattern = "env_.*_Chelsa.env$", full.names = TRUE)

for (e_file in new_env_files) {
  print(e_file)
  file_name <- paste0(paste0(basename(sub("\\.csv$", "", e_file)), ".env")); print(file_name)
  tmp <- read.delim(e_file, header = F, sep = "\t", dec = ".")
  write.env(tmp, file.path(go_main, "env", file_name))

  X_pred <- read.env(file.path(go_main, "env", file_name))
  e_name <- basename(sub("\\.env$", "", file_name))
  message(">>> Processing environement: ", e_name)
  
  for (s_name in names(scenarios)) {
    message("Data:", s_name)
    dat <- scenarios[[s_name]]
    
    # 1. All SNPs
    message("Full offset")
    run_and_save_go(dat$geno, X, X_pred, NULL, s_name, paste0(e_name, "_All"), info_ref)
    
    # 2. FDR SNPs
    message("FDR offset")
    fdr_idx <- dat$meta$LFMM_pos[dat$meta$FDR_signif]
    run_and_save_go(dat$geno, X, X_pred, fdr_idx, s_name, paste0(e_name, "_FDR"), info_ref)
    
    # 3. Bonferroni SNPs
    message("Bonferroni offset")
    bonf_idx <- dat$meta$LFMM_pos[dat$meta$Bonf_signif]
    run_and_save_go(dat$geno, X, X_pred, bonf_idx, s_name, paste0(e_name, "_Bonf"), info_ref)
  }
}