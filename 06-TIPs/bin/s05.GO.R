#!/usr/bin/env Rscript
rm(list = ls())
source("./bin/fonction.R")
source("./bin/config.R")

# Project <- "dom_admixed_rm"
# K=5
# output_folder <- paste0("./out/", Project, "/s03.GO/")
# computational_folder <- paste0("./out/", Project, "/s01.computation/"
# lfmm_folder <- paste0("./out/", Project, "/s02.lfmm/")
# 

message("Folders creation...")
create_directory(go_folder)
work.directory <- paste0(go_folder, "file/")
fig_output <- paste0(go_folder, "figure/")
env_output <- paste0(go_folder, 'env/')
create_directory(work.directory)
create_directory(fig_output)
create_directory(env_output)
message("...Done")


file_name <- basename(sub("\\.lfmm$", "", lfmm))

## Read data
info_ref_all_SYL = read.table(ref, header = TRUE, sep = '\t')

if (file.exists(paste0(computational_folder,"file/", file_name, "_filt.lfmm_imputed.geno"))){
  Y <- read.geno(paste0(computational_folder,"file/", file_name, "_filt.lfmm_imputed.geno"))
} else {
  Y <- read.geno(paste0(computational_folder,"file/", file_name, "_filt.geno"))
}

env_name <- basename(sub("\\.csv$", "", env.1))
X <- read.env(paste0(lfmm_folder, "env/", env_name,".env"))

# Define names
env_name <- basename(sub("\\.csv$", "", env.2))
out_file <- file.path(env_output, paste0(env_name,".env"))

# Check if file already exists
next_env <- read.delim(env.2, header = FALSE, sep = '\t', dec = '.')
env.format.2 <- scale(next_env[c(1:length(next_env))], center = FALSE, scale = FALSE)
write.env(env.format.2, out_file)
message(paste("Created:", env_name, ".env"))

X.pred <- read.env(paste0(env_output, env_name, ".env"))

# next_env <- read.delim(env.2, header = TRUE, sep = ',', dec = '.')
# new_env.format.2 <- scale(next_env[c(2:length(next_env))], center = FALSE, scale = FALSE)
# write.env(new_env.format.2, paste0(input_folder, env.2, '.env'))

####----- Prediction -----
message("Prediction with all SNPs...")

# Compute genomic offset with all SNPs
g.gap.scaled <- genetic.gap(input = Y, env = X, pred.env = X.pred, scale = TRUE, K = K)

# Compute squared environmental distances
Delta <- X[, 1:19] - X.pred[, 1:19]
squared.env.dist <- rowSums(Delta^2)

# Prepare data frame for plotting
data_test <- cbind(
  squared.env.dist, 
  g.gap.scaled$offset, 
  g.gap.scaled$distance
)
colnames(data_test) <- c("squared.env.dist", "offset", "offset_Relstad")
write.table(x=data_test, paste0(work.directory, file_name, "_all_SNPs_offset_prediction.tsv"), row.names=F)


# Sort sample info
sorted_info <- info_ref_all_SYL[order(info_ref_all_SYL$Order_vcf), ]

# Ensure figure output folder exists
if(!dir.exists(fig_output)) dir.create(fig_output, recursive = TRUE)

# --- Plots for all SNPs ---
library(ggplot2)
library(ggrepel)

for (metric in c("offset", "offset_Relstad")) {
  
  p <- ggplot(data_test, aes(
    x = log(squared.env.dist),
    y = .data[[metric]],
    color = sorted_info$Country.code..ISO.3166.,
    label = sorted_info$Site
  )) +
    geom_point(size = 5) +
    geom_text_repel(size = 3, max.overlaps = 40) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    labs(x = "Log(Environmental distance)", y = metric, color = "Country") +
    theme_bw()
  
  # Save PNG
  ggsave(
    filename = paste0(fig_output, file_name, "_Genomic_offset_", metric, ".png"),
    plot = p,
    width = 10, height = 10, units = "in", dpi = 300
  )
  
  # Save SVG
  ggsave(
    filename = paste0(fig_output, file_name, "_Genomic_offset_", metric, ".svg"),
    plot = p,
    width = 10, height = 10, units = "in"
  )
}

####----- Prediction with candidate SNPs -----
message("Prediction with candidate SNPs...")

# Load candidate SNPs
GEA_snp <- read.csv(paste0(lfmm_folder, "file/", file_name, "_SNPs_associated_with_all_bio.csv"), header = TRUE)
list_GEA_snp <- GEA_snp$x

# Compute genomic offset with candidate SNPs
g.gap.scaled.bis <- genetic.gap(
  input = Y, env = X, pred.env = X.pred,
  scale = TRUE, K = K, candidate.loci = list_GEA_snp
)

# Prepare data frame for plotting
data_test_bis <- cbind(
  squared.env.dist,
  g.gap.scaled.bis$offset,
  g.gap.scaled.bis$distance
)
colnames(data_test_bis) <- c("squared.env.dist", "offset", "offset_Relstad")
write.table(x=data_test_bis, paste0(work.directory, file_name, "_GEA_SNPs_offset_prediction.tsv"), row.names=F)

# --- Plots for candidate SNPs ---
for (metric in c("offset", "offset_Relstad")) {
  
  p <- ggplot(data_test_bis, aes(
    x = log(squared.env.dist),
    y = .data[[metric]],
    color = sorted_info$Country.code..ISO.3166.,
    label = sorted_info$Site
  )) +
    geom_point(size = 5) +
    geom_text_repel(size = 3, max.overlaps = 40) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    labs(x = "Log(Environmental distance)", y = metric, color = "Country") +
    theme_bw()
  
  # Save PNG
  ggsave(
    filename = paste0(fig_output, file_name, "_Genomic_offset_GEA_SNP_", metric, ".png"),
    plot = p,
    width = 10, height = 10, units = "in", dpi = 300
  )
  
  # Save SVG
  ggsave(
    filename = paste0(fig_output, file_name, "_Genomic_offset_GEA_SNP_", metric, ".svg"),
    plot = p,
    width = 10, height = 10, units = "in"
  )
}

message("...Done.")
