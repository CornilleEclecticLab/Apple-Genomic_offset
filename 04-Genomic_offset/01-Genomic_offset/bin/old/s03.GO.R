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


## Paths
work.directory = paste0(output_folder, "file/")
fig_output = paste0(output_folder, "figure/")
file_name <- basename(sub("\\.lfmm$", "", paste0(input_folder, lfmm)))

## Read data
info_ref_all_SYL = read.table(paste0(input_folder, ref), header = TRUE, sep = '\t')
Y <- read.geno(paste0(work.directory, file_name, "_filt.lfmm_imputed.geno"))
X <- read.env(paste0(input_folder, env.1, ".env"))
X.pred <- read.env(paste0(input_folder, env.2, ".env"))

next_env <- read.delim(env.2, header = TRUE, sep = ',', dec = '.')
new_env.format.1 <- scale(next_env[c(2:length(next_env))], center = FALSE, scale = FALSE)
write.env(new_env.format.1, paste0(input_folder, env.2, '.env'))

####----- Prediction -----
message("Prediction with all SNPs...")
g.gap.scaled <- genetic.gap(input = Y, env = X, pred.env = X.pred, scale = TRUE, K = K)
Delta = X[, 1:19] - X.pred[, 1:19]
squared.env.dist = rowSums(Delta^2)

data_test = cbind(squared.env.dist, g.gap.scaled$offset, g.gap.scaled$distance)
colnames(data_test) = c("squared.env.dist", "offset", "offset_Relstad")

sorted_info = info_ref_all_SYL[order(info_ref_all_SYL$ALL_SYL_wo_admixed), ]

## Plots
library(ggplot2)
library(ggrepel)

for (metric in c("offset", "offset_Relstad")) {
  png(paste0(fig_output, file_name, "_Genomic_offset_", metric, ".png"),
      width = 10, height = 10, units = "in", res = 300)
  ggplot(data_test, aes(x = log(squared.env.dist), y = data_test[, metric],
                        color = sorted_info$Country.code..ISO.3166., label = sorted_info$Site)) +
    geom_point(size = 5) +
    geom_text_repel(size = 3, max.overlaps = 40) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    labs(x = "Log(Environmental distance)", y = metric, color = "Country") +
    theme_bw()
  dev.off()
}

####----- Prediction with candidate SNPs -----
GEA_snp <- read.csv(paste0(work.directory, file_name, "_SNPs_associated_with_all_bio.csv"), header = TRUE)
list_GEA_snp <- GEA_snp$x

g.gap.scaled.bis <- genetic.gap(input = Y, env = X, pred.env = X.pred,
                                scale = TRUE, K = K, candidate.loci = list_GEA_snp)

data_test_bis = cbind(squared.env.dist, g.gap.scaled.bis$offset, g.gap.scaled.bis$distance)
colnames(data_test_bis) = c("squared.env.dist", "offset", "offset_Relstad")

for (metric in c("offset", "offset_Relstad")) {
  png(paste0(fig_output, file_name, "_Genomic_offset_GEA_SNP_", metric, ".png"),
      width = 10, height = 10, units = "in", res = 300)
  ggplot(data_test_bis, aes(x = log(squared.env.dist), y = data_test_bis[, metric],
                            color = sorted_info$Country.code..ISO.3166., label = sorted_info$Site)) +
    geom_point(size = 5) +
    geom_text_repel(size = 3, max.overlaps = 40) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    labs(x = "Log(Environmental distance)", y = metric, color = "Country") +
    theme_bw()
  dev.off()
}
message("...Done.")
