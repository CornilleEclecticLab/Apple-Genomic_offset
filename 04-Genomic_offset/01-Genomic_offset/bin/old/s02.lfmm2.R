#!/usr/bin/env Rscript
rm(list = ls())
source("./bin/fonction.R")
source("./bin/config.R")

# output_folder <- paste0("./out/", Project, "/s02.lfmm/")
# input_folder <- "./in/"
# lfmm <- "chr_concat.geno20_maf005.anno_filtered_dom_admixed_rm.lfmm"
# env.1 <- "Bioclim_orchard_chelsa_current.csv"
# vcf <- "chr_concat.geno20_maf005.anno_filtered.vcf.gz"
# computation_folder <- paste0("./out/", Project,"/s01.computation/file/")
# Y from s01.computation
# K=5
# Y <- read.geno(paste0(computation_folder, file_name, "_filt.lfmm_imputed.geno"))


## Input checking
message("Input file checking...")
check_input_files() # Need to check if the file exist, not needed in the input 
message("...Done.")

## Output folder creation
message("Folders creation...")
create_directory(lfmm_folder)
work.directory <- paste0(lfmm_folder, "file/")
fig_output <- paste0(lfmm_folder, "figure/")
env_output <- paste0(lfmm_folder, 'env/')
create_directory(work.directory)
create_directory(fig_output)
create_directory(env_output)
message("...Done")

## Paths
work.directory = paste0(lfmm_folder, "file/")
fig_output = paste0(lfmm_folder, "figure/")
file_name <- basename(sub("\\.lfmm$", "", lfmm))

## Read SNPs and geno data
snp_names_filt <- data.frame(extract_SNPs_names(vcf))
Y <- read.geno(paste0(computational_folder, file_name, "_filt.lfmm_imputed.geno"))

####----- Make env file -----
message("Creating env files...") # make a loop and create it if it doesn't exist

# Define names
env_name <- basename(sub("\\.csv$", "", env.1))
out_file <- file.path(env_output, paste0(env_name, ".env"))

# Check if file already exists
if (file.exists(out_file)) {
  message(paste("Skipping:", env_name, "— .env file already exists."))
} else {
  env <- read.delim(env.1, header = TRUE, sep = '\t', dec = ".")
  env.format.1 <- scale(env[c(3:length(env))], center = FALSE, scale = FALSE)
  write.env(env.format.1, out_file)
  message(paste("Created:", env_name, ".env"))
}
message("...Done.")


####----- LFMM2 -----
message("Running LFMM2...")
X <- read.env(paste0(env_output, env_name, '.env'))
mod.lfmm2 <- lfmm2(input = Y, env = X, K = K)

png(paste0(fig_output, file_name, '_plot_mod.lfmm2.png'), width = 8, height = 5, units = "in", res = 300)
plot(mod.lfmm2@U, col = "grey", pch = 19, xlab = "Factor 1", ylab = "Factor 2")
dev.off()

# Manhattan plot
pv <- lfmm2.test(object = mod.lfmm2, input = Y, env = X, full = TRUE)
png(paste0(fig_output, file_name, '_manhattan_plot.lfmm2.png'), width = 8, height = 5, units = "in", res = 300)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19,
     xlab = "SNP Index", ylab = "-log10(pvalue)",
     main = "Manhattan plot of log10 p-value of each SNP for all climate variables")
abline(h = -log10(0.1/length(pv$pvalues)), lty = 1, col = "orange", lwd = 2)
dev.off()

df_pv <- data.frame(
  SNPs = snp_names_filt,
  pvalues = as.numeric(pv$pvalues),
  fscores = as.numeric(pv$fscores),
  adj.r.squared = as.numeric(pv$adj.r.squared)
)
write.csv(df_pv, file = paste0(work.directory, file_name, "_GEAs_full.csv"), row.names = FALSE)

GEA_snp_all <- which(-log10(pv$pvalues) > -log10(0.1 / length(pv$pvalues)))
write.csv(GEA_snp_all, file = paste0(work.directory, file_name, "_SNPs_associated_with_all_bio.csv"), row.names = FALSE)

####----- Variable-specific plots -----
pv.2 <- lfmm2.test(object = mod.lfmm2, input = Y, env = X, full = FALSE)
df_all <- data.frame()
name_plot <- colnames(env.format.1)

png(paste0(fig_output, file_name, '_GEAs_plot_all_bio.png'), width = 8, height = 5, units = "in", res = 300)
par(mfrow = c(4, 5), mar = c(2, 2, 2, 2))

for (i in 1:length(name_plot)) {
  df <- data.frame(
    SNPs = snp_names_filt,
    pval = as.numeric(pv.2$pvalues[i, ]),
    qval = as.numeric(qvalue(pchisq((pv.2$zscores[i, ] ^ 2) / pv.2$gif,
                                    df = 1, lower = FALSE))$qvalue),
    zsc = as.numeric(pv.2$zscores[i, ])
  )
  names(df)[2:4] <- c(paste0("pvalues_", name_plot[i]),
                      paste0("qvalue_", name_plot[i]),
                      paste0("zscores_", name_plot[i]))
  if (nrow(df_all) == 0) df_all <- df else df_all <- merge(df_all, df, by = "SNPs")
  
  plot(-log10(pv.2$pvalues[i, ]), col = "grey", cex = .6, pch = 19,
       main = name_plot[i], xlab = "SNP index", ylab = "-log(pvalue)")
  abline(h = -log10(0.1 / length(pv.2$pvalues[i, ])), lty = 2, col = "orange")
}
dev.off()

write.csv(df_all, file = paste0(work.directory, file_name, "_GEAs_bio_sep.csv"), row.names = FALSE)
message("...Done.")
