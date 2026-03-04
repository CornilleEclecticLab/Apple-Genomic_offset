#!/usr/bin/env Rscript
rm(list = ls())
source("./bin/fonction.R")
source("./bin/config.R")

##------------------------------------------------------------
## Input checking
##------------------------------------------------------------
message("Input file checking...")
#check_input_files()
message("...Done.")

##------------------------------------------------------------
## Folder setup
##------------------------------------------------------------
message("Folders creation...")
create_directory(lfmm_folder)
work.directory <- file.path(lfmm_folder, "file/")
fig_output <- file.path(lfmm_folder, "figure/")
env_output <- file.path(lfmm_folder, "env/")
create_directory(work.directory)
create_directory(fig_output)
create_directory(env_output)
message("...Done.")

##------------------------------------------------------------
## Environment file creation
##------------------------------------------------------------
message("Creating env files...")
env_name <- basename(sub("\\.csv$", "", env.1))
out_file <- file.path(env_output, paste0(env_name, ".env"))

info_ref_all_SYL = read.table(ref, header = TRUE, sep = '\t')
sorted_info <- info_ref_all_SYL[order(info_ref_all_SYL$Order_vcf_dom_admixed_rm), ]

env <- read.delim(env.1, header = TRUE, sep = '\t', dec = ".")

env_ordered <- env[match(sorted_info$Seq_ID, env$Seq_ID), ]

env.format.1 <- as.data.frame(scale(env_ordered[, grep("^Bio", colnames(env))], center = FALSE, scale = FALSE))
write.env(env.format.1, out_file)

message(paste("Created:", env_name, ".env"))

message("...Done.")

##------------------------------------------------------------
## Prepare environment matrix
##------------------------------------------------------------
X <- read.env(out_file)

##------------------------------------------------------------
## Find all LFMM files (e.g. per chromosome)
##------------------------------------------------------------
lfmm_files <- list.files(paste0("../00-File_preparation/out/dom_admixed_rm/s01.lfmm_file/"), pattern = "\\.lfmm$", full.names = TRUE)
if (length(lfmm_files) == 0) stop("No .lfmm files found in input folder.")

##------------------------------------------------------------
## Initialize combined results for Manhattan plot
##------------------------------------------------------------
combined_pvalues <- data.frame()
chrom_labels <- c()

##------------------------------------------------------------
## Process each LFMM file (loop over chromosomes)
##------------------------------------------------------------
for (lfmm in lfmm_files) {
  message("Processing: ", basename(lfmm))

  file_name <- basename(sub("\\.lfmm$", "", lfmm))
  snp_names_filt <- read.csv(paste0(computational_folder,"file/", file_name, "_snp_names_filtered.csv"))[,1]
  print(dim(snp_names_filt))
  Y <- read.geno(paste0(computational_folder, "file/", file_name, "_filt.lfmm_imputed.geno"))

  ##--------------------------------------------------------
  ## Run LFMM2
  ##--------------------------------------------------------
  mod.lfmm2 <- lfmm2(input = Y, env = X, K = K)

  ## Plot factors
  png(file.path(fig_output, paste0(file_name, "_plot_mod.lfmm2.png")),
      width = 8, height = 5, units = "in", res = 300)
  plot(mod.lfmm2@U, col = "grey", pch = 19, xlab = "Factor 1", ylab = "Factor 2")
  dev.off()

  ##--------------------------------------------------------
  ## Full test (all variables)
  ##--------------------------------------------------------
  pv <- lfmm2.test(object = mod.lfmm2, input = Y, env = X, full = TRUE)

  # Save results
  df_pv <- data.frame(
    SNPs = snp_names_filt,
    pvalues = as.numeric(pv$pvalues),
    fscores = as.numeric(pv$fscores),
    adj.r.squared = as.numeric(pv$adj.r.squared)
  )
  out_csv <- file.path(work.directory, paste0(file_name, "_GEAs_full.csv"))
  write.csv(df_pv, out_csv, row.names = FALSE)

  # Store pvalues for global Manhattan
  df_pv$chromosome <- gsub("chr", "", file_name)  # assuming file_name like "chr1"
  combined_pvalues <- rbind(combined_pvalues, df_pv)
  chrom_labels <- c(chrom_labels, unique(df_pv$chromosome))

  # Identify significant SNPs
  GEA_snp_all <- which(-log10(pv$pvalues) > -log10(0.1 / length(pv$pvalues)))
  write.csv(GEA_snp_all,
            file = file.path(work.directory, paste0(file_name, "_SNPs_associated_with_all_bio.csv")),
            row.names = FALSE)

  ##--------------------------------------------------------
  ## Variable-specific tests
  ##--------------------------------------------------------
  pv.2 <- lfmm2.test(object = mod.lfmm2, input = Y, env = X, full = FALSE)
  df_all <- data.frame()
  name_plot <- colnames(env.format.1)

  png(file.path(fig_output, paste0(file_name, "_GEAs_plot_all_bio.png")),
      width = 8, height = 5, units = "in", res = 300)
  par(mfrow = c(4, 5), mar = c(2, 2, 2, 2))

  for (i in seq_along(name_plot)) {
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
    #print(head(df))
    if (nrow(df_all) == 0) df_all <- df else df_all <- merge(df_all, df, by = "SNPs")

    plot(-log10(pv.2$pvalues[i, ]), col = "grey", cex = .6, pch = 19,
         main = name_plot[i], xlab = "SNP index", ylab = "-log(pvalue)")
    abline(h = -log10(0.1 / length(pv.2$pvalues[i, ])), lty = 2, col = "orange")
  }
  dev.off()

  write.csv(df_all,
            file = file.path(work.directory, paste0(file_name, "_GEAs_bio_sep.csv")),
            row.names = FALSE)
}

##------------------------------------------------------------
## Global Manhattan plot across chromosomes
##------------------------------------------------------------
message("Creating global Manhattan plot...")
combined_pvalues$logp <- -log10(combined_pvalues$pvalues)
combined_pvalues$chromosome <- as.factor(combined_pvalues$chromosome)

# Create index for continuous x-axis across chromosomes
combined_pvalues$index <- seq_len(nrow(combined_pvalues))

# Compute genome-wide threshold
thresh <- -log10(0.1 / nrow(combined_pvalues))

# Plot
png(file.path(fig_output, "Global_Manhattan_plot_LFMM2.png"),
    width = 12, height = 6, units = "in", res = 300)

plot(combined_pvalues$index, combined_pvalues$logp,
     col = "grey", pch = 19, cex = 0.4,
     xlab = "Chromosome",
     ylab = "-log10(p-value)",
     main = "Global Manhattan plot across all chromosomes",
     xaxt = "n")
axis(1, at = tapply(combined_pvalues$index, combined_pvalues$chromosome, median),
     labels = levels(combined_pvalues$chromosome))
abline(h = thresh, col = "orange", lwd = 2, lty = 2)
dev.off()
message("...Done.")
