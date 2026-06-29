#!/usr/bin/env Rscript

# -----------------------------
# ---- Clean environment ------
# -----------------------------

rm(list = ls())
source("./bin/fonction.R")
#source("./bin/config.R")

# -----------------------------
# ---- Input arguments --------
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 5){
  stop("Usage: Rscript bin/s02.postprocess_LFMM.R <project_name> <env_csv_file> <position_data> <geno_LDpruned>")
}
PROJECT <- args[1]
env_csv <- args[2]  # e.g., "../00-File_preparation/out/dom_admixed_rm/env.csv"
position.data <- args[3]
K <- as.numeric(args[4])
geno_LDpruned <- args[5]

# PROJECT="dom_admixed_rm"
# env_csv="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/dom_admixed_rm/Bioclim_dom_admixed_rm_chelsa_current.csv"
# position.data="./in/data_tree_dom_admixed_rm_V2.txt"
# K=5
# geno_LDpruned <- "./out/dom_admixed_rm/s01.computation/chr_concat.geno20_maf005.anno_dom_admixed_rm.LDpruned/file/chr_concat.geno20_maf005.anno_dom_admixed_rm.LDpruned_filt.lfmm_imputed.geno"

PROJECT="all_wild"
env_csv="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/all_wild/climate_extract_current_all_wild.csv"
position.data="./in/data_tree_all_wild.txt"
K=6
geno_LDpruned <- "./out/all_wild/s01.computation/chr_concat.geno20_maf005.anno_all_wild.LDpruned/file/chr_concat.geno20_maf005.anno_all_wild.LDpruned_filt.lfmm_imputed.geno"

# -----------------------------
# ---- Define main folders ----
# -----------------------------

comp_main <- paste0("./out/", PROJECT, "/s01.computation/")
chrom_folders <- list.dirs(comp_main, recursive = FALSE) %>% 
  grep(paste0("chr_chr.*anno_",PROJECT), ., value = TRUE) ; print(chrom_folders)

lfmm2_main <- paste0("./out/", PROJECT, "/s02.lfmm2/")
file_dir <- file.path(lfmm2_main, "file")
fig_dir  <- file.path(lfmm2_main, "figure")
env_dir  <- file.path(lfmm2_main, "env")

create_directory(file_dir)
create_directory(fig_dir)
create_directory(env_dir)

# -----------------------------
# ---- Convert CSV to .env ----
# -----------------------------

env_name <- basename(sub("\\.csv$", "", env_csv))
out_env_file <- file.path(env_dir, paste0(env_name, ".env"))

env_data <- read.delim(env_csv, header = TRUE, sep = ",", dec = "."); dim(env_data)
list_accession = read.delim(file = position.data, header = T, sep='\t', dec = ",")
# optional: order rows if needed, e.g., matching Seq_ID with reference
env_ordered <- env_data[match(list_accession$Seq_ID, env_data$Seq_ID), ]
env_ordered_format <- env_ordered[-1]

# Convert to .env format
write.env(env_ordered_format, out_env_file)
message("Created env file: ", out_env_file)

# -----------------------------
# -- Read environment matrix --
# -----------------------------

X <- read.env(out_env_file)

# -----------------------------
# -- Run first lfmm2 command --
# ---- with LD pruned data ----
# -----------------------------

Y_LDpruned <- read.geno(geno_LDpruned)
# Run LFMM2 with the LD pruned dataset, here Y must be Y_LDpruned can be charged outside
mod.lfmm2 <- lfmm2(input = Y_LDpruned, env = X, K = K)

png(file.path(fig_dir, "Plot_mod.lfmm2.png"),
      width = 8, height = 5, units = "in", res = 300)
  plot(mod.lfmm2@U, col = "grey", pch = 19, xlab = "Factor 1", ylab = "Factor 2")
dev.off()

# -----------------------------
#  Initialize combined results 
# -----------------------------

combined_pvalues <- data.frame()
combined_pvalues_rm <- data.frame()
combined_pvalues_bio_sep <- data.frame()
full_geno <- NULL

# -----------------------------
# --- Loop over chromosomes ---
# -----------------------------

for (chr_folder in chrom_folders) {
  chr_name <- basename(chr_folder)
  message("Processing chromosome folder: ", chr_name)

  # Create per-chromosome folders in s02.lfmm2
  chr_file_dir <- file.path(file_dir, chr_name)
  chr_fig_dir  <- file.path(fig_dir, chr_name)
  create_directory(chr_file_dir)
  create_directory(chr_fig_dir)

  # Read filtered SNP names
  snp_csv <- list.files(file.path(chr_folder, "file"), pattern = "_snp_names_filtered.csv$", full.names = TRUE)
  if(length(snp_csv) != 1) stop("Expected 1 SNP CSV in ", chr_folder)
  snp_names <- read.csv(snp_csv[1])[,1] # Get SNPs position in the genome

  # Read imputed geno
  geno_file <- list.files(file.path(chr_folder, "file"), pattern = "_filt.lfmm_imputed.geno$", full.names = TRUE)
  if(length(geno_file) != 1) stop("Expected 1 imputed geno in ", chr_folder)
  Y_chr <- read.geno(geno_file[1])
  
  # Full test
  pv <- lfmm2.test(object = mod.lfmm2, input = Y_chr, env = X, full = TRUE)
  df_pv <- data.frame(SNPs = snp_names,
                      pvalues = as.numeric(pv$pvalues),
                      fscores = as.numeric(pv$fscores),
                      adj.r.squared = as.numeric(pv$adj.r.squared))

  write.csv(df_pv, file = file.path(chr_file_dir, paste0(chr_name, "_GEAs_full.csv")), row.names = FALSE)

  # SNPs position in LFMM
  GEA_snp_chr <- which(-log10(df_pv$pvalues) > -log10(0.05/length(df_pv$pvalues))); length(GEA_snp_chr)
  write.csv(GEA_snp_chr, file.path(chr_file_dir, paste0(chr_name, "_SNPs_associated_with_all_bio_LFMM_Bonferroni.csv")), row.names = F)

  qv_lfmm2 <- qvalue::qvalue(df_pv$pvalues, fdr.level = 0.05)
  candidates <- which(qv_lfmm2$significant)
  write.csv(candidates, file.path(chr_file_dir, paste0(chr_name, "_SNPs_associated_with_all_bio_LFMM_FDR.csv")), row.names = F)

  # SNPs position in vcf
  df_sig <- df_pv[-log10(df_pv$pvalues) > -log10(0.05 / length(df_pv$pvalues)), ]
  write.csv(df_sig, file.path(chr_file_dir, paste0(chr_name, "_SNPs_associated_with_all_bio_VCF_Bonferroni.csv")), row.names = F)

  df_fdr <- df_pv[which(qv_lfmm2$significant),]
  write.csv(df_fdr, file.path(chr_file_dir, paste0(chr_name, "_SNPs_associated_with_all_bio_VCF_FDR.csv")), row.names = F)
  
  # Store for global Manhattan
  df_pv$chromosome  <- sub(".*(chr[0-9]+).*", "\\1", chr_name)
  combined_pvalues <- rbind(combined_pvalues, df_pv)

  # Plot all bio per chromosomes
  plot_density( 
    df = df_pv,
    FIG_TITLE = paste0("P-value distribution - ", sub(".*(chr[0-9]+).*", "\\1", chr_name)),
    OUT_FILE = file.path(chr_fig_dir, paste0(chr_name, "_pvalues_density.png")))

  plot_histo(
    df = df_pv,
    FIG_TITLE = paste0("P-value distribution - ", sub(".*(chr[0-9]+).*", "\\1", chr_name)),
    OUT_FILE = file.path(chr_fig_dir, paste0(chr_name, "_pvalues_histogram.png")))
  # ----------------------
  # -- pvalue filtering --
  # ----------------------

  # pvalue >= 0.9 removing
  SNPs_rm <- df_pv[df_pv$pvalues >= 0.9,]
  indices_to_removes <- as.numeric(rownames(SNPs_rm))
  Y_chr_filtered <- Y_chr[, -indices_to_removes]
  snp_names_filt <- snp_names[!(snp_names %in% SNPs_rm$SNPs)]

  # Initialize or append
  if (is.null(full_geno)) {
    full_geno <- Y_chr_filtered
  } else {
    # Check that rows match (same individuals)
    if (nrow(Y_chr_filtered) != nrow(full_geno)) {
      stop(paste("❌ Row mismatch in", chr_name, ":", 
                 nrow(Y_chr_filtered), "vs", nrow(full_geno)))
    }
    full_geno <- cbind(full_geno, Y_chr_filtered)
  }

  # New test
  pv.2 <- lfmm2.test(object = mod.lfmm2, input = Y_chr_filtered, env = X, full = TRUE)
  df.pv.2 <- data.frame(SNPs = snp_names_filt,
                      pvalues = as.numeric(pv.2$pvalues),
                      fscores = as.numeric(pv.2$fscores),
                      adj.r.squared = as.numeric(pv.2$adj.r.squared))

  # Raw storing
  df.pv.2$chromosome  <- sub(".*(chr[0-9]+).*", "\\1", chr_name)
  combined_pvalues_rm <- rbind(combined_pvalues_rm, df.pv.2)

  plot_density( 
    df = df.pv.2,
    FIG_TITLE = paste0("P-value distribution - ", sub(".*(chr[0-9]+).*", "\\1", chr_name)),
    OUT_FILE = file.path(chr_fig_dir, paste0(chr_name, "_pvalues_filtered_density.png")))

  plot_histo(
    df = df.pv.2,
    FIG_TITLE = paste0("P-value distribution - ", sub(".*(chr[0-9]+).*", "\\1", chr_name)),
    OUT_FILE = file.path(chr_fig_dir, paste0(chr_name, "_pvalues_filtered_histogram.png")))

  # SNPs position in LFMM
  # GEA_snp_chr <- which(-log10(df.pv.2$pvalues) > -log10(0.05/length(df.pv.2$pvalues))); length(GEA_snp_chr)
  # write.csv(GEA_snp_chr, file.path(chr_file_dir, paste0(chr_name, "_SNPs_associated_with_all_bio_LFMM_Bonferroni_pval_rm.csv")), row.names = F)

  # qv_lfmm2 <- qvalue::qvalue(df.pv.2$pvalues, fdr.level = 0.05)
  # candidates <- which(qv_lfmm2$significant)
  # write.csv(candidates, file.path(chr_file_dir, paste0(chr_name, "_SNPs_associated_with_all_bio_LFMM_FDR_pval_rm.csv")), row.names = F)

  # # SNPs position in vcf
  # df_sig <- df.pv.2[-log10(df.pv.2$pvalues) > -log10(0.05 / length(df.pv.2$pvalues)), ]
  # write.csv(df_sig, file.path(chr_file_dir, paste0(chr_name, "_SNPs_associated_with_all_bio_VCF_Bonferroni_pval_rm.csv")), row.names = F)

  # df_fdr <- df.pv.2[which(qv_lfmm2$significant),]
  # write.csv(df_fdr, file.path(chr_file_dir, paste0(chr_name, "_SNPs_associated_with_all_bio_VCF_FDR_pval_rm.csv")), row.names = F)

  ## ---- LFMM2 per Bio variable

  message("LFMM2 per Bio variable")

  pv.3 <- lfmm2.test(object = mod.lfmm2, input = Y_chr_filtered, env = X, full = FALSE)
  df_all <- data.frame()
  name_plot <- colnames(env_data)[2:length(colnames(env_data))]
  
  png(file.path(chr_fig_dir, paste0(chr_name,'_GEAs_plot_bio_sep.png')), width = 8, height = 5, units = "in", res = 300)
  par(mfrow = c(4, 5), mar = c(2, 2, 2, 2))

  for (i in 1:length(name_plot)) {
    
    bio_dir <- file.path(chr_file_dir, name_plot[i])
    create_directory(bio_dir)

    df <- data.frame(
      SNPs = snp_names_filt,
      pvalues = as.numeric(pv.3$pvalues[i, ]),
      qvalues = as.numeric(qvalue(pchisq((pv.3$zscores[i, ] ^ 2) / pv.3$gif,
                                      df = 1, lower = FALSE))$qvalue),
      zsc = as.numeric(pv.3$zscores[i, ])
    )
    message(name_plot[i])

    plot_density(
      df = df,
      FIG_TITLE = paste0("P-value distribution - ", sub(".*(chr[0-9]+).*", "\\1", chr_name), " - " , name_plot[i]),
      OUT_FILE = file.path(chr_fig_dir, paste0(chr_name, "_pvalues_density_",name_plot[i],".png")))

    print(name_plot[i])
    print(dim(df))
    names(df)[2:4] <- c(paste0("pvalues_", name_plot[i]),
                        paste0("qvalue_", name_plot[i]),
                        paste0("zscores_", name_plot[i]))
    
    df_bio_sep <- df[-log10(df$pvalues) > -log10(0.05 / length(df$pvalues)), ]; dim(df_bio_sep)
    write.csv(df_bio_sep, file = file.path(bio_dir, paste0(chr_name, "_GEAs_",name_plot[i],"_position_VCF_Bonferroni.csv")), row.names = FALSE)

    qv_lfmm2_bio_sep <- qvalue::qvalue(df$pvalues, fdr.level = 0.05)
    df_fdr_bio_sep <- df[which(qv_lfmm2_bio_sep$significant),]
    write.csv(df_fdr_bio_sep, file.path(bio_dir, paste0(chr_name, "_GEAs_",name_plot[i],"_position_VCF_FDR.csv")), row.names = F)

    if (nrow(df_all) == 0) df_all <- df else df_all <- merge(df_all, df, by = "SNPs")
    
  }

  dev.off()

  head(df_all)
  write.csv(df_all, file = file.path(chr_file_dir, paste0(chr_name, "_GEAs_bio_sep.csv")), row.names = FALSE)

  df_all$chromosome  <- sub(".*(chr[0-9]+).*", "\\1", chr_name)
  combined_pvalues_bio_sep <- rbind(combined_pvalues_bio_sep, df_all)

}

# -----------------------------
# - SNP significant selection -
# -----------------------------

GEA_snp_all_chr <- which(-log10(combined_pvalues$pvalues) > -log10(0.05/length(combined_pvalues$pvalues))); length(GEA_snp_all_chr) # Bonferroni
write.csv(GEA_snp_all_chr, file.path(file_dir, "All_chr_SNPs_associated_with_all_bio_LFMM2_position_Bonferroni.csv"), row.names = F)

qv_lfmm2_all <- qvalue::qvalue(combined_pvalues$pvalues, fdr.level = 0.05)
candidates_all <- which(qv_lfmm2_all$significant)
write.csv(candidates_all, file.path(file_dir, "All_chr_SNPs_associated_with_all_bio_LFMM2_position_FDR.csv"), row.names = F)

# SNPs position in vcf
df_sig_all <- combined_pvalues[-log10(combined_pvalues$pvalues) > -log10(0.05/ length(combined_pvalues$pvalues)), ]; dim(df_sig_all)
write.csv(df_sig_all, file.path(file_dir, "All_chr_SNPs_associated_with_all_bio_VCF_position_Bonferroni.csv"), row.names = F)

df_fdr_all <- combined_pvalues[which(qv_lfmm2_all$significant),]; dim(df_fdr_all)
write.csv(df_fdr_all, file.path(file_dir, "All_chr_SNPs_associated_with_all_bio_VCF_position_FDR.csv"), row.names = F)

length(intersect(df_fdr_all$SNPs, df_sig_all$SNPs))

combined_pvalues_final <- combined_pvalues %>%
dplyr::mutate(
  Bonf_signif = pvalues < 0.05 / length(pvalues),
  FDR_signif  = qv_lfmm2_all$significant
)

write.csv(combined_pvalues_final,file.path(file_dir, "All_chr_SNPs.csv"), row.names = F)

# --------------------------
# -- with pvalue filtered --
# --------------------------

GEA_snp_all_chr <- which(-log10(combined_pvalues_rm$pvalues) > -log10(0.05/length(combined_pvalues_rm$pvalues))); length(GEA_snp_all_chr) # Bonferroni
write.csv(GEA_snp_all_chr, file.path(file_dir, "All_chr_SNPs_associated_with_all_bio_LFMM2_position_Bonferroni_pvalue_filtered.csv"), row.names = F)

qv_lfmm2_all <- qvalue::qvalue(combined_pvalues_rm$pvalues, fdr.level = 0.05, lambda = seq(0.05, 0.89, 0.05))
# we specified the lambda because of the truncation done before for pvalues >= 0.9
message(paste0("Here the pi0 estimated (% of SNPs not significant): ", 
               round(qv_lfmm2_all$pi0, 4), 
               ".\nIt's used to see if the truncation doesn't overestimate significant SNPs."))
candidates_all <- which(qv_lfmm2_all$significant)
write.csv(candidates_all, file.path(file_dir, "All_chr_SNPs_associated_with_all_bio_LFMM2_position_FDR_pvalue_filtered.csv"), row.names = F)

# SNPs position in vcf
df_sig_all <- combined_pvalues_rm[-log10(combined_pvalues_rm$pvalues) > -log10(0.05/ length(combined_pvalues_rm$pvalues)), ]; dim(df_sig_all)
write.csv(df_sig_all, file.path(file_dir, "All_chr_SNPs_associated_with_all_bio_VCF_position_Bonferroni_pvalue_filtered.csv"), row.names = F)

df_fdr_all <- combined_pvalues_rm[which(qv_lfmm2_all$significant),]; dim(df_fdr_all)
write.csv(df_fdr_all, file.path(file_dir, "All_chr_SNPs_associated_with_all_bio_VCF_position_FDR_pvalue_filtered.csv"), row.names = F)

length(intersect(df_fdr_all$SNPs, df_sig_all$SNPs))

combined_pvalues_rm_final <- combined_pvalues_rm %>%
dplyr::mutate(
  Bonf_signif = pvalues < 0.05 / length(pvalues),
  FDR_signif  = qv_lfmm2_all$significant
)

write.csv(combined_pvalues_rm_final,file.path(file_dir, "All_chr_SNPs_pvalue_filtered.csv"), row.names = F)

# --------------------------------
# - write full genome Y_filtered -
# --------------------------------

write.geno(full_geno, file.path(file_dir,paste0("Full_chr.geno20_maf005.anno_",PROJECT,".pvalue_filt.geno")))
# -----------------------------
# ---- Global pvalues plot ----
# -----------------------------

plot_density( 
  df = combined_pvalues,
  FIG_TITLE = "P-value distribution all chromosome",
  OUT_FILE = file.path(fig_dir, "P-value_distribution_density_LFMM2.png"))

plot_histo(
  df = combined_pvalues,
  FIG_TITLE = "P-value distribution all chromosome",
  OUT_FILE = file.path(fig_dir, "P-value_distribution_histogram_LFMM2.png"))

# -----------------------------
# --- Global Manhattan plot ---
# -----------------------------

# combined_pvalues$logp <- -log10(combined_pvalues$pvalues)
# combined_pvalues$chromosome <- as.factor(combined_pvalues$chromosome)
# combined_pvalues$index <- seq_len(nrow(combined_pvalues))
# thresh <- -log10(0.1 / nrow(combined_pvalues))

# # Ensure chromosome is a factor with proper order
# combined_pvalues$chromosome <- factor(combined_pvalues$chromosome, 
#                                       levels = unique(combined_pvalues$chromosome))

# # Define colors (alternating for readability)
# chromosomes <- levels(combined_pvalues$chromosome)
# chrom_colors <- rep(c("grey30", "black"), length.out = length(chromosomes))

# # Plot
# p <- ggplot(combined_pvalues, aes(x = index, y = logp, color = chromosome)) +
#   geom_point(size = 0.8) +
#   scale_color_manual(values = chrom_colors) +
#   geom_hline(yintercept = thresh, color = "orange", linetype = "dashed", linewidth = 1) +
#   labs(
#     x = "Chromosome",
#     y = expression(-log[10](p-value)),
#     title = "Global Manhattan plot across all chromosomes"
#   ) +
#   scale_x_continuous(
#     breaks = tapply(combined_pvalues$index, combined_pvalues$chromosome, median),
#     labels = chromosomes
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.background = element_rect(fill = "grey95", color = NA),   # Background of the whole plot
#     panel.background = element_rect(fill = "white", color = "grey70"),  # Background of panel area
#     panel.grid.major = element_line(color = "grey85"),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "none"
#   )

# # Save plot
# ggsave(
#   filename = file.path(fig_dir, "Global_Manhattan_plot_LFMM2.png"),
#   plot = p,
#   width = 15, height = 6, units = "in", dpi = 300
# )

# -----------------------------------------
# Manhattan colored for FDR and Bonferroni
# -----------------------------------------

manhattan_plot_bonf_FDR(
  df = combined_pvalues_final,
  TITLE = "Global Manhattan plot across all chromosomes", 
  OUT_FILE = file.path(fig_dir, "Global_Manhattan_plot_LFMM2_log10p.png"))

manhattan_plot_bonf_FDR(
  df = combined_pvalues_rm_final,
  TITLE = "Global Manhattan plot across all chromosomes", 
  OUT_FILE = file.path(fig_dir, "Global_Manhattan_plot_LFMM2_log10p_pvalue_filtered.png"))



# -----------------------------
# ---- Bio separated plot -----
# -----------------------------

for (bio in name_plot) {

  cat("\n====================\n")
  cat("Results for", bio, "\n")
  cat("====================\n")

  # -------------------------------
  # Data preparation
  # -------------------------------

  tmp <- combined_pvalues_bio_sep[, c(
    "SNPs",
    "chromosome",
    grep(paste0("_", bio, "$"),
         colnames(combined_pvalues_bio_sep),
         value = TRUE)
  )]

  pcol <- paste0("pvalues_", bio)

  tmp$chromosome <- factor(tmp$chromosome)
  tmp$logp <- -log10(tmp[[pcol]])

  # Sécurité pvalues == 0
  tmp$logp[is.infinite(tmp$logp)] <-
    max(tmp$logp[is.finite(tmp$logp)]) + 1

  thresh <- -log10(0.05 / nrow(tmp))

  # -------------------------------
  # Manhattan index WITH GAP (METHOD 1)
  # -------------------------------

  gap_size <- 50

  tmp <- tmp %>%
    dplyr::arrange(chromosome, SNPs) %>%
    dplyr::group_by(chromosome) %>%
    dplyr::mutate(chr_index = row_number()) %>%
    dplyr::ungroup()

  chr_sizes <- tmp %>%
    dplyr::group_by(chromosome) %>%
    dplyr::summarise(chr_len = max(chr_index)) %>%
    dplyr::mutate(chr_start = cumsum(lag(chr_len + gap_size, default = 0)))

  tmp_final <- tmp %>%
    dplyr::left_join(chr_sizes, by = "chromosome") %>%
    dplyr::mutate(index = chr_index + chr_start)

  # -------------------------------
  # Significance
  # -------------------------------

  qv_lfmm2_bio <- qvalue::qvalue(tmp_final[[pcol]], fdr.level = 0.05)

  tmp_final <- tmp_final %>%
    dplyr::mutate(
      Bonf_signif = .data[[pcol]] < 0.05 / nrow(tmp_final),
      FDR_signif  = qv_lfmm2_bio$significant
    )

  # -------------------------------
  # Chromosome labels & colors
  # -------------------------------

  chromosomes <- levels(tmp_final$chromosome)

  chrom_centers <- chr_sizes %>%
    dplyr::mutate(center = chr_start + chr_len / 2)

  chrom_colors <- setNames(
    rep(c("grey30", "black"), length.out = length(chromosomes)),
    chromosomes
  )

  # -------------------------------
  # Manhattan plot (simple)
  # -------------------------------

  p <- ggplot(tmp_final, aes(x = index, y = logp, color = chromosome)) +
    geom_point(size = 0.8) +
    scale_color_manual(values = chrom_colors) +
    geom_hline(yintercept = thresh,
               color = "orange",
               linetype = "dashed",
               linewidth = 1) +
    scale_x_continuous(
      breaks = chrom_centers$center,
      labels = chrom_centers$chromosome
    ) +
    labs(
      x = "Chromosome",
      y = expression(-log[10](p-value)),
      title = paste0(bio, " – Manhattan plot (LFMM2)")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.background  = element_rect(fill = "grey95", color = NA),
      panel.background = element_rect(fill = "white", color = "grey70"),
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "none"
    )

  ggsave(
    filename = file.path(fig_dir, paste0(bio, "_Manhattan_plot_LFMM2.png")),
    plot = p,
    width = 15,
    height = 6,
    units = "in",
    dpi = 300
  )

  # -------------------------------
  # Manhattan plot Bonf + FDR
  # -------------------------------

  chrom_colors2 <- setNames(
    rep(c("grey30", "grey60"), length.out = length(chromosomes)),
    chromosomes
  )

  p <- ggplot(tmp_final, aes(x = index, y = logp)) +

    geom_point(
      data = subset(tmp_final, !Bonf_signif & !FDR_signif),
      aes(color = chromosome),
      size = 0.8
    ) +

    geom_point(
      data = subset(tmp_final, Bonf_signif & !FDR_signif),
      color = "blue",
      size = 1.2
    ) +

    geom_point(
      data = subset(tmp_final, !Bonf_signif & FDR_signif),
      color = "red",
      size = 1.2
    ) +

    geom_point(
      data = subset(tmp_final, Bonf_signif & FDR_signif),
      color = "yellow",
      size = 1.4
    ) +

    geom_hline(
      yintercept = thresh,
      color = "orange",
      linetype = "dashed",
      linewidth = 1
    ) +

    scale_color_manual(values = chrom_colors2) +
    scale_x_continuous(
      breaks = chrom_centers$center,
      labels = chrom_centers$chromosome
    ) +

    labs(
      x = "Chromosome",
      y = expression(-log[10](p-value)),
      title = paste0(bio, " – Manhattan plot (Bonferroni + FDR)")
    ) +

    theme_minimal(base_size = 14) +
    theme(
      plot.background  = element_rect(fill = "grey95", color = NA),
      panel.background = element_rect(fill = "white", color = "grey70"),
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "none"
    )

  ggsave(
    filename = file.path(fig_dir,
                         paste0(bio, "_Manhattan_plot_LFMM2_Bonf_FDR.png")),
    plot = p,
    width = 20,
    height = 10,
    units = "in",
    dpi = 300
  )
}


message("✅ LFMM2 post-processing complete in ", lfmm2_main)
