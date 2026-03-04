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

# args <- commandArgs(trailingOnly = TRUE)
# if(length(args) < 4){
#   stop("Usage: Rscript bin/s02.postprocess_LFMM.R <project_name> <env_current_csv_file> <position_data>")
# }
# PROJECT <- args[1]
# env_current_csv <- args[2]  # e.g., "../00-File_preparation/out/dom_admixed_rm/env.csv"
# env_future_csv  <- args[3]
# position.data <- args[4]
# K <- as.numeric(args[5])

PROJECT <- "dom_admixed_rm"
env_current_csv <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/dom_admixed_rm/Bioclim_dom_admixed_rm_chelsa_current.csv"
env_future_csv <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/dom_admixed_rm/env_Saclay_Chelsa"
position.data <- "./in/data_tree_dom_admixed_rm_V2.txt"
K <- 5

# -----------------------------
# ---- Define main folders ----
# -----------------------------

comp_main <- paste0("./out/", PROJECT, "/s01.computation/")
chrom_folders <- list.dirs(comp_main, recursive = FALSE) %>% 
  grep("chr_chr.*anno_dom_admixed_rm", ., value = TRUE)


GO_spectrum <- paste0("./out/", PROJECT, "/s04.2.GO_spectrum/")
file_dir <- file.path(GO_spectrum, "file")
fig_dir  <- file.path(GO_spectrum, "figure")
env_dir  <- file.path(GO_spectrum, "env")
top_snps <- file.path(GO_spectrum, "top_snps_subsets")

create_directory(GO_spectrum)
create_directory(file_dir)
create_directory(fig_dir)
create_directory(env_dir)
create_directory(top_snps)

# -----------------------------
# ---- Convert CSV to .env ----
# -----------------------------

env_name <- basename(sub("\\.csv$", "", env_current_csv))
out_env_file <- file.path(env_dir, paste0(env_name, ".env"))

env_data <- read.delim(env_current_csv, header = TRUE, sep = "\t", dec = "."); dim(env_data)
list_accession = read.delim(file = position.data, header = T, sep='\t', dec = ",")
# optional: order rows if needed, e.g., matching Seq_ID with reference
env_ordered <- env_data[match(list_accession$Seq_ID, env_data$Seq_ID), ]
env_ordered_format <- env_ordered[-1]

# Convert to .env format
write.env(env_ordered_format, out_env_file)

# ---- Convert environments ----
env_future_name  <- basename(sub("\\.csv$", "", env_future_csv))
env_future  <- read.delim(env_future_csv, header = FALSE, sep = "\t", dec = ".")
write.env(env_future, file.path(env_dir, paste0(env_future_name, ".env"))) # there isn't any Seq_ID column here so no [-1]

message("Created env file: ", out_env_file)

# -----------------------------
# -- Read environment matrix --
# -----------------------------

X <- read.env(out_env_file)
X.pred <- read.env(file.path(env_dir, paste0(env_future_name, ".env")))

# -----------------------------
#  Initialize combined results 
# -----------------------------

geno_LDpruned <- "./out/dom_admixed_rm/s01.computation/chr_concat.geno20_maf005.anno_dom_admixed_rm.LDpruned/file/chr_concat.geno20_maf005.anno_dom_admixed_rm.LDpruned_filt.lfmm_imputed.geno"
Y_LDpruned <- read.geno(geno_LDpruned)
# Run LFMM2 with the LD pruned dataset, here Y must be Y_LDpruned can be charged outside
mod.lfmm2 <- lfmm2(input = Y_LDpruned, env = X, K = K)

combined_pvalues <- data.frame()
combined_pvalues_rm <- data.frame()
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
  df.pv <- data.frame(SNPs = snp_names,
                      pvalues = as.numeric(pv$pvalues),
                      fscores = as.numeric(pv$fscores),
                      adj.r.squared = as.numeric(pv$adj.r.squared))

  # Raw storing
  df.pv$chromosome  <- sub(".*(chr[0-9]+).*", "\\1", chr_name)
  combined_pvalues <- rbind(combined_pvalues, df.pv)

  # pvalue >= 0.9 removing
  SNPs_rm <- df.pv[df.pv$pvalues >= 0.9,]
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
}

histo_pvalues(
  df = combined_pvalues,
  out_file = paste0(fig_dir,"/Raw_P-value_distribution_histogram_LFMM2.png"))

histo_pvalues( 
  df = combined_pvalues_rm,
  out_file = paste0(fig_dir,"/Filt_P-value_distribution_histogram_LFMM2.png"))

######### Raw Manhattan plot
combined_pvalues$logp <- -log10(combined_pvalues$pvalues)
combined_pvalues$chromosome <- as.factor(combined_pvalues$chromosome)
combined_pvalues$index <- seq_len(nrow(combined_pvalues))
thresh <- -log10(0.1 / nrow(combined_pvalues))

# Ensure chromosome is a factor with proper order
combined_pvalues$chromosome <- factor(combined_pvalues$chromosome, 
                                      levels = unique(combined_pvalues$chromosome))

# Define colors (alternating for readability)
chromosomes <- levels(combined_pvalues$chromosome)
chrom_colors <- rep(c("grey30", "black"), length.out = length(chromosomes))

# Plot
p <- ggplot(combined_pvalues, aes(x = index, y = logp, color = chromosome)) +
  geom_point(size = 0.8) +
  scale_color_manual(values = chrom_colors) +
  geom_hline(yintercept = thresh, color = "orange", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](p-value)),
    title = "Global Manhattan plot across all chromosomes"
  ) +
  scale_x_continuous(
    breaks = tapply(combined_pvalues$index, combined_pvalues$chromosome, median),
    labels = chromosomes
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "grey95", color = NA),   # Background of the whole plot
    panel.background = element_rect(fill = "white", color = "grey70"),  # Background of panel area
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Save plot
ggsave(
  filename = file.path(fig_dir, "Raw_Global_Manhattan_plot_LFMM2.png"),
  plot = p,
  width = 15, height = 6, units = "in", dpi = 300
)


######### Filtered Manhattan plot

combined_pvalues_rm$logp <- -log10(combined_pvalues_rm$pvalues)
combined_pvalues_rm$chromosome <- as.factor(combined_pvalues_rm$chromosome)
combined_pvalues_rm$index <- seq_len(nrow(combined_pvalues_rm))
thresh <- -log10(0.1 / nrow(combined_pvalues_rm))

# Ensure chromosome is a factor with proper order
combined_pvalues_rm$chromosome <- factor(combined_pvalues_rm$chromosome, 
                                      levels = unique(combined_pvalues_rm$chromosome))

# Define colors (alternating for readability)
chromosomes <- levels(combined_pvalues_rm$chromosome)
chrom_colors <- rep(c("grey30", "black"), length.out = length(chromosomes))

# Plot
p <- ggplot(combined_pvalues_rm, aes(x = index, y = logp, color = chromosome)) +
  geom_point(size = 0.8) +
  scale_color_manual(values = chrom_colors) +
  geom_hline(yintercept = thresh, color = "orange", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](p-value)),
    title = "Global Manhattan plot across all chromosomes"
  ) +
  scale_x_continuous(
    breaks = tapply(combined_pvalues_rm$index, combined_pvalues_rm$chromosome, median),
    labels = chromosomes
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "grey95", color = NA),   # Background of the whole plot
    panel.background = element_rect(fill = "white", color = "grey70"),  # Background of panel area
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Save plot
ggsave(
  filename = file.path(fig_dir, "Filt_Global_Manhattan_plot_LFMM2.png"),
  plot = p,
  width = 15, height = 6, units = "in", dpi = 300
)

#combined_pvalues_rm_sorted <- combined_pvalues_rm %>% arrange(pvalues)

# 3. Define the two sequences and combine them
# Part 1: 10k to 100k (by 10k)
# Part 2: 200k to 700k (by 100k)
steps <- c(seq(10000, 100000, by = 10000), 
           seq(200000, 700000, by = 100000))

# 4. Execute the loop
for (n in steps) {
  # Safety check: ensure 'n' does not exceed the total number of SNPs in your data
  if (n > nrow(combined_pvalues_rm)) {
    cat("Warning: Requested", n, "SNPs, but only", nrow(combined_pvalues_rm), "available. Skipping.\n")
    next
  }

  # 1. Get the row indices of the 10,000 smallest p-values
  # This gives you the "original row names" as a numeric vector
  top_indices <- order(combined_pvalues_rm$pvalues)[1:n]

  # 2. If you want to see the data WITHOUT changing the order of the original object:
  top_snps_data <- combined_pvalues_rm[top_indices, ]
  # Save the file
  file_name <- paste0(top_snps,"/top_", n, "_snps.csv")
  write_csv(top_snps_data, file_name)
  
  # cat("Successfully extracted the top", n, "SNPs.\n")
  
  # type <- sub(".*_(.*)\\.csv$", "\\1", basename(file))
  # GEA_snp <- read.csv(file)
  
  candidate.loci <- top_indices
  g.gap.gea <- genetic.gap(input = full_geno, env = X, pred.env = X.pred, scale = TRUE, K = K, candidate.loci = candidate.loci)

  Delta <- X - X.pred
  squared.env.dist = rowSums(Delta^2)
  data_test_gea = cbind(squared.env.dist, g.gap.gea$offset, g.gap.gea$distance)

  colnames(data_test_gea) = c("squared.env.dist", "offset", "offset_Relstad")
  write.csv(data_test_gea, file.path(file_dir, paste0("All_chr_offset_GEA_",n,"_SNPs.csv")), row.names = FALSE)

  plot_go_env(
    data_test_gea,
    list_accession,
    type=paste0(n,"_SNPs"),
    fig_dir,
    "GEA_SNPs"
  )
}