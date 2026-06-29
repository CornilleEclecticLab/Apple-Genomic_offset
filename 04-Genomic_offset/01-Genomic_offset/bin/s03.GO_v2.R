#!/usr/bin/env Rscript

# ============================================================
# Script: s03.Genomic_offset.R
# Purpose: Run genomic offset (genetic.gap) per chromosome,
#          then global prediction and plotting
# ============================================================

# ---- Clean environment ----
rm(list = ls())
source("./bin/fonction.R")
#source("./bin/config.R")

# ---- Input arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript bin/s03.Genomic_offset.R <project_name> <env_current.csv> <env_future.csv> <reference_info> <K_value>")
}

PROJECT <- args[1]
env_current_csv <- args[2]
env_future_csv  <- args[3]
ref_info        <- args[4]
K <- as.numeric(args[5])

# PROJECT="dom_admixed_rm"
# env_current_csv="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/dom_admixed_rm/climate_extract_current_dom_admixed_rm.csv"
# env_future_csv="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s03.orchard_generalization/dom_admixed_rm/env_Saclay_Chelsa"
# ref_info="./in/data_tree_dom_admixed_rm_V2.txt"
# K=5


message("🔧 Running Genomic Offset for project: ", PROJECT)

# ---- Define folders ----
comp_main <- paste0("./out/", PROJECT, "/s01.computation/")
lfmm2_main <- paste0("./out/", PROJECT, "/s02.lfmm2/")
go_main <- paste0("./out/", PROJECT, "/s03.GO/")

file_dir <- file.path(go_main, "file")
fig_dir  <- file.path(go_main, "figure")
env_dir  <- file.path(go_main, "env")

create_directory(go_main)
create_directory(file_dir)
create_directory(fig_dir)
create_directory(env_dir)

# ---- Reference info ----
info_ref_all_SYL <- read.table(ref_info, header = TRUE, sep = "\t")
sorted_info = info_ref_all_SYL[order(info_ref_all_SYL$Order_vcf), ]

# ---- Convert environments ----
env_current_name <- basename(sub("\\.csv$", "", env_current_csv))
env_future_name  <- basename(sub("\\.csv$", "", env_future_csv))

# Read and convert CSVs to .env
env_current <- read.delim(env_current_csv, header = TRUE, sep = ",", dec = ".")
env_current <- env_current[match(sorted_info$Seq_ID, env_current$Seq_ID), ]
env_future  <- read.delim(env_future_csv, header = FALSE, sep = "\t", dec = ".")

# print(dim(env_current))
# print(dim(env_future))

# Remove Seq_ID
write.env(env_current[-1], file.path(env_dir, paste0(env_current_name, ".env")))
write.env(env_future, file.path(env_dir, paste0(env_future_name, ".env"))) # there isn't any Seq_ID column here so no [-1]

X <- read.env(file.path(env_dir, paste0(env_current_name, ".env")))
X.pred <- read.env(file.path(env_dir, paste0(env_future_name, ".env")))

# ---- Initialize storage ----
chrom_folders <- list.dirs(comp_main, recursive = FALSE) %>% 
  grep(paste0("chr_chr.*anno_",PROJECT), ., value = TRUE)
message("Found ", length(chrom_folders), " chromosome folders.")

# ---- Run genetic.gap for the full genome ----
full_geno_pvalue_filt <- read.geno(paste0(lfmm2_main, "file/Full_chr.geno20_maf005.anno_",PROJECT,".pvalue_filt.geno"))

full_geno <- NULL

for (chr_folder in chrom_folders) {
  chr_name <- basename(chr_folder)
  message("🧬 Grouping chromosome: ", chr_name)

  # Find and read the geno file
  geno_file <- list.files(file.path(chr_folder, "file"), pattern = "_filt.lfmm_imputed.geno$", full.names = TRUE)
  gene_data <- read.geno(geno_file)

  # Initialize or append
  if (is.null(full_geno)) {
    full_geno <- gene_data
  } else {
    # Check that rows match (same individuals)
    if (nrow(gene_data) != nrow(full_geno)) {
      stop(paste("❌ Row mismatch in", chr_name, ":", 
                 nrow(gene_data), "vs", nrow(full_geno)))
    }
    full_geno <- cbind(full_geno, gene_data)
  }

  message("✅ Done grouping: ", chr_name)
}

# message("Start full genome ...")

# ---- GEA SNPs ----

## new file type All_chr_SNPs_associated_with_all_bio_LFMM2_position_[Bonferroni/FDR].csv
## position in geno file

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
    geom_text_repel(size = 2.5, max.overlaps = 100) +
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


snp_file  <- list.files(paste0(lfmm2_main,"file"), pattern = "All_chr_SNPs_associated_with_all_bio_LFMM2_position_", full.names = TRUE)

for (file in snp_file) {
  message("Working on: ", file)
  
  # 1. Déterminer si c'est filtré et extraire le type
  is_filtered <- grepl("pvalue_filtered", file)
  
  if (is_filtered) {
    type <- sub(".*_([^_]+)_pvalue_filtered\\.csv$", "\\1", basename(file))
    input_geno <- full_geno_pvalue_filt
    suffix_out <- paste0("All_chr_offset_GEA_SNP_pvalue_filtered_", type, ".csv")
    plot_prefix <- "All_chr_Genomic_offset_pvalue_filtered_"
  } else {
    type <- sub(".*_(.*)\\.csv$", "\\1", basename(file))
    input_geno <- full_geno
    suffix_out <- paste0("All_chr_offset_GEA_SNP_", type, ".csv")
    plot_prefix <- "All_chr_Genomic_offset_"
  }
  
  # 2. Traitement commun
  GEA_snp <- read.csv(file)
  candidate.loci <- GEA_snp$x
  
  g.gap.gea <- genetic.gap(input = input_geno, env = X, pred.env = X.pred, 
                           scale = TRUE, K = K, candidate.loci = candidate.loci)
  
  Delta <- X - X.pred
  squared.env.dist <- rowSums(Delta^2)
  data_test_gea <- cbind(squared.env.dist, g.gap.gea$offset, g.gap.gea$distance)
  colnames(data_test_gea) <- c("squared.env.dist", "offset", "offset_Relstad")
  
  # 3. Sauvegarde et Plot
  write.csv(data_test_gea, file.path(file_dir, suffix_out), row.names = FALSE)
  
  plot_go_env(data_test_gea, sorted_info, type = type, fig_dir, 
              "GEA_SNPs", plot_prefix)
}

# ---- All SNPs ----

g.gap.full <- genetic.gap(input = full_geno, env = X, pred.env = X.pred, scale = TRUE, K = K)

Delta <- X - X.pred
squared.env.dist = rowSums(Delta^2)
data_test_full = cbind(squared.env.dist, g.gap.full$offset, g.gap.full$distance)
colnames(data_test_full) = c("squared.env.dist", "offset", "offset_Relstad")
write.csv(data_test_full, file.path(file_dir, "All_chr_offset_all_SNP.csv"), row.names = FALSE)

plot_go_env(
  data_test_full,
  sorted_info,
  type=NULL,
  fig_dir,
  "Full_SNPs",
  "All_chr_Genomic_offset_"
)

g.gap.full.pvalue_filtered <- genetic.gap(input = full_geno_pvalue_filt, env = X, pred.env = X.pred, scale = TRUE, K = K)

Delta <- X - X.pred
squared.env.dist = rowSums(Delta^2)
data_test_full = cbind(squared.env.dist, g.gap.full.pvalue_filtered$offset, g.gap.full.pvalue_filtered$distance)
colnames(data_test_full) = c("squared.env.dist", "offset", "offset_Relstad")
write.csv(data_test_full, file.path(file_dir, "All_chr_offset_all_SNP_pvalue_filtered.csv"), row.names = FALSE)

plot_go_env(
  data_test_full,
  sorted_info,
  type=NULL,
  fig_dir,
  "Full_SNPs",
  "All_chr_Genomic_offset_pvalue_filtered_"
)


message("... End full genome.")
