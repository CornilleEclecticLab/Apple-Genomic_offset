#!/usr/bin/env Rscript

# ============================================================
# Script: s04.Future_GO.R
# Purpose: Run genomic offset (genetic.gap) for future climate
#          conditions
# ============================================================

# ---- Clean environment ----
rm(list = ls())
source("./bin/fonction.R")
#source("./bin/config.R")

# ---- Input arguments ----
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) < 5) {
#   stop("Usage: Rscript bin/s03.Genomic_offset.R <project_name> <env_current.csv> <env_future.csv> <reference_info> <K_value>")
# }

# PROJECT <- args[1]
# env_current_csv <- args[2]
# env_future_csv  <- args[3]
# ref_info        <- args[4]
# K <- as.numeric(args[5])

PROJECT <- "dom_admixed_rm"
env_folder <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/dom_admixed_rm/"
ref_info <- "./in/data_tree_dom_admixed_rm_V2.txt"
K <- 5

message("🔧 Running Genomic Offset for project: ", PROJECT)

# ---- Define folders ----
comp_main <- paste0("./out/", PROJECT, "/s01.computation/")
lfmm2_main <- paste0("./out/", PROJECT, "/s02.lfmm2/")
go_main <- paste0("./out/", PROJECT, "/s04.GO_future/")

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

env_file_present <- list.files(env_folder, pattern="climate_extract_current_*", full.names=T ); print(env_file_present)
env_file_future <- list.files(env_folder, pattern="climate_extract_future_*", full.names=T ); print(env_file_future)
env_file_past <- list.files(env_folder, pattern="climate_extract_LGM_*", full.names=T ); print(env_file_past)


make_env_file <- function(file_list, output_folder, sorted_info){
    for (file in file_list){
        env_name <- basename(sub("\\.csv$", "", file))
        env_data <- read.delim(file, header = TRUE, sep = ",", dec = ".")
        env_data <- env_data[match(sorted_info$Seq_ID, env_data$Seq_ID), ];print(head(env_data))
        write.env(env_data[-1], file.path(output_folder, paste0(env_name, ".env")))
    }
}

make_env_file(env_file_present, env_dir, sorted_info)
make_env_file(env_file_future, env_dir, sorted_info)
make_env_file(env_file_past, env_dir, sorted_info)

# ---- Convert environments ----
env_current_name <- basename(sub("\\.csv$", "", env_current_csv))
env_future_name  <- basename(sub("\\.csv$", "", env_future_csv))

# Read and convert CSVs to .env
env_current <- read.delim(env_current_csv, header = TRUE, sep = "\t", dec = ".")
env_current <- env_current[match(sorted_info$Seq_ID, env_current$Seq_ID), ]
env_future  <- read.delim(env_future_csv, header = FALSE, sep = "\t", dec = ".")

# print(dim(env_current))
# print(dim(env_future))

# Remove Seq_ID
write.env(env_current[-1], file.path(env_dir, paste0(env_current_name, ".env")))
write.env(env_future, file.path(env_dir, paste0(env_future_name, ".env"))) # there isn't any Seq_ID column here so no [-1]

X <- read.env(file.path(env_dir, paste0(env_current_name, ".env")))
X.pred <- read.env(file.path(env_dir, paste0(env_future_name, ".env")))










full_geno_pvalue_filt <- read.geno(paste0(lfmm2_main, "file/Full_chr.geno20_maf005.anno_dom_admixed_rm.pvalue_filt.geno"))


ull_geno <- NULL

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








## genetic gap

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