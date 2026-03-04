#!/usr/bin/env Rscript
# full pipeline: merge offsets, compute growth, correlations, and plots for many offset files
rm(list = ls())
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)

source("./bin/fonction.R")
source("./bin/config.R")

# Output root
OUT_FOLD <- paste0("./out/",Project,"/s04.Validation/")
create_directory(OUT_FOLD)

# Input constants you already used
# SNP_data <- 'GEA_SNPs'
# message("SNP dataset: ", SNP_data)

# Reference files (leave as in your original script)
info_ref_all_SYL <- read.table(ref, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
sorted_info <- info_ref_all_SYL[order(info_ref_all_SYL$Order_vcf), ]

ind_info <- read.delim("./in/Individuals_mother_information_Saclay.tsv", header = TRUE, stringsAsFactors = FALSE)
pheno_trait <- read.table("./in/Phenotypic_traits_information_Saclay.tsv", header = TRUE, sep = '\t', stringsAsFactors = FALSE)

# Variables to test correlations against (target variables)
target_vars <- c(
  "RGR", "AGR",
  "Summer_Trunk_Measurement", "Winter_Trunk_Measurement",
  "Avg_Chl", "Avg_Flav", "Avg_NBI"
)

# We'll calculate correlations vs both offset columns
offset_vars <- c("offset", "offset_Relstad")

# Helper: safe cor.test wrapper
safe_cor_test <- function(x, y) {
  ok <- complete.cases(x, y)
  n <- sum(ok)
  if (n < 3) {
    return(list(cor = NA_real_, p.value = NA_real_, n = n))
  }
  ct <- tryCatch(cor.test(x[ok], y[ok], method = "pearson"), error = function(e) NULL)
  if (is.null(ct)) return(list(cor = NA_real_, p.value = NA_real_, n = n))
  list(cor = ct$estimate[[1]], p.value = ct$p.value, n = n)
}

# Provide vector of offset prediction files (replace with your actual 18 paths)
# Examples: paste0(go_folder, "file/offset_file1.csv"), ...
offset_files <- c(
  paste0(go_folder, "file/All_chr_offset_all_SNP.csv"), # keep at least one example
  paste0(go_folder, "file/All_chr_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr01.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr01.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr02.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr02.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr03.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr03.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr04.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr04.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr05.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr05.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr06.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr06.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr07.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr07.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr08.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr08.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr09.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr09.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr10.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr10.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr11.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr11.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr12.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr12.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr13.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr13.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr14.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr14.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr15.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr15.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr16.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr16.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv"),
  paste0(go_folder, "file/chr_chr17.geno20_maf005.anno_dom_admixed_rm_offset_GEA_SNP.csv"),
  paste0(go_folder, "file/chr_chr17.geno20_maf005.anno_dom_admixed_rm_offset_allSNP.csv")
)


# Run for all offset files (loop)
results_list <- list()
for (f in offset_files) {
  if (!file.exists(f)) {
    warning("Offset file not found, skipping: ", f)
    next
  }
  results_list[[basename(f)]] <- process_one_file(f)
}

# Optionally: write a combined correlations file across all offset files
all_cors <- bind_rows(lapply(results_list, function(x) x$cor_summary))
if (nrow(all_cors) > 0) {
  write_csv(all_cors, file.path(OUT_FOLD, "all_offset_files_correlations_combined.csv"))
}


# ==== Load correlation table ====
#corr_df <- read.csv(paste0(OUT_FOLD, "all_offset_files_correlations_combined.csv"), header = TRUE)
corr_df <- all_cors
# ------------------------------------------------------------------
#   Filename cleaning function
# ------------------------------------------------------------------
clean_filename <- function(x) {

  x <- gsub("^chr_", "", x)  # remove prefix "chr_"
  
  # remove long middle block
  x <- gsub("\\.geno20_maf005\\.anno_dom_admixed_rm_offset_", "_", x)
  
  # remove trailing underscore if created
  x <- gsub("_$", "", x)

  return(x)
}

# Apply it directly in the dataframe
corr_df$file_clean <- clean_filename(corr_df$file)

# ==== Output dir ====
global_summary_dir <- paste0(OUT_FOLD,"GLOBAL_SUMMARY")
dir.create(global_summary_dir, showWarnings = FALSE, recursive = TRUE)

# ==== Unique target variables and years ====
target_vars <- unique(corr_df$target_var)
year_types <- unique(corr_df$Year)

# ======================================================================
#   BARPLOTS PER TARGET VARIABLE, SEPARATED BY OFFSET AND BY YEAR
# ======================================================================

for (tv in target_vars) {
  for (yr in year_types) {

    df_sub <- corr_df %>% 
      dplyr::filter(target_var == tv, Year == yr)

    if (nrow(df_sub) == 0) next

    # Reorder files inside offset_var panels by correlation value
    df_sub <- df_sub %>%
      group_by(offset_var) %>%
      mutate(file_clean = factor(file_clean, levels = unique(file_clean[order(cor)]))) %>%
      ungroup()

    p_bar <- ggplot(df_sub, aes(x = file_clean, y = cor, fill = offset_var)) +
      geom_col(alpha = 0.75) +
      geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
      facet_wrap(~ offset_var, scales = "free_x") +
      theme_bw(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none"
      ) +
      labs(
        title = paste0("Correlation per file for ", tv, " (Year: ", yr, ")"),
        x = "File",
        y = "Correlation value"
      )

    ggsave(
      file.path(global_summary_dir, 
                paste0("CORR_bar_", tv, "_Year_", yr, "_by_offset.png")),
      p_bar,
      width = 20, height = 20, dpi = 300
    )

    # ---------- DISTRIBUTION PLOT ----------
    p_dist <- ggplot(df_sub, aes(x = cor, fill = offset_var)) +
      geom_histogram(position = "identity", bins = 30, alpha = 0.7) +
      facet_wrap(~ offset_var, scales = "free_y") +
      theme_bw(base_size = 16) +
      theme(
        legend.position = "bottom"
      ) +
      labs(
        title = paste0("Distribution of correlation per file for ", tv, " (Year: ", yr, ")"),
        x = "Correlation value",
        y = "Count",
        fill = "Offset type"
      )

    ggsave(
      file.path(global_summary_dir, paste0("CORR_dist_", tv, "_Year_", yr, "_by_offset.png")),
      p_dist,
      width = 15, height = 10, dpi = 300
    )
  }
}
