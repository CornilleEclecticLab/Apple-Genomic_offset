#!/usr/bin/env Rscript
rm(list = ls())
source("./bin/fonction.R")
source("./bin/config.R")


# Project <- "dom_admixed_rm"
# ref <- "dom_admixed_rm_info_sample.txt"
# vcf <- "chr_concat.geno20_maf005.anno_filtered.vcf.gz"
# lfmm <- "chr_concat.geno20_maf005.anno_filtered_dom_admixed_rm.lfmm"
# On peut mettre les liens complets pour les inputs mais faut retirer la verif des inputs et les paste0 avec le input folder ??
# output_folder <- paste0("./out/", Project,"/s01.computation/")
# input_folder <- "./in/"
# K=5
# faire des runs bash pour chaque s0 ...


## Input checking
message("Input file checking...")
check_input_files() # Need to check if the file exist, not needed in the input 
message("...Done.")

## Output folder creation
message("Folders creation...")
create_directory(computational_folder)
work.directory = paste0(computational_folder, "file/")
fig_output = paste0(computational_folder, "figure/")
create_directory(work.directory)
create_directory(fig_output)
message("...Done")

## SNPs name extraction
snp_names <- data.frame(extract_SNPs_names(vcf))

## File name
file_name <- basename(sub("\\.lfmm$", "", lfmm))

## lfmm to geno
message("Converting lfmm to geno...")
lfmm2geno(lfmm, paste0(work.directory, file_name, ".geno"))
geno.file = paste0(work.directory, file_name, ".geno")
geno.data = read.geno(geno.file)
message("...Done.")

####----- FILTERING -----
message("FILTERING ...")
all_equal <- apply(geno.data, 2, function(col) length(unique(col)) == 1)
one_different <- apply(geno.data, 2, check_one_different)
column_to_remove <- which(all_equal | one_different)
print(paste("There is", length(column_to_remove), "column(s) to remove"))

geno.data_filt <- geno.data[, -column_to_remove]
snp_names_filt <- snp_names[-column_to_remove, ]

write.geno(geno.data_filt, paste0(work.directory, file_name, "_filt.geno"))
geno2lfmm(paste0(work.directory, file_name, "_filt.geno"), paste0(work.directory, file_name, "_filt.lfmm"))
message("...Done.")

####----- IMPUTATION -----
message("Imputation...")
project.missing = snmf(
  paste0(work.directory, file_name, "_filt.lfmm"),
  K = K,
  entropy = TRUE,
  repetitions = nb_rep,
  project = "new"
)

best_run = which.min(cross.entropy(project.missing, K = K))
impute(
  project.missing,
  paste0(work.directory, file_name, "_filt.lfmm"),
  method = 'mode',
  K = K,
  run = best_run
)

dat.imp = read.lfmm(paste0(work.directory, file_name, "_filt.lfmm_imputed.lfmm"))
table(dat.imp)

lfmm2geno(
  paste0(work.directory, file_name, "_filt.lfmm_imputed.lfmm"),
  paste0(work.directory, file_name, "_filt.lfmm_imputed.geno")
)
message("...Done.")
