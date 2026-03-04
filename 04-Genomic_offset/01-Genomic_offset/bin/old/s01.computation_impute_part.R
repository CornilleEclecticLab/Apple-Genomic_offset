#!/usr/bin/env Rscript
rm(list = ls())
source("./bin/fonction.R")
source("./bin/config.R")

create_directory(computational_folder)
work.directory = paste0(computational_folder, "file/")
fig_output = paste0(computational_folder, "figure/")

file_name <- basename(sub("\\.lfmm$", "", lfmm))

project.missing = load.snmfProject("./out/dom_admixed_rm/s01.Computation/file/chr_concat.geno20_maf005.anno_filtered_dom_admixed_rm_filt.snmfProject")

best_run = which.min(cross.entropy(project.missing, K = K)) 
# from the out file of the full computation script, but use here to be faster

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