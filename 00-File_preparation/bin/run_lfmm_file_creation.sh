#!/usr/bin/env bash

mkdir ./out/s00.sbatch
mkdir ./out/s01.lfmm_file

sbatch bin/s01_lfmm_file_creation.sh /scratch/mc10739/01-SNPs_calling_filtering/SNPs_filtering/out/gatk/s08.pop_gen_filter/ *.vcf.gz ./out/s01.lfmm_file

