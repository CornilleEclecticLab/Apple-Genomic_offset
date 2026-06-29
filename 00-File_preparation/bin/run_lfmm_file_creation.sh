#!/usr/bin/env bash

# Project="dom_admixed_rm"
# Project="All_data_wo_orientalis"
Project="all_wild"
mkdir -p ./out/s00.sbatch
mkdir -p ./out/$Project/s01.lfmm_file

# sbatch bin/s01.lfmm_file_creation.sh   *.vcf.gz ./out/$Project/s01.lfmm_file
# sbatch bin/s01.lfmm_file_creation.sh /scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/01-Population_structure/03-membership_coefficient_investigating/out/All_data_wo_orientalis/s04.filter_vcf/ *.vcf.gz ./out/$Project/s01.lfmm_file
sbatch bin/s01.lfmm_file_creation.sh /scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/01-Population_structure/03-membership_coefficient_investigating/out/All_data_wo_orientalis/s04.5.filter_crops_ind/ *.vcf.gz ./out/$Project/s01.lfmm_file