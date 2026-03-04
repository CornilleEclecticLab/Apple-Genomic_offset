#!/usr/bin/env bash

Project="dom_admixed_rm"
mkdir -p ./out/s00.sbatch
mkdir -p ./out/$Project/s02.lfmm_file

sbatch bin/s02.lfmm_file_creation.sh /scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/04-Genomic_offset/00-File_preparation/out/$Project/s01.split_vcf_by_chr/ *.vcf.gz ./out/$Project/s02.lfmm_file


## to change need to do it for 17 files (chromosomes)