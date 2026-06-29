#!/usr/bin/env bash
#SBATCH --job-name=s02.lfmm2
#SBATCH --output=./out/s00.sbatch/s02.lfmm2_dom_admixed_rm_v3.out
#SBATCH --error=./out/s00.sbatch/s02.lfmm2_dom_admixed_rm_v3.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=100G

module load miniconda

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-lea-env

PROJECT="dom_admixed_rm"
env_file="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/dom_admixed_rm/climate_extract_current_dom_admixed_rm.csv"
position_data="./in/data_tree_dom_admixed_rm_V2.txt"
K_value=5
geno_LDpruned="./out/dom_admixed_rm/s01.computation/chr_concat.geno20_maf005.anno_dom_admixed_rm.LDpruned/file/chr_concat.geno20_maf005.anno_dom_admixed_rm.LDpruned_filt.lfmm_imputed.geno"

# PROJECT="all_wild"
# env_file="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/all_wild/climate_extract_current_all_wild.csv"
# position_data="./in/data_tree_all_wild.txt"
# K_value=6
# geno_LDpruned="./out/all_wild/s01.computation/chr_concat.geno20_maf005.anno_all_wild.LDpruned/file/chr_concat.geno20_maf005.anno_all_wild.LDpruned_filt.lfmm_imputed.geno"


# Rscript bin/s02.lfmm2_by_chr_v2.R $PROJECT $env_file $position_data $K_value $geno_LDpruned

Rscript bin/s02.lfmm2_by_chr_v3.R $PROJECT $env_file $position_data $K_value $geno_LDpruned