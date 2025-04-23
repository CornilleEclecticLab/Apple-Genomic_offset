#!/usr/bin/env bash
#SBATCH --job-name=snmf_structure
#SBATCH --output=./out/s00.sbatch/snmf_structure.out
#SBATCH --error=./out/s00.sbatch/snmf_structure.err
#SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4


module load all gencore/2  
module load r-base/4.2.3

LFMM_FILE="../../00-File_preparation/out/s01.lfmm_file/chr_concat.geno20_maf005.anno.syno.LDpruned.lfmm"
POS_FILE="in/info_position_ind.txt"
OUT_FOLD="out/s01.snmf"
OUT_FIG="out/s01.snmf/fig"

mkdir -p out/s00.sbatch
mkdir -p $OUT_FOLD
mkdir -p $OUT_FIG

Rscript bin/snmf_structure.R $LFMM_FILE $POS_FILE $OUT_FOLD $OUT_FIG

