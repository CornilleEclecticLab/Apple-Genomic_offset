#!/usr/bin/env bash

#SBATCH -J s01.plink_pca
#SBATCH -o out/s00.sbatch/s01.plink_pca.out
#SBATCH -e out/s00.sbatch/s01.plink_pca.err
#SBATCH -c 4
#SBATCH --mem=50G

module load all gencore/3
module load plink/1.9

PROJECT="dom_admixed_rm"
OUTPUT_DIR="./out/${PROJECT}/s01.plink_pca"
INPUT_FILE="../00-File_preparation/out/${PROJECT}/s01.filter_vcf/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered_${PROJECT}"


mkdir -p out/s00.sbatch
mkdir -p ${OUTPUT_DIR}

plink --bfile $INPUT_FILE --pca --chr-set 17 no-xy --out $OUTPUT_DIR/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered_$PROJECT
