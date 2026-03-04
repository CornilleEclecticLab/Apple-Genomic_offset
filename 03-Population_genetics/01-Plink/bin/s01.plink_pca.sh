#!/usr/bin/env bash

#SBATCH -J s01.plink_pca
#SBATCH -o out/s00.sbatch/s01.plink_pca.out
#SBATCH -e out/s00.sbatch/s01.plink_pca.err
#SBATCH -c 4
#SBATCH --mem=50G

module load all gencore/3
module load plink2

PROJECT="dom_admixed_rm"
# PROJECT="All_data_wo_orientalis"

OUTPUT_DIR="./out/${PROJECT}/s01.plink_pca"

INPUT_FILE="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/01-Population_structure/03-membership_coefficient_investigating/out/All_data_wo_orientalis/s04.filter_vcf/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm"
# INPUT_FILE="/scratch/mc10739/01-SNPs_calling_filtering/SNPs_filtering/Crops-wild_dataset/out/gatk/s08.pop_gen_filter/chr_concat.geno20_maf005.anno.vcf.gz"

mkdir -p out/s00.sbatch
mkdir -p ${OUTPUT_DIR}

# plink2 --vcf $INPUT_FILE --make-bed --pca --chr-set 17 no-xy --out $OUTPUT_DIR/chr_concat.geno20_maf005.anno.syno.LDpruned_$PROJECT

plink2 --bfile $INPUT_FILE --pca --chr-set 17 no-xy --out $OUTPUT_DIR/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm
