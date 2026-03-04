#!/bin/bash

#SBATCH --job-name=s01.stacks
#SBATCH --output=./out/s00.sbatch/s01.stacks.out
#SBATCH --error=./out/s00.sbatch/s01.stacks.err
#SBATCH --time=5:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1


module load all gencore/3
module load stacks

PROJECT="dom_admixed_rm_stack_by_ind"
OUTPUT_DIR="./out/${PROJECT}/s01.stacks"
Pop_file="./in/group_geo_${PROJECT}.txt"
vcf_file="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/01-Population_structure/03-membership_coefficient_investigating/out/All_data_wo_orientalis/s04.filter_vcf/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.vcf.gz"

mkdir -p "$OUTPUT_DIR"

populations -V "${vcf_file}" -M "${Pop_file}" -O "${OUTPUT_DIR}" --fstats