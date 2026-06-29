#!/bin/bash

#SBATCH --job-name=s02.run_vcf_to_lfmm
#SBATCH --output=./out/s00.sbatch/s02.run_vcf_to_lfmm.out
#SBATCH --error=./out/s00.sbatch/s02.run_vcf_to_lfmm.err
#SBATCH --time=5:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1


module load all gencore/3
module load python

PROJECT="tips_dom_admixed_rm_V3"
OUTPUT_DIR="./out/${PROJECT}/s02.vcf_to_lfmm" # Directory for outputs
vcf="./out/${PROJECT}/s01.filter_vcf/GraffiTE.merged.genotypes_${PROJECT}.vcf.gz"

mkdir -p $OUTPUT_DIR


python3 bin/s02.vcf_to_lfmm.py $vcf "${OUTPUT_DIR}/GraffiTE.merged.genotypes_${PROJECT}"