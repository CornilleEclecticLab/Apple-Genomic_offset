#!/bin/bash

#SBATCH --job-name=s01.filter_vcf
#SBATCH --output=./out/s00.sbatch/s01.filter_vcf.out
#SBATCH --error=./out/s00.sbatch/s01.filter_vcf.err
#SBATCH --time=15:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1

module load all gencore/3
module load bcftools
module load plink/1.9
module load samtools

# Inputs
# Project name for output files
#PROJECT="tips_dom_admixed_rm_V2"
PROJECT="tips_dom_admixed_rm_V3"
#VCF_LIST="./in/vcf_list_${PROJECT}.txt" # File listing VCFs (one per line)
VCF_LIST="./in/vcf_list_tips_dom_admixed_rm.txt"
#REMOVE_LIST="./in/${PROJECT}_list.txt" # A file with one sample ID per line to remove
REMOVE_LIST="./in/tips_dom_admixed_rm_list_V2.txt"
OUTPUT_DIR="./out/${PROJECT}/s01.filter_vcf" # Directory for outputs

mkdir -p "$OUTPUT_DIR"

while read -r VCF; do
    if [[ ! -f "$VCF" ]]; then
        echo "Warning: File not found - $VCF"
        continue
    fi

    BASENAME=$(basename "$VCF" .vcf.gz)
    FILTERED_VCF="${OUTPUT_DIR}/${BASENAME}_${PROJECT}.vcf.gz"
    PLINK_PREFIX="${OUTPUT_DIR}/${BASENAME}_${PROJECT}"

    echo "Processing $VCF..."

    # Remove individuals
    # bcftools view --samples-file ^"$REMOVE_LIST" "$VCF" --force-samples -Oz -o "$FILTERED_VCF"
    # tabix -f -p vcf "$FILTERED_VCF"
    bcftools view --samples-file ^"$REMOVE_LIST" "$VCF" --force-samples \
        | bcftools filter -e 'AC==0 || AC==AN' -Oz -o "$FILTERED_VCF"
    tabix -f -p vcf "$FILTERED_VCF"
    
    # Optional: skip if output is empty
    if [[ $(bcftools view -H "$FILTERED_VCF" | wc -l) -eq 0 ]]; then
        echo "Warning: No variants left in $FILTERED_VCF after filtering. Skipping..."
        continue
    fi

    # Convert to PLINK .bed format
    plink --vcf "$FILTERED_VCF" \
        --make-bed \
        --const-fid \
        --out "$PLINK_PREFIX" \
        --set-missing-var-ids @:\# \
        --keep-allele-order \
        --allow-extra-chr

done < "$VCF_LIST"

echo "All VCFs processed."
