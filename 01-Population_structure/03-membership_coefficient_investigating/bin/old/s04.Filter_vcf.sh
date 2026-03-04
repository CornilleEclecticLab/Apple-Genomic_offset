#!/bin/bash

#SBATCH --job-name=s02.filter_vcf
#SBATCH --output=./out/s00.sbatch/s02.filter_vcf.out
#SBATCH --error=./out/s00.sbatch/s02.filter_vcf.err
#SBATCH --time=4:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1

module load all gencore/3
module load bcftools
module load plink/1.9
module load samtools

# Inputs
VCF_LIST="./in/vcf_file_to_filtered_admix.txt" # File listing VCFs (one per line)
REMOVE_LIST="./out/fastStructure_admix_list_K5_0.2.txt" # A file with one sample ID per line to remove
Threshold="0_2" # Threshold for naming used before
OUTPUT_DIR="./out/s02.filter_vcf" # Directory for outputs

mkdir -p "$OUTPUT_DIR"

while read -r VCF; do
    if [[ ! -f "$VCF" ]]; then
        echo "Warning: File not found - $VCF"
        continue
    fi

    BASENAME=$(basename "$VCF" .vcf.gz)
    FILTERED_VCF="${OUTPUT_DIR}/${BASENAME}_admix_dom_rm_${Threshold}.vcf.gz"
    PLINK_PREFIX="${OUTPUT_DIR}/${BASENAME}_admix_dom_rm_${Threshold}"

    echo "Processing $VCF..."

    # Remove individuals
    bcftools view --samples-file ^$REMOVE_LIST -Oz -o "$FILTERED_VCF" "$VCF"
    tabix -p vcf "$FILTERED_VCF"

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
