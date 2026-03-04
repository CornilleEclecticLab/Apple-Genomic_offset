#!/bin/bash
#SBATCH --job-name=LD_prunning
#SBATCH --output=./out/s00.sbatch/LD_prunning.out
#SBATCH --error=./out/s00.sbatch/LD_prunning.err
#SBATCH --time=4:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4

# Load necessary modules
module load all gencore/3
module load plink/1.9
module load samtools
module load python
module load vcftools
# module load htslib  # Uncomment if bgzip/tabix are not in your path

# Variables
#PROJECT="dom_admixed_rm"
PROJECT="all_wild"
PREFIX="chr_concat.geno20_maf005.anno_all_wild"
IN_VCF="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/01-Population_structure/03-membership_coefficient_investigating/out/All_data_wo_orientalis/s04.5.filter_crops_ind/${PREFIX}.vcf.gz"
OUT_DIR="./out/${PROJECT}/s03.LD_prunning"
extension="*.vcf.gz"

# Ensure directories exist
mkdir -p "${OUT_DIR}"

## 1. Identify SNPs in LD
# This creates ${PREFIX}.prune.in (list of SNPs to keep)
plink --vcf "${IN_VCF}" \
    --make-bed \
    --indep-pairwise 50 5 0.2 \
    --const-fid \
    --out "${OUT_DIR}/${PREFIX}" \
    --set-missing-var-ids @:# \
    --keep-allele-order \
    --allow-extra-chr \
    --threads "${SLURM_CPUS_PER_TASK}"

## 2. Extract unlinked SNPs using your custom Python script
# Ensure the path to your python script is correct relative to where you submit the job
python3 ./bin/filter_by_po.py \
 "${IN_VCF}" \
 "${OUT_DIR}/${PREFIX}.prune.in" \
 "${OUT_DIR}/${PREFIX}.LDpruned.vcf"

## 3. Compress and Index
bgzip -f "${OUT_DIR}/${PREFIX}.LDpruned.vcf"
tabix -p vcf "${OUT_DIR}/${PREFIX}.LDpruned.vcf.gz"

## 4. LFMM creation
# Define a clean base name without extensions for output files
BASE_NAME="${PREFIX}.LDpruned"

echo "Processing ${BASE_NAME}..."

# 4.1 Convert VCF to 012 format
# We use --gzvcf because we know we just bgzipped it in Step 3
vcftools --gzvcf "${OUT_DIR}/${BASE_NAME}.vcf.gz" \
    --012 \
    --out "${OUT_DIR}/${BASE_NAME}"

# 4.2 Reformat for LFMM
# LFMM requires a matrix where -1 (missing) is converted to 9, 
# and the first column (index) is removed.
# VCFtools outputs the positions in a separate file ending in .012.pos
sed "s:-1:9:g" "${OUT_DIR}/${BASE_NAME}.012" | cut -f 2- > "${OUT_DIR}/${BASE_NAME}.lfmm"

# 4.3 Save coordinates/loci for reference
# VCFtools already creates this file as .012.pos
cp "${OUT_DIR}/${BASE_NAME}.012.pos" "${OUT_DIR}/coordinates_loci_${BASE_NAME}.txt"

echo "LFMM conversion complete: ${OUT_DIR}/${BASE_NAME}.lfmm"