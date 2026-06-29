#!/bin/sh

#SBATCH -J s01_bis.pixy_for_small_dataset
#SBATCH -o job.%J.out
#SBATCH -e job.%J.er
#SBATCH -c 10
#SBATCH --mem=20G

# Set up environment
module load all gencore/3
module load pixy
module load bcftools

mkdir -p ./out/dom_admixed_rm_54K/s01_bis.pixy_for_small_dataset_pi_by_ind/

#--stats pi fst dxy
# pixy --stats fst dxy \
#     --vcf /scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/01-Population_structure/03-membership_coefficient_investigating/out/All_data_wo_orientalis/s04.filter_vcf/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.vcf.gz \
#     --populations ./in/group_geo_dom_admixed_rm.txt \
#     --window_size 60000 \
#     --n_cores 10 \
#     --output_folder /scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/03-Population_genetics/02-Pixy/out/dom_admixed_rm_54K/s01_bis.pixy_for_small_dataset/ \
#     --output_prefix s01_bis.pixy_for_small_dataset \
#     --bypass_invariant_check

pixy --stats pi \
    --vcf /scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/01-Population_structure/03-membership_coefficient_investigating/out/All_data_wo_orientalis/s04.filter_vcf/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.vcf.gz \
    --populations ./in/group_geo_invariant_variant_filter_dom_admixed_rm_pi_by_ind.txt \
    --window_size 60000 \
    --n_cores 10 \
    --output_folder /scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/03-Population_genetics/02-Pixy/out/dom_admixed_rm_54K/s01_bis.pixy_for_small_dataset_pi_by_ind/ \
    --output_prefix s01_bis.pixy_for_small_dataset \
    --bypass_invariant_check
