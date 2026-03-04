#!/bin/bash
#SBATCH -J s01.make_kinship_analysis
#SBATCH -o ./out/dom_admixed_rm/s00.sbatch/s01.make_kinship_analysis.out
#SBATCH -e ./out/dom_admixed_rm/s00.sbatch/s01.make_kinship_analysis.err

module load all gencore/3
module load plink2

plink2 \
  --vcf /scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/01-Population_structure/03-membership_coefficient_investigating/out/All_data_wo_orientalis/s04.filter_vcf/chr_concat.geno20_maf005.anno_dom_admixed_rm.vcf.gz \
  --make-king-table \
  --king-cutoff 0.354 \
  --make-bed \
  --out ./out/dom_admixed_rm/s01.kin_table/chr_concat.geno20_maf005.anno_dom_admixed_rm
