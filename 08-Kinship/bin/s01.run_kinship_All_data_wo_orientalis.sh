#!/bin/bash
#SBATCH -J s01.make_kinship_analysis
#SBATCH -o ./out/All_data_wo_orientalis/s00.sbatch/s01.make_kinship_analysis.out
#SBATCH -e ./out/All_data_wo_orientalis/s00.sbatch/s01.make_kinship_analysis.err

module load all gencore/3
module load plink2

plink2 \
  --vcf /scratch/mc10739/01-SNPs_calling_filtering/SNPs_filtering/Crops-wild_dataset/out/gatk/s08.pop_gen_filter/chr_concat.geno20_maf005.anno.vcf.gz \
  --make-king-table \
  --king-cutoff 0.354 \
  --make-bed \
  --out ./out/All_data_wo_orientalis/s01.kin_table/chr_concat.geno20_maf005.anno
