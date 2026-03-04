#!/usr/bin/env bash
#SBATCH --job-name=s01.Comp_chr.geno20_maf005.anno_all_wild.LD_prunned
#SBATCH --output=./out/all_wild/s00.sbatch/s01.Comp_chr.geno20_maf005.anno_all_wild.LDpruned.out
#SBATCH --error=./out/all_wild/s00.sbatch/s01.Comp_chr.geno20_maf005.anno_all_wild.LDpruned.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=150G

module load miniconda

# Activate conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-lea-env

module load all gencore/3 
module load samtools

# Run R script
# Rscript bin/s01.computation_by_chr_v2.R ../00-File_preparation/out/dom_admixed_rm/s03.LD_prunning/chr_concat.geno20_maf005.anno_dom_admixed_rm.LDpruned.vcf.gz dom_admixed_rm 5 ../00-File_preparation/out/dom_admixed_rm/s03.LD_prunning/chr_concat.geno20_maf005.anno_dom_admixed_rm.LDpruned.lfmm
Rscript bin/s01.computation_by_chr_v2.R ../00-File_preparation/out/all_wild/s03.LD_prunning/chr_concat.geno20_maf005.anno_all_wild.LDpruned.vcf.gz all_wild 6 ../00-File_preparation/out/all_wild/s03.LD_prunning/chr_concat.geno20_maf005.anno_all_wild.LDpruned.lfmm ./in/data_tree_all_wild.txt