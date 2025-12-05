#!/usr/bin/env bash
#SBATCH --job-name=s06.Plink_distance
#SBATCH --output=./out/s00.sbatch/s06.Plink_distance.out
#SBATCH --error=./out/s00.sbatch/s06.Plink_distance.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=50G

module load all gencore/3
module load plink/1.9

PROJECT="tips_dom_admixed_rm"
OUT_FOLD="./out/$PROJECT/s06.plink_distance"

mkdir -p $OUT_FOLD

plink --bfile ./out/$PROJECT/s01.filter_vcf/GraffiTE.merged.genotypes_$PROJECT \
      --pca --distance square --chr-set 17 no-xy \
      --out $OUT_FOLD/GraffiTE.merged.genotypes_$PROJECT
