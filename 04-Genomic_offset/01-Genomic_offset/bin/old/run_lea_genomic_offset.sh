#!/usr/bin/env bash
#SBATCH --job-name=lea_genomic_offset
#SBATCH --output=./out/s00.sbatch/lea_genomic_offset.out
#SBATCH --error=./out/s00.sbatch/lea_genomic_offset.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=100G

module load miniconda

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-lea-env

Rscript bin/lea_genomic_offset.R