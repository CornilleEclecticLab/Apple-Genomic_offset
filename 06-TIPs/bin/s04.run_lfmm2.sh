#!/usr/bin/env bash
#SBATCH --job-name=s04.lfmm2
#SBATCH --output=./out/s00.sbatch/s04.lfmm2.out
#SBATCH --error=./out/s00.sbatch/s04.lfmm2.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=50G

module load miniconda

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-lea-env

Rscript bin/s04.lfmm2.R