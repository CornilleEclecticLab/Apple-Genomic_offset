#!/usr/bin/env bash
#SBATCH --job-name=s05.GO
#SBATCH --output=./out/s00.sbatch/s05.GO.out
#SBATCH --error=./out/s00.sbatch/s05.GO.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=50G

module load miniconda

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-lea-env

Rscript bin/s05.GO.R