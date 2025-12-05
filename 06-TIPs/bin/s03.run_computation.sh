#!/usr/bin/env bash
#SBATCH --job-name=s03.computation
#SBATCH --output=./out/s00.sbatch/s03.computation.out
#SBATCH --error=./out/s00.sbatch/s03.computation.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=100G

module load miniconda

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-lea-env

Rscript bin/s03.computation.R