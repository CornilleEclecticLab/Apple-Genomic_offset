#!/usr/bin/env bash
#SBATCH --job-name=s01.Computation
#SBATCH --output=./out/s00.sbatch/s01.Computation.out
#SBATCH --error=./out/s00.sbatch/s01.Computation.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=150G

module load all gencore/3
module load miniconda
module load samtools

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-lea-env

#Rscript bin/s01.Computation.R
Rscript bin/s01.computation_by_chr.R