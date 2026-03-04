#!/usr/bin/env bash
#SBATCH --job-name=s04.2.GO_spectrum
#SBATCH --output=./out/s00.sbatch/s04.2.GO_spectrum.out
#SBATCH --error=./out/s00.sbatch/s04.2.GO_spectrum.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=80G

module load miniconda

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-lea-env

Rscript ./bin/s04.2.GO_spectrum.R