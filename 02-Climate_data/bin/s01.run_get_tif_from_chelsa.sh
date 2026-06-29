#!/bin/bash

#SBATCH -J extract_climate_data.sh
#SBATCH --output ./out/s00.sbatch/s01.get_tif_from_chelsa_LGM.out
#SBATCH --error ./out/s00.sbatch/s01.get_tif_from_chelsa_LGM.err
#SBATCH --time=48:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1

module load miniconda

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-chelsa-env

Rscript bin/s01.get_tif_from_chelsa.R