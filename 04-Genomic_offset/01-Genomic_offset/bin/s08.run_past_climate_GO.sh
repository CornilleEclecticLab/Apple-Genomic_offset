#!/usr/bin/env bash
#SBATCH --job-name=s08.Past_climate_GO
#SBATCH --output=./out/s00.sbatch/s08.Past_climate_GO_dom_admixed_rm.out
#SBATCH --error=./out/s00.sbatch/s08.Past_climate_GO_dom_admixed_rm.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=150G

module load miniconda

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-lea-env

PROJECT="dom_admixed_rm"
env_folder="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/dom_admixed_rm/"
ref_info="./in/data_tree_dom_admixed_rm_V2.txt"
K=5

Rscript bin/s08.Past_climate_GO.R $PROJECT $env_folder $ref_info $K