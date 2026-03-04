#!/usr/bin/env bash
#SBATCH --job-name=s03.GO
#SBATCH --output=./out/s00.sbatch/s03.GO.out
#SBATCH --error=./out/s00.sbatch/s03.GO.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=150G

module load miniconda

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-lea-env
PROJECT="dom_admixed_rm"
env_file="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/dom_admixed_rm/Bioclim_dom_admixed_rm_chelsa_current.csv"
new_env_file="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/dom_admixed_rm/env_Saclay_Chelsa"
position_data="./in/data_tree_dom_admixed_rm_V2.txt"
K_value=5

Rscript bin/s03.GO_v2.R $PROJECT $env_file $new_env_file $position_data $K_value