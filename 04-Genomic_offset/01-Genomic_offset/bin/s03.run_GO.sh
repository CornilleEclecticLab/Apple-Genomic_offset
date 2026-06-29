#!/usr/bin/env bash
#SBATCH --job-name=s03.GO_V3
#SBATCH --output=./out/s00.sbatch/s03.GO_all_wild_V3_new_env.out
#SBATCH --error=./out/s00.sbatch/s03.GO_all_wild_V3_new_env.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=150G

module load miniconda

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-lea-env

# PROJECT="dom_admixed_rm"
# env_file="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/dom_admixed_rm/climate_extract_current_dom_admixed_rm.csv"
# new_env_file="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s03.orchard_generalization/dom_admixed_rm/env_Saclay_Chelsa"
# position_data="./in/data_tree_dom_admixed_rm_V2.txt"
# K_value=5

# PROJECT="all_wild"
# env_file="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/all_wild/climate_extract_current_all_wild.csv"
# new_env_file="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s03.orchard_generalization/all_wild/env_Saclay_Chelsa"
# position_data="./in/data_tree_all_wild.txt"
# K_value=6

# PROJECT="dom_admixed_rm"
# ENV_CURR="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/dom_admixed_rm/climate_extract_current_dom_admixed_rm.csv"
# # NEW_ENV_D="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s03.orchard_generalization/dom_admixed_rm" 
# NEW_ENV_D="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s05.make_bio_var/dom_admixed_rm"
# REF_INFO="./in/data_tree_dom_admixed_rm_V2.txt"
# K=5

PROJECT="all_wild"
ENV_CURR="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/all_wild/climate_extract_current_all_wild.csv"
# NEW_ENV_D="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s03.orchard_generalization/all_wild" 
NEW_ENV_D="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s05.make_bio_var/all_wild"
REF_INFO="./in/data_tree_all_wild.txt"
K=6


# Rscript bin/s03.GO_v2.R $PROJECT $env_file $new_env_file $position_data $K_value
Rscript bin/s03.GO_V3.R $PROJECT $ENV_CURR $NEW_ENV_D $REF_INFO $K