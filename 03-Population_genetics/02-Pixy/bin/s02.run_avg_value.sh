#!/bin/bash

#SBATCH -J avg_value.sh
#SBATCH --output ./out/full_genome_dom_admixed_rm_no_filter_pi_by_ind/s00.sbatch/avg_value.out
#SBATCH --error ./out/full_genome_dom_admixed_rm_no_filter_pi_by_ind/s00.sbatch/avg_value.err

module load miniconda
#activate conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate avg_value

python3 ./bin/avg_value.py
