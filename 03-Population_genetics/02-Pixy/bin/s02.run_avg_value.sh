#!/bin/bash

#SBATCH -J avg_value.sh
#SBATCH --output ./out/test/s00.sbatch/avg_value.out
#SBATCH --error ./out/test/s00.sbatch/avg_value.err


#activate conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate avg_value

python3 ./bin/avg_value.py
