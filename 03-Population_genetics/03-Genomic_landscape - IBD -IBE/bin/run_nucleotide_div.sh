#!/usr/bin/env bash
#SBATCH --job-name=Nucleotide_diversity
#SBATCH --output=./out/s00.sbatch/Nucleotide_diversity.out
#SBATCH --error=./out/s00.sbatch/Nucleotide_diversity.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=150G

module load miniconda
# Active conda env temporairement (pas nécessaire si déjà dans le bon env)
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-spatial-env

which Rscript

# Définit R_HOME **temporairement juste pour ce script**
#export R_HOME=$(R RHOME)

# Exécute le script R normalement
Rscript ./bin/Nucleotide_diversity.R 