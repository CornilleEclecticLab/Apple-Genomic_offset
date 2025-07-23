#!/bin/bash

# Active conda env temporairement (pas nécessaire si déjà dans le bon env)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate env_IBD

which Rscript

# Définit R_HOME **temporairement juste pour ce script**
#export R_HOME=$(R RHOME)

# Exécute le script R normalement
Rscript ./bin/IBD.R 
