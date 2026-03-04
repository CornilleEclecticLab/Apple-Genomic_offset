#!/bin/bash

#SBATCH -J s02.max_membership_coefficient.sh
#SBATCH --output ./out/dom_admixed_rm/s00.sbatch/s02.max_membership_coefficient.out
#SBATCH --error ./out/dom_admixed_rm/s00.sbatch/s02.max_membership_coefficient.err

module load all gencore/2
module load r-base/4.2.3

Rscript ./bin/s02.max_membership_coefficient.R