#!/bin/bash

#SBATCH -J s05.max_membership_coefficient.sh
#SBATCH --output ./out/dom_admixed_rm/s00.sbatch/s05.max_membership_coefficient.out
#SBATCH --error ./out/dom_admixed_rm/s00.sbatch/s05.max_membership_coefficient.err

module load all gencore/2
module load r-base/4.2.3

Rscript ./bin/s05.max_membership_coefficient.R