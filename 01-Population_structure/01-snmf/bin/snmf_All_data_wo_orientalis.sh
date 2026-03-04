#!/usr/bin/env bash

#SBATCH --job-name=All_data_wo_orientalis_snmf_structure
#SBATCH --output=./out/All_data_wo_orientalis/s00.sbatch/snmf_structure.out
#SBATCH --error=./out/All_data_wo_orientalis/s00.sbatch/snmf_structure.err
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4

module load all gencore/2  
module load r-base/4.2.3

Rscript bin/snmf_structure.R ../../00-File_preparation/out/All_data_wo_orientalis/s01.lfmm_file/chr_concat.geno20_maf005.anno.syno.LDpruned.lfmm ./in/data_tree_wo_orientalis_V2.txt ./out/All_data_wo_orientalis/s01.snmf ./out/All_data_wo_orientalis/s01.snmf/fig
