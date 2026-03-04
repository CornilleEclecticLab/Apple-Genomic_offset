#!/bin/bash

#SBATCH -J s04.run_plot_structure.sh
#SBATCH --output ./out/All_data_wo_admixed/s00.sbatch/s04.run_plot_structure.out
#SBATCH --error ./out/All_data_wo_admixed/s00.sbatch/s04.run_plot_structure.err

module load all gencore/2
module load r-base/4.2.3

Project_name="dom_admixed_rm"
INFO_FILE="./in/data_tree_dom_admixed_rm_V2.txt"

# Project_name="All_data_wo_orientalis"
# INFO_FILE="./in/data_tree_wo_orientalis_V2.txt"

OUT_FOLDER=./out/$Project_name
Q_FILES_DIR=$OUT_FOLDER/s03.pophelper


Rscript ./bin/plot_population_structure_V2.R $Project_name $OUT_FOLDER $Q_FILES_DIR $INFO_FILE
