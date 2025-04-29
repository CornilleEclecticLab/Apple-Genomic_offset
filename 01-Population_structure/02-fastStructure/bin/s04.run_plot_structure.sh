#!/bin/bash

module load all gencore/2
module load r-base/4.2.3

Project_name="All_data_wo_orientalis"
OUT_FOLDER=./out/$Project_name
Q_FILES_DIR=$OUT_FOLDER/s03.pophelper
INFO_FILE="./in/data_tree_wo_orientalis.txt"

Rscript ./bin/plot_population_structure.R $Project_name $OUT_FOLDER $Q_FILES_DIR $INFO_FILE