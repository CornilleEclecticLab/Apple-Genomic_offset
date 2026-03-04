#!/bin/bash

# num_clumpak_folder=1762167588 #Depending of clumpak results
Project_name="dom_admixed_rm"
# Project_name="All_data_wo_orientalis"
# ZIP_file="./out/"$Project_name"/s02.clumpak/"$num_clumpak_folder".zip"

# unzip -q $ZIP_file -d ./out/$Project_name/s02-5.clumpak_2/

# mv ./out/$Project_name/s02-5.clumpak_2/$num_clumpak_folder/* \
#    ./out/$Project_name/s02-5.clumpak_2/

# rm -rf ./out/$Project_name/s02-5.clumpak_2/$num_clumpak_folder/

python3 ./bin/get_clumpak_result_to_pophelper_V2.py ./out/$Project_name/s02.clumpak/ 2 10  ./out/$Project_name/s03.pophelper/
