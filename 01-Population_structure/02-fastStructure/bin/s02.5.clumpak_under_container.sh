#!/bin/bash

module load all gencore/3
module load apptainer

#apptainer build clumpak.sif docker://aurelim/clumpak:latest

PROJECT=dom_admixed_rm
# PROJECT=All_data_wo_orientalis

faststructure_dir=/out/${PROJECT}

apptainer run \
  --bind /scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/01-Population_structure/02-fastStructure:/mnt \
  clumpak.sif \
  --id /mnt/${PROJECT}_run \
  --dir /mnt/out/${PROJECT}/s02.clumpak/ \
  --file /mnt/$faststructure_dir/input_clumpak.zip \
  --inputtype admixture

  