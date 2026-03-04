#!/usr/bin/env bash

# Project_name="All_data_wo_orientalis"
# LFMM_FILE="../../00-File_preparation/out/$Project_name/s01.lfmm_file/chr_concat.geno20_maf005.anno.syno.LDpruned.lfmm"
# POS_FILE="./in/data_tree_wo_orientalis_V2.txt"

# Project_name="dom_admixed_rm"
# LFMM_FILE="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/00-File_preparation/out/dom_admixed_rm/s01.lfmm_file/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.lfmm"
# POS_FILE="./in/data_tree_dom_admixed_rm_V2.txt"

Project_name="all_wild"
LFMM_FILE="/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/00-File_preparation/out/all_wild/s01.lfmm_file/chr_concat.geno20_maf005.anno.syno.LDpruned_all_wild.lfmm"
POS_FILE="./in/data_tree_all_wild.txt"

OUT_FOLD="./out/$Project_name/s01.snmf"
OUT_FIG="./out/$Project_name/s01.snmf/fig"
SCRIPT="./bin/snmf_${Project_name}.sh"

mkdir -p out/$Project_name/s00.sbatch
mkdir -p "$OUT_FOLD"
mkdir -p "$OUT_FIG"


cat <<EOF > "$SCRIPT"
#!/usr/bin/env bash

#SBATCH --job-name=${Project_name}_snmf_structure
#SBATCH --output=./out/$Project_name/s00.sbatch/snmf_structure.out
#SBATCH --error=./out/$Project_name/s00.sbatch/snmf_structure.err
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-snmf-env

Rscript bin/snmf_structure.R $LFMM_FILE $POS_FILE $OUT_FOLD $OUT_FIG
EOF

sbatch "$SCRIPT"

