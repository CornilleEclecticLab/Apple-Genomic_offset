#!/usr/bin/env bash

Project_name="All_data_wo_orientalis"
LFMM_FILE="../../00-File_preparation/out/s01.lfmm_file/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered.lfmm"
POS_FILE="./in/data_tree_wo_orientalis.txt"
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

module load all gencore/2  
module load r-base/4.2.3

Rscript bin/snmf_structure.R $LFMM_FILE $POS_FILE $OUT_FOLD $OUT_FIG
EOF

sbatch "$SCRIPT"

