import os

# ==== USER PARAMETERS ====
# PROJECT = "dom_admixed_rm"
PROJECT = "all_wild"
VCF_DIR = f"../00-File_preparation/out/{PROJECT}/s01.split_vcf_by_chr"
LFMM_DIR = f"../00-File_preparation/out/{PROJECT}/s02.lfmm_file"
OUTPUT_DIR = f"./bin/s01.computation_by_chr_{PROJECT}"
SBATCH_DIR = f"./out/{PROJECT}/s00.sbatch"
R_SCRIPT = "bin/s01.computation_by_chr_v2.R"
# K_value = 5
K_value = 6
BIN_DIR = "./bin"
position_data = "./in/data_tree_all_wild.txt"

# ==== CREATE DIRECTORIES ====
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(SBATCH_DIR, exist_ok=True)

# ==== LOOP OVER VCF FILES ====
for vcf_file in sorted(os.listdir(VCF_DIR)):
    if vcf_file.endswith(".vcf.gz"):
        chr_name = os.path.splitext(os.path.splitext(vcf_file)[0])[0]
        bash_file = os.path.join(OUTPUT_DIR, f"s01.computation_{chr_name}.sh")

        lfmm_file = os.path.join(LFMM_DIR, f"{chr_name}.lfmm")

        with open(bash_file, "w") as f:
            f.write(f"""#!/usr/bin/env bash
#SBATCH --job-name=s01.Comp_{chr_name}
#SBATCH --output={SBATCH_DIR}/s01.Comp_{chr_name}.out
#SBATCH --error={SBATCH_DIR}/s01.Comp_{chr_name}.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH -c 1
#SBATCH --mem=150G

module load miniconda

# Activate conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate r-lea-env

module load all gencore/3
module load samtools

# Run R script
Rscript {R_SCRIPT} {os.path.join(VCF_DIR, vcf_file)} {PROJECT} {K_value} {lfmm_file} {position_data}
""")

print("✅ All computation scripts created successfully.")

command_file = os.path.join(BIN_DIR, f"s01.computation_by_chr_{PROJECT}.command_to_run.txt")
with open(command_file, "w") as f:
    subdir = f"s01.computation_by_chr_{PROJECT}"
    for file in os.listdir(os.path.join(BIN_DIR, subdir)):
        if file.startswith("s01") and file.endswith(".sh"):
            path = os.path.join(BIN_DIR, subdir, file)
            f.write(f"sbatch {path}\n")
