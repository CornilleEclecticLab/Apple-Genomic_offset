import os
import pandas as pd

# ==== User parameters ====
PROJECT = "All_data_wo_orientalis"
VCF_LIST_FILE = "./in/vcf_file_to_filtered_admix.txt"
OUT_FILE_EXTENTION = "dom_admixed_rm"
REMOVE_FOLD = f"./out/{PROJECT}/s01.filter_crops_wild"
OUTPUT_DIR = f"./out/{PROJECT}/s04.filter_vcf"
BIN_FOLDER = "./bin"
SBATCH_OUTPUT_DIR = f"./out/{PROJECT}/s00.sbatch"

# Create directories if they don't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(BIN_FOLDER,"s04.filter_vcf"), exist_ok=True)
os.makedirs(SBATCH_OUTPUT_DIR, exist_ok=True)

# ==== Generate final remove list ====
remove_files = [os.path.join(REMOVE_FOLD, f) for f in os.listdir(REMOVE_FOLD) if os.path.isfile(os.path.join(REMOVE_FOLD, f))]
all_samples = []
for f in remove_files:
    df = pd.read_csv(f, header=0)  # or sep="\t" if tab-delimited
    all_samples.extend(df.iloc[:, 0].astype(str).tolist())
final_remove_list = sorted(set(all_samples))
print(f"Total unique individuals to remove: {len(final_remove_list)}")

# Write the final remove list to a temporary file
REMOVE_LIST_FILE = os.path.join(BIN_FOLDER, "final_remove_list.txt")
with open(REMOVE_LIST_FILE, "w") as f:
    for sample in final_remove_list:
        f.write(sample + "\n")

# ==== Read VCF list ====
with open(VCF_LIST_FILE, "r") as f:
    vcf_files = [line.strip() for line in f if line.strip()]

# ==== Generate bash scripts ====
for vcf in vcf_files:
    if not os.path.isfile(vcf):
        print(f"Warning: File not found - {vcf}")
        continue

    basename = os.path.basename(vcf).replace(".vcf.gz", "")
    bash_file = os.path.join(BIN_FOLDER, "s04.filter_vcf", f"s04_filter_{basename}.sh")
    filtered_vcf = os.path.join(OUTPUT_DIR, f"{basename}_{OUT_FILE_EXTENTION}.vcf.gz")
    plink_prefix = os.path.join(OUTPUT_DIR, f"{basename}_{OUT_FILE_EXTENTION}")

    with open(bash_file, "w") as bf:
        bf.write(f"""#!/bin/bash

#SBATCH --job-name=s04.filter_vcf_{basename}
#SBATCH --output={SBATCH_OUTPUT_DIR}/s04.filter_vcf_{basename}.out
#SBATCH --error={SBATCH_OUTPUT_DIR}/s04.filter_vcf_{basename}.err
#SBATCH --time=10:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1

module load all gencore/3
module load bcftools
module load plink/1.9
module load samtools

# Remove individuals and filter
# Remove individuals and filter
bcftools view -S ^"{REMOVE_LIST_FILE}" "{vcf}" | \\
bcftools filter -e 'AC==0 || AC==AN' -Oz -o "{filtered_vcf}"

tabix -f -p vcf {filtered_vcf}

# Optional: skip if output is empty
if [[ $(bcftools view -H {filtered_vcf} | wc -l) -eq 0 ]]; then
    echo "Warning: No variants left in {filtered_vcf} after filtering. Skipping..."
else
    # Convert to PLINK .bed format
    plink --vcf {filtered_vcf} \\
          --make-bed \\
          --const-fid \\
          --out {plink_prefix} \\
          --set-missing-var-ids @:\\# \\
          --keep-allele-order \\
          --allow-extra-chr
fi

echo "Processing {vcf} finished."
""")
    
    print(f"Bash file created: {bash_file}")

# ==== Generate master command file ====
command_file = os.path.join(BIN_FOLDER, "s04.filter_vcf.command_to_run.txt")
with open(command_file, "w") as f:
    for file in os.listdir(os.path.join(BIN_FOLDER, "s04.filter_vcf")):
        if file.startswith("s04") and file.endswith(".sh"):
            f.write(f"sbatch {os.path.join(BIN_FOLDER, 's04.filter_vcf', file)}\n")
