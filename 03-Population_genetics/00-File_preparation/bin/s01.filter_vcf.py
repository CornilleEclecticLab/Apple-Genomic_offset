import os

# ==== User parameters ====
VCF_LIST = "./in/vcf_list_full_genome_dom_admixed_rm.txt"
REMOVE_LIST = "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/01-Population_structure/03-membership_coefficient_investigating/bin/final_remove_list.txt"
PROJECT = "full_genome_dom_admixed_rm_no_filter"
OUTPUT_DIR = f"./out/{PROJECT}/s01.filter_vcf"
SBATCH_DIR = f"./bin/s01.filter_vcf_{PROJECT}"
SBATCH_OUTPUT_DIR = f"./out/s00.sbatch"
BIN_FOLDER = "./bin"

# ==== Create output directories ====
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(SBATCH_DIR, exist_ok=True)
os.makedirs(SBATCH_OUTPUT_DIR, exist_ok=True)

# ==== Template for SLURM script ====
SCRIPT_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={sbatch_output}/{basename}.out
#SBATCH --error={sbatch_output}/{basename}.err
#SBATCH --time=20:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1

module load all gencore/3
module load bcftools
module load plink/1.9
module load samtools

VCF="{vcf}"
REMOVE_LIST="{remove_list}"
PROJECT="{project}"
OUTPUT_DIR="{output_dir}"

BASENAME=$(basename "$VCF" .vcf.gz)
FILTERED_VCF="${{OUTPUT_DIR}}/${{BASENAME}}_${{PROJECT}}.vcf.gz"
PLINK_PREFIX="${{OUTPUT_DIR}}/${{BASENAME}}_${{PROJECT}}"

echo "Processing $VCF..."

# Remove individuals
bcftools view --samples-file ^"$REMOVE_LIST" "$VCF" | \\
# bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.05' | \\ #remove for the variant/invariant no filtering test
 bcftools filter --exclude 'AC==AN || AC==0' -Oz -o "$FILTERED_VCF"
tabix -f -p vcf "$FILTERED_VCF"

# Optional: skip if output is empty
if [[ $(bcftools view -H "$FILTERED_VCF" | wc -l) -eq 0 ]]; then
    echo "Warning: No variants left in $FILTERED_VCF after filtering. Skipping..."
    exit 0
fi

# Convert to PLINK .bed format
plink --vcf "$FILTERED_VCF" \\
    --make-bed \\
    --const-fid \\
    --out "$PLINK_PREFIX" \\
    --set-missing-var-ids @:\\# \\
    --keep-allele-order \\
    --allow-extra-chr

echo "Finished $VCF."
"""

# ==== Generate one script per VCF ====
with open(VCF_LIST) as f:
    vcfs = [line.strip() for line in f if line.strip()]

for vcf in vcfs:
    if not os.path.exists(vcf):
        print(f"⚠️ Warning: file not found -> {vcf}")
        continue

    basename = os.path.basename(vcf).replace(".vcf.gz", "")
    job_name = f"s01.filter_{basename}"
    script_path = os.path.join(SBATCH_DIR, f"{job_name}.sh")

    with open(script_path, "w") as out:
        out.write(SCRIPT_TEMPLATE.format(
            job_name=job_name,
            sbatch_output=SBATCH_OUTPUT_DIR,
            basename=basename,
            vcf=vcf,
            remove_list=REMOVE_LIST,
            project=PROJECT,
            output_dir=OUTPUT_DIR
        ))

    print(f"✅ Created: {script_path}")

print("\nAll SLURM scripts created.")

command_file = os.path.join(BIN_FOLDER, f"s01.filter_vcf_{PROJECT}.command_to_run.txt")
with open(command_file, "w") as f:
    subdir = f"s01.filter_vcf_{PROJECT}"
    for file in os.listdir(os.path.join(BIN_FOLDER, subdir)):
        if file.startswith("s01") and file.endswith(".sh"):
            path = os.path.join(BIN_FOLDER, subdir, file)
            f.write(f"sbatch {path}\n")


print(f"To submit all: for f in {SBATCH_DIR}/*.sh; do sbatch $f; done \n\
Or slurm_parallel_ja_submit.sh {BIN_FOLDER}/s01.filter_vcf_{PROJECT}.command_to_run.txt")