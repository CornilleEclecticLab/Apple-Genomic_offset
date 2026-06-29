import datetime
import sys
import textwrap
import os

# Configuration

script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
Project = "invariant_variant_filter_dom_admixed_rm_54K"
output_dir = os.path.join(work_dir, 'out', Project, "s05.merged_invariant_variant_by_chr")
sub_script_dir = os.path.join(bin_dir, f"{script_basename}_{Project}")
sbatch_out_dir = "out/s00.sbatch"

variant_in_dir = "out/invariant_variant_filter_dom_admixed_rm_54K/s04.split_vcf_variant_filtered_by_chr"
invariant_in_dir = "out/invariant_variant_filter_dom_admixed_rm/s03.filter_invariant"

# Create necessary directories
os.makedirs(sub_script_dir, exist_ok=True)
os.makedirs(sbatch_out_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# Template for the bash script
script_template = """#!/usr/bin/env bash

#SBATCH -J s05.merge_chr{chr_num}
#SBATCH -o {sbatch_out}/s05.merge_chr{chr_num}.out
#SBATCH -e {sbatch_out}/s05.merge_chr{chr_num}.err
#SBATCH -c 4
#SBATCH --mem=8G

module load all gencore/3
module load bcftools

# Define paths
INV_FILE="{inv_dir}/chr{chr_num}.filtered_invariant_invariant_variant_filter_dom_admixed_rm.vcf.gz"
# VAR_FILE="{var_dir}/chr_chr{chr_pad}.chr_concat.geno20_maf005.anno_dom_admixed_rm.vcf.gz"
VAR_FILE="{var_dir}/chr_chr{chr_pad}.chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.vcf.gz"
MERGED_FILE="{out_dir}/chr{chr_num}.merged_variant_invariant.vcf.gz"

echo "Processing Chromosome {chr_num}..."

# Concatenate variant and invariant sites
bcftools concat -a "$INV_FILE" "$VAR_FILE" -O z -o "$MERGED_FILE"

# Index the output
bcftools index -t "$MERGED_FILE"

echo "Chromosome {chr_num} complete."
"""

for i in range(1, 18):
    chr_num = str(i)
    chr_pad = f"{i:02d}"  # Pads 1 to 01
    
    # Fill the template
    content = script_template.format(
        chr_num=chr_num,
        chr_pad=chr_pad,
        sbatch_out=sbatch_out_dir,
        inv_dir=invariant_in_dir,
        var_dir=variant_in_dir,
        out_dir=output_dir
    )
    
    # Write the file
    file_path = os.path.join(sub_script_dir, f"merge_chr{chr_pad}.sh")
    with open(file_path, "w") as f:
        f.write(content)

# List all sub-scripts
files = [f for f in os.listdir(sub_script_dir) if os.path.isfile(os.path.join(sub_script_dir, f))]
output_file = os.path.join("./bin", f"s05.merged_invariant_variant_by_chr_{Project}.command_to_run.txt")

with open(output_file, "w") as f:
    for file in files:
        f.write(f"sbatch {os.path.join(sub_script_dir, file)}\n")


print(f"Success: 17 scripts generated in '{sub_script_dir}'")