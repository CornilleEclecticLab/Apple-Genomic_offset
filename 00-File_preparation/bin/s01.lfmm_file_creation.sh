#!/usr/bin/env bash
#SBATCH --job-name=s01.lfmm_file_creation
#SBATCH --output=./out/s00.sbatch/s01.lfmm_file_creation.out
#SBATCH --error=./out/s00.sbatch/s01.lfmm_file_creation.err
#SBATCH --cpus-per-task=4
#SBATCH -c 10
#SBATCH --mem=100G

# Set up environment
module load all gencore/3
module load vcftools


# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <vcf_folder> <extension> <output_folder>"
    exit 1
fi

# Assign input arguments to variables
vcf_folder=$1
extension=$2
output_folder=$3

mkdir -p $output_folder
# Find files matching the specified extension
files=( $(find "$vcf_folder" -maxdepth 1 -type f -name "$extension") )

extension="${extension/\*/}"
# Process each VCF file
for file in "${files[@]}"; do
    base_name=$(basename "$file" "$extension")
    echo "Processing $base_name..."

    # If the extension is *.vcf.gz, use --gzvcf, otherwise use --vcf
    if [[ "$extension" == *.vcf.gz ]]; then
        vcftools --gzvcf "$file" --012 --out "${output_folder}/${base_name}"
    else
        vcftools --vcf "$file" --012 --out "${output_folder}/${base_name}"
    fi

    # Extract coordinates and loci
    cut -f1,2 "${vcf_folder}/${base_name}${extension}" | grep -v "#" > "${output_folder}/coordinates_loci_${base_name}.txt"

    # Convert -1 to 9 in the .012 file and reformat for LFMM
    sed "s:-1:9:g" "${output_folder}/${base_name}.012" | cut -f 2- > "${output_folder}/2_${base_name}.012"
    mv "${output_folder}/2_${base_name}.012" "${output_folder}/${base_name}.lfmm"
done

echo "Processing complete."

