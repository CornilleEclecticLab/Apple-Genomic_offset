import os
import random
import argparse

# By Xilong CHEN, adapted to Python by Maxime Criado
# Create Date (Perl): 4/11/2020
# Python Adaptation: 2025/04/28

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate fastStructure SLURM jobs and helper scripts.")
    parser.add_argument('-b', '--bedfile', required=True, help="BED file prefix (without .bed extension)")
    parser.add_argument('-o', '--output_dir', default='.', help="Output directory")
    parser.add_argument('-minK', type=int, default=2, help="Minimum K value (default: 2)")
    parser.add_argument('-maxK', type=int, default=19, help="Maximum K value (default: 19)")
    parser.add_argument('-r', '--repeat', type=int, default=2, help="Repeat number for each K (default: 2)")
    parser.add_argument('-f', '--output_in_file', default='.', help="Output directory in file")
    return parser.parse_args()

def create_faststructure_jobs(bedfile, output_dir, minK, maxK, repeat, output_in_file):
    faststructure_py = "structure.py"

    for k in range(minK, maxK + 1):
        k_dir = os.path.join(output_dir, f"K{k}")
        os.makedirs(k_dir, exist_ok=True)

        for run_idx in range(1, repeat + 1):
            seed = random.randint(0, 100000)
            job_name = f"K{k}.run{run_idx}.sh"
            script_path = os.path.join(k_dir, f"run{run_idx}.sh")
            with open(script_path, 'w') as f:
                f.write(f"""#!/bin/bash

#SBATCH -J {job_name}
#SBATCH --output {output_in_file}/s00.sbatch/{job_name}.out
#SBATCH --error {output_in_file}/s00.sbatch/{job_name}.err

#Activating conda
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate faststructure

{faststructure_py} -K {k} --cv 10 --full --seed={seed} --input={bedfile} --output={output_in_file}/s01.faststructure/run{run_idx}

""")

def create_launcher_script(output_dir, minK, maxK):
    launcher_script = os.path.join(output_dir, "00.group_submits.sh")
    with open(launcher_script, 'w') as f:
        f.write(f"""#!/bin/bash

for rep in $(seq 2 10)
do
    ls {output_dir}/*/*.sh | awk '{{print "sbatch -c 1 --mem=1G " $0}}' > {output_dir}/faststructure_job_list.txt
    
done
""")

def create_zip_script(output_dir, minK, maxK, output_in_file):
    zip_script = os.path.join("bin","s02.zip_for_clumpak.sh")
    with open(zip_script, 'w') as f:
        f.write(f"""#!/bin/bash

cd {output_in_file}/s01.faststructure
for i in $(seq 2 10)
do
    zip -q K${{i}}.zip *.${{i}}.meanQ
done

zip -q ../input_clumpak.zip K*.zip
rm -rf K*.zip

""")

def main():
    args = parse_arguments()
    create_faststructure_jobs(args.bedfile, args.output_dir, args.minK, args.maxK, args.repeat, args.output_in_file)
    create_launcher_script(args.output_dir, args.minK, args.maxK)
    create_zip_script(args.output_dir, args.minK, args.maxK, args.output_in_file)
    print("✔️ Job scripts and helper scripts generated successfully.")

if __name__ == "__main__":
    main()
