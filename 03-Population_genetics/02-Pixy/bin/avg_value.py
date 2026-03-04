import os
import pandas as pd

# Configuration
# project = "dom_admixed_rm"
# project = "full_genome_dom_admixed_rm_no_filter"
project = "full_genome_dom_admixed_rm_no_filter_pi_by_ind"

input_dir = f"./out/{project}/s01.Pixy_by_chromosome"
output_dir = f"./out/{project}/s02.avg_value"
os.makedirs(output_dir, exist_ok=True)

# Function to read and group data by population(s)
def read_and_group(file_path, file_type):
    df = pd.read_csv(file_path, sep='\t')
    pop_windows = {}

    if file_type == 'pi':
        df = df.sort_values(by=['pop', 'chromosome'])
        for _, row in df.iterrows():
            key = row['pop']
            pop_windows.setdefault(key, []).append({
                'count_diffs': row['count_diffs'],
                'count_comparisons': row['count_comparisons']
            })

    elif file_type == 'dxy':
        df = df.sort_values(by=['pop1', 'pop2', 'chromosome'])
        for _, row in df.iterrows():
            key = f"{row['pop1']}-{row['pop2']}"
            pop_windows.setdefault(key, []).append({
                'count_diffs': row['count_diffs'],
                'count_comparisons': row['count_comparisons']
            })

    elif file_type == 'fst':
        df = df.sort_values(by=['pop1', 'pop2', 'chromosome'])
        for _, row in df.iterrows():
            key = f"{row['pop1']}-{row['pop2']}"
            pop_windows.setdefault(key, []).append({
                'fst_wind': row['avg_wc_fst']
            })

    else:
        raise ValueError(f"Unknown file type: {file_type}")

    return pop_windows

# Function to compute average values and write to output
def compute_and_write(pop_windows, output_file, file_type):
    with open(output_file, 'w') as f:
        if file_type in ['fst', 'dxy']:
            metric_name = "avg_fst" if file_type == "fst" else "avg_dxy"
            f.write("pop1,pop2," + metric_name + "\n")
            for key, values in pop_windows.items():
                metrics = []
                for v in values:
                    if file_type == "fst":
                        val = v['fst_wind']
                    else:  # dxy
                        if pd.notnull(v['count_comparisons']) and v['count_comparisons'] != 0:
                            val = v['count_diffs'] / v['count_comparisons']
                        else:
                            val = None
                    if val is not None and pd.notnull(val):
                        metrics.append(val)
                avg_val = sum(metrics) / len(metrics) if metrics else float('nan')
                pop1, pop2 = key.split("-")
                f.write(f"{pop1},{pop2},{avg_val}\n")

        elif file_type == 'pi':
            f.write("pop,avg_pi\n")
            for key, values in pop_windows.items():
                total_diffs = sum(v['count_diffs'] for v in values if pd.notnull(v['count_diffs']))
                total_comparisons = sum(v['count_comparisons'] for v in values if pd.notnull(v['count_comparisons']))
                avg_val = total_diffs / total_comparisons if total_comparisons else float('nan')
                f.write(f"{key},{avg_val}\n")

# Loop through all files in the input directory
for path, _, files in os.walk(input_dir):
    for file in files:
        if file.endswith(('pi.txt', 'dxy.txt', 'fst.txt')):
            print(f"Processing file: {file}")
            file_type = file.split(".")[2].split("_")[1]  # 'pi', 'dxy', or 'fst'
            file_path = os.path.join(path, file)

            # Read and group data
            pop_windows = read_and_group(file_path, file_type)

            # Extract chromosome from filename
            chrom = file.split(".")[2].split("_")[0]

            # Output file naming
            suffix = {
                "pi": "_avg_pi.csv",
                "dxy": "_avg_dxy.csv",
                "fst": "_avg_fst.csv"
            }[file_type]

            output_file = os.path.join(output_dir, f"{chrom}{suffix}")

            # Compute and write output
            compute_and_write(pop_windows, output_file, file_type)
