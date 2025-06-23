#!/usr/bin/env python3

import sys
import glob
import os
import re

if len(sys.argv) != 5:
    print(f"Usage: python {sys.argv[0]} <clumpak_result_folder> <min_K> <max_K> <output_folder>")
    sys.exit(1)

re_folder = sys.argv[1]
min_k = int(sys.argv[2])
max_k = int(sys.argv[3])
output_folder = sys.argv[4]

os.makedirs(output_folder, exist_ok=True)

for k_f in range(min_k, max_k + 1):
    cluster_dirs = glob.glob(os.path.join(re_folder, f"K={k_f}", "*Cluster*"))

    for cluster_path in cluster_dirs:
        cluster_name = os.path.basename(cluster_path)
        clumpp_file = os.path.join(cluster_path, "CLUMPP.files", "ClumppIndFile.output")

        if os.path.isfile(clumpp_file):
            print(f"Found: {clumpp_file}")
            k_f2 = f"{k_f:02d}"
            output_path = os.path.join(output_folder, f"K.{k_f2}.{cluster_name}.Q")

            with open(clumpp_file, 'r') as infile, open(output_path, 'w') as outfile:
                for line in infile:
                    line = re.sub(r'^.+?:\s*', '', line)
                    outfile.write(line)
        else:
            print(f"Warning: Missing file - {clumpp_file}")

