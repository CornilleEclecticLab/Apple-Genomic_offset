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
    pattern = os.path.join(re_folder, f"K={k_f}", "MajorCluster", "CLUMPP.files", "ClumppIndFile.output")
    qfiles = glob.glob(pattern)

    for f in qfiles:
        print(f)
        # ✅ Correct the extraction: get only the 'MajorCluster' part
        cluster_name = os.path.basename(os.path.dirname(os.path.dirname(f)))  # 'MajorCluster'

        print(f"{k_f}/K.{k_f}.{cluster_name}.Q")
        k_f2 = f"{k_f:02d}"
        output_path = os.path.join(output_folder, f"K.{k_f2}.{cluster_name}.Q")

        # ✅ Now no broken path anymore
        with open(f, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                line = re.sub(r'^.+?:\s*', '', line)
                outfile.write(line)
