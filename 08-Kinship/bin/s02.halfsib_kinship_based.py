#!/usr/bin/env python3
import networkx as nx
import os

# ----------------------------
# Input parameters
# ----------------------------
PROJECT = "All_data_wo_orientalis"
KIN_OUT = f"./out/{PROJECT}/s01.kin_table"
KIN_FILE = os.path.join(KIN_OUT, "chr_concat.geno20_maf005.anno.kin0")
HALF_SIB_THRESHOLD = 0.076 #0.125
OUT_DIR= f"./out/{PROJECT}/s02.halfsib_kinship_based/"
OUTPUT_FILE = os.path.join(OUT_DIR, f"halfsib_groups_{str(HALF_SIB_THRESHOLD)}.txt")
os.makedirs(OUT_DIR, exist_ok=True)

# ----------------------------
# Build graph
# ----------------------------
G = nx.Graph()

with open(KIN_FILE) as f:
    lines = f.read().strip().split("\n")

header = lines[0].lstrip("#").split()
id1_idx = header.index("IID1")
id2_idx = header.index("IID2")
kin_idx = header.index("KINSHIP")

for row in lines[1:]:
    cols = row.split()
    ind1, ind2 = cols[id1_idx], cols[id2_idx]
    kin = float(cols[kin_idx])
    if kin >= HALF_SIB_THRESHOLD:
        G.add_edge(ind1, ind2)

# Add isolated individuals (no half-sibling connection)
all_individuals = set()
for row in lines[1:]:
    all_individuals.update([row.split()[id1_idx], row.split()[id2_idx]])

for ind in all_individuals:
    if ind not in G:
        G.add_node(ind)

# ----------------------------
# Extract groups
# ----------------------------
groups = list(nx.connected_components(G))

# ----------------------------
# Print and save groups
# ----------------------------
with open(OUTPUT_FILE, "w") as out_f:
    print(f"Groups of individuals with kinship >= {HALF_SIB_THRESHOLD} (half-sibling or closer):\n")
    out_f.write(f"# Groups of individuals with kinship >= {HALF_SIB_THRESHOLD}\n\n")
    
    for i, grp in enumerate(groups, 1):
        grp_sorted = sorted(grp)
        line = f"Group {i} ({len(grp_sorted)} individuals): {', '.join(grp_sorted)}"
        print(line)
        out_f.write(line + "\n")

print(f"\nGroups saved to: {OUTPUT_FILE}")
