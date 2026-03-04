import os
import networkx as nx
import matplotlib.pyplot as plt

CUTOFF = 0.354
KEPT_COLOR = "#1EA01E"
REMOVED_COLOR = "#666666"

KIN_OUT = "./out/All_data_wo_orientalis/s01.kin_table"
FIG_OUT = "./out/All_data_wo_orientalis/s01.kin_fig"
VCF_BASE = "chr_concat.geno20_maf005.anno"

# Input files from PLINK
kept_file = os.path.join(KIN_OUT, VCF_BASE + ".king.cutoff.in.id")
removed_file = os.path.join(KIN_OUT, VCF_BASE + ".king.cutoff.out.id")
kin_file = os.path.join(KIN_OUT, VCF_BASE + ".kin0")

G = nx.Graph()
node_color = {}
legend_labels = {}

# ----------------- Kept samples -----------------
if os.path.exists(kept_file):
    with open(kept_file) as f_in:
        kept_ls = [line.strip() for line in f_in if line.strip() and not line.startswith("#")]
    for s in kept_ls:
        G.add_node(s)
        node_color[s] = KEPT_COLOR
else:
    kept_ls = []

legend_labels[KEPT_COLOR] = f"Kept: {len(kept_ls)} genotypes"

# ----------------- Removed samples -----------------
if os.path.exists(removed_file):
    with open(removed_file) as f_out:
        removed_ls = [line.strip() for line in f_out if line.strip() and not line.startswith("#")]
    for s in removed_ls:
        G.add_node(s)
        node_color[s] = REMOVED_COLOR
else:
    removed_ls = []

legend_labels[REMOVED_COLOR] = f"Removed (KING kinship > {CUTOFF}): {len(removed_ls)} genotypes"

# ----------------- Kinship edges -----------------
edge_weights = []
if os.path.exists(kin_file):
    with open(kin_file) as f_kin:
        lines = f_kin.read().strip().split("\n")
    if lines:
        header = lines[0].lstrip("#").split()
        id1_idx = header.index("IID1")
        id2_idx = header.index("IID2")
        kin_idx = header.index("KINSHIP")
        for row in lines[1:]:
            cols = row.split()
            ind1, ind2 = cols[id1_idx], cols[id2_idx]
            kin = float(cols[kin_idx])
            if kin >= CUTOFF and ind1 in node_color and ind2 in node_color:
                G.add_edge(ind1, ind2, weight=kin)
                edge_weights.append(kin)

# ----------------- Plot -----------------
plt.figure(figsize=(7.5, 7.5))
pos = nx.spring_layout(G, seed=4)
nx.draw(
    G,
    pos,
    node_color=[node_color[n] for n in G.nodes()],
    node_size=12,
    edge_color=edge_weights if edge_weights else "lightgray",
    edge_cmap=plt.cm.spring_r,
    width=1,
    with_labels=False
)

for color, text in legend_labels.items():
    plt.scatter([], [], c=[color], label=text)

plt.legend(title="Samples", loc="upper right", bbox_to_anchor=(1.35, 1.05))

if edge_weights:
    sm = plt.cm.ScalarMappable(
        cmap=plt.cm.spring_r,
        norm=plt.Normalize(vmin=min(edge_weights), vmax=max(edge_weights))
    )
    cax = plt.axes([1, 0.15, 0.02, 0.6])
    plt.colorbar(sm, cax=cax, label=f"KING kinship (cutoff={CUTOFF})")

plt.savefig(os.path.join(FIG_OUT, VCF_BASE + ".kinship_network.png"), dpi=300, bbox_inches="tight")
plt.show()
