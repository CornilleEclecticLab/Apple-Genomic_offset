import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ----------------------------
# Parameters
# ----------------------------
# PROJECT = "All_data_wo_orientalis"
PROJECT = "dom_admixed_rm"
KIN_OUT = f"./out/{PROJECT}/s01.kin_table"
# KIN_FILE = os.path.join(KIN_OUT, "chr_concat.geno20_maf005.anno.kin0")
KIN_FILE = os.path.join(KIN_OUT, "chr_concat.geno20_maf005.anno_dom_admixed_rm.kin0")

# META_FILE = "./in/data_tree_wo_orientalis_V2.txt"
META_FILE = "./in/data_tree_dom_admixed_rm_V3.txt"

FIG_OUT = f"./out/{PROJECT}/s03.kinship_heatmap/"
os.makedirs(FIG_OUT, exist_ok=True)

FIG_FILE = os.path.join(FIG_OUT, "kinship_heatmap.png")

# --------------------------------------------------
# 1) READ METADATA AND ORDER BY Order_plot
# --------------------------------------------------
print("Loading metadata…")
meta = pd.read_csv(META_FILE, sep="\t")

required_cols = {"Seq_ID", "Order_plot"}
if not required_cols.issubset(meta.columns):
    raise ValueError(f"Metadata must contain columns: {required_cols}")

meta_sorted = meta.sort_values("Order_plot")
ordered_ids = meta_sorted["Seq_ID"].tolist()

print(f"Ordering heatmap by Order_plot, {len(ordered_ids)} individuals.")

# --------------------------------------------------
# 2) READ KINSHIP TABLE (robust)
# --------------------------------------------------
print("Loading kinship table…")
df = pd.read_csv(KIN_FILE, sep="\t", comment=None)

df.columns = [c.strip() for c in df.columns]
print("Columns detected:", df.columns.tolist())

# Normalize the IID1 column name
if "#IID1" in df.columns:
    df.rename(columns={"#IID1": "IID1"}, inplace=True)

if "IID1" not in df.columns:
    raise ValueError("Could not find column IID1 or #IID1 in kinship file.")

required = {"IID1", "IID2", "KINSHIP"}
if not required.issubset(df.columns):
    raise ValueError(f"Kinship file must contain columns: {required}")

df = df[["IID1", "IID2", "KINSHIP"]].copy()

# --------------------------------------------------
# 3) BUILD SYMMETRIC MATRIX
# --------------------------------------------------
kin_matrix = df.pivot(index="IID1", columns="IID2", values="KINSHIP")

all_individuals = sorted(set(df["IID1"]).union(df["IID2"]))
kin_matrix = kin_matrix.reindex(index=all_individuals, columns=all_individuals)

# fill symmetry
kin_matrix = kin_matrix.combine_first(kin_matrix.T)

# --------------------------------------------------
# 4) APPLY Order_plot ORDER
# --------------------------------------------------
print("Reordering matrix according to Order_plot…")

ordered_ids_final = [x for x in ordered_ids if x in kin_matrix.index]

missing = [x for x in ordered_ids if x not in kin_matrix.index]
if missing:
    print(f"Warning: {len(missing)} samples in metadata not found in kinship matrix.")

kin_matrix = kin_matrix.loc[ordered_ids_final, ordered_ids_final]

print(f"Final ordered matrix size: {kin_matrix.shape[0]} x {kin_matrix.shape[1]}")

# --------------------------------------------------
# 5) PLOT
# --------------------------------------------------
plt.figure(figsize=(12, 10))
sns.heatmap(
    kin_matrix,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    square=True,
    linewidths=0.2,
    linecolor="gray",
    cbar_kws={"label": "KINSHIP coefficient"}
)

plt.title("Pairwise KING kinship coefficients")
plt.xlabel("Individuals")
plt.ylabel("Individuals")
plt.xticks(rotation=90, fontsize=7)
plt.yticks(rotation=0, fontsize=7)

plt.tight_layout()
plt.savefig(FIG_FILE, dpi=300)
plt.show()

plt.figure(figsize=(12, 10))
sns.heatmap(
    kin_matrix,
    cmap="coolwarm",
    vmin=-0.5,
    vmax=0.5,
    square=True,
    linewidths=0.2,
    linecolor="gray",
    cbar_kws={"label": "KINSHIP coefficient"}
)

plt.title("Pairwise KING kinship coefficients")
plt.xlabel("Individuals")
plt.ylabel("Individuals")
plt.xticks(rotation=90, fontsize=7)
plt.yticks(rotation=0, fontsize=7)

plt.tight_layout()

# Save as SVG
plt.savefig(FIG_FILE.replace(".png", ".svg"), format="svg")

plt.show()

print(f"Heatmap saved to: {FIG_FILE}")
