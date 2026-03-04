import os
import pandas as pd
import re

def dom_admixed_filtering(Project, model, in_fold, K_value, domestica_group, Threshold=0.1):
    """
    Generic function to filter admixed individuals based on membership coefficients.
    
    Parameters
    ----------
    Project : str
        Name of the project (used in paths)
    model : str
        Model name ('FastStructure', 'SNMF', etc.)
    in_fold : str
        Input folder containing the membership coefficient files
    K_value : int
        Number of clusters (K)
    domestica_group : str
        Column name corresponding to the domestica ancestry group (e.g. 'V2' or 'V3')
    Threshold : float
        Minimum membership coefficient to consider an individual admixed
    """

    out_fold = f"./out/{Project}/s01.filter_crops_wild/"
    os.makedirs(out_fold, exist_ok=True)

    # --- Check input folder ---
    if not os.path.exists(in_fold):
        raise FileNotFoundError(f"The folder {in_fold} does not exist.")

    # --- Find file corresponding to K_value ---
    file_to_filter = next((f for f in os.listdir(in_fold) if re.search(str(K_value), f)), None)
    if file_to_filter is None:
        raise FileNotFoundError(f"No file found in {in_fold} matching K={K_value}")

    file_path = os.path.join(in_fold, file_to_filter)
    print(f"[{model}] Found file: {file_path}")

    # --- Read the data ---
    # Try comma or tab separator automatically
    try:
        file_df = pd.read_csv(file_path, sep=",", header=0)
    except pd.errors.ParserError:
        file_df = pd.read_csv(file_path, sep="\t", header=0)

    # Check columns exist
    required_cols = ["Seq_ID", "Country", "Site", domestica_group]
    for col in required_cols:
        if col not in file_df.columns:
            raise KeyError(f"Column '{col}' not found in {file_path}")

    # --- Filter admixed individuals ---
    admix_df = file_df[file_df[domestica_group] > Threshold][required_cols]
    print(f"[{model}] {admix_df.shape[0]} admixed individuals found (>{Threshold}).")

    # --- Save result ---
    out_file = os.path.join(out_fold, f"{model}_admixed_K{K_value}_{domestica_group}.csv")
    admix_df['Seq_ID'].to_csv(out_file, index=False)
    print(f"[{model}] Filtered admixed individuals saved to: {out_file}\n")

    return admix_df


# =========================================================
# Example usage for both models
# =========================================================

Project = "All_data_wo_orientalis"
K_value = 6
Threshold = 0.1

# --- FastStructure ---
fast_in = f"../02-fastStructure/out/{Project}/s05.membership_coefficient/"
admix_fast = dom_admixed_filtering(
    Project=Project,
    model="FastStructure",
    in_fold=fast_in,
    K_value=K_value,
    domestica_group="V1",
    Threshold=Threshold
)

# --- SNMF ---
snmf_in = f"../01-snmf/out/{Project}/s01.snmf/K_value"
admix_snmf = dom_admixed_filtering(
    Project=Project,
    model="SNMF",
    in_fold=snmf_in,
    K_value=K_value,
    domestica_group="V1",
    Threshold=Threshold
)
