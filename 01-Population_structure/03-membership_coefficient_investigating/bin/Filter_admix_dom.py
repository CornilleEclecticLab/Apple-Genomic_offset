import os
import pandas as pd
import re


K_value=5
threshold=0.2

# For fastStructure
model="fastStructure"
input_folder="../02-fastStructure/out/All_data_wo_orientalis/s05.membership_coefficient"
ouput_folder="./out"
domestica_group="V2"

# For snmf
#model="snmf"
#input_folder="../01-snmf/out/All_data_wo_orientalis/s01.snmf/K_value"
#ouput_folder="./out"
#domestica_group="V1" # K=5
#domestica_group="V3" # K=6


# Check if the input folder exists and contains the expected files
if not os.path.exists(input_folder):
    raise FileNotFoundError(f"The folder {input_folder} does not exist.")
    
file_to_filter = next((f for f in os.listdir(input_folder) if re.search(str(K_value), f)), None)
if file_to_filter is None:
    raise FileNotFoundError(f"No file found in {input_folder} matching K={K_value}")

# Read the file
input_file = os.path.join(input_folder, file_to_filter)
file_df = pd.read_csv(input_file, sep=",", header=0)
print(file_df.head())

# if V2 is higher than 0.01 then it's admix
admix_df = file_df[file_df[domestica_group] >= threshold]
admix_df = admix_df[["Seq_ID", "Country", "Site", domestica_group]]
print(admix_df)

# Save the list of the Seq_ID of admixed individuals
admix_list = admix_df["Seq_ID"].tolist()
admix_list = [str(i) for i in admix_list]
admix_list = list(set(admix_list))
admix_list.sort()
print(len(admix_list))
# Save the list to a file
output_file = os.path.join(ouput_folder, f"{model}_admix_list_K{K_value}_{threshold}.txt")
with open(output_file, "w") as f:
    for item in admix_list:
        f.write("%s\n" % item)


