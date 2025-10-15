import os
import pandas as pd
import re

K_value=6
#model="FastStructure"
#input_folder="../02-fastStructure/out/All_data_wo_orientalis/s05.membership_coefficient"
#ouput_folder="./out/s01.membership_coefficient_mean_by_site"
#
##model="snmf"
##input_folder="../01-snmf/out/admix_dom_removed_0.2/s01.snmf/K_value"
##ouput_folder="./out/s01.membership_coefficient_mean_by_site"
##threshold=0.2
#
#
#
## Check if the input folder exists and contains the expected files
#if not os.path.exists(input_folder):
#    raise FileNotFoundError(f"The folder {input_folder} does not exist.")
#    
#file_to_filter = next((f for f in os.listdir(input_folder) if re.search(str(K_value), f)), None)
#if file_to_filter is None:
#    raise FileNotFoundError(f"No file found in {input_folder} matching K={K_value}")
#
## Read the file
#input_file = os.path.join(input_folder, file_to_filter)
#file_df = pd.read_csv(input_file, sep=",", header=0)
#
## mean value of V1, V2, ... by site and also mean the longitude and latitude
#file_df = file_df.loc[:, ["Site", "Latitude", "Longitude"] + [col for col in file_df.columns if col.startswith("V")]]
#
## Nan value in site put Domestica
#file_df["Site"] = file_df["Site"].fillna("Domestica")
#
## Nan value in Longitude and Latitude put Longitude -7.075 and latitude 48.005
#file_df["Longitude"] = file_df["Longitude"].fillna(-7.075)
#file_df["Latitude"] = file_df["Latitude"].fillna(48.005)
#
## change the , in . for the Longitude and Latitude
#file_df["Longitude"] = file_df["Longitude"].astype(str).str.replace(",", ".")
#file_df["Latitude"] = file_df["Latitude"].astype(str).str.replace(",", ".")
#
## and also mean the longitude and latitude
#file_df["Longitude"] = pd.to_numeric(file_df["Longitude"], errors='coerce')
#file_df["Latitude"] = pd.to_numeric(file_df["Latitude"], errors='coerce')
#
## mean by site
#file_df = file_df.groupby("Site").mean()
#file_df = file_df.reset_index()
#file_df = file_df.loc[:, ["Site", "Latitude", "Longitude"] + [col for col in file_df.columns if col.startswith("V")]]
#
#file_df = file_df.rename(columns={"Longitude": "Longitude_mean", "Latitude": "Latitude_mean"})
#file_df = file_df.rename(columns={col: f"{col}_mean" for col in file_df.columns if col.startswith("V")})
#
## save the file
#if not os.path.exists(ouput_folder):
#    os.makedirs(ouput_folder)
#output_file = os.path.join(ouput_folder, f"K{K_value}_{model}__membership_coefficient_mean_by_site.csv")
#file_df.to_csv(output_file, sep=",", index=False)
#print(f"File saved in {output_file}")




model="snmf"
#input_file="./in/Table_ms - ST3_ membership table with fastStructure K6.csv"
#input_file="./in/Structure de pop snmf _ FS - K6_Minor_FS.csv"
input_file="./in/Table_ms - SNMF K6.csv"
#model="snmf"
#input_file="./in/Table_ms - ST4_ membership table with snmf k6.csv"
#input_file="./in/Structure de pop snmf _ FS - K9_snmf.csv"
#input_folder="./in"
ouput_folder="./out/s01.membership_coefficient_mean_by_site"

file_df = pd.read_csv(input_file, sep=",", header=0, decimal=".")
print(file_df)

# mean value of V1, V2, ... by Sub_group and also mean the longitude and latitude
file_df = file_df.loc[:, ["Sub_group", "Latitude", "Longitude"] + [col for col in file_df.columns if col.startswith("V")]]



# Nan value in Sub_group put Domestica
#file_df["Sub_group"] = file_df["Sub_group"].fillna("Domestica")

# Nan value in Longitude and Latitude put Longitude -7.075 and latitude 48.005
file_df["Longitude"] = file_df["Longitude"].fillna(-7.075)
file_df["Latitude"] = file_df["Latitude"].fillna(48.005)

# change the , in . for the Longitude and Latitude
file_df["Longitude"] = file_df["Longitude"].astype(str).str.replace(",", ".")
file_df["Latitude"] = file_df["Latitude"].astype(str).str.replace(",", ".")

# and also mean the longitude and latitude
file_df["Longitude"] = pd.to_numeric(file_df["Longitude"], errors='coerce')
file_df["Latitude"] = pd.to_numeric(file_df["Latitude"], errors='coerce')
file_df.iloc[:, 3:] = file_df.iloc[:, 3:].apply(pd.to_numeric, errors='coerce')
print(file_df)

# mean by Sub_group
file_df = file_df.groupby("Sub_group").mean(numeric_only=True)
print(file_df)
file_df = file_df.reset_index()
file_df = file_df.loc[:, ["Sub_group", "Latitude", "Longitude"] + [col for col in file_df.columns if col.startswith("V")]]

file_df = file_df.rename(columns={"Longitude": "Longitude_mean", "Latitude": "Latitude_mean"})
file_df = file_df.rename(columns={col: f"{col}_mean" for col in file_df.columns if col.startswith("V")})


# save the file
if not os.path.exists(ouput_folder):
    os.makedirs(ouput_folder)
output_file = os.path.join(ouput_folder, f"K{K_value}_Major_{model}_membership_coefficient_mean_by_sub_group.csv")
file_df.to_csv(output_file, sep=",", index=False)
print(f"File saved in {output_file}")
