## config file
library(raster)
library(rgdal)
library(dplyr)
library(stringr)


# Project = "dom_admixed_rm"
Project = "all_wild_samples"
output_folder = paste0("./out/",Project,"/")

## get tif and convert to ascii

db_link = ""
folder_name = ""
bio.nb_tif = 1:19 # Based on the number of bio climatic variable : DO NOT CHANGE IT 


## extract climate data from coordinate

#coordinate_file = "./in/coordinate_dom_admixed_rm.txt"
# coordinate_file = "./in/group_coordinate_dom_admixed_rm_V2.txt"
# coordinate_file = "./in/coordinate_dom_admixed_rm_V2.txt"
coordinate_file = "./in/coordinate_all_wild_samples.txt"

ascii_folder = "./in/CHELSA_current_ascii/"
# ascii_folder = "./in/CHELSA_LGM_ascii/"
#ascii_folder = "~/Bureau/Maxime_PhD/02-Analysis/Climatic_data/CHELSA_future_ascii/"
bio.nb = 1:19 ## Based on the number of bio climatic variable : DO NOT CHANGE IT
bio_num = c(1:20) #there is a ID column + 19 bio variables
out_name = paste0("Bioclim_",Project,"_chelsa_current")
# out_name = paste0("Bioclim_orchard_chelsa_current")