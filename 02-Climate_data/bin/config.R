## config file
library(raster)
library(rgdal)
library(dplyr)
library(stringr)

output_folder = "./out/"

## get tif and convert to ascii

db_link = ""
folder_name = ""
bio.nb_tif = 1:19 # Based on the number of bio climatic variable : DO NOT CHANGE IT 


## extract climate data from coordinate

coordinate_file = "./in/info_samples.txt"
ascii_folder = "~/Bureau/Maxime_PhD/02-Analysis/Climatic_data/CHELSA_ascii/"
bio.nb = 1:19 ## Based on the number of bio climatic variable : DO NOT CHANGE IT
bio_num = c(1:20) #there is a ID column + 19 bio variables
out_name = "Bioclim_sample_chelsa"
