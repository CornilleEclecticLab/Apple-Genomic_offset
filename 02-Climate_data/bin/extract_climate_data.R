## Function to extract climate data based on coordinate
## Rscript extract_climate_data.R -param 1 = file with coordinate (3 columns ID_site,X,Y) -param 2= folder with ascii file

## coord_file : folder with the coordinate
## ascii_folder : folder of ascii file from which you want to extract data 
## bio.num : list of bio variable wanted
## out_name : name of the output file

source("./bin/config.R")
source("./bin/function.R")

create_directory(output_folder)

bio.nb <- 1:19 ## we are going to get the 19 value and filtered depending on the bio wanted
varnames <- list.files(ascii_folder)
print(varnames)

data_coord <- read.delim(coordinate_file, sep = "\t", dec = ".", header = T)
data_coord$Latitude <- as.numeric(as.character(data_coord$Latitude))
data_coord$Longitude <- as.numeric(as.character(data_coord$Longitude))
print(str(data_coord))
print(data_coord$Latitude)
print(data_coord$Longitude)

bioclim = stack()
for (i in 1:length(bio.nb)) {
  bioclim <-
    addLayer(bioclim, raster(paste(ascii_folder, varnames[i], sep = "")))
}

env_data <-
  extract(
    bioclim,
    data_coord[c("Latitude", "Longitude")],
    method = 'simple',
    buffer = 30000 ,
    fun = mean,
    df = TRUE,
    na.rm = TRUE
  )

reduce_bio_name <- c(
  "ID",
  "Bio1",
  "Bio2",
  "Bio3",
  "Bio4",
  "Bio5",
  "Bio6",
  "Bio7",
  "Bio8",
  "Bio9",
  "Bio10",
  "Bio11",
  "Bio12",
  "Bio13",
  "Bio14",
  "Bio15",
  "Bio16",
  "Bio17",
  "Bio18",
  "Bio19"
)

colnames(env_data) <- reduce_bio_name
env_data_filtered <- env_data[,bio_num]

env_data_filtered <- cbind(data_coord[1], env_data_filtered[,2:length(env_data_filtered)])

write.table(env_data_filtered, paste0(output_folder, out_name,".csv"), quote=FALSE, sep="\t", row.names = F)
