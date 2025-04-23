## code to extract bio variable from climate database and convert the tif file into ascii
## Rscript get_tif_and_conv_ascii.R -param 1 = link -param 2 = folder_name
library(raster)
library(rgdal)
library(dplyr)
library(stringr)

source("./bin/config.R")
source("./bin/function.R")

out = "./out"
create_directory(out)
output_tif = paste0(out, folder_name, "_tif")
output_ascii = paste0(out, folder_name, "_ascii")

bio.nb1 <- 1:19 # define the bioclim numbering
bio.nb2 <- sprintf("%02d", 1:19)

path1 <- "../Climatic_data/"
path2 <-
  "https://os.zhdk.cloud.switch.ch/chelsav1/pmip3/IPSL_IPSL-CM5A-LR/"
path3 <- "CHELSA_LGM/" # to create if not created
create_directory(path3)
path4 <- "CHELSA_LGM_ascii/" # to create if not created
create_directory(path4)
#### Download environmental data ####
# Current climate -- download rasters from CHELSA
# browseURL(
#   "https://chelsa-climate.org/downloads/",
#   browser = getOption("browser"),
#   encodeIfNeeded = F
# )
getOption("timeout") # check the time out setting
options(timeout = 1000) # set a longer time if needed
for (i in 19:length(bio.nb1)) {
  if (i<10){
    download.file(
      paste(path2, "CHELSA_PMIP_IPSL-CM5A-LR_BIO_0", bio.nb1[i], ".tif", sep = ""),
      paste(
        path1,
        path3,
        "CHELSA_PMIP_IPSL-CM5A-LR_BIO_",
        bio.nb2[i],
        ".tif",
        sep = ""
      )
    )
  } else {
    download.file(
      paste(path2, "CHELSA_PMIP_IPSL-CM5A-LR_BIO_", bio.nb1[i], ".tif", sep = ""),
      paste(
        path1,
        path3,
        "CHELSA_PMIP_IPSL-CM5A-LR_BIO_",
        bio.nb2[i],
        ".tif",
        sep = ""
      )
    )
  }
}

varnames <- list.files(path3)
bioclim = stack()
for (i in 1:length(bio.nb1)) {
  bioclim <-
    addLayer(bioclim, raster(paste(path3, varnames[i], sep = "")))
}


#### Convertion Tif to Acsii #####
for (i in 1:length(bio.nb1)) {
  read_file <- raster(paste(path3, varnames[i], sep = ""))
  writeRaster(read_file,
              paste(path4, varnames[i], sep = ""),
              format = "ascii",
              overwrite = TRUE)
}
