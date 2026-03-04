# module load miniconda
# conda activate r-chelsa-env

# ===================
# = current climate =
# ===================

# https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio01/1981-2010/CHELSA_bio01_1981-2010_V.2.1.tif
rm(list=ls())

bio_nb <- 1:19
bio_padded <- sprintf("%02d", bio_nb)

base_url <- "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim"
dest_root <- "./out/s01.tif/current_climate"

current_dir <- dest_root
  if (!dir.exists(current_dir)) {
        dir.create(current_dir, recursive = TRUE)
      }

options(timeout = 100000) # set a longer time if needed
for (i in bio_nb) {
  # Format bio number (e.g., bio01)
   bio_str <- paste0("bio", sprintf("%02d", i))
  
  file_name <- paste0("CHELSA_", bio_str, "_1981-2010_V.2.1.tif")

  url_path <- paste(base_url, bio_str, "1981-2010", file_name, sep = "/")

  # 3. Construct Destination Path
  dest_path <- file.path(current_dir, file_name)

  # 4. Download
  if (!file.exists(dest_path)) {
    message(paste("Accessing:", url_path))
    tryCatch({
      # Added mode = "wb" for binary download (crucial for .tif files on Windows)
      download.file(url_path, destfile = dest_path, mode = "wb", quiet = FALSE)
    }, error = function(e) {
      message(paste("Skipping: Could not find", bio_str, "for", ssp, "in", time))
    })
  }
}


# =======================
# == Future scenarios ===
# =======================

# rm(list=ls())
# Setup Parameters
# bio_nb <- 1:19
# bio_padded <- sprintf("%02d", bio_nb)

# base_url <- "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim"
# ssp_models <- c("ssp126", "ssp370", "ssp585")
# time_periods <- c("2011-2040", "2041-2070", "2071-2100")
# gcm_model <- "IPSL-CM6A-LR"
# version <- "V.2.1"
# dest_root <- "./out/s01.tif/"

# if (!dir.exists(dest_root)) {
#         dir.create(dest_root, recursive = TRUE)
#       }

# #getOption("timeout") # check the time out setting
# options(timeout = 100000) # set a longer time if needed
# # Main Loop
# for (i in bio_nb) {
#   # Format bio number (e.g., bio01)
#   bio_str <- paste0("bio", sprintf("%02d", i))
  
#   for (time in time_periods) {
#     for (ssp in ssp_models) {
      
#       # 1. Create directory structure
#       current_dir <- file.path(dest_root,"future_climate",time, ssp)
#       if (!dir.exists(current_dir)) {
#         dir.create(current_dir, recursive = TRUE)
#       }
      
#       # 2. Construct the URL components carefully
#       # The filename requires the GCM name in lowercase
#       file_name <- paste0("CHELSA_", tolower(gcm_model), "_", ssp, "_", bio_str, "_", time, "_", version, ".tif")
      
#       # Construct URL: base/bio01/2011-2040/IPSL-CM6A-LR/ssp126/filename
#       url_path <- paste(base_url, bio_str, time, gcm_model, ssp, file_name, sep = "/")
      
#       # 3. Construct Destination Path
#       dest_path <- file.path(current_dir, file_name)
      
#       # 4. Download
#       if (!file.exists(dest_path)) {
#         message(paste("Accessing:", url_path))
#         tryCatch({
#           # Added mode = "wb" for binary download (crucial for .tif files on Windows)
#           download.file(url_path, destfile = dest_path, mode = "wb", quiet = FALSE)
#         }, error = function(e) {
#           message(paste("Skipping: Could not find", bio_str, "for", ssp, "in", time))
#         })
#       } 
#     }
#   }
# }


# ================
# = Past climate =
# ================
# rm(list=ls())
# https://os.zhdk.cloud.switch.ch/chelsav1/pmip3/IPSL_IPSL-CM5A-LR/CHELSA_PMIP_IPSL-CM5A-LR_BIO_01.tif

# bio_nb <- 1:19
# bio_padded <- sprintf("%02d", bio_nb)

# base_url <- "https://os.zhdk.cloud.switch.ch/chelsav1/pmip3/IPSL_IPSL-CM5A-LR"
# dest_root <- "./out/s01.tif/LGM_climate"

# current_dir <- dest_root
#   if (!dir.exists(current_dir)) {
#         dir.create(current_dir, recursive = TRUE)
#       }

# options(timeout = 100000) # set a longer time if needed
# for (i in bio_nb) {
#   # Format bio number (e.g., bio01)
#   bio_str <- sprintf("%02d", i)
  
#   file_name <- paste0("CHELSA_PMIP_IPSL-CM5A-LR_BIO_", bio_str, ".tif")

#   url_path <- paste(base_url, file_name, sep = "/")

#   # 3. Construct Destination Path
#   dest_path <- file.path(current_dir, file_name)

#   # 4. Download
#   if (!file.exists(dest_path)) {
#     message(paste("Accessing:", url_path))
#     tryCatch({
#       # Added mode = "wb" for binary download (crucial for .tif files on Windows)
#       download.file(url_path, destfile = dest_path, mode = "wb", quiet = FALSE)
#     }, error = function(e) {
#       message(paste("Skipping: Could not find", bio_str, "for", ssp, "in", time))
#     })
#   }
# }