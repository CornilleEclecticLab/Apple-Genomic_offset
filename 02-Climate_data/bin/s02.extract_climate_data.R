rm(list=ls())

library(terra)
library(dplyr)

project = "dom_admixed_rm"
coords <- read.table("./in/coordinate_dom_admixed_rm_V2.txt", header = T, sep="\t", dec=",")
# Convert to spatial points (WGS84)
pts <- vect(coords, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

output_dir <- paste0("./out/s02.climate_data/",project)
if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }

# ==========================
# current climate extraction
# ==========================
message("Processing current climate")

current_root <- "./out/s01.tif/current_climate"
tif_files <- list.files(current_root, pattern = "\\.tif$", full.names = TRUE)
climate_stack <- rast(tif_files)
      
# Extract values
extracted <- extract(climate_stack, pts, ID = FALSE)
# Clean up column names (removes the long CHELSA file prefix)
# e.g., "CHELSA_ipsl..._bio01_..." becomes "bio01"
names(extracted) <- gsub(".*_(bio\\d+)_.*", "\\1", names(extracted))

# Add metadata columns
combined <- cbind(Seq_ID = coords$Seq_ID, extracted)

file_out_name <- paste0("climate_extract_current_", project, ".csv")
full_out_path <- file.path(output_dir, file_out_name)

# 5. Save individual file
write.csv(combined, full_out_path, row.names = FALSE)

# ========================
# future climate extraction
# ========================
message("Processing future climate")

future_root <- "./out/s01.tif/future_climate"
# 2. Setup paths (matching your download script)
ssp_models <- c("ssp126", "ssp370", "ssp585")
time_periods <- c("2011-2040", "2041-2070", "2071-2100")

# 3. Storage list for results
all_extractions <- list()

# 4. Loop through folders to extract
for (time in time_periods) {
  for (ssp in ssp_models) {
    
    folder_path <- file.path(future_root, time, ssp)
    tif_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)
    
    if (length(tif_files) > 0) {
      message(paste("Processing:", time, "|", ssp))
      
      # Load all 19 variables as a stack
      climate_stack <- rast(tif_files)
      
      # Extract values
      extracted <- extract(climate_stack, pts, ID = FALSE)
      
      # Clean up column names (removes the long CHELSA file prefix)
      # e.g., "CHELSA_ipsl..._bio01_..." becomes "bio01"
      names(extracted) <- gsub(".*_(bio\\d+)_.*", "\\1", names(extracted))
      
      # Add metadata columns
      combined <- cbind(Seq_ID = coords$Seq_ID, extracted)
      
      file_out_name <- paste0("climate_extract_future_", project ,"_",time, "_", ssp, ".csv")

      full_out_path <- file.path(output_dir, file_out_name)
      
      # 5. Save individual file
      write.csv(combined, full_out_path, row.names = FALSE)

    }
  }
}

# ==========================
# LGM climate extraction
# ==========================
message("Processing LGM climate")

past_root <- "./out/s01.tif/LGM_climate"
tif_files <- list.files(past_root, pattern = "\\.tif$", full.names = TRUE)
climate_stack <- rast(tif_files)
      
# Extract values
extracted <- extract(climate_stack, pts, ID = FALSE)
# Clean up column names (removes the long CHELSA file prefix)
# e.g., "CHELSA_ipsl..._bio01_..." becomes "bio01"
# 1. Extract and convert to lowercase
names(extracted) <- gsub(".*_(BIO)_(\\d+).*", "bio\\2", names(extracted))

# Add metadata columns
combined <- cbind(Seq_ID = coords$Seq_ID, extracted)

file_out_name <- paste0("climate_extract_LGM_", project, ".csv")
full_out_path <- file.path(output_dir, file_out_name)

# 5. Save individual file
write.csv(combined, full_out_path, row.names = FALSE)