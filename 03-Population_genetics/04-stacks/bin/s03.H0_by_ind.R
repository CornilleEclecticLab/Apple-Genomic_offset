rm(list=ls())

library(tmap)
library(raster)
library(RColorBrewer)
library(sf)      # For spatial data handling
library(dplyr)
library(purrr)
library(maps)
library(gstat)
library(ggplot2)
library(ggrepel)

# ==== Parameters ====
Project="dom_admixed_rm_stack_by_ind"
input_file <- paste0("./out/",Project,"/s01.stacks/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.p.sumstats_summary.tsv")        # replace with your file path
ref_file <- "./in/data_tree_dom_admixed_rm_V2.txt"
out_folder=paste0("./out/",Project,"/s02.fig_stacks/")
variant_file <- paste0(out_folder,"/variant_positions.txt")
all_positions_file <- paste0(out_folder,"/all_positions.txt")
dir.create(out_folder, recursive = TRUE)

# ==== Read the file ====
lines <- readLines(input_file)

# ==== Find the split line ====
split_index <- grep("^# All positions", lines)

# ==== Extract variant positions ====
variant_lines <- lines[1:(split_index-1)]
writeLines(variant_lines[-1], variant_file)

# ==== Extract all positions ====
all_positions_lines <- lines[split_index:length(lines)]
writeLines(all_positions_lines[-1], all_positions_file)

# ==== Optional: read them as data frames ====
variant_df <- read.delim(variant_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(variant_df)[1] <- "Seq_ID"

all_positions_df <- read.delim(all_positions_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(all_positions_df)[1] <- "Seq_ID"

## all_positions_df

info_ref_all_SYL <- read.table(ref_file, header = TRUE, sep = "\t", dec=',')
sorted_info = info_ref_all_SYL[order(info_ref_all_SYL$Order_vcf), ]

merge_position_ref <- merge(all_positions_df, sorted_info, by="Seq_ID")

###

data_sf <- st_as_sf(merge_position_ref, coords = c("Longitude", "Latitude"), crs = 4326)

# Load European boundary map
#"https://public.opendatasoft.com/explore/dataset/world-administrative-boundaries/export/?q=Europe&refine.continent=Europe&location=2,37.48899,0.00845&basemap=jawg.light&dataChart=eyJxdWVyaWVzIjpbeyJjb25maWciOnsiZGF0YXNldCI6IndvcmxkLWFkbWluaXN0cmF0aXZlLWJvdW5kYXJpZXMiLCJvcHRpb25zIjp7fX0sImNoYXJ0cyI6W3siYWxpZ25Nb250aCI6dHJ1ZSwidHlwZSI6ImNvbHVtbiIsImZ1bmMiOiJDT1VOVCIsInNjaWVudGlmaWNEaXNwbGF5Ijp0cnVlLCJjb2xvciI6IiNGRjUxNUEifV0sInhBeGlzIjoic3RhdHVzIiwibWF4cG9pbnRzIjo1MCwic29ydCI6IiJ9XSwidGltZXNjYWxlIjoiIiwiZGlzcGxheUxlZ2VuZCI6dHJ1ZSwiYWxpZ25Nb250aCI6dHJ1ZX0%3D"
europe_boundaries <- st_read("./in/world-administrative-boundaries_EUR.geojson")
print(names(europe_boundaries))
countries_to_keep <- c("Austria", "Belgium", "Bulgaria","Croatia", "Cyprus", "Czech Republic",
                       "Denmark", "Estonia", "France", "Germany", "Greece","Hungary",
                       "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands",
                       "Poland", "Romania", "Slovakia", "Slovenia", "Switzerland","Albania", 
                       "Andorra", "Belarus", "Bosnia and Herzegovina", "Georgia", "Kosovo",
                       "Liechtenstein", "Macedonia", "Moldova", "Monaco", "Montenegro", 
                       "San Marino","Ukraine", "Vatican", "Spain", "Moldova, Republic of") # List of countries you want to keep
europe_boundaries <- europe_boundaries %>%
  filter(name %in% countries_to_keep)

# Ensure CRS of europe_boundaries matches that of data_sf (4326)
europe_boundaries <- st_transform(europe_boundaries, crs = 4326)

bbox <- st_bbox(europe_boundaries)

# Create a regular grid within the bounding box
grd <- expand.grid(
  X = seq(bbox["xmin"], bbox["xmax"], length.out = 500),
  Y = seq(bbox["ymin"], bbox["ymax"], length.out = 500)
)
coordinates(grd) <- c("X", "Y")
gridded(grd) <- TRUE
fullgrid(grd) <- TRUE
crs(grd) <- st_crs(europe_boundaries)$proj4string

# Perform interpolation using IDW
idw_model <- gstat(id = "H0", formula = Obs_Het ~ 1, data = data_sf)
P.idw <- predict(idw_model, newdata = grd)

# Convert to raster object
r <- raster(P.idw)
crs(r) <- st_crs(europe_boundaries)$proj4string

# Clip raster to the European boundaries
r.m <- mask(r, st_as_sf(europe_boundaries))

# ============================
# Save as PNG
# ============================
data_sf_jitter <- st_jitter(data_sf, amount = 0.1)

png(paste0(out_folder, "Heterozigoty_obs_across_europe_", Project, "_site_grouped_level.png"),
    width = 10, height = 10, units = "in", res = 300)

tm_shape(r.m) +
  tm_raster(
    n = 12,
    palette = colorRampPalette(c("#1D60B4","#43CBF3", "#ADDA90","#F4E54E","#FFA43C", "#FF462B"))(100),
    title = "Heterozygotie observed (H0)"
  ) +
  tm_shape(data_sf_jitter) +
  tm_dots(col = "#37323E", size = 0.2) +
  tm_layout(legend.outside = TRUE)

dev.off()

svg(paste0(out_folder, "Heterozigoty_obs_across_europe_", Project, "_site_grouped_level.svg"),
    width = 10, height = 10)

tm_shape(r.m) +
  tm_raster(
    n = 12,
    palette = colorRampPalette(c("#1D60B4","#43CBF3", "#ADDA90","#F4E54E","#FFA43C", "#FF462B"))(100),
    title = "Heterozygotie observed (H0)"
  ) +
  tm_shape(data_sf_jitter) +
  tm_dots(col = "#37323E", size = 0.2) +
  tm_layout(legend.outside = TRUE)

dev.off()