rm(list = ls())

# install.packages("purrr")  # For functional programming and easy list handling
# install.packages(c("ggplot2", "akima", "sf", "ggspatial"))
# install.packages("maps")
# install.packages("tmap")
library(tmap)
library(raster)
library(RColorBrewer)
library(sf)      # For spatial data handling
library(dplyr)
library(purrr)
library(maps)
library(gstat)

#Project = "dom_admixed_rm" #"dom_admixed_rm_pi_by_ind"
# Project = "full_genome_dom_admixed_rm_no_filter" ##variant/invariant dataset
Project = "full_genome_dom_admixed_rm_no_filter_pi_by_ind"
out_folder = paste0("./out/",Project,"/s01.Metrics/")
dir.create(out_folder, recursive = T)

#pi_folder = "output/avg_pi_EC_syl.Chr.light.reformat/" #site level
#pi_folder = "output/avg_pi_Ind_level_EC_syl.light.reformat/" #ind level
# pi_folder = "output/avg_pi_Site_grouped_level_EC_syl.light.reformat/"
pi_folder = paste0("../02-Pixy/out/",Project,"/s02.avg_value/")
mean_coordinate = read.delim(paste0("./in/group_coordinate_",Project,".txt"), sep = "\t", dec = ",")
#mean_coordinate = read.delim("./in/group_coordinate_dom_admixed_rm_pi_by_ind.txt", sep = "\t", dec = ",")
summary(mean_coordinate)
files = list.files(pi_folder, full.names = T, pattern = "pi")

# Step 2: Read all CSV files and combine them into one dataset
combined_data <- files %>%
  map_df(~ read.csv(.))  # Use purrr::map_df to read and combine all files

# Step 3: Inspect the combined data
head(combined_data)

# Step 4: Group by Site and calculate the mean avg_pi for each site
group_means <- combined_data %>%
  group_by(pop) %>%
  summarise(mean_avg_pi = mean(avg_pi, na.rm = TRUE)) %>%
  rename(Group = pop)

data_site <- merge(group_means, mean_coordinate, by = "Group")

# Assuming 'data_site' is your data frame with site data
data_sf <- st_as_sf(data_site, coords = c("Longitude", "Latitude"), crs = 4326)

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

# Define a gradient color palette (using RColorBrewer)
palette_gradient <- colorRampPalette(c("#1D60B4","#43CBF3", "#ADDA90","#F4E54E","#FFA43C", "#FF462B"))(100)

# ============================
# Plot using tmap (site-level dots)
# ============================
tm_shape(europe_boundaries) + 
  tm_polygons() +
  tm_shape(data_sf) +
  tm_dots(
    col = "mean_avg_pi", 
    palette = colorRampPalette(c("#1D60B4","#43CBF3", "#ADDA90","#F4E54E","#FFA43C", "#FF462B"))(100),  # FIX: use palette= function
    title = "Sampled precipitation \n(in inches)", 
    size = 0.7
  ) +
  tm_text("Group", just = "left", xmod = 0.5, size = 0.7) +
  tm_layout(legend.outside = TRUE)

# Define the bounding box of the European boundaries
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
idw_model <- gstat(id = "avg_pi", formula = mean_avg_pi ~ 1, data = data_sf)
P.idw <- predict(idw_model, newdata = grd)

# Convert to raster object
r <- raster(P.idw)
crs(r) <- st_crs(europe_boundaries)$proj4string

# Clip raster to the European boundaries
r.m <- mask(r, st_as_sf(europe_boundaries))

# ============================
# Save as PNG
# ============================
png(paste0(out_folder, "Nucleotide_diversity_across_europe_", Project, "_site_grouped_level.png"),
    width = 10, height = 10, units = "in", res = 300)

tm_shape(r.m) +
  tm_raster(
    n = 12, 
    palette = colorRampPalette(c("#1D60B4","#43CBF3", "#ADDA90","#F4E54E","#FFA43C", "#FF462B"))(100),  # FIX: use palette= function
    title = "Nucleotide diversity (Pi)"
  ) +
  tm_shape(data_sf) +
  tm_dots(col = "#37323E", size = 0.2) +
  #tm_text("Group", just = "left", ymod = 1, size = 0.7) +
  tm_layout(legend.outside = TRUE)

dev.off()

# ============================
# Save as SVG
# ============================
svg(paste0(out_folder, "Nucleotide_diversity_across_europe_", Project, "_site_grouped_level.svg"),
    width = 10, height = 10)

tm_shape(r.m) +
  tm_raster(
    n = 12, 
    palette = colorRampPalette(c("#1D60B4","#43CBF3", "#ADDA90","#F4E54E","#FFA43C", "#FF462B"))(100),  # FIX: use palette= function
    title = "Nucleotide diversity (Pi)"
  ) +
  tm_shape(data_sf) +
  tm_dots(col = "#37323E", size = 0.2) +
  #tm_text("Group", just = "left", ymod = 1, size = 0.7) +
  tm_layout(legend.outside = TRUE)

dev.off()



