# ==============================================================================
# PROJECT: Bio19 Climatic Extrapolation Maps
# OBJECTIVE: Map Bio19 values for Current + 9 Future timelines with sample points
# ==============================================================================

# ---- 1. Clean environment & Load Libraries ----
rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(sf)
  library(raster)
  library(gstat)
  library(rnaturalearth)
  library(patchwork)
  library(scales)
  library(purrr)
  library(tidyr)
})

# ---- 2. Paths & Config ----
PROJECT   <- "dom_admixed_rm"
Bio_col   <- 8
env_dir   <- file.path("./out", PROJECT, "s09.Bio_Map/env")
fig_dir   <- file.path("./out", PROJECT, "s09.Bio_Map/figure")
ref_info  <- "./in/data_tree_dom_admixed_rm_V2.txt"

dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 3. Load & Prepare Data ----
info_ref <- read.table(ref_info, header = TRUE, sep = "\t")

# A. Reference for interpolation (matches .env rows)
sorted_info <- info_ref[order(info_ref$Order_vcf), ] %>%
  mutate(
    Longitude = as.numeric(gsub(",", ".", Longitude)),
    Latitude  = as.numeric(gsub(",", ".", Latitude))
  )

# B. Dataframe for the point layer (geom_point needs raw columns, not sf)
pts <- sorted_info %>%
  filter(!is.na(Longitude), !is.na(Latitude)) %>%
  group_by(Site, Longitude, Latitude) %>%
  summarise(n_trees = n(), .groups = "drop")

# ---- 4. Spatial Setup ----
# Define extent for Europe
grid_extent <- extent(-15, 40, 35, 65)
r <- raster(grid_extent, res = 0.4) # Resolution: 0.1 for high quality, 0.4 for speed
europe_sf <- ne_countries(scale = "medium", returnclass = "sf", continent = "Europe")
europe_sp <- as_Spatial(st_geometry(europe_sf))
projection(r) <- crs(europe_sf)

# ---- 5. Mapping Function ----

generate_bio_map <- function(env_path, pts_df, bio_idx, global_range, sorted_info, raster_template, mask_poly) {
  
  # Load .env file
  env_data <- read.table(env_path, header = FALSE) 
  
  if (ncol(env_data) < bio_idx) {
    stop(paste("File", basename(env_path), "has only", ncol(env_data), "cols."))
  }
  
  # Prepare interpolation data
  interp_df <- sorted_info %>%
    mutate(val = env_data[[bio_idx]]) %>%
    filter(!is.na(Longitude), !is.na(Latitude), !is.na(val)) %>%
    group_by(Longitude, Latitude) %>%
    summarise(mean_val = mean(val, na.rm = TRUE), .groups = "drop")
  
  # Convert to spatial for gstat
  interp_sf <- st_as_sf(interp_df, coords = c("Longitude", "Latitude"), crs = 4326)
  
  # Title Logic
  fname <- basename(env_path)
  title_lab <- if(grepl("future", fname)) {
    paste(str_extract(fname, "\\d{4}-\\d{4}"), str_extract(fname, "ssp\\d+"))
  } else {
    "Current"
  }

  # Interpolation (Inverse Distance Weighting)
  idw_mod    <- gstat(formula = mean_val ~ 1, data = as(interp_sf, "Spatial"))
  idw_raster <- interpolate(raster_template, idw_mod)
  idw_masked <- mask(idw_raster, mask_poly)
  
  # Convert to DF for ggplot
  idw_df <- as.data.frame(idw_masked, xy = TRUE)
  colnames(idw_df) <- c("long", "lat", "bio_val")

  # Plotting
  p <- ggplot() +
    geom_sf(data = europe_sf, fill = "gray98", color = "gray85", size = 0.1) +
    geom_raster(data = idw_df %>% filter(!is.na(bio_val)), 
                aes(x = long, y = lat, fill = bio_val), interpolate = TRUE) +
    # Use 'pts_df' which is the raw dataframe with Longitude/Latitude columns
    geom_point(data = pts_df, 
               aes(x = Longitude, y = Latitude, size = n_trees), 
               shape = 21, color = "black", fill = NA, stroke = 0.6, alpha = 0.5) +
    scale_fill_gradientn(
      colors = c("#4575b4", "#abd9e9", "#fee090", "#f46d43", "#d73027"),
      name = paste0("Bio", bio_idx),
      limits = global_range,
      oob = scales::squish
    ) +
    scale_size_continuous(range = c(0.5, 4), guide = "none") +
    coord_sf(xlim = c(-12, 35), ylim = c(36, 64), expand = FALSE) +
    theme_minimal(base_size = 8) +
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank()) +
    labs(title = title_lab)

  return(p)
}

# ---- 6. Execution ----
env_files <- list.files(env_dir, pattern = "\\.env$", full.names = TRUE)

# Calculate Global Range
message("Calculating global range...")
ranges_list <- map(env_files, ~{
  d <- read.table(.x, header = FALSE)
  range(d[[Bio_col]], na.rm = TRUE)
})
global_range <- range(unlist(ranges_list), na.rm = TRUE)

# Generate all plots into a named list
message("Generating maps...")
map_list <- map(env_files, ~{
  p <- generate_bio_map(.x, pts, Bio_col, global_range, sorted_info, r, europe_sp)
  # Name the plot by the filename to help sorting
  attr(p, "filename") <- basename(.x)
  return(p)
})
names(map_list) <- map_chr(map_list, ~attr(.x, "filename"))

# ---- 7. Separate and Sort Plots ----

# 1. Identify the "Current" plot
current_key <- names(map_list)[grepl("current", names(map_list))]
p_current <- map_list[[current_key]] + 
  labs(subtitle = "Baseline Reference") +
  theme(plot.title = element_text(size = 14, face = "bold"))

# 2. Identify and Sort "Future" plots
future_keys <- names(map_list)[grepl("future", names(map_list))]
future_plots <- map_list[future_keys]

# Sort logic: 
# We want rows = Timeline (2041, 2061, 2081) and columns = SSP (ssp126, ssp245, ssp585)
sort_df <- data.frame(
  key = future_keys,
  timeline = str_extract(future_keys, "\\d{4}"),
  ssp = str_extract(future_keys, "ssp\\d+")
) %>%
  arrange(timeline, ssp) # This ensures Row 1 is 2041-2060, Col 1 is ssp126, etc.

ordered_future_plots <- future_plots[sort_df$key]

# ---- 8. Create the Complex Layout ----

# Combine the 9 future plots into a 3x3 grid
future_grid <- wrap_plots(ordered_future_plots, ncol = 3)

# Define the final layout: 
# 'A' is the Current map (taking up more height/width)
# 'B' is the 3x3 grid
# We use the design layout where 'A' is one big block and 'B' is the rest
final_plot <- (p_current | future_grid) + 
  plot_layout(widths = c(1, 2.5)) + 
  plot_annotation(
    title = paste("Bio", Bio_col, "Evolution: Current vs. Future Scenarios"),
    subtitle = "Left: Current Baseline | Right Grid: Rows (Time) x Columns (SSP)",
    theme = theme(plot.title = element_text(size = 20, face = "bold"))
  ) +
  plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom", 
        legend.key.width = unit(3, "cm"), # Make legend longer to see colors clearly
        legend.title = element_text(face = "bold"))

# ---- 9. Save ----
output_file <- file.path(fig_dir, paste0("Bio", Bio_col, "_Custom_Layout.png"))
ggsave(output_file, final_plot, width = 22, height = 12, dpi = 300)

message("✅ Success! Custom layout saved to: ", output_file)