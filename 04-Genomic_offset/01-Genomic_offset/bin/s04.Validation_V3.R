rm(list = ls())

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(broom.mixed)
library(lme4)
library(ggcorrplot)

# for map ploting
library(ggplot2)
library(sf)
library(raster)
library(gstat)
library(rnaturalearth)
library(rnaturalearthdata)


# --- 0. Setup & Config ---
PROJECT <- "dom_admixed_rm"
go_folder <- paste0("./out/", PROJECT, "/s03.GO/")
val_folder <- paste0("./out/", PROJECT, "/s04.Validation_V3/") # New version folder

traits_base <- c("Avg_Chl", "Avg_Flav", "Avg_NBI", "Bud_Break")
offsets <- c("offset", "offset_Relstad")

offset_files <- c(
  paste0(go_folder, "file/All_chr_offset_all_SNP.csv"),
  paste0(go_folder, "file/All_chr_offset_GEA_SNP_Bonferroni.csv"),
  paste0(go_folder, "file/All_chr_offset_GEA_SNP_FDR.csv"),
  paste0(go_folder, "file/All_chr_offset_GEA_SNP_pvalue_filtered_Bonferroni.csv"),
  paste0(go_folder, "file/All_chr_offset_GEA_SNP_pvalue_filtered_FDR.csv")
)

# --- 1. Load Data ---
pheno_raw <- read.table("./in/Data_all_orchard_pheno_25022026.tsv", header=TRUE, sep="\t", dec=".")
ind.info <- read.table("./in/data_tree_dom_admixed_rm_V2.txt", header=TRUE, sep="\t", dec=",")

# --- 2. Main Loop ---
for (file_path in offset_files) {
  name <- sub(".csv", "", basename(file_path))
  message(">>> Processing: ", name)
  
  # Directory Setup
  current_out_dir <- paste0(val_folder, name, "/")
  dir.create(paste0(current_out_dir, "Models/"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(current_out_dir, "Plots/"), recursive = TRUE, showWarnings = FALSE)

  # Load and Merge Offset
  go_data <- read.csv(file_path)
  ind.go <- cbind(ind.info, go_data)
  
  # Extrapolation Logic: Calculate site means for missing offsets
  site_means <- ind.go %>%
    group_by(Site) %>%
    summarise(across(all_of(offsets), \(x) mean(x, na.rm = TRUE), .names = "mean_{.col}"))

  # Final Integrated Dataset for Saclay
  df_saclay <- pheno_raw %>%
    filter(Location == "Saclay") %>%
    rename(individual_ID = Secondary_ID) %>%
    filter(!grepl("Msyl-Ref", individual_ID)) %>%
    left_join(ind.go, by = "individual_ID") %>%
    left_join(site_means, by = c("Site_Name" = "Site")) %>%
    mutate(
      offset = coalesce(offset, mean_offset),
      offset_Relstad = coalesce(offset_Relstad, mean_offset_Relstad),
      across(contains("Date"), as.Date)
    )

  # --- 3. Growth Rate Calculations (AGR/RGR) ---
  df_saclay <- df_saclay %>%
    arrange(individual_ID, Year) %>%
    group_by(individual_ID) %>%
    mutate(
    # Summer-to-Summer Growth (Long term)
    days_ss = as.numeric(Summer_Trunk_Date - lag(Summer_Trunk_Date)),
    AGR_Annual_ss = (Summer_Trunk_Measurement - lag(Summer_Trunk_Measurement)) / (days_ss / 365),
    RGR_Annual_ss = (log(Summer_Trunk_Measurement) - log(lag(Summer_Trunk_Measurement))) / (days_ss / 365),

    # Winter-to-Winter Growth (Long term)
    days_ww = as.numeric(Winter_Trunk_Date - lag(Winter_Trunk_Date)),
    AGR_Annual_ww = (Winter_Trunk_Measurement - lag(Winter_Trunk_Measurement)) / (days_ww / 365),
    RGR_Annual_ww = (log(Winter_Trunk_Measurement) - log(lag(Winter_Trunk_Measurement))) / (days_ww / 365),

    # Seasonal Growth (Winter to Summer)
    days_ws = as.numeric(Summer_Trunk_Date - Winter_Trunk_Date),
    AGR_Seasonal = (Summer_Trunk_Measurement - Winter_Trunk_Measurement) / (days_ws / 30), # Growth per month
    RGR_Seasonal = (log(Summer_Trunk_Measurement) - log(Winter_Trunk_Measurement)) / (days_ws / 30)
  ) %>%
    ungroup()

  # Final Trait List
  all_traits <- c(traits_base, "AGR_Annual_ss", "RGR_Annual_ss", "AGR_Annual_ww", "RGR_Annual_ww","AGR_Seasonal", "RGR_Seasonal")

  # --- 4. The Mixed-Model Loop (The Selection Test) ---
  model_results <- list()
  
  for (trait in all_traits) {
    # Clean data for model
    df_mod <- df_saclay %>% filter(!is.na(.data[[trait]]), !is.na(offset), is.finite(.data[[trait]]))
    
    if(nrow(df_mod) < 20) next # Skip if not enough data
    
    # Logic: Log transform everything except Bud_Break and Rates
    if(trait %in% c("Bud_Break", "AGR_Annual_ss", "RGR_Annual_ss", "AGR_Annual_ww", "RGR_Annual_ww","AGR_Seasonal", "RGR_Seasonal")) {
      formula_str <- paste0(trait, " ~ I(offset^2) + (1 | Year) + (1 | individual_ID)")
    } else {
      formula_str <- paste0("log(", trait, ") ~ I(offset^2) + (1 | Year) + (1 | individual_ID)")
    }
    
    try({
      m <- lmer(as.formula(formula_str), data = df_mod)
      res <- tidy(m, conf.int = TRUE) %>%
        filter(term == "I(offset^2)") %>%
        mutate(trait = trait, n_obs = nrow(df_mod))
      model_results[[trait]] <- res
    }, silent = TRUE)
  }

  summary_df <- bind_rows(model_results)
  write.csv(summary_df, paste0(current_out_dir, "Models/Selection_Gradients_Summary.csv"), row.names = FALSE)

  # --- 5. Visualizations ---
  # A. Forest Plot of Selection Gradients
  p_forest <- ggplot(summary_df, aes(x = reorder(trait, estimate), y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
    coord_flip() + theme_minimal() +
    labs(title = "Selection Gradients Summary", subtitle = name, x = "Trait", y = "Estimate (Beta)")
  
  ggsave(paste0(current_out_dir, "Plots/Forest_Plot_Gradients.pdf"), p_forest, width = 8, height = 6)

  # B. Correlation Matrix (using individual means)
  mean_ind <- df_saclay %>%
    group_by(individual_ID) %>%
    summarise(across(all_of(c(all_traits, offsets)), \(x) mean(x, na.rm = TRUE)))
  
  M <- cor(mean_ind %>% dplyr::select(-individual_ID), use = "pairwise.complete.obs")
  p_corr <- ggcorrplot(M[offsets, all_traits], lab = TRUE, colors = c("#E46726", "white", "#6D9EC1"))
  ggsave(paste0(current_out_dir, "Plots/Correlation_Matrix.pdf"), p_corr, width = 10, height = 5)


  ## representation for bud breaks
  library(ggplot2)
  library(viridis)

  # 1. Filter and clean data specifically for this plot
  plot_df <- df_saclay %>%
    filter(!is.na(Bud_Break), !is.na(offset))

  # 2. Create the plot
  p_budbreak <- ggplot(plot_df, aes(x = offset, y = Bud_Break)) +
    # Add points colored by Year to show the 2 years of data
    geom_point(aes(color = as.factor(Year)), alpha = 0.5, size = 2.5) +
    
    # Add the quadratic trend line (matches the I(offset^2) model)
    # We use formula y ~ I(x^2) to show the effect of the squared offset
    geom_smooth(method = "lm", formula = y ~ I(x^2), 
                color = "black", size = 1.5, se = TRUE) +
    
    # Formatting
    theme_minimal(base_size = 14) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    labs(
      title = "Impact of Genomic Offset on Phenology",
      subtitle = "Maladapted individuals (high offset) exhibit significantly earlier Bud Break",
      x = "Genomic Offset (Distance to Optimum)",
      y = "Bud Break (Julian Days)",
      color = "Year"
    ) +
    
    # Add an annotation with your model stats
    annotate("label", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
            label = "Effect size (ß): -6.88\nt-value: -6.33 (p < 0.001)",
            fill = "white", alpha = 0.8)

  # 3. Save the plot
  ggsave("Bud_Break_Offset_Relationship.pdf", p_budbreak, width = 8, height = 6)

  # Display the plot
  print(p_budbreak)



  ## GO representation

  # library(ggplot2)
  # library(sf)
  # library(rnaturalearth)
  # library(rnaturalearthdata)
  # library(ggspatial)
  # library(dplyr)

  # # 1. Prepare Spatial Data (Average offset per source site)
  map_data <- df_saclay %>%
    group_by(Site, Longitude, Latitude, Country) %>%
    summarise(
      mean_offset = mean(offset, na.rm = TRUE),
      n_trees = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(Longitude), !is.na(Latitude))

  # # 2. Get a base map of Europe (assuming your trees are from Europe)
  # world <- ne_countries(scale = "medium", returnclass = "sf")

  # # 3. Create the Map
  # p_map <- ggplot(data = world) +
  #   # Draw the background map
  #   geom_sf(fill = "#f9f9f9", color = "gray80") +
    
  #   # Zoom into Europe (adjust coordinates to fit your specific sites)
  #   coord_sf(xlim = c(-15, 30), ylim = c(35, 65), expand = FALSE) +
    
  #   # Plot the sites. Color = Offset Value, Size = Number of trees
  #   geom_point(data = map_data, 
  #             aes(x = Longitude, y = Latitude, color = mean_offset, size = n_trees), 
  #             alpha = 0.8) +
    
  #   # Add a star for Saclay (the common garden/planting site)
  #   # Replace with actual Saclay coordinates
  #   annotate("text", x = 2.16, y = 48.72, label = "★ Saclay", 
  #           color = "red", fontface = "bold", vjust = -1) +
    
  #   # Aesthetics
  #   scale_color_viridis_c(option = "magma", direction = -1, name = "Genomic Offset") +
  #   theme_minimal() +
  #   labs(
  #     title = "Geographic Origin of Genomic Maladaptation",
  #     subtitle = "Higher offset indicates greater genetic distance from the Saclay environment",
  #     size = "Sample Size (n)"
  #   ) +
  #   annotation_scale(location = "bl", width_hint = 0.2) +
  #   annotation_north_arrow(location = "bl", which_north = "true", 
  #                         pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
  #                         style = north_arrow_fancy_orienteering)

  # # 4. Save
  # ggsave("Map_Genomic_Offset.pdf", p_map, width = 10, height = 8)

  # # Display
  # print(p_map)

# 1. Prepare Point Data
# Use the individual means you already calculated
  pts <- df_saclay %>%
    group_by(Site, Longitude, Latitude, Country) %>%
    summarise(
      mean_offset = mean(offset, na.rm = TRUE),
      n_trees = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(Longitude), !is.na(Latitude))

  # Convert to a spatial object (sf)
  pts_sf <- st_as_sf(pts, coords = c("Longitude", "Latitude"), crs = 4326)

  # 2. Create an empty Grid for the Europe area
  # Define the extent (adjust based on your data)
  grid_extent <- extent(-15, 40, 35, 65)
  r <- raster(grid_extent, res = 0.2) # 0.2 degree resolution for smoothness
  crs(r) <- crs(pts_sf)

  # 3. Perform Interpolation (IDW - Inverse Distance Weighting)
  # This calculates the value for every pixel based on nearby points
  idw_mod <- gstat(formula = mean_offset ~ 1, data = as(pts_sf, "Spatial"))
  idw_raster <- interpolate(r, idw_mod)

  # 4. Mask the Raster to the Europe landmass
  # This prevents the "heat" from showing up over the ocean
  # europe_mask <- ne_countries(scale = "medium", returnclass = "sp", continent = "Europe")
  # idw_masked <- mask(idw_raster, europe_mask)

  europe_sf <- ne_countries(scale = "medium", returnclass = "sf", continent = "Europe")

  # 2. Convert sf to Spatial object for the 'raster' package functions
  # This solves the 'arg should be one of sf, sv' error
  europe_sp <- as_Spatial(st_geometry(europe_sf))

  # 3. Perform the mask using the converted object
  # idw_raster is the result of your interpolate() function
  idw_masked <- mask(idw_raster, europe_sp)

  # 4. Optional: Crop to the specific extent of Europe to remove empty space
  idw_final <- crop(idw_masked, extent(europe_sp))

  # 5. Convert to Dataframe for ggplot
  idw_df <- as.data.frame(idw_masked, xy = TRUE)
  colnames(idw_df) <- c("long", "lat", "offset_val")

  # 6. Plot the Extrapolation Map
  p_extrap <- ggplot() +
    # Draw the background map
    geom_sf(data = st_as_sf(europe_sf), fill = "white", color = "gray90") +
    
    # Add the interpolated heat map layer
    geom_raster(data = idw_df %>% filter(!is.na(offset_val)), 
                aes(x = long, y = lat, fill = offset_val), interpolate = TRUE) +
    
    # Overlay the original sampling points (optional)
    # geom_point(data = pts, aes(x = Longitude, y = Latitude), 
    #           color = "black", size = 0.5, alpha = 0.3) +
    geom_point(data = pts, 
              aes(x = Longitude, y = Latitude, color = mean_offset, size = n_trees), 
              alpha = 0.8) +
    
    # Add a star for Saclay (the common garden/planting site)
    # Replace with actual Saclay coordinates
    annotate("text", x = 2.16, y = 48.72, label = "★ Saclay", 
            color = "red", fontface = "bold", vjust = -1) +
    
    # Match the color scale to your image (Blue -> Yellow -> Red)
    scale_fill_gradientn(colors = c("#9ecae1", "#ffeda0", "#feb24c", "#f03b20"), 
                        name = "Offset", na.value = "transparent") +
    
    coord_sf(xlim = c(-15, 40), ylim = c(35, 65), expand = FALSE) +
    theme_minimal() +
    labs(title = "Spatial Extrapolation of Genomic Offset",
        x = "Longitude", y = "Latitude")

  # Save
  ggsave(paste0(current_out_dir, "Plots/Extrapolation_Map_Offset.png"), p_extrap, width = 10, height = 7, dpi = 300)
}