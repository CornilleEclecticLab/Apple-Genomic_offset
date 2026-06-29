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
library(svglite)


# --- 0. Setup & Config ---
PROJECT <- "dom_admixed_rm"
# PROJECT <- "all_wild"
go_folder <- paste0("./out/", PROJECT, "/s03.GO_V3/")
val_folder <- paste0("./out/", PROJECT, "/s04.5.Annual_validation/") # New version folder

traits_base <- c("Avg_Chl", "Avg_Flav", "Avg_NBI", "Bud_Break")
offsets <- c("offset", "offset_Relstad")

orchard <- "Saclay"

annee_pattern <- "2021|2022|2023|2024|2025"

# Lister uniquement les fichiers qui correspondent à ce motif
offset_files <- list.files(
  path = paste0(go_folder, "file/"), 
  pattern = annee_pattern, 
  full.names = TRUE # Mettez TRUE si vous voulez le chemin complet
)

# offset_files <- c(
# paste0(go_folder, "file/GO_Filt_env_",orchard,"_Chelsa_All.csv"),
# paste0(go_folder, "file/GO_Filt_env_",orchard,"_Chelsa_Bonf.csv"),
# paste0(go_folder, "file/GO_Filt_env_",orchard,"_Chelsa_FDR.csv"),
# paste0(go_folder, "file/GO_Full_env_",orchard,"_Chelsa_All.csv"),
# paste0(go_folder, "file/GO_Full_env_",orchard,"_Chelsa_Bonf.csv"),
# paste0(go_folder, "file/GO_Full_env_",orchard,"_Chelsa_FDR.csv")
# )

# --- 1. Load Data ---
pheno_raw <- read.table("./in/Data_all_orchard_pheno_25022026.tsv", header=TRUE, sep="\t", dec=".")
ind.info <- read.table("./in/data_tree_dom_admixed_rm_V2.txt", header=TRUE, sep="\t", dec=",")
# ind.info <- read.table("./in/data_tree_all_wild.txt", header=TRUE, sep="\t", dec=",")

# --- 2. Main Loop ---
for (file_path in offset_files) {
  name <- sub(".csv", "", basename(file_path))
  year <- str_extract(name, "\\d{4}")

  message(">>> Processing: ", name)
  
  # Directory Setup
  current_out_dir <- paste0(val_folder, name, "/")
  dir.create(paste0(current_out_dir, "Models/"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(current_out_dir, "Plots/"), recursive = TRUE, showWarnings = FALSE)

  # Load and Merge Offset
  go_data <- read.csv(file_path)[c("squared.env.dist", "offset", "offset_Relstad")]
  ind.go <- cbind(ind.info, go_data)
  
  # Extrapolation Logic: Calculate site means for missing offsets
  site_means <- ind.go %>%
    group_by(Site) %>%
    summarise(across(all_of(offsets), \(x) mean(x, na.rm = TRUE), .names = "mean_{.col}"))

  df_saclay <- pheno_raw %>%
    filter(Location == orchard) %>%
    filter(Year == year) %>%
    rename(individual_ID = Secondary_ID) %>%
    filter(!grepl("Msyl-Ref", individual_ID)) %>%
    # --- NEW FILTER START ---
    # Remove rows where Month.Year_Death has any value (is not empty)
    filter(is.na(Month.Year_Death) | Month.Year_Death == "") %>%
    # --- NEW FILTER END ---
    left_join(ind.go, by = "individual_ID") %>%
    left_join(site_means, by = c("Site_Name" = "Site")) %>%
    mutate(
      offset = coalesce(offset, mean_offset),
      offset_Relstad = coalesce(offset_Relstad, mean_offset_Relstad),
      across(contains("Date"), as.Date)
    )

  df_saclay <- df_saclay %>% 
  arrange(individual_ID, Year) %>% 
  group_by(individual_ID) %>% 
  mutate( 
    # # --- 1. Summer-to-Summer Growth (Yearly) ---
    # # Verification: Current Summer must be > Previous Summer
    # days_ss = as.numeric(Summer_Trunk_Date - lag(Summer_Trunk_Date)), 
    # AGR_Annual_ss = if_else(Summer_Trunk_Measurement > lag(Summer_Trunk_Measurement) & days_ss > 0,
    #                         (Summer_Trunk_Measurement - lag(Summer_Trunk_Measurement)) / (days_ss / 365), 
    #                         NA_real_), 
    # RGR_Annual_ss = if_else(!is.na(AGR_Annual_ss),
    #                         (log(Summer_Trunk_Measurement) - log(lag(Summer_Trunk_Measurement))) / (days_ss / 365), 
    #                         NA_real_), 

    # # --- 2. Winter-to-Winter Growth (Yearly) ---
    # # Verification: Current Winter must be > Previous Winter
    # days_ww = as.numeric(Winter_Trunk_Date - lag(Winter_Trunk_Date)), 
    # AGR_Annual_ww = if_else(Winter_Trunk_Measurement > lag(Winter_Trunk_Measurement) & days_ww > 0,
    #                         (Winter_Trunk_Measurement - lag(Winter_Trunk_Measurement)) / (days_ww / 365), 
    #                         NA_real_), 
    # RGR_Annual_ww = if_else(!is.na(AGR_Annual_ww),
    #                         (log(Winter_Trunk_Measurement) - log(lag(Winter_Trunk_Measurement))) / (days_ww / 365), 
    #                         NA_real_), 

    # --- 3. Seasonal Growth (Within Year: Summer to Winter) ---
    # Verification: Winter measurement MUST be higher than Summer measurement
    days_sw = as.numeric(Winter_Trunk_Date - Summer_Trunk_Date), 
    AGR_Seasonal = if_else(Winter_Trunk_Measurement > Summer_Trunk_Measurement & days_sw > 0,
                           (Winter_Trunk_Measurement - Summer_Trunk_Measurement) / (days_sw / 30), 
                           NA_real_), 
    RGR_Seasonal = if_else(!is.na(AGR_Seasonal),
                           (log(Winter_Trunk_Measurement) - log(Summer_Trunk_Measurement)) / (days_sw / 30), 
                           NA_real_)
  ) %>% 
  ungroup()

  # Final Trait List
  all_traits <- c(traits_base, "AGR_Seasonal", "RGR_Seasonal")

  library(performance)
  library(broom.mixed)

  model_results_list <- list()

  for (trait in all_traits) {
    
    # 1. Data Cleaning (As per your previous logic)
    df_mod <- df_saclay %>% 
      filter(!is.na(.data[[trait]]), !is.na(offset), !is.na(Site_Name), !is.na(Origin)) %>%
      group_by(Year) %>%
      mutate(W_rel = .data[[trait]] / mean(.data[[trait]], na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(log_W = log(W_rel)) %>%
      filter(is.finite(log_W))
    
    if(nrow(df_mod) < 30) next # Increased slightly as we add more parameters

    # 2. Determine Dependent Variable
    dep_var <- if(trait %in% c("Bud_Break", "AGR_Seasonal")) "log_W" else "W_rel"

    # 3. Define the three model formulas
    formulas <- list(
      Base   = paste0(dep_var, " ~ offset + I(offset^2)"),
      w_Site = paste0(dep_var, " ~ offset + I(offset^2) + Site_Name"),
      Full   = paste0(dep_var, " ~ offset + I(offset^2) + Site_Name + Origin")
    )

    # 4. Run and compare
    for (mod_name in names(formulas)) {
      try({
        m <- lm(as.formula(formulas[[mod_name]]), data = df_mod)
        
        # Extract Performance Metrics
        perf <- model_performance(m)
        
        # Extract specific "offset" coefficient (Linear term)
        res <- tidy(m, conf.int = TRUE) %>%
          filter(term == "offset") %>%
          mutate(
            trait = trait,
            model_type = mod_name,
            n_obs = nrow(df_mod),
            R2 = perf$R2,
            R2_adj = perf$R2_adjusted,
            # R2_marginal = perf$R2_marginal,
            # R2_conditional = perf$R2_conditional,
            AIC = perf$AIC
          )
        
        # Store using a unique key
        model_results_list[[paste(trait, mod_name, sep="_")]] <- res
      }, silent = FALSE)
    }

    plot_df <- df_mod

    info_model <- model_results_list[[paste(trait, "Base", sep="_" )]]
    print(info_model)

    p_budbreak_1 <- ggplot(plot_df, aes(x = offset, y = .data[[trait]])) +
        geom_point(aes(color = as.factor(Origin)), alpha = 0.5, size = 2.5) +
        geom_smooth(method = "lm", formula = y ~ x, color = "black", size = 1.5, se = TRUE) +
        theme_minimal(base_size = 14) +
        #scale_color_manual(values = c("#E69F00", "#56B4E9")) +
        labs(title = "Impact of genomic offset on phenology",
        x = "Genomic offset",
        y = trait, #"Bud Break (Julian Days)",
        color = "Origin") #+
        # annotate("label", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
        # label = paste0("estimate:",round(info_model$estimate, digits=3),
        # "\nt-value:",round(info_model$statistic, digits=3),
        # "\nR2: ",round(info_model$R2_marginal, digits=3)),
        # fill = "white", alpha = 0.5)
    ggsave(paste0(current_out_dir,"Raw_",trait,"_Offset_Relationship_Origin.png"), p_budbreak_1, width = 8, height = 6, dpi=300)
    ggsave(paste0(current_out_dir,"Raw_",trait,"_Offset_Relationship_Origin.svg"), p_budbreak_1, width = 8, height = 6)


    p_budbreak <- ggplot(plot_df, aes(x = offset, y = .data[[dep_var]])) +
      geom_point(aes(color = as.factor(Origin)), alpha = 0.5, size = 2.5) +
      # Note: Changed 'size' to 'linewidth' in geom_smooth to avoid modern ggplot2 warnings
      geom_smooth(method = "lm", formula = y ~ x, 
                  color = "black", linewidth = 1.5, se = TRUE) +
      theme_minimal(base_size = 14) +
      labs(
        title = "Impact of genomic offset on phenology",
        x = "Genomic offset",
        y = paste0(trait, " (Relative fitness)"), 
        color = "Origin"
      ) +
      # 1. Turn off clipping so the label can be drawn outside the plot panel
      coord_cartesian(clip = "off") +
      # 2. Add extra padding to the right margin to make room for the label
      # margin(top, right, bottom, left)
      theme(plot.margin = margin(10, 150, 10, 10)) + 
      # 3. Adjust hjust and vjust to push the label into the right margin
      annotate("label", 
              x = Inf, y = Inf, 
              hjust = -0.15, vjust = 1, 
              label = paste0("estimate: ", round(info_model$estimate, digits=3),
                              "\nt-value: ", round(info_model$statistic, digits=3),
                              "\nR2: ", round(info_model$R2, digits=3)),
              fill = "white", alpha = 0.5)
    # 3. Save the plot
    ggsave(paste0(current_out_dir,trait,"_Offset_Relationship_Origin.png"), p_budbreak, width = 8, height = 6, dpi=300)
    ggsave(paste0(current_out_dir,trait,"_Offset_Relationship_Origin.svg"), p_budbreak, width = 8, height = 6)
  }

  # Combine into one master summary
  summary_df <- bind_rows(model_results_list)
  write.csv(summary_df, paste0(current_out_dir, "Models/Selection_Gradients_Summary_Test.csv"), row.names = FALSE)
# View a comparison for one trait
# summary_df %>% filter(trait == "AGR_Seasonal") %>% select(model_type, R2_marginal, AIC)

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
    summarise(across(all_of(c(all_traits, offsets)), \(x) log(mean(x, na.rm = TRUE))))
  
  M <- cor(mean_ind %>% dplyr::select(-individual_ID), use = "pairwise.complete.obs")
  p_corr <- ggcorrplot(M[all_traits, offsets], lab = TRUE, colors = c("#E46726", "white", "#6D9EC1"))
  #ggsave(paste0(current_out_dir, "Plots/Correlation_Matrix.pdf"), p_corr, width = 10, height = 5)
  ggsave(paste0(current_out_dir,"Plots/test_log_Correlation_Matrix.png"), p_corr, width = 10, height = 5, dpi=300)
  ggsave(paste0(current_out_dir,"Plots/test_log_Correlation_Matrix.svg"), p_corr, width = 10, height = 5)


  # Calculate correlations per year
  corr_year <- df_saclay %>%
    group_by(Year) %>%
    summarise(across(all_of(all_traits), 
                    ~cor(.x, offset, use = "pairwise.complete.obs"))) %>%
    pivot_longer(-Year, names_to = "trait", values_to = "correlation")

  # Plot
  p_year <- ggplot(corr_year, aes(x = as.factor(Year), y = trait, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = round(correlation, 2))) +
    scale_fill_gradient2(low = "#E46726", mid = "white", high = "#6D9EC1", midpoint = 0) +
    labs(title = "Correlation: Traits vs Offset (per Year)", x = "Year", y = "Trait") +
    theme_minimal()

  ggsave(paste0(current_out_dir,"Correlation_per_Year.png"), p_year, width = 8, height = 6, dpi=300)


  # Calculate correlations per population
  corr_pop <- df_saclay[!df_saclay$Origin=="Spain",] %>%
    group_by(Origin) %>%
    summarise(across(all_of(all_traits), 
                    ~cor(.x, offset, use = "pairwise.complete.obs"))) %>%
    pivot_longer(-Origin, names_to = "trait", values_to = "correlation")

  # Plot
  p_pop <- ggplot(corr_pop, aes(x = Origin, y = trait, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = round(correlation, 2))) +
    scale_fill_gradient2(low = "#E46726", mid = "white", high = "#6D9EC1", midpoint = 0) +
    labs(title = "Correlation: Traits vs Offset (per Population)", x = "Population", y = "Trait") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(paste0(current_out_dir,"Correlation_per_Pop.png"), p_pop, width = 10, height = 6, dpi=300)


  corr_pop_rel <- df_saclay %>%
  filter(Origin != "Spain") %>%
  group_by(Origin, Year) %>% # Grouping by year ensures relative fitness is local
  mutate(across(all_of(all_traits), 
                ~ .x / mean(.x, na.rm = TRUE), 
                .names = "{.col}_rel")) %>%
  ungroup() %>%
  group_by(Origin) %>%
  summarise(across(ends_with("_rel"), 
                   ~cor(.x, offset, use = "pairwise.complete.obs"))) %>%
  pivot_longer(-Origin, names_to = "trait", values_to = "correlation")

  # Clean up trait names for the plot (remove the "_rel" suffix for beauty)
  corr_pop_rel$trait <- gsub("_rel", "", corr_pop_rel$trait)

  p_pop_rel <- ggplot(corr_pop_rel, aes(x = Origin, y = trait, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
  scale_fill_gradient2(low = "#E46726", mid = "white", high = "#6D9EC1", 
                       midpoint = 0, limit = c(-1, 1)) +
  labs(title = "Correlation: Relative Traits vs Offset",
       subtitle = "Relative Trait = Trait / Mean(Trait per Origin & Year)",
       x = "Population (Origin)", y = "Relative Trait") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Save the figure
  ggsave(paste0(current_out_dir,"Correlation_Relative_Traits_Pop.png"), p_pop_rel, width = 10, height = 7, dpi=300)
  ggsave(paste0(current_out_dir,"Correlation_Relative_Traits_Pop.svg"), p_pop_rel, width = 10, height = 7)

  # Calculate correlations per year AND population
  corr_year_pop <- df_saclay[!df_saclay$Origin=="Spain",] %>%
    group_by(Year, Origin) %>%
    summarise(across(all_of(all_traits), 
                    ~cor(.x, offset, use = "pairwise.complete.obs")),
              .groups = 'drop') %>%
    pivot_longer(cols = all_of(all_traits), names_to = "trait", values_to = "correlation")

  # Plot using Facets
  p_year_pop <- ggplot(corr_year_pop, aes(x = Origin, y = trait, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = round(correlation, 2))) +
    facet_wrap(~Year, ncol = 1) + # Each year gets its own matrix
    scale_fill_gradient2(low = "#E46726", mid = "white", high = "#6D9EC1", midpoint = 0) +
    labs(title = "Trait-Offset Correlation per Year and Population", 
        x = "Population", y = "Trait") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(paste0(current_out_dir,"Correlation_per_Year_Pop.png"), p_year_pop, width = 12, height = 10, dpi=300)

  ## representation for bud breaks
  library(ggplot2)
  library(viridis)

  # # 1. Prepare Spatial Data (Average offset per source site)
  map_data <- df_saclay %>%
    group_by(Site, Longitude, Latitude, Country) %>%
    summarise(
      mean_offset = mean(offset, na.rm = TRUE),
      n_trees = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(Longitude), !is.na(Latitude))

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
    geom_point(data = pts, 
           aes(x = Longitude, y = Latitude, size = n_trees), 
           shape = 21,            # Circle with a border
           color = "black",       # Fixed color for the border
           fill = NA,             # Empty center
           stroke = 0.8,          # Adjust border thickness
           alpha = 0.6) +
    
    # Add Saclay as a geom_point with aes() mapping to trigger a legend
    geom_point(data = data.frame(lon = 2.16, lat = 48.72, Site = "Saclay"), 
               aes(x = lon, y = lat, color = Site), 
               size = 3,          # Adjust size as needed
               shape = 16) +      # Solid circle
    
    # Match the fill color scale to your image (Blue -> Yellow -> Red)
    scale_fill_gradientn(colors = c("#9ecae1", "#ffeda0", "#feb24c", "#f03b20"), 
                        name = "Offset", na.value = "transparent") +
    
    # Add the manual color scale to create the legend for Saclay
    scale_color_manual(name = "Orchard", values = c("Saclay" = "red")) +
    
    coord_sf(xlim = c(-15, 40), ylim = c(35, 65), expand = FALSE) +
    theme_minimal() +
    labs(title = "Spatial Extrapolation of Genomic Offset",
        x = "Longitude", y = "Latitude")

  # Save
  ggsave(paste0(current_out_dir, "Plots/Extrapolation_Map_Offset.png"), p_extrap, width = 10, height = 7, dpi = 300)
  ggsave(paste0(current_out_dir, "Plots/Extrapolation_Map_Offset.svg"), p_extrap, width = 10, height = 7)




  # 1. Prepare Bud Break Data
  # Assuming 'bud_break' is the column name in your df_saclay
  map_data_bb <- df_saclay %>%
    group_by(Site, Longitude, Latitude, Country) %>%
    summarise(
      mean_bb = mean(Bud_Break, na.rm = TRUE), # Use your specific column name here
      n_trees = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(Longitude), !is.na(Latitude), !is.na(mean_bb))

  # 2. Get Europe Base Map
  europe_sf <- ne_countries(scale = "medium", returnclass = "sf", continent = "Europe")

  # 3. Create the Point Map
  p_bud_break <- ggplot() +
    # Draw the background map
    geom_sf(data = europe_sf, fill = "gray95", color = "gray80", size = 0.2) +
    
    # Add the sampling points: Color and Size mapped to Bud Break value
    geom_point(data = map_data_bb, 
              aes(x = Longitude, y = Latitude, 
                  fill = mean_bb, 
                  size = mean_bb), 
              shape = 21,      # Circle with border
              color = "black", # Border color
              stroke = 0.5, 
              alpha = 0.8) +
    
    # Highlight Saclay Orchard
    geom_point(data = data.frame(lon = 2.16, lat = 48.72, Site = "Saclay"), 
              aes(x = lon, y = lat, color = Site), 
              size = 4, 
              shape = 18) + # Diamond shape for the orchard
    
    # Color Scale: Green (Early/Low) to Red (Late/High)
    scale_fill_distiller(palette = "RdYlGn", 
                        direction = 1, # 1 for Green at low, Red at high
                        name = "Bud Break Value\n(Lower = Earlier)") +
    
    # Size Scale: Helps visualize the "lateness"
    scale_size_continuous(range = c(2, 7), name = "Bud Break Value") +
    
    # Manual color for Saclay
    scale_color_manual(name = "Common Garden", values = c("Saclay" = "blue")) +
    
    # Set Map Extent
    coord_sf(xlim = c(-12, 35), ylim = c(35, 65), expand = FALSE) +
    
    # Clean Theme
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "aliceblue", color = NA), # Light blue sea
      legend.position = "right",
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = "Spatial Distribution of Bud Break Phenology",
      subtitle = "Points represents site means from origin populations",
      x = "Longitude", 
      y = "Latitude"
    )

  # 4. Save
  ggsave(paste0(current_out_dir, "Plots/Bud_Break_Distribution_Map.png"), p_bud_break, width = 10, height = 7, dpi = 300)
  ggsave(paste0(current_out_dir, "Plots/Bud_Break_Distribution_Map.svg"), p_bud_break, width = 10, height = 7)
}