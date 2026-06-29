# Clear environment
rm(list=ls())

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggrepel) 
library(stringr)

## ==========================================
## 1. Définition des chemins et fichiers
## ==========================================
PROJECT <- "dom_admixed_rm"
ENV_DIR <- paste0("/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/02-Climate_data/out/s02.climate_data/",PROJECT,"/")
OUT_DIR <- paste0("./out/",PROJECT,"/s11.biovar_temp_comparison/")
GO_FUTURE_F <- paste0("./out/", PROJECT,"/s07.GO_future_V2/file/")
GO_PAST_F <- paste0("./out/",PROJECT,"/s08.Past_climate_GO/file/")

dir.create(OUT_DIR, recursive=T, showWarnings=FALSE)

# Charger les informations des individus
INFO_IND <- read.delim("./in/data_tree_dom_admixed_rm_V2.txt", header=T, sep="\t", dec=",")

# Charger les données du climat actuel
current_file <- paste0(ENV_DIR, "/climate_extract_current_dom_admixed_rm.csv")
current_data <- read.csv(current_file, header = TRUE)

# Lister tous les autres fichiers de scénarios
all_files <- list.files(ENV_DIR, pattern = "climate_extract_", full.names = TRUE)
scenario_files <- all_files[ !all_files %in% current_file ]

## ==========================================
## 2. Préparation de la PCA de référence (Current)
## ==========================================
message("Calcul de la PCA de référence sur le climat actuel...")

bio_columns <- paste0("bio", sprintf("%02d", 1:19))
pca_current <- prcomp(current_data[, bio_columns], center = TRUE, scale. = TRUE)

# --- Extraction et mise à l'échelle des flèches des variables ---
arrow_data <- as.data.frame(pca_current$rotation[, 1:2]) %>% 
  mutate(Variable = rownames(.))

arrow_scale_factor <- 5 
arrow_data <- arrow_data %>% 
  mutate(
    x_end = PC1 * arrow_scale_factor,
    y_end = PC2 * arrow_scale_factor
  )

# Coordonnées du climat actuel + Métadonnées
df_current_pc <- as.data.frame(pca_current$x) %>% 
  mutate(Seq_ID = current_data$Seq_ID) %>% 
  left_join(INFO_IND %>% select(Seq_ID, Country, Site, Latitude, Longitude), by = "Seq_ID")

# --- INITIALISATION DES TABLEAUX RÉSUMÉS ---
global_summary <- data.frame()
correlation_summary <- data.frame() # <- AJOUT : Stockage des corrélations

## ==========================================
## 3. Boucle de comparaison, distances et graphique
## ==========================================
for (scen_path in scenario_files) {
  
  scen_name <- basename(scen_path) %>% 
    gsub("climate_extract_", "", .) %>% 
    gsub("_dom_admixed_rm.csv", "", .) %>% 
    gsub(".csv", "", .)
  
  if (length(grep("future", scen_name)) != 0){
    short_name <- str_extract(scen_name, "[0-9]{4}-[0-9]{4}_ssp[0-9]{3}")
    go_file <- list.files(GO_FUTURE_F, pattern=paste0("*FDR_filtered_",short_name))[1]
    go_data <- read.csv(paste0(GO_FUTURE_F, go_file))
  } else if (length(grep("LGM", scen_name)) != 0) {
    go_file <- list.files(GO_PAST_F, pattern=paste0("*FDR_filtered_",scen_name))[1]
    go_data <- read.csv(paste0(GO_PAST_F, go_file))
  } else { 
    message("No corresponding GO file !\n Passing to the next comparison.")
    next
  }
  
  message(paste("Traitement du scénario :", scen_name))
  scen_data <- read.csv(scen_path, header = TRUE)
  
  # A. PROJECTION PASSIVE du scénario
  scen_scaled <- scale(scen_data[, bio_columns], 
                       center = pca_current$center, 
                       scale = pca_current$scale)
  scen_projected <- as.data.frame(scen_scaled %*% pca_current$rotation) %>% 
    mutate(Seq_ID = scen_data$Seq_ID) %>% 
    left_join(INFO_IND %>% select(Seq_ID, Country, Site, Latitude, Longitude), by = "Seq_ID")
  
  # B. CALCUL DE LA DISTANCE EUCLIDIENNE
  matched_current <- df_current_pc %>% filter(Seq_ID %in% scen_projected$Seq_ID) %>% arrange(Seq_ID)
  matched_scenario <- scen_projected %>% filter(Seq_ID %in% df_current_pc$Seq_ID) %>% arrange(Seq_ID)
  
  mat_current <- as.matrix(matched_current %>% select(starts_with("PC")))
  mat_scenario <- as.matrix(matched_scenario %>% select(starts_with("PC")))
  euclidian_distances <- sqrt(rowSums((mat_current - mat_scenario)^2))
  
  df_distances <- data.frame(
    Seq_ID = matched_current$Seq_ID,
    Country = matched_current$Country,
    Site = matched_current$Site,
    Scenario = scen_name,
    Euclidian_Distance = euclidian_distances
  )
  
  dist_file_name <- paste0(OUT_DIR, "Environmental_distances_Current_vs_", scen_name, ".txt")
  write.table(df_distances, file = dist_file_name, sep = "\t", row.names = FALSE, quote = FALSE)

  # C. CALCUL ET STOCKAGE DE LA CORRÉLATION (Distance vs Genomic Offset)
  merge_dist_go <- merge(df_distances, go_data, by=c("Seq_ID","Country", "Site"))
  
  # Calcul du coefficient de corrélation de Pearson
  cor_score <- cor(merge_dist_go$offset, merge_dist_go$Euclidian_Distance, use = "complete.obs")
  
  # Stockage dans le dataframe dédié
  correlation_summary <- rbind(
    correlation_summary, 
    data.frame(Scenario = scen_name, Genomic_Offset_File = go_file, Pearson_r = cor_score)
  )
  
  # Statistiques basiques de distance
  mean_dist <- mean(euclidian_distances)
  max_dist <- max(euclidian_distances)
  global_summary <- rbind(global_summary, data.frame(Scenario = scen_name, Mean_Distance = mean_dist, Max_Distance = max_dist))
  
  # D. PRÉPARATION DU GRAPHIQUE
  df_current_plot <- matched_current %>% mutate(Env_Time = "Current")
  df_scen_plot <- matched_scenario %>% mutate(Env_Time = scen_name)
  df_combined <- rbind(df_current_plot, df_scen_plot)
  df_combined$Env_Time <- factor(df_combined$Env_Time, levels = c("Current", scen_name))
  
  var_explained <- round(100 * (pca_current$sdev^2 / sum(pca_current$sdev^2)), 1)
  
  # E. GÉNÉRATION DU GRAPHIQUE
  p <- ggplot() +
    # 1. Tracé des flèches des variables
    geom_segment(
      data = arrow_data,
      aes(x = 0, y = 0, xend = x_end, yend = y_end),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "grey", 
      alpha = 0.4, 
      linewidth = 0.5
    ) +
    geom_text(
      data = arrow_data,
      aes(x = x_end * 1.05, y = y_end * 1.05, label = Variable),
      color = "grey", 
      fontface = "bold", 
      size = 3.5
    ) +
    
    # 2. Lignes pointillées
    geom_line(data = df_combined, aes(x = PC1, y = PC2, group = Seq_ID), color = "grey75", alpha = 0.5, linetype = "dashed") +
    
    # 3. Points
    geom_point(data = df_combined, aes(x = PC1, y = PC2, color = Country, shape = Env_Time), alpha = 0.8, size = 2.5) + 
    scale_shape_manual(name = "Climate Period", values = c("Current" = 16, setNames(2, scen_name))) +
    
    # 4. Étiquettes de texte
    geom_text_repel(
      data = df_combined,
      aes(x = PC1, y = PC2, label = Seq_ID, color = Country),
      size = 2.5,                 
      max.overlaps = 200,        
      box.padding = 0.3,        
      point.padding = 0.2,      
      show.legend = FALSE       
    ) +
    
    labs(
      title = paste("PCA Climate Shift & Variables: Current vs", scen_name),
      subtitle = paste0("Grey arrows show environmental variables (scaled x", arrow_scale_factor, "). Mean Distance: ", round(mean_dist, 2), " | r(Dist, GO) = ", round(cor_score, 3)),
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)")
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title = element_text(face = "bold", size = 14)
    )
  
  # Sauvegarde du graphique
  out_name <- paste0(OUT_DIR, "PCA_with_arrows_Current_vs_", scen_name, ".png")
  ggsave(out_name, plot = p, width = 11, height = 9, dpi = 300) 
  
  cat("Graphique biplot complet généré pour :", scen_name, " (r =", round(cor_score, 3), ")\n")
}

## ==========================================
## 4. Sauvegarde des résumés globaux
## ==========================================
# Sauvegarde du résumé des distances
write.table(global_summary, file = paste0(OUT_DIR, "Global_scenarios_distance_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Sauvegarde du résumé des corrélations (Distance vs Genomic Offset)
write.table(correlation_summary, file = paste0(OUT_DIR, "Global_GO_correlation_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

message("Toutes les analyses biplots et les résumés de corrélations sont terminés !")