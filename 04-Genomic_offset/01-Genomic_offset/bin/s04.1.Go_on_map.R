rm(list=ls())

library(fields)
library(maps)

source("./bin/fonction.R")
source("./bin/config.R")

# -------------------------------
# 1. Load metadata & offsets
# -------------------------------
info_ref_all_SYL <- read.table("./in/data_tree_cluster_dom_admixed_rm.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
sorted_info <- info_ref_all_SYL[order(info_ref_all_SYL$Order_vcf), ]

go_file <- "./out/dom_admixed_rm/s03.GO/file/All_chr_offset_GEA_SNP.csv"
data_offset <- read.table(go_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)

sorted_info_offset <- cbind(sorted_info, data_offset)

# -------------------------------
# 2. Fix coordinates (comma -> dot)
# -------------------------------
sorted_info_offset$Longitude <- as.numeric(gsub(",", ".", sorted_info_offset$Longitude))
sorted_info_offset$Latitude  <- as.numeric(gsub(",", ".", sorted_info_offset$Latitude))

valid_pts <- sorted_info_offset[complete.cases(sorted_info_offset[, c("Longitude", "Latitude")]), ]

# Extract offset
offset <- as.numeric(valid_pts$offset)
valid_pts <- valid_pts[!is.na(offset), ]
offset <- offset[!is.na(offset)]

coords <- cbind(valid_pts$Longitude, valid_pts$Latitude)

# ================================
# 3. FIRST PLOT: sample coordinates
# ================================
# Définir un dégradé de couleurs
cols <- colorRampPalette(c("blue", "red"))(100)  # 100 couleurs

# Normaliser les valeurs d'offset entre 1 et 100 pour l'indexation
offset_scaled <- as.numeric(cut(valid_pts$offset, breaks = 100))

png("./out/test_map.png", width = 1200, height = 900)

plot(
  x = valid_pts$Longitude,
  y = valid_pts$Latitude,
  cex = 0.6,
  pch = 19,
  col = cols[offset_scaled],  # couleur selon offset
  xlab = "Longitude",
  ylab = "Latitude",
  main = "Sample coordinates colored by offset",
  las = 1
)

maps::map(add = TRUE, interior = FALSE, col = "grey40")

dev.off()



# # ================================
# # 4. SECOND PLOT: interpolated offset map
# # ================================

# # TPS interpolation model
# tps_model <- fields::Tps(coords, offset)

# # Prediction grid
# lon_seq <- seq(min(valid_pts$Longitude), max(valid_pts$Longitude), length = 200)
# lat_seq <- seq(min(valid_pts$Latitude),  max(valid_pts$Latitude), length = 200)
# grid <- expand.grid(Longitude = lon_seq, Latitude = lat_seq)

# # Predict
# z_pred <- predict(tps_model, newdata = grid)
# z_matrix <- matrix(z_pred, nrow = length(lon_seq), ncol = length(lat_seq))

# # Plot
# png("./out/test_map_go.png", width = 1200, height = 900)

# my.colors <- colorRampPalette(c("lightblue3", "orange2", "red3"))(100)

# fields::image.plot(
#   x = lon_seq,
#   y = lat_seq,
#   z = z_matrix,
#   col = my.colors,
#   xlab = "Longitude",
#   ylab = "Latitude",
#   las = 1
# )

# maps::map(add = TRUE, interior = FALSE, col = "grey40")

# points(
#   x = valid_pts$Longitude,
#   y = valid_pts$Latitude,
#   col = "grey40",
#   cex = 0.4
# )

# dev.off()





# -------
library(ggplot2)
library(plotly)
library(sf)
library(rnaturalearth)

palblind <- colorRampPalette(c("lightblue3", "orange2", "red3"))(100)
# Données spatiales Afrique
# africa <- ne_countries(continent = "Europe", scale = "medium", returnclass = "sf")
europe <- ne_countries(continent = "Europe", scale="medium", returnclass="sf")

#passeport_490points = read.csv("accessions_bbox.csv")
# Créer un data.frame avec coordonnées et décalage génomique
# df_gap_env_Canbbox <- data.frame(
#                      Long = X_pred_bbox$lon,
#                      Lat = X_pred_bbox$lat,
#                      gap = g_gap_candidate_Canbbox$distance)

# df_gap_GO_Canbbox <- data.frame(
#                      Long = X_pred_bbox$lon,
#                      Lat = X_pred_bbox$lat,
#                      gap = g_gap_candidate_Canbbox$offset)

# #df_gap_GO_Canbbox_outliers <- df_gap_GO_Canbbox[df_gap_GO_Canbbox$gap <6, ]
# #df_gap_env_Canbbox_outliers <- df_gap_env_Canbbox[df_gap_env_Canbbox$gap <1.2, ]
# plot(df_gap_GO_Canbbox$gap)
# abline(h = 1.6, col = "red", lty = 2)

# Carte soit GO soit distance genomique
# p_Canbbox = ggplot() +
#   geom_sf(data = europe, fill = "gray95", color = "black") +
#   geom_point(data = valid_pts, aes(x = Longitude, y = Latitude, color = offset), size = 2.5) +
#   scale_color_gradientn(colours = palblind, name = "Genomic offset") +
#   theme_bw() +
#   coord_sf(xlim = c(-15, 35), ylim = c(35, 65)) +
#   labs(title = "Genomic offset map",
#        x = "Longitude", y = "Latitude")

# Ajout d'un léger décalage pour éviter la superposition (jitter)
valid_pts$Longitude_j <- jitter(valid_pts$Longitude, amount = 0.2)
valid_pts$Latitude_j  <- jitter(valid_pts$Latitude, amount = 0.6)

valid_pts <- valid_pts %>%
  mutate(
    cluster_FS = ifelse(cluster_FS %in% c("wild-wild admixed"),
                        cluster_FS,
                        "pure wild")
  )

# Code du plot
p_Canbbox = ggplot() +
  geom_sf(data = europe, fill = "gray95", color = "black") +
  geom_point(
    data = valid_pts,
    aes(
      x = Longitude_j,
      y = Latitude_j,
      color = offset,
      shape = cluster_FS      # <- forme selon statut
    ),
    size = 1
  ) +
  scale_color_gradientn(colours = palblind, name = "Genomic offset") +
  scale_shape_manual(
    values = c(
      "wild-wild admixed" = 17,   # triangle
      "pure wild" = 16       # cercle
    )
  ) +
  theme_bw() +
  coord_sf(xlim = c(-15, 35), ylim = c(35, 65)) +
  labs(title = "Genomic offset map",
       x = "Longitude", y = "Latitude")

ggsave(filename="test.png", plot=p_Canbbox)

# png("offset_map_interpolated.png", width = 1600, height = 1200, res = 200, type = "cairo")

# fields::image.plot(valid_pts$Longitute, valid_pts$Latitude, valid_pts$offset, 
#                    col = palblind,
#                    las = 1,
#                    xlab = "Longitude",
#                    ylab = "Latitude")
# maps::map(add = TRUE, interior = FALSE, col = "grey40")
# points(valid_pts$Longitude, valid_pts$Latitude, col = "grey40", cex = .3)
# dev.off()


## --- 


library(fields)
library(maps)

# Clean data
valid_clean <- valid_pts[
  is.finite(valid_pts$Longitude) &
  is.finite(valid_pts$Latitude) &
  is.finite(valid_pts$offset), ]

x <- valid_clean$Longitude
y <- valid_clean$Latitude
z <- valid_clean$offset

# Cap extreme values
z_cap <- z
z_cap[z_cap > 0.11] <- 0.11

# Fit TPS model
tps_model <- Tps(cbind(x, y), z)

# Define grid explicitly as a list of ranges
grid_pred <- predictSurface(tps_model, grid.list = list(
  x = seq(-15, 35, length.out = 150),
  y = seq(35, 65, length.out = 150)
))

# Extract matrices
long.mat <- grid_pred$x
lat.mat  <- grid_pred$y
z.mat    <- grid_pred$z

# Save figure
library(ragg) 
png("./offset_map_interpolated.png", width = 1600, height = 1200, res = 200)
my.colors <- colorRampPalette(c("lightblue3", "orange2", "red3"))(100)

fields::image.plot(long.mat, lat.mat, z.mat,
                   col = my.colors,
                   xlab = "Longitude",
                   ylab = "Latitude",
                   las = 1)

maps::map(add = TRUE, interior = FALSE, col = "black")
points(x, y, col = "grey40", cex = 0.3, pch = 19)
dev.off()
