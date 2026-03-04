rm(list=ls())

library(dplyr)
library(purrr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(elevatr)
library(svglite)
library(ragg)


#function
create_folder <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Folder created: ", path)
  } else {
    message("Folder already exists: ", path)
  }
}

save_plot_png_svg <- function(plot, filename, width = 10, height = 6, dpi = 300) {
  # Save as PNG
  ggsave(paste0(filename, ".png"), plot = plot, width = width, height = height, dpi = dpi, units = "in")
  
  # Save as SVG
  ggsave(paste0(filename, ".svg"), plot = plot, width = width, height = height, units = "in")
  
  message("Plot saved as: ", filename, ".png and .svg")
}


#input
Project="full_genome_dom_admixed_rm_no_filter"
# Project="dom_admixed_rm"
fst_folder = paste0("../02-Pixy/out/",Project,"/s02.avg_value/")
coordinate_folder = paste0("./in/group_coordinate_",Project,".txt")
out_folder=paste0("out/",Project)
create_folder(out_folder)
out_metrics=paste0(out_folder,"/s01.Metrics")
create_folder(out_metrics)
out_IBD=paste0(out_folder,"/s02.IBD")
create_folder(out_IBD)

#read
coordinate = read.delim(coordinate_folder, sep = "\t", dec = ",")
coord_order_long <- coordinate[order(coordinate$Longitude),]


##-------- FST
fst_files = list.files(fst_folder, full.names = T, pattern = 'fst')
combined_data <- fst_files %>%
  map_df(~ read.csv(.)) %>% # Use purrr::map_df to read and combine all files 
  mutate(avg_fst = ifelse(avg_fst < 0, 0, avg_fst))   # set negatives to 0

#Group by pop and calculate the mean avg_pi for each pop
site_means <- combined_data %>%
  group_by(pop1, pop2) %>%
  summarise(mean_avg_fst = mean(avg_fst, na.rm = TRUE))

data_sym <- site_means

## it create a symetric table
for(i in 1:nrow(site_means)) {
  data_sym <- rbind(data_sym, data.frame(
    pop1 = site_means$pop2[i],
    pop2 = site_means$pop1[i],
    mean_avg_fst = site_means$mean_avg_fst[i]
  ))
}

# Reshape data into wide format
data_wide <- dcast(data_sym, pop1 ~ pop2, value.var = "mean_avg_fst", fill = 0)

# Convert to matrix, rownames as population names
rownames(data_wide) <- data_wide$pop1
data_wide$pop1 <- NULL

data_matrix <- as.matrix(data_wide)
diag(data_matrix) <- 0

# Create a matrix with only the upper triangle
data_matrix_ordered <- data_matrix[coord_order_long$Group, coord_order_long$Group]
upper_triangle_matrix <- data_matrix_ordered
upper_triangle_matrix[lower.tri(upper_triangle_matrix, diag = T)] <- NA

# Melt the data frame for ggplot
data_melt <- melt(as.matrix(upper_triangle_matrix), na.rm = F)


# Create heatmap
plot1 <- ggplot(data_melt, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#E8EC67", high = "#D7263D", na.value = "white", name = "Fst Value") +
  geom_text(aes(label = round(value, 3)), color = "black", size = 3) +
  theme_bw() +
  labs(x = "Population 1", y = "Population 2", title = "Pairwise Mean Average Fst Heatmap") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

save_plot_png_svg(plot = plot1, filename = paste0(out_metrics, "/Pairwise_Fst_",Project))

###--------- Dxy

dxy_files = list.files(fst_folder, full.names = T, pattern = 'dxy')
combined_data <- dxy_files %>%
  map_df(~ read.csv(.))  # Use purrr::map_df to read and combine all files

#Group by pop and calculate the mean avg_pi for each pop
site_means <- combined_data %>%
  group_by(pop1, pop2) %>%
  summarise(mean_avg_dxy = mean(avg_dxy, na.rm = TRUE))

data_sym <- site_means

## it create a symetric table
for(i in 1:nrow(site_means)) {
  data_sym <- rbind(data_sym, data.frame(
    pop1 = site_means$pop2[i],
    pop2 = site_means$pop1[i],
    mean_avg_dxy = site_means$mean_avg_dxy[i]
  ))
}

# Reshape data into wide format
data_wide <- dcast(data_sym, pop1 ~ pop2, value.var = "mean_avg_dxy", fill = 0)

# Convert to matrix, rownames as population names
rownames(data_wide) <- data_wide$pop1
data_wide$pop1 <- NULL

data_matrix <- as.matrix(data_wide)
diag(data_matrix) <- 0

# Create a matrix with only the upper triangle
data_matrix_ordered <- data_matrix[coord_order_long$Group, coord_order_long$Group]
lower_triangle_matrix <- data_matrix_ordered
lower_triangle_matrix[upper.tri(upper_triangle_matrix, diag = T)] <- NA

# Melt the data frame for ggplot
data_melt <- melt(as.matrix(lower_triangle_matrix), na.rm = F)


# Create heatmap
plot4 <- ggplot(data_melt, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#D5FFF3", high = "#820B8A", na.value = "white", name = "Dxy Value") +
  geom_text(aes(label = round(value, 3)), color = "black", size = 3) +
  theme_bw() +
  labs(x = "Population 1", y = "Population 2", title = "Pairwise Mean Average Dxy Heatmap") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

save_plot_png_svg(plot = plot4, filename = paste0(out_metrics, "/Pairwise_Dxy_",Project))


###------ geograpghic distance
library(geosphere)
#Calculate the distances between the pairs of sites

dist_matrix <- distm(
  matrix(c(coordinate$Longitude, coordinate$Latitude), ncol = 2),
  fun = distHaversine
)

# Convert the distance matrix from meters to kilometers
dist_matrix_km <- dist_matrix / 1000

# Assign row and column names for clarity
rownames(dist_matrix_km) <- colnames(dist_matrix_km) <- coordinate$Group

# Print the distance matrix
dist_matrix_km

####--- elevation
coord_for_elevation <- data.frame(x = coordinate$Longitude, y = coordinate$Latitude)
data_elevation <- get_elev_point(coord_for_elevation, prj = 4326, src = "aws")
elevation <- cbind(coordinate, data_elevation$elevation)

# Calculate the difference matrix
diff_elevation_matrix <- outer(elevation$`data_elevation$elevation`, elevation$`data_elevation$elevation`, FUN = "-")

# Add row and column names
rownames(diff_elevation_matrix) <- colnames(diff_elevation_matrix) <- elevation$Group

####--- correlation
data_fst <- data_matrix

lower_tri_indices <- lower.tri(dist_matrix_km, diag = FALSE)
dist_vector <- dist_matrix_km[lower_tri_indices]
fst_vector <- data_fst[lower_tri_indices]
elev_vector <- abs(diff_elevation_matrix[lower_tri_indices])

table_fst <- as.data.frame(as.table(data_fst))
write.csv(data_fst, paste0(out_metrics,"/matrix_fst.csv"), row.names = F)

table_dist <- as.data.frame(as.table(dist_matrix_km))
write.csv(dist_matrix_km, paste0(out_metrics,"/matrix_dist.csv"), row.names = F)

table_elev <- as.data.frame(as.table(diff_elevation_matrix))
write.csv(dist_matrix_km, paste0(out_metrics,"/matrix_elevation.csv"), row.names = F)


table_fst<- table_fst[as.vector(upper.tri(data_fst)), ]; dim(table_fst)
colnames(table_fst) <- c("pop1", "pop2", "Fst")
table_dist<- table_dist[as.vector(upper.tri(dist_matrix_km)), ]; dim(table_dist)
colnames(table_dist) <- c("pop1" , "pop2", "distance")
table_elev<- table_elev[as.vector(upper.tri(diff_elevation_matrix)), ]; dim(table_elev)
colnames(table_elev) <- c("pop1" , "pop2", "elevation")

data_fst_dist <- cbind(table_fst, log(table_dist[,3]), abs(table_elev$elevation))
colnames(data_fst_dist) <- c("pop1" , "pop2", "Fst" ,"distance", "elevation")

# Create a data frame to hold the vectors
data <- data.frame(
  Distance = log(dist_vector),
  Fst = fst_vector
)

data_moins_fst_dist <- data_fst_dist %>% 
  mutate(moins_fst = Fst/(1-Fst)) %>%
  mutate(label = paste0(pop1,"_vs_",pop2))


# Compute the correlation matrix
cor.test(data_moins_fst_dist$moins_fst, data_moins_fst_dist$elevation) # correlation with elevation
correlation_matrix <- cor(data, method = "pearson") # or "spearman" for rank correlation

# Print the correlation matrix
print(correlation_matrix)

##### Test de mentel
library(vegan)

Fst_geo <- mantel(data_fst, dist_matrix_km, method = "spearman", permutations = 9999, na.rm = T)
print(Fst_geo)


## Fst
plot2 <- ggplot(data_moins_fst_dist, aes(x = distance, y = Fst)) + 
  geom_point(size = 5, alpha = 0.75, colour = "black",shape = 21, aes(fill = elevation)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) +
  geom_text(aes(label = label), size = 2.5)+
  labs(x = "Distance (km)", y = "Fst", fill = "Elevation difference (m)") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 11, face = "bold")) +
  scale_fill_continuous(high = "navy", low = "skyblue") +
  annotate("text", x = -Inf, y = Inf, label = " M.sylvestris", vjust = 2, hjust = -0.1, size = 5, color = "navy") +
  annotate("text", x = -Inf, y = Inf, label = "Mantel r = 0.6088", vjust = 3.5, hjust = -0.1, size = 5, color = "navy") +
  annotate("text", x = -Inf, y = Inf, label = " (P = 0.0038)", vjust = 5, hjust = -0.1, size = 5, color = "navy")

save_plot_png_svg(plot = plot2, filename = paste0(out_IBD, "/Figure_IBD_fst_",Project))

  # Fst/1-Fst
plot3 <- ggplot(data_moins_fst_dist, aes(x = distance, y = moins_fst)) + 
  geom_point(size = 5, alpha = 0.75, colour = "black",shape = 21, aes(fill = elevation)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) + 
  geom_text(aes(label = label), size = 2.5)+
  labs(x = "Distance (km)", y = "Fst/1-Fst", fill = "Elevation difference (m)") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 11, face = "bold")) +
  scale_fill_continuous(high = "navy", low = "skyblue") +
  annotate("text", x = -Inf, y = Inf, label = " M.sylvestris", vjust = 2, hjust = -0.1, size = 5, color = "navy") +
  annotate("text", x = -Inf, y = Inf, label = "Mantel r = 0.6088", vjust = 3.5, hjust = -0.1, size = 5, color = "navy") +
  annotate("text", x = -Inf, y = Inf, label = " (P = 0.0038)", vjust = 5, hjust = -0.1, size = 5, color = "navy")

save_plot_png_svg(plot = plot3, filename = paste0(out_IBD, "/Figure_IBD_moins_fst_",Project))

plot4 <- ggplot(data_moins_fst_dist, aes(x = elevation, y = moins_fst)) + 
  geom_point(size = 5, alpha = 0.75, colour = "black",shape = 21, aes(fill = distance)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) + 
  geom_text(aes(label = label), size = 2.5)+
  labs(x = "Elevation difference (m)", y = "Fst/1-Fst", fill = "Distance (km)") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 11, face = "bold")) +
  scale_fill_continuous(high = "navy", low = "skyblue") +
  annotate("text", x = -Inf, y = Inf, label = " M.sylvestris", vjust = 2, hjust = -0.1, size = 5, color = "navy") +
  annotate("text", x = -Inf, y = Inf, label = "Mantel r = 0.6088", vjust = 3.5, hjust = -0.1, size = 5, color = "navy") +
  annotate("text", x = -Inf, y = Inf, label = " (P = 0.0038)", vjust = 5, hjust = -0.1, size = 5, color = "navy")

save_plot_png_svg(plot = plot4, filename = paste0(out_IBD, "/Figure_IBD_elevation_moins_fst_",Project))

plot5 <- ggplot(data_moins_fst_dist, aes(x = elevation, y = Fst)) + 
  geom_point(size = 5, alpha = 0.75, colour = "black",shape = 21, aes(fill = distance)) + 
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) + 
  geom_text(aes(label = label), size = 2.5)+
  labs(x = "Elevation difference (m)", y = "Fst", fill = "Distance (km)") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 11, face = "bold")) +
  scale_fill_continuous(high = "navy", low = "skyblue") +
  annotate("text", x = -Inf, y = Inf, label = " M.sylvestris", vjust = 2, hjust = -0.1, size = 5, color = "navy") +
  annotate("text", x = -Inf, y = Inf, label = "Mantel r = 0.6088", vjust = 3.5, hjust = -0.1, size = 5, color = "navy") +
  annotate("text", x = -Inf, y = Inf, label = " (P = 0.0038)", vjust = 5, hjust = -0.1, size = 5, color = "navy")

save_plot_png_svg(plot = plot4, filename = paste0(out_IBD, "/Figure_IBD_elevation_fst_",Project))