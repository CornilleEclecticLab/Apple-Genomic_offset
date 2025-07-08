library(dplyr)
library(purrr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(elevatr)

fst_folder = "output/avg_fst_Site_grouped_level_EC_syl.light.reformat/"
coordinate = read.delim("input/mean_coordinate_site_grouped", sep = "\t", dec = ",")

coord_order_long <- coordinate[order(coordinate$Longitude),]

fst_files = list.files(fst_folder, full.names = T)

combined_data <- fst_files %>%
  map_df(~ read.csv(.))  # Use purrr::map_df to read and combine all files

# Step 3: Inspect the combined data
head(combined_data)
#combined_data <- combined_data[match(combined_data$pop1,coord_order_long$Site),]
combined_data$pop2

# Step 4: Group by Site and calculate the mean avg_pi for each site
site_means <- combined_data %>%
  group_by(pop1, pop2) %>%
  summarise(mean_avg_fst = mean(avg_fst, na.rm = TRUE))

#site_means <- site_means[match(site_means$pop2,coord_order_long$Site),]

# pop1_un <- unique(site_means$pop1)
# pop2_un <- unique(site_means$pop2)
# 
# pop_un <- unique(c(pop1_un,pop2_un))
# 
# new_form <- data.frame()
# for (i in 1:nrow(site_means)){
#   print(site_means$pop1[i])
# }

data_sym <- site_means


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
data_matrix_ordered <- data_matrix[coord_order_long$Site, coord_order_long$Site]
upper_triangle_matrix <- data_matrix_ordered
upper_triangle_matrix[lower.tri(upper_triangle_matrix, diag = F)] <- NA


# Create the heatmap
# heatmap(upper_triangle_matrix, 
#         col = heat.colors(256), 
#         scale = "none", 
#         margins = c(8, 8),
#         xlab = "Population 2", 
#         ylab = "Population 1", 
#         main = "Pairwise Mean Average Fst Heatmap")

# Melt the data frame for ggplot
data_melt <- melt(as.matrix(upper_triangle_matrix), na.rm = F)

# data_melt$Var1 <- factor(data_melt$Var1, levels = coord_order_long$Site)
# 
# # Sort the dataframe by the reordered 'Var2'
# data_melt <- data_melt[order(data_melt$Var1), ]
# 
# data_melt$Var2 <- factor(data_melt$Var2, levels = coord_order_long$Site)
# 
# # Sort the dataframe by the reordered 'Var2'
# data_melt <- data_melt[order(data_melt$Var2), ]

# Create heatmap
ggplot(data_melt, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#1D60B4", high = "#FF462B", na.value = "white", name = "Fst Value") +
  geom_text(aes(label = round(value, 3)), color = "black", size = 3) +
  theme_minimal() +
  labs(x = "Population 1", y = "Population 2", title = "Pairwise Mean Average Fst Heatmap") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

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
rownames(dist_matrix_km) <- colnames(dist_matrix_km) <- coordinate$Site

# Print the distance matrix
dist_matrix_km

####--- elevation
coord_for_elevation <- data.frame(x = coordinate[,3], y = coordinate[,2])
data_elevation <- get_elev_point(coord_for_elevation, prj = 4326, src = "aws")
elevation <- cbind(coordinate, data_elevation$elevation)

# Calculate the difference matrix
diff_elevation_matrix <- outer(elevation$`data_elevation$elevation`, elevation$`data_elevation$elevation`, FUN = "-")

# Add row and column names
rownames(diff_elevation_matrix) <- colnames(diff_elevation_matrix) <- elevation$Site

####--- correlation
data_fst <- data_matrix

lower_tri_indices <- lower.tri(dist_matrix_km, diag = FALSE)
dist_vector <- dist_matrix_km[lower_tri_indices]
fst_vector <- data_fst[lower_tri_indices]
elev_vector <- abs(diff_elevation_matrix[lower_tri_indices])

table_fst <- as.data.frame(as.table(data_fst))
write.csv(data_fst, "./output/matrix_fst_site_grouped.csv", row.names = F)
table_fst<- table_fst[as.vector(upper.tri(data_fst)), ]; dim(table_fst)
colnames(table_fst) <- c("pop1", "pop2", "Fst")

table_dist <- as.data.frame(as.table(dist_matrix_km))
write.csv(dist_matrix_km, "./output/matrix_dist_site_grouped.csv", row.names = F)
table_dist<- table_dist[as.vector(upper.tri(dist_matrix_km)), ]; dim(table_dist)
colnames(table_dist) <- c("pop1" , "pop2", "distance")

table_elev <- as.data.frame(as.table(diff_elevation_matrix))
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
correlation_matrix <- cor(data, method = "pearson") # or "spearman" for rank correlation

# Print the correlation matrix
print(correlation_matrix)

##### Test de mentel
library(vegan)

Fst_geo <- mantel(data_fst, dist_matrix_km, method = "spearman", permutations = 9999, na.rm = T)
print(Fst_geo)


## Fst
ggplot(data_moins_fst_dist, aes(x = distance, y = Fst)) + 
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
  annotate("text", x = -Inf, y = Inf, label = "r = 0.8292", vjust = 2, hjust = -0.1, size = 5, fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, label = "p = 0.0004", vjust = 4, hjust = -0.1, size = 5, fontface = "bold")

  # Fst/1-Fst
ggplot(data_moins_fst_dist, aes(x = distance, y = moins_fst)) + 
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
  annotate("text", x = -Inf, y = Inf, label = "Mentel r = 0.8292", vjust = 3.5, hjust = -0.1, size = 5, color = "navy") +
  annotate("text", x = -Inf, y = Inf, label = " (P = 4e-4)", vjust = 5, hjust = -0.1, size = 5, color = "navy")


## sans gradient

# ggplot(data_moins_fst_dist, aes(x = distance, y = moins_fst)) + 
#   geom_point(size = 5, alpha = 0.75, colour = "black",shape = 21, fill = "grey") + 
#   geom_smooth(method = "lm", colour = "black", alpha = 0.2) + 
#   geom_text(aes(label = label), size = 3)+
#   labs(x = "log Distance", y = "Fst/1-Fst") + 
#   theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
#          axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
#          axis.title= element_text(face = "bold", size = 14, colour = "black"), 
#          panel.background = element_blank(), 
#          panel.border = element_rect(fill = NA, colour = "black"),
#          legend.position = "top",
#          legend.text = element_text(size = 10, face = "bold"),
#          legend.title = element_text(size = 11, face = "bold")) +
#   #scale_fill_continuous(high = "navy", low = "skyblue") +
#   annotate("text", x = -Inf, y = Inf, label = "r = 0.8292", vjust = 2, hjust = -0.1, size = 5, fontface = "bold") +
#   annotate("text", x = -Inf, y = Inf, label = "p = 0.0004", vjust = 4, hjust = -0.1, size = 5, fontface = "bold")
# 