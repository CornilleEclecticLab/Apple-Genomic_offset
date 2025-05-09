rm(list = ls())

# Set the library - JUBAIL
#library(dplyr, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
#library(ggplot2, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
#library(gridExtra, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
#library(rgl, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
#library(plotly, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
#library(scatterplot3d, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')

# Set the library - local
library(dplyr)
library(ggplot2)
library(gridExtra)
library(rgl)
library(plotly)
library(scatterplot3d)

# Set input parameters:
# args <- commandArgs(trailingOnly = TRUE)
# 
# Project_name <- args[1]
# OUT_FOLDER <- args[2]
# # Directory containing the Q files
# Q_FILES_DIR <- args[3]
#Information file
#INFO_FILE <- args[4]
#info <- read.table(INFO_FILE, header=T, sep='\t', dec='.')

Project_name="All_data_wo_orientalis"
OUT_FOLDER="./out/All_data_wo_orientalis/s06.max_membership_coefficient"
Q_FILES_DIR="out/All_data_wo_orientalis/s05.membership_coefficient"


if (!dir.exists(OUT_FOLDER)) {
  dir.create(OUT_FOLDER, recursive = TRUE)
}


# FUNCTION :

extract_k_number <- function(file_path) {
  k_number <- sub(".*K(\\d+).*", "\\1", file_path)
  return(as.numeric(k_number))
}

files = list.files(path = Q_FILES_DIR, pattern = ".csv$", full.names = T)

# filtration and order based on the K value
files <- files[!grepl("K1(?!0)", files, perl = TRUE)]
files_ordered <- files[order(sapply(files, extract_k_number))]

# Plot the max membership coefficient
plots_list <- list()
for (file in files_ordered) {
  print(file)
  struct_file <- read.csv(file, header = TRUE, sep = ",", dec = ".")
  #struct_file <- struct_file[, 1:(ncol(struct_file) - 6)] # 6 is the number of useless column to plot the max membership coeff
  #sub.struct <- struct_file[,1:length(struct_file)]
  #sub.struct <- data.frame(sub.struct)
  struct_file <- struct_file[, grepl("V", names(struct_file))]
  max.pos <- max.col(struct_file, ties.method = c("random", "first", "last"))
  print(max.pos)
  
  max_cluster <- NA
  for (i in 1:nrow(struct_file)) {
    max_cluster[i] <- struct_file[i, max.pos[i] ]
  }
  max_cluster <- data.frame(max_cluster)
  
  max_cluster$ID <- seq(1,nrow(struct_file),1)
  max_cluster$grp <- max.pos
  max_cluster$individuals <- struct_file$Seq_ID
  
  fig <-
    ggplot(max_cluster, aes(x = max_cluster)) + geom_density(aes(y = ..count..), fill = "red3") +
    scale_x_continuous(name = 'Proportion of assignation',
                       breaks = round(seq(
                         min(max_cluster$max_cluster),
                         max(max_cluster$max_cluster),
                         by = 0.05
                       ), 1)) +
    scale_y_continuous(name = 'Density') + labs(title = paste0("K ", ncol(struct_file))) + geom_vline(
      xintercept = 0.82,
      linetype = "dashed",
      color = "blue",
      size = 1
    )  + 
    # geom_vline(
    #   xintercept = 0.85,
    #   linetype = "dashed",
    #   color = "black",
    #   size = 1
    # )  
  theme_bw()
  print(fig)
  plots_list[[file]] <- fig
}
combined_plots <- grid.arrange(grobs = plots_list, ncol = 3)
#plot(combined_plots[1,2]) #K=3
ggsave(paste0(OUT_FOLDER,"/Max_membership_distribution.png"), combined_plots, device = "png",dpi = 300)



# threshold <- 0.82
# step <- 0
# for (file in files) {
#   print(file)
#   K <- sub(".*K(\\d+).*", "\\1", file)
#   print(K)
#   # Lire le fichier CSV
#   struct_file <- read.delim(file, header = TRUE, sep = ",", dec = ".")
#   
#   # Retirer les 7 dernières colonnes
#   struct_file <- struct_file[, 1:(ncol(struct_file) - 5)]
#   # Initialiser le dictionnaire pour les seuils par population
#   define_pop_struct <- {}
#   
#   # Parcourir chaque ligne pour trouver les valeurs maximales par population
#   for (i in 1:nrow(struct_file)) {
#     row_values <- struct_file[i, 1:ncol(struct_file)-1]
#     max_value <- max(row_values)
#     max_col <- names(struct_file)[which.max(row_values)]
#     
#     if (max_value >= threshold){
#       define_pop_struct[[struct_file[i,ncol(struct_file)]]] <- max_col
#     } else {
#       define_pop_struct[[struct_file[i,ncol(struct_file)]]] <- "admix"
#     }
#   }
#   
#   define_pop_struct_df <- data.frame(
#     Seq_ID = names(define_pop_struct),
#     K = unlist(define_pop_struct),
#     row.names = NULL
#   )
#   #print(define_pop_struct_df)
#   colnames(define_pop_struct_df)[2] <- paste0('K', K)
#   if (step == 0){
#     all_df <-  define_pop_struct_df
#     step <- 1
#   }else {
#     all_df <- merge(all_df, define_pop_struct_df, by = "Seq_ID")
#   }
#   
# }
# print(all_df)

## K=3

# k3 <- all_df[,c(1,4)]
# #k3_noadmix$Seq_ID <- k3[k3$K3 != "admix",]
# admix <- k3[k3$K3 == "admix",]
# 
# 
# 
# pca = read.delim("plink_pca/EC_syl_for_pop.eigenvec", header = F, sep = " ")
# colnames(pca) <- c("Seq_ID", "Pop", paste0("PCA",1:20))
# 
# #ref_struct = read.delim("all_syl_wo_ITA/FS_wo_admix/EC_syl/structure_accession_for_all_K_thresh_0.95.txt", header = T, sep = ",")
# ref_struct = all_df
# 
# pca_ref_struct <- merge(pca, ref_struct, by="Seq_ID")
# pca_ref_struct <- pca_ref_struct %>%
#   mutate(country = case_when(
#            K3 == "V1" ~ "Austria",
#            K3 == "V2" ~ "West Europe",
#            K3 == "V3" ~ "Romania",
#            K3 == "admix" ~ "admix"
#          ))
# 
# coordinate <- read.delim("coordinate_seq_id", sep = ";", header = T, dec = ",")
# coor_struct <- merge(pca_ref_struct, coordinate, by = "Seq_ID")
# write.csv(coor_struct, "./K3_coordinate_cluster.csv", row.names = F)
# 
# color_palette_struct <- c( "#CC3333", "#0000ff","grey","green")
# #struct_colors <- setNames(color_palette_struct, c("V1","V2","V3","admix"))
# struct_colors <- setNames(color_palette_struct, c("Austria","West Europe","Romania","admix"))
# #pca_ref_struct$Color  <- struct_colors[pca_ref_struct$K3]
# #struct_colors <- setNames(color_palette_struct, c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","admix"))
# 
# # pca_ref_struct$Color  <- struct_colors[pca_ref_struct$K3]
# 
# pca_ref_struct$K3_numeric <- as.numeric(factor(pca_ref_struct$country, levels = unique(pca_ref_struct$country)))
# 
# png(paste0('./3D_PCA_population_structure_K3.png'), width = 8, height = 5, units = "in", res = 300)
# scatterplot3d(pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3, color = color_palette_struct[pca_ref_struct$K3_numeric], pch = 19, 
#               main = "3D PCA",
#               xlab = paste('PC1 (5.65% explained)'),
#               ylab = paste('PC2 (2.9% explained)'),
#               zlab = paste('PC3 (2.06% explained)'),
#               angle = 70)
# legend("topleft", legend = unique(pca_ref_struct$country),
#        col = color_palette_struct,
#        pch = 19, cex = 0.8)
# dev.off()
