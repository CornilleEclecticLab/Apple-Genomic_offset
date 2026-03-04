rm(list = ls())

# Set the library - JUBAIL
library(dplyr, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(ggplot2, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(gridExtra, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(rgl, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(plotly, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')
library(scatterplot3d, lib='/home/mc10739/R/x86_64-conda-linux-gnu-library/4.4')

# Set the library - local
#library(dplyr)
#library(ggplot2)
#library(gridExtra)
#library(rgl)
#library(plotly)
#library(scatterplot3d)

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

Project_name="dom_admixed_rm"
# Project_name="All_data_wo_orientalis"
OUT_FOLDER = paste0("./out/",Project_name,"/s06.max_membership_coefficient")
Q_FILES_DIR = paste0("./out/",Project_name,"/s05.membership_coefficient")


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

  clean_name <- sub("_value_sorted\\.csv$", "", sub(".*/", "", file))
  clean_name

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
      xintercept = 0.99,
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

  ggsave(paste0(OUT_FOLDER,"/Max_membership_distribution_",clean_name,".png"), fig, device = "png",dpi = 300, width = 12, height = 8 )
  plots_list[[file]] <- fig
}
combined_plots <- grid.arrange(grobs = plots_list, ncol = 5)
#plot(combined_plots[1,2]) #K=3
ggsave(paste0(OUT_FOLDER,"/Max_membership_distribution.png"), combined_plots, device = "png",dpi = 300, width = 12, height = 8 )

