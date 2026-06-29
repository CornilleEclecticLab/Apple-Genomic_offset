rm(list=ls())

library(raster)
library(rgdal)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(ggrepel)

Project <- "dom_admixed_rm"
CLIM_FOLD <- "./out/s02.climate_data/"
OUT_FOLD <- "./out/s04.PCA/"
PCA_FOLD <- paste0(OUT_FOLD,Project, "/PCA/")
CORR_FOLD <- paste0(OUT_FOLD,Project, "/Corr/")

if (!dir.exists(OUT_FOLD)) {
        dir.create(OUT_FOLD, recursive = TRUE)
        dir.create(PCA_FOLD, recursive = TRUE)
        dir.create(CORR_FOLD, recursive = TRUE)
      }



clim_file <- paste0(CLIM_FOLD, Project ,"/climate_extract_current_", Project,".csv")
ind_file <- "./in/coordinate_dom_admixed_rm_V2.txt"


climate_data <- read.csv(clim_file)
ind_data <- read.table(ind_file, sep="\t", header = T, dec=",")

all_data <- merge(ind_data, climate_data, by="Seq_ID")

idx <- (length(colnames(all_data)) - 18):length(colnames(all_data))
# colnames(all_data)[idx] <- reduce_bio_name

all_data <- all_data %>%
  mutate(Name = Seq_ID)

# This returns the names "bio1", "bio2", etc. instead of TRUE/FALSE
bio_names <- grep("bio", colnames(all_data), value = TRUE)

# Then run your PCA using those names
pca_all_data <- PCA(all_data[, bio_names], 
                    scale.unit = TRUE, 
                    ncp = 4, 
                    graph = FALSE)

coord_PCA <- cbind(all_data["Name"], pca_all_data$ind$coord)
coord_PCA_ind <- coord_PCA[1:length(all_data$Name),]

col_vals <- c(
  "Romania"     = "#1b9e77",
  "Denmark"     = "#d95f02",
  "France"     = "#7570b3",
  "Austria"     = "#e7298a",
  "Belgique"     = "#66a61e",
  "Italy"     = "#e6ab02"
  #"Orchard" = "#666666"
)

all_data$Country <- factor(
all_data$Country,
levels = names(col_vals)
)

shape_vals <- c(16, 17, 15, 3, 7, 8, 18)

# Make site names unique
unique_sites <- make.unique(as.character(all_data$Site))

# Assign them as rownames
rownames(pca_all_data$ind$coord) <- unique_sites

png(paste0(PCA_FOLD,"/PCA_climate_vars.png"), width = 2000, height = 1500, res = 300)

factoextra::fviz_pca_biplot(
  pca_all_data,
  geom.ind   = c("point", "text"),  # points AND labels
  label.ind  = "all",               # use rownames of PCA object
  repel      = TRUE,
  labelsize  = 2,                   # smaller text
  col.ind    = all_data$Country,
  col.var    = "grey40",
  palette    = col_vals,
  pointsize  = 2,
  mean.point = FALSE,
  legend.title = "Legend",
  title       = "PCA"
) +
  theme_minimal() +
  scale_shape_manual(values = shape_vals)

dev.off()


# 1. Subset and calculate correlation using the sorted names
cor_matrix_ordered <- cor(all_data[, bio_names], use = "complete.obs")

# 2. Generate the plot
png(paste0(CORR_FOLD, "/Corr_climate_vars_ordered.png"), width = 1500, height = 1500, res = 300)

corrplot(cor_matrix_ordered, 
         method = "color",       
         type = "upper",         
         order = "original",     # <--- This keeps the bio01 to bio19 sequence
         addCoef.col = "black",  
         number.cex = 0.4,       # Small numbers to fit 19 columns
         tl.cex = 0.6,           # Small labels to fit 19 columns
         tl.col = "black",       
         tl.srt = 45,            
         diag = FALSE)           

dev.off()