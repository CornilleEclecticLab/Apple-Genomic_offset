# load function
library(rgl)
library(plotly)
library(scatterplot3d)
library(ggplot2)

### -------- input K6 Major FS

# METHOD <- "FS"
# K_value <- "K6"
# PROJECT <- "All_data_wo_orientalis"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/cluster_assignation_FS_K6_Major")
# OUTPUT_NAME_3D="/PCA_3D_K6_Major_"
# OUTPUT_NAME_2D="/PCA_2D_K6_Major_"
# desired_levels <- c(
#   "V1-P.West_Eur",
#   "V2-P.Dom",
#   "V3-P.Aus",
#   "V4-P.Ita",
#   "V5-P.Rom",
#   "V6-P.Lor",
#   "crops-wild admixed",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#56b4e9ff",        # V1 - P.West_Eur
#   "#d55e00ff",  # V2 - Pdom
#   "#009e73ff",  # V3 - P.Aus
#   "#e69f00ff",       # V4 - P.Est_Fr
#   "#0072b2ff",  # V5 - P.Rom
#   "#cc79a7ff",  # V6 - P.Lor
#   "black",  # crops-wild admixed
#   "grey"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23)

### --------- input K6 FS Minor
# METHOD <- "FS"
# K_value <- "K6"
# PROJECT <- "All_data_wo_orientalis"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/cluster_assignation_FS_K6_Minor")
# OUTPUT_NAME_3D="/PCA_3D_K6_Minor_"
# OUTPUT_NAME_2D="/PCA_2D_K6_Minor_"
# desired_levels <- c(
#   "V1-P.West_Eur",
#   "V2-P.Dom",
#   "V3-P.Aus",
#   "V4-P.Ita",
#   "V5-P.Rom",
#   "crops-wild admixed"
# )
# color_palette_struct <- c(
#   "#56b4e9ff",        # V1 - P.West_Eur
#   "#d55e00ff",  # V2 - Pdom
#   "#009e73ff",  # V3 - P.Aus
#   "#e69f00ff",       # V4 - P.Est_Fr
#   "#0072b2ff",  # V5 - P.Rom
#   "#d7ad68" # crops-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21)

## --------- input K5 FS
# METHOD <- "FS"
# K_value <- "K5"
# PROJECT <- "All_data_wo_orientalis"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/cluster_assignation_FS_K5_Major")
# OUTPUT_NAME_3D="/PCA_3D_K5_Major"
# OUTPUT_NAME_2D="/PCA_2D_K5_Major"
# desired_levels <- c(
#   "V1-P.West_Eur",
#   "V2-P.Dom",
#   "V3-P.Aus",
#   "V4-P.Ita",
#   "V5-P.Rom",
#   "crops-wild admixed"
# )
# color_palette_struct <- c(
#   "#56b4e9ff",  # V1 - P.West_Eur
#   "#d55e00ff",  # V2 - Pdom
#   "#009e73ff",  # V3 - P.Aus
#   "#e69f00ff",  # V4 - P.Est_Fr
#   "#0072b2ff",  # V5 - P.Rom
#   "black" # crops-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21)

### --------- input K6 snmf
# METHOD <- "snmf"
# K_value <- "K6"
# PROJECT <- "All_data_wo_orientalis"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/cluster_assignation_snmf_K6")
# OUTPUT_NAME_3D="/PCA_3D_K6_"
# OUTPUT_NAME_2D="/PCA_2D_K6_"
# desired_levels <- c(
#   "V1-P.Ita",
#   "V2-P.Aus",
#   "V3-P.Dom",
#   "V4-P.West_Eur",
#   "V5-P.Lor",
#   "V6-P.Rom",
#   "crops-wild admixed",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#e69f00ff",  # V1-P.Ita
#   "#009e73ff",  # V2-P.Aus
#   "#d55e00ff",  # V3-P.Dom
#   "#56b4e9ff",  # V4-P.West_Eur
#   "#cc79a7ff",  # V5-P.Lor
#   "#0072b2ff",  # V6-P.Rom
#   "black",  # crops-wild admixed
#   "grey"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25)

### ------- input K5 snmf
METHOD <- "snmf"
K_value <- "K5"
PROJECT <- "All_data_wo_orientalis"
PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered.eigenvec")
OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
ref_struct_file <- paste0("./in/cluster_assignation_snmf_K5")
OUTPUT_NAME_3D="/PCA_3D_K5_"
OUTPUT_NAME_2D="/PCA_2D_K5_"
desired_levels <- c(
  "V1-P.Dom",
  "V2-P.Rom",
  "V3-P.West_Eur",
  "V4-P.Aus",
  "V5-P.Ita",
  "crops-wild admixed",
  "wild-wild admixed"
)
color_palette_struct <- c(
  "#d55e00ff",  # V1 - Pdom
  "#0072b2ff",  # V2 - P.Rom
  "#56b4e9ff",  # V3 - P.West_Eur
  "#009e73ff",  # V4 - P.Aus
  "#e69f00ff",  # V5 - P.Ita
  "black",  # crops-wild admixed
  "grey"   # wild-wild admixed
)
shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24)