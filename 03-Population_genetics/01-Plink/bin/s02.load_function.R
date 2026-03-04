# load function
#library(rgl)
library(plotly)
library(scatterplot3d)
library(ggplot2)
library(svglite)

## All_data_wo_orientalis
## K6 snmf
# METHOD <- "snmf"
# K_value <- "K6"
# PROJECT <- "All_data_wo_orientalis"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_All_data_wo_orientalis.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K6_snmf_0.70_cluster_assignation.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "V6",
#   "crops-wild_admixed_snmf",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#FF7B13",
#   "#FF54B3",
#   "#0072B2",
#   "#56B4E9",
#   "#F0E442",
#   "#59EC13",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)

## K6 Major FastStructure
# METHOD <- "FS"
# K_value <- "K6"
# PROJECT <- "All_data_wo_orientalis"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_All_data_wo_orientalis.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K6_Major_FastStructure_0.80_cluster_assignation.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "V6",
#   "crops-wild_admixed_FastStructure",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#FF7B13", 
#   "#FF54B3", 
#   "#0072B2", 
#   "#59EC13", 
#   "#56B4E9", 
#   "#F0E442",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)


## Dom_admixed_rm

## ## K5 Major FastStructure
# METHOD <- "FS"
# K_value <- "K5_Major"
# PROJECT <- "dom_admixed_rm"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K5_Major_FastStructure_0.99_cluster_assignation_11012026.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "crops-wild_admixed_FastStructure",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#56B4E9", 
#   "#59EC13", 
#   "#F0E442", 
#   "#0072B2", 
#   "#FF54B3",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)

## ## K5 MinorCluster1 FastStructure
# METHOD <- "FS"
# K_value <- "K5_MinorCluster1"
# PROJECT <- "dom_admixed_rm"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K5_MinorCluster1_FastStructure_0.99_cluster_assignation_11012026.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "crops-wild_admixed_FastStructure",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#56B4E9", 
#   "#59EC13", 
#   "#F0E442", 
#   "#0072B2", 
#   "#FF54B3",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)

## K5 snmf
# METHOD <- "snmf"
# K_value <- "K5"
# PROJECT <- "dom_admixed_rm"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K5_snmf_0.7_cluster_assignation_11012026.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "crops-wild_admixed_FastStructure",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#F0E442", 
#   "#59EC13", 
#   "#0072B2", 
#   "#56B4E9", 
#   "#FF54B3",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)






## ## K6 Major FastStructure
# METHOD <- "FS"
# K_value <- "K6_Major"
# PROJECT <- "dom_admixed_rm"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K6_Major_FastStructure_0.99_cluster_assignation_11012026.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "V6",
#   "crops-wild_admixed_FastStructure",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#56B4E9", 
#   "#59EC13", 
#   "#F0E442", 
#   "#0072B2", 
#   "#FF54B3",
#   "#ec0e0eff",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)

## ## K6 MinorCluster1 FastStructure
# METHOD <- "FS"
# K_value <- "K6_MinorCluster1"
# PROJECT <- "dom_admixed_rm"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K6_MinorCluster1_FastStructure_0.99_cluster_assignation_11012026.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "V6",
#   "crops-wild_admixed_FastStructure",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#56B4E9", 
#   "#59EC13", 
#   "#F0E442", 
#   "#0072B2", 
#   "#FF54B3",
#   "#ec0e0eff",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)

## K6 snmf
# METHOD <- "snmf"
# K_value <- "K6"
# PROJECT <- "dom_admixed_rm"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K6_snmf_0.73_cluster_assignation_11012026.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "V6",
#   "crops-wild_admixed_FastStructure",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#F0E442", 
#   "#59EC13", 
#   "#0072B2", 
#   "#ec0e0eff", 
#   "#FF54B3", 
#   "#56B4E9",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)




# ## K7 Major FastStructure
# METHOD <- "FS"
# K_value <- "K7_Major"
# PROJECT <- "dom_admixed_rm"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K7_Major_FastStructure_0.99_cluster_assignation_11012026.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "V6",
#   "V7",
#   "crops-wild_admixed_FastStructure",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#56B4E9", 
#   "#59EC13", 
#   "#F0E442", 
#   "#0072B2", 
#   "#FF54B3",
#   "#ec0e0eff",
#   "brown",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)

## ## K7 MinorCluster1 FastStructure
# METHOD <- "FS"
# K_value <- "K7_MinorCluster1"
# PROJECT <- "dom_admixed_rm"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K7_MinorCluster1_FastStructure_0.99_cluster_assignation_11012026.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "V6",
#   "V7",
#   "crops-wild_admixed_FastStructure",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#56B4E9", 
#   "#59EC13", 
#   "#F0E442", 
#   "#0072B2", 
#   "#FF54B3",
#   "#ec0e0eff",
#   "brown",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)

## K7 snmf
# METHOD <- "snmf"
# K_value <- "K7"
# PROJECT <- "dom_admixed_rm"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K7_snmf_0.72_cluster_assignation_11012026.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "V6",
#   "V7",
#   "crops-wild_admixed_FastStructure",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#56B4E9", 
#   "#59EC13", 
#   "#FF54B3", 
#   "#0072B2", 
#   "#ec0e0eff", 
#   "#F0E442",
#   "brown",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)




# ## K8 Major FastStructure
# METHOD <- "FS"
# K_value <- "K8_Major"
# PROJECT <- "dom_admixed_rm"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K8_Major_FastStructure_0.99_cluster_assignation_11012026.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "V6",
#   "V7",
#   "V8",
#   "crops-wild_admixed_FastStructure",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#56B4E9", 
#   "#59EC13", 
#   "#F0E442", 
#   "#0072B2", 
#   "#FF54B3", 
#   "#FF7B13",
#   "#ec0e0eff",
#   "brown",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
# )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)

## ## K8 MinorCluster1 FastStructure
# METHOD <- "FS"
# K_value <- "K8_MinorCluster1"
# PROJECT <- "dom_admixed_rm"
# PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.eigenvec")
# OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
# ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K8_MinorCluster1_FastStructure_0.99_cluster_assignation_11012026.csv")
# OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
# OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
# desired_levels <- c(
#   "V1",
#   "V2",
#   "V3",
#   "V4",
#   "V5",
#   "V6",
#   "V7",
#   "V8",
#   "crops-wild_admixed_FastStructure",
#   "wild-wild admixed"
# )
# color_palette_struct <- c(
#   "#56B4E9", 
#   "#59EC13", 
#   "#F0E442", 
#   "#0072B2", 
#   "#FF54B3", 
#   "#FF7B13",
#   "#ec0e0eff",
#   "brown",
#   "#888888",  # crops-wild admixed
#   "#101010"   # wild-wild admixed
#  )
# shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)

## K7 snmf
METHOD <- "snmf"
K_value <- "K8"
PROJECT <- "dom_admixed_rm"
PCA_file <- paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.eigenvec")
OUTPUT_FOLDER <- paste0('./out/',PROJECT,"/s02.PCA_figure")
ref_struct_file <- paste0("./in/",PROJECT,"_11012026","/K8_snmf_0.56_cluster_assignation_11012026.csv")
OUTPUT_NAME_3D=paste0("/PCA_3D_",K_value,"_")
OUTPUT_NAME_2D=paste0("/PCA_2D_",K_value,"_")
desired_levels <- c(
  "V1",
  "V2",
  "V3",
  "V4",
  "V5",
  "V6",
  "V7",
  "V8",
  "crops-wild_admixed_FastStructure",
  "wild-wild admixed"
)
color_palette_struct <- c(
  "brown", 
  "#56B4E9", 
  "#FF54B3", 
  "#F0E442", 
  "#FF7B13", 
  "#0072B2",
  "#ec0e0eff",
  "#59EC13",
  "#888888",  # crops-wild admixed
  "#101010"   # wild-wild admixed
)
shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21)