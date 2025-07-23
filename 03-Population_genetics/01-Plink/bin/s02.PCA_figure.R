library(rgl)
library(plotly)
library(scatterplot3d)

PROJECT="dom_admixed_rm"

pca = read.delim(paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered_",PROJECT,".eigenvec"), header = F, sep = " ")
colnames(pca) <- c("Pop", "Seq_ID", paste0("PCA",1:20))
dim(pca)

ref_struct = read.delim(paste0("in/population_structure_assignation_group_",PROJECT,".txt"))
ref_struct <- subset(
  ref_struct,
  !(Cluster_FS == "V2 - Pdom" | Cluster_snmf == "V3 - P.dom")
); dim(ref_struct)

pca_ref_struct <- merge(pca, ref_struct, by="Seq_ID")



#color_palette_struct_snmf_dom_rm <- c("#d55e00ff","#d55e00ff","#d55e00ff", "#009e73ff", "#e69f00ff","#0072b2ff", "#56b4e9ff","#cc79a7ff","grey")
#color_palette_struct_FS_dom_rm <- c("#d55e00ff","#d55e00ff","#d55e00ff", "#56b4e9ff", "#009e73ff", "#e69f00ff","#0072b2ff", "#cc79a7ff","grey")

color_palette_struct_snmf_dom_admixed_rm <- c("#009e73ff", "#e69f00ff","#0072b2ff", "#56b4e9ff","#cc79a7ff","grey")
color_palette_struct_FS_dom_admixed_rm <- c("#56b4e9ff", "#009e73ff", "#e69f00ff","#0072b2ff", "#cc79a7ff","grey")

plotly_df_struct <- pca_ref_struct
plotly_df_struct$Cluster_FS <- factor(plotly_df_struct$Cluster_FS)
plotly_df_struct$Cluster_snmf <- factor(plotly_df_struct$Cluster_snmf)

fig.2 <- plot_ly(
  data = plotly_df_struct,
  x = ~PCA1,
  y = ~PCA2,
  z = ~PCA3,
  color = ~Cluster_snmf,
  colors = ~color_palette_struct_snmf_dom_admixed_rm,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 5)
)

# Ajouter des titres et une légende
fig.2 <- fig.2 %>%
  layout(
    scene = list(
      xaxis = list(title = paste0('PC1 (6.66% explained)')),
      yaxis = list(title = 'PC2 (4.52% explained)'),
      zaxis = list(title = 'PC3 (3.05% explained)')
    ),
    title = '3D PCA Plot snmf with Color by pop struct',
    legend = list(title = 'Country Code')
  )

# # Afficher le graphique interactif
fig.2



fig.2 <- plot_ly(
  data = plotly_df_struct,
  x = ~PCA1,
  y = ~PCA2,
  z = ~PCA3,
  color = ~Cluster_FS,
  colors = ~color_palette_struct_FS_dom_admixed_rm,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 5)
)

# Ajouter des titres et une légende
fig.2 <- fig.2 %>%
  layout(
    scene = list(
      xaxis = list(title = paste0('PC1 (6.66% explained)')),
      yaxis = list(title = 'PC2 (4.52% explained)'),
      zaxis = list(title = 'PC3 (3.05% explained)')
    ),
    title = '3D PCA Plot FS with Color by pop struct',
    legend = list(title = 'Country Code')
  )

# # Afficher le graphique interactif
fig.2


# png_file <- paste0("all_syl_wo_ITA/faststructure_figure/PCA_K=7_EC_SYL.png")
# orca(fig.2, filename = png_file, width = 6, height = 6, units = "in", dpi = 300)
# 
# 
pca_ref_struct$Cluster_FS <- as.numeric(factor(pca_ref_struct$Cluster_FS, levels = unique(pca_ref_struct$Cluster_FS)))
# 
# png_file <- "all_syl_wo_ITA/faststructure_figure/PCA_K7_EC_SYL.png"
# png(png_file, width = 1600, height = 1600, res = 300)
# pdf_file <- "all_syl_wo_ITA/faststructure_figure/PCA_K7_EC_SYL.pdf"
# pdf(pdf_file, width = 6, height = 6) 
# 
scatterplot3d(pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3, color = color_palette_struct[pca_ref_struct$Cluster_FS], pch = 19,
              main = "3D PCA",
              xlab = paste('PC1 (6.66% explained)'),
              ylab = paste('PC2 (4.52% explained)'),
              zlab = paste('PC3 (3.05% explained)'),
              angle = 240)
legend("topright", legend = unique(pca_ref_struct$Cluster_FS),
       col = color_palette_struct,
       pch = 19)


pca_ref_struct$Cluster_snmf <- as.numeric(factor(pca_ref_struct$Cluster_snmf, levels = unique(pca_ref_struct$Cluster_snmf)))

scatterplot3d(pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3, color = color_palette_struct[pca_ref_struct$Cluster_snmf], pch = 19,
              main = "3D PCA",
              xlab = paste('PC1 (6.66% explained)'),
              ylab = paste('PC2 (4.52% explained)'),
              zlab = paste('PC3 (3.05% explained)'),
              angle = 280)
legend("topright", legend = unique(pca_ref_struct$Cluster_snmf),
       col = color_palette_struct,
       pch = 19)
# 
# dev.off()




#-----------
##### PCA All_data_wo_orientalis
#-----------
rm(list=ls())

library(rgl)
library(plotly)
library(scatterplot3d)
#--------
## PCA 3D K6
#--------

PROJECT="All_data_wo_orientalis"

OUTPUT_FOLDER=paste0('./out/',PROJECT,"/s02.PCA_figure")
dir.create(OUTPUT_FOLDER,recursive = T)

pca = read.delim(paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered.eigenvec"), header = F, sep = " ")
colnames(pca) <- c("Pop", "Seq_ID", paste0("PCA",1:20))
dim(pca)

ref_struct = read.delim(paste0("./in/population_structure_assignation_group_all_data_wo_orientalis_K6"))
dim(ref_struct)
# ref_struct <- subset(
#   ref_struct,
#   !(Cluster_FS == "V2 - Pdom" | Cluster_snmf == "V3 - P.dom")
# ); dim(ref_struct)

pca_ref_struct <- merge(pca, ref_struct, by="Seq_ID")

# color_palette_struct_snmf_all_data <- c("#d55e00ff","#d55e00ff","#d55e00ff", "#009e73ff", "#e69f00ff","#0072b2ff", "#56b4e9ff","#cc79a7ff","grey")
color_palette_struct_FS_all_data <- c("#d55e00ff","#d55e00ff","#d55e00ff", "#56b4e9ff", "red", "#009e73ff", "#e69f00ff","#0072b2ff", "#cc79a7ff","grey")

#color_palette_struct_snmf_dom_admixed_rm <- c("#009e73ff", "#e69f00ff","#0072b2ff", "#56b4e9ff","#cc79a7ff","grey")
#color_palette_struct_FS_dom_admixed_rm <- c("#56b4e9ff", "#009e73ff", "#e69f00ff","#0072b2ff", "#cc79a7ff","grey")

plotly_df_struct <- pca_ref_struct
plotly_df_struct$Cluster_FS <- factor(plotly_df_struct$Cluster_FS)
plotly_df_struct$Cluster_snmf <- factor(plotly_df_struct$Cluster_snmf)

fig.2 <- plot_ly(
  data = plotly_df_struct,
  x = ~PCA1,
  y = ~PCA2,
  z = ~PCA3,
  color = ~Cluster_FS,
  colors = ~color_palette_struct_FS_all_data,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 5)
)

# Ajouter des titres et une légende
fig.2 <- fig.2 %>%
  layout(
    scene = list(
      xaxis = list(title = ('PC1 (9.80% explained)')),
      yaxis = list(title = 'PC2 (8.21% explained)'),
      zaxis = list(title = 'PC3 (4.75% explained)')
    ),
    title = '3D PCA Plot snmf with Color by pop struct',
    legend = list(title = 'Country Code')
  )

# # Afficher le graphique interactif
fig.2

# Set the desired order of levels manually (modify this list as needed)
desired_levels <- c(
  "V1 - P.West_Eur",
  "V2 - Pdom",
  "V3 - P.Aus",
  "V4 - P.Est_Fr",
  "V5 - P.Rom",
  "V6 - P.Lor",
  "crops-wild admixed",
  "crops-wild admixed FS",
  "crops-wild admixed snmf",
  "wild-wild admixed"
)

# Color palette (must be in the same order as 'desired_levels')
color_palette_struct_FS_all_data_scatterplot <- c(
  "#56b4e9ff",        # V1 - P.West_Eur
  "#d55e00ff",  # V2 - Pdom
  "#009e73ff",  # V3 - P.Aus
  "#e69f00ff",       # V4 - P.Est_Fr
  "#0072b2ff",  # V5 - P.Rom
  "#cc79a7ff",  # V6 - P.Lor
  "#d7ad68",  # crops-wild admixed
  "#d7ad68",  # crops-wild admixed FS
  "#d7ad68",  # crops-wild admixed snmf
  "grey"   # wild-wild admixed
)

# Reorder factor levels
pca_ref_struct$Cluster_FS <- factor(pca_ref_struct$Cluster_FS, levels = desired_levels)

# Numeric index for color mapping
pca_ref_struct$Cluster_FS_num <- as.numeric(pca_ref_struct$Cluster_FS)

# Output file paths
png_file <- paste0(OUTPUT_FOLDER, "/PCA_3D_K6_", PROJECT, ".png")
pdf_file <- paste0(OUTPUT_FOLDER, "/PCA_3D_K6_", PROJECT, ".pdf")

# PNG output
png(png_file, width = 1600, height = 1600, res = 300)
scatterplot3d(
  pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3,
  color = color_palette_struct_FS_all_data_scatterplot[pca_ref_struct$Cluster_FS_num],
  pch = 19,
  main = "",
  xlab = 'PC1 (9.80% explained)',
  ylab = 'PC2 (8.21% explained)',
  zlab = 'PC3 (4.75% explained)',
  angle = 40
)
legend("topright",
       legend = desired_levels,
       col = color_palette_struct_FS_all_data_scatterplot,
       cex = 0.5,
       pch = 19,
       bg = "white")
dev.off()

# PDF output
pdf(pdf_file, width = 6, height = 6)
scatterplot3d(
  pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3,
  color = color_palette_struct_FS_all_data_scatterplot[pca_ref_struct$Cluster_FS_num],
  pch = 19,
  main = "3D PCA",
  xlab = 'PC1 (9.80% explained)',
  ylab = 'PC2 (8.21% explained)',
  zlab = 'PC3 (4.75% explained)',
  angle = 40
)
legend("topright",
       legend = desired_levels,
       col = color_palette_struct_FS_all_data_scatterplot,
       cex = 0.6,
       pch = 19,
       bg = "white")
dev.off()


library(ggplot2)

shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25)

# Ensure the group column is a factor with desired levels
pca_ref_struct$Cluster_FS <- factor(pca_ref_struct$Cluster_FS, levels = desired_levels)

# Name color and shape vectors
names(color_palette_struct_FS_all_data_scatterplot) <- desired_levels
names(shape_palette) <- desired_levels

# Plot
ggplot(pca_ref_struct, aes(x = PCA1, y = PCA2,
                           color = Cluster_FS,
                           shape = Cluster_FS,
                           fill = Cluster_FS)) +
  geom_point(size = 3, stroke = 0.5) +
  scale_color_manual(values = rep("black", length(desired_levels))) +
  scale_fill_manual(values = color_palette_struct_FS_all_data_scatterplot) +
  scale_shape_manual(values = shape_palette) +
  labs(
    title = "PCA_2D_FS_K6",
    x = "PC1 (9.80% explained)",
    y = "PC2 (8.21% explained)",
    color = "Group",
    shape = "Group",
    fill = "Group"
  ) +
  theme_bw() +
  theme(legend.position = "right")


#----------
# PCA 2D K5 #
#----------

rm(list=ls())

PROJECT="All_data_wo_orientalis"

OUTPUT_FOLDER=paste0('./out/',PROJECT,"/s02.PCA_figure")
dir.create(OUTPUT_FOLDER,recursive = T)

pca = read.delim(paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered.eigenvec"), header = F, sep = " ")
colnames(pca) <- c("Pop", "Seq_ID", paste0("PCA",1:20))
dim(pca)

ref_struct = read.delim(paste0("./in/population_structure_assignation_group_all_data_wo_orientalis_K5"))
dim(ref_struct)

pca_ref_struct <- merge(pca, ref_struct, by="Seq_ID")
# Set the desired order of levels manually (modify this list as needed)
desired_levels <- c(
  "V1 - P.West_Eur",
  "V2 - P.Dom",
  "V3 - P.Aus",
  "V4 - P.Est_Fr",
  "V5 - P.Rom",
  "crops-wild admixed",
  "crops-wild admixed FS",
  "crops-wild admixed snmf"
) # FS

# desired_levels <- c(
#   "V1 - P.West_Eur",
#   "V2 - Pdom",
#   "V3 - P.Aus",
#   "V4 - P.Est_Fr",
#   "V5 - P.Rom",
#   "V6 - P.Lor",
#   "crops-wild admixed",
#   "crops-wild admixed FS",
#   "crops-wild admixed snmf",
#   "wild-wild admixed"
# ) # snmf

# Color palette (must be in the same order as 'desired_levels')
color_palette_struct_FS_all_data_scatterplot <- c(
  "#56b4e9ff",        # V1 - P.West_Eur
  "#d55e00ff",  # V2 - Pdom
  "#009e73ff",  # V3 - P.Aus
  "#e69f00ff",       # V4 - P.Est_Fr
  "#0072b2ff",  # V5 - P.Rom
  "#d7ad68",  # crops-wild admixed
  "#d7ad68",  # crops-wild admixed FS
  "#d7ad68"  # crops-wild admixed snmf
)

# Reorder factor levels
pca_ref_struct$Cluster_FS <- factor(pca_ref_struct$Cluster_FS, levels = desired_levels)

# Numeric index for color mapping
pca_ref_struct$Cluster_FS_num <- as.numeric(pca_ref_struct$Cluster_FS)

# Output file paths
png_file <- paste0(OUTPUT_FOLDER, "/PCA_3D_K5_", PROJECT, ".png")
pdf_file <- paste0(OUTPUT_FOLDER, "/PCA_3D_K5_", PROJECT, ".pdf")
svg_file <- paste0(OUTPUT_FOLDER, "/PCA_3D_K5_", PROJECT, ".svg")
# PNG output
png(png_file, width = 1600, height = 1600, res = 300)
scatterplot3d(
  pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3,
  color = color_palette_struct_FS_all_data_scatterplot[pca_ref_struct$Cluster_FS_num],
  pch = 19,
  main = "",
  xlab = 'PC1 (9.80% explained)',
  ylab = 'PC2 (8.21% explained)',
  zlab = 'PC3 (4.75% explained)',
  angle = 40
)
legend("topright",
       legend = desired_levels,
       col = color_palette_struct_FS_all_data_scatterplot,
       cex = 0.5,
       pch = 19,
       bg = "white")
dev.off()

# PDF output
pdf(pdf_file, width = 6, height = 6)
scatterplot3d(
  pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3,
  color = color_palette_struct_FS_all_data_scatterplot[pca_ref_struct$Cluster_FS_num],
  pch = 19,
  main = "3D PCA",
  xlab = 'PC1 (9.80% explained)',
  ylab = 'PC2 (8.21% explained)',
  zlab = 'PC3 (4.75% explained)',
  angle = 40
)
legend("topright",
       legend = desired_levels,
       col = color_palette_struct_FS_all_data_scatterplot,
       cex = 0.6,
       pch = 19,
       bg = "white")
dev.off()

library(ggplot2)

shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23)

# Ensure the group column is a factor with desired levels
pca_ref_struct$Cluster_FS <- factor(pca_ref_struct$Cluster_FS, levels = desired_levels)

# Name color and shape vectors
names(color_palette_struct_FS_all_data_scatterplot) <- desired_levels
names(shape_palette) <- desired_levels

# Plot
p <- ggplot(pca_ref_struct, aes(x = PCA1, y = PCA2,
                                color = Cluster_FS,
                                shape = Cluster_FS,
                                fill = Cluster_FS)) +
  geom_point(size = 3, stroke = 0.5) +
  scale_color_manual(values = rep("black", length(desired_levels))) +
  scale_fill_manual(values = color_palette_struct_FS_all_data_scatterplot) +
  scale_shape_manual(values = shape_palette) +
  labs(
    title = "PCA_2D_FS_K5",
    x = "PC1",
    y = "PC2",
    color = "Group",
    shape = "Group",
    fill = "Group"
  ) +
  theme_bw() +
  theme(legend.position = "right")

# Save as SVG
ggsave(paste0(OUTPUT_FOLDER, "/PCA_2D_FS_K5_", PROJECT, ".svg"), plot = p, width = 6, height = 6)

# Save as PNG (300 dpi for publication-quality)
ggsave(paste0(OUTPUT_FOLDER, "/PCA_2D_FS_K5_", PROJECT, ".png"), plot = p, width = 6, height = 6, dpi = 300)


#----------
# PCA 2D K6 #
#----------

rm(list=ls())

PROJECT="All_data_wo_orientalis"

OUTPUT_FOLDER=paste0('./out/',PROJECT,"/s02.PCA_figure")
dir.create(OUTPUT_FOLDER,recursive = T)

pca = read.delim(paste0("./out/",PROJECT,"/s01.plink_pca/chr_concat.geno20_maf005.anno.syno.LDpruned_filtered.eigenvec"), header = F, sep = " ")
colnames(pca) <- c("Pop", "Seq_ID", paste0("PCA",1:20))
dim(pca)

ref_struct = read.delim(paste0("./in/population_structure_assignation_group_all_data_wo_orientalis_K6"))
dim(ref_struct)

pca_ref_struct <- merge(pca, ref_struct, by="Seq_ID")
# Set the desired order of levels manually (modify this list as needed)
desired_levels <- c(
  "V1 - P.West_Eur",
  "V2 - P.Dom",
  "V3 - P.Aus",
  "V4 - P.Est_Fr",
  "V5 - P.Rom",
  "V6 - P.Lor",
  "crops-wild admixed",
  "crops-wild admixed FS",
  "crops-wild admixed snmf",
  "wild-wild admixed"
) # FS

# desired_levels <- c(
#   "V1 - P.West_Eur",
#   "V2 - Pdom",
#   "V3 - P.Aus",
#   "V4 - P.Est_Fr",
#   "V5 - P.Rom",
#   "V6 - P.Lor",
#   "crops-wild admixed",
#   "crops-wild admixed FS",
#   "crops-wild admixed snmf",
#   "wild-wild admixed"
# ) # snmf

# Color palette (must be in the same order as 'desired_levels')
color_palette_struct_FS_all_data_scatterplot <- c(
  "#56b4e9ff",        # V1 - P.West_Eur
  "#d55e00ff",  # V2 - Pdom
  "#009e73ff",  # V3 - P.Aus
  "#e69f00ff",       # V4 - P.Est_Fr
  "#0072b2ff",  # V5 - P.Rom
  "#cc79a7ff",  # V6 - P.Lor
  "#d7ad68",  # crops-wild admixed
  "#d7ad68",  # crops-wild admixed FS
  "#d7ad68",  # crops-wild admixed snmf
  "grey"   # wild-wild admixed
)

# Reorder factor levels
pca_ref_struct$Cluster_FS <- factor(pca_ref_struct$Cluster_FS, levels = desired_levels)

# Numeric index for color mapping
pca_ref_struct$Cluster_FS_num <- as.numeric(pca_ref_struct$Cluster_FS)

# Output file paths
png_file <- paste0(OUTPUT_FOLDER, "/PCA_3D_K6_", PROJECT, ".png")
pdf_file <- paste0(OUTPUT_FOLDER, "/PCA_3D_K6_", PROJECT, ".pdf")
svg_file <- paste0(OUTPUT_FOLDER, "/PCA_3D_K6_", PROJECT, ".svg")
# PNG output
png(png_file, width = 1600, height = 1600, res = 300)
scatterplot3d(
  pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3,
  color = color_palette_struct_FS_all_data_scatterplot[pca_ref_struct$Cluster_FS_num],
  pch = 19,
  main = "",
  xlab = 'PC1 (9.80% explained)',
  ylab = 'PC2 (8.21% explained)',
  zlab = 'PC3 (4.75% explained)',
  angle = 40
)
legend("topright",
       legend = desired_levels,
       col = color_palette_struct_FS_all_data_scatterplot,
       cex = 0.5,
       pch = 19,
       bg = "white")
dev.off()

# PDF output
pdf(pdf_file, width = 6, height = 6)
scatterplot3d(
  pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3,
  color = color_palette_struct_FS_all_data_scatterplot[pca_ref_struct$Cluster_FS_num],
  pch = 19,
  main = "3D PCA",
  xlab = 'PC1 (9.80% explained)',
  ylab = 'PC2 (8.21% explained)',
  zlab = 'PC3 (4.75% explained)',
  angle = 40
)
legend("topright",
       legend = desired_levels,
       col = color_palette_struct_FS_all_data_scatterplot,
       cex = 0.6,
       pch = 19,
       bg = "white")
dev.off()

library(ggplot2)

shape_palette <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25)

# Ensure the group column is a factor with desired levels
pca_ref_struct$Cluster_FS <- factor(pca_ref_struct$Cluster_FS, levels = desired_levels)

# Name color and shape vectors
names(color_palette_struct_FS_all_data_scatterplot) <- desired_levels
names(shape_palette) <- desired_levels

# Plot
p <- ggplot(pca_ref_struct, aes(x = PCA1, y = PCA2,
                                color = Cluster_FS,
                                shape = Cluster_FS,
                                fill = Cluster_FS)) +
  geom_point(size = 3, stroke = 0.5) +
  scale_color_manual(values = rep("black", length(desired_levels))) +
  scale_fill_manual(values = color_palette_struct_FS_all_data_scatterplot) +
  scale_shape_manual(values = shape_palette) +
  labs(
    title = "PCA_2D_FS_K5",
    x = "PC1",
    y = "PC2",
    color = "Group",
    shape = "Group",
    fill = "Group"
  ) +
  theme_bw() +
  theme(legend.position = "right")

# Save as SVG
ggsave(paste0(OUTPUT_FOLDER, "/PCA_2D_FS_K6_", PROJECT, ".svg"), plot = p, width = 6, height = 6)

# Save as PNG (300 dpi for publication-quality)
ggsave(paste0(OUTPUT_FOLDER, "/PCA_2D_FS_K6_", PROJECT, ".png"), plot = p, width = 6, height = 6, dpi = 300)

