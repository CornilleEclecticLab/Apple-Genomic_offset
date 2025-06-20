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
      xaxis = list(title = paste0('PC1 (5.65% explained)')),
      yaxis = list(title = 'PC2 (2.9% explained)'),
      zaxis = list(title = 'PC3 (2.06% explained)')
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
      xaxis = list(title = paste0('PC1 (5.65% explained)')),
      yaxis = list(title = 'PC2 (2.9% explained)'),
      zaxis = list(title = 'PC3 (2.06% explained)')
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
              xlab = paste('PC1 (5.65% explained)'),
              ylab = paste('PC2 (2.9% explained)'),
              zlab = paste('PC3 (2.06% explained)'),
              angle = 280)
legend("topright", legend = unique(pca_ref_struct$Cluster_FS),
       col = color_palette_struct,
       pch = 19)


pca_ref_struct$Cluster_snmf <- as.numeric(factor(pca_ref_struct$Cluster_snmf, levels = unique(pca_ref_struct$Cluster_snmf)))

scatterplot3d(pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3, color = color_palette_struct[pca_ref_struct$Cluster_snmf], pch = 19,
              main = "3D PCA",
              xlab = paste('PC1 (5.65% explained)'),
              ylab = paste('PC2 (2.9% explained)'),
              zlab = paste('PC3 (2.06% explained)'),
              angle = 280)
legend("topright", legend = unique(pca_ref_struct$Cluster_snmf),
       col = color_palette_struct,
       pch = 19)
# 
# dev.off()



# ref = read.delim("all_syl_wo_ITA/ref_sample_all_syl_wo_ITA_plus_d-info_wo_admix", header = T, sep = "\t", dec = ",")
# ref_struct = read.delim("all_syl_wo_ITA/FS_wo_admix/EC_syl/structure_accession_for_all_K_thresh_0.95.txt", header = T, sep = ",")
# 
# pca_ref <- merge(pca, ref, by = "Seq_ID")
# pca_ref_struct <- merge(pca, ref_struct, by="Seq_ID")
# 
# color_palette <- c("green","yellow","blue","red")
# country_colors <- setNames(color_palette, c("FRA", "ROU", "DNK", "AUT"))
# 
# color_palette_struct <- c( "#0099FF", "#CC3333","grey","purple")#,"pink","purple","brown","black","lightblue","lightgreen","grey")
# struct_colors <- setNames(color_palette_struct, c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","admix"))
# 
# pca_ref$Color <- country_colors[pca_ref$Country.code..ISO.3166.]
# ## change the K depending of the plot
# pca_ref_struct$Color  <- struct_colors[pca_ref_struct$K7]
# 
# #scene <- plot3d(pca_ref[, c("PCA1", "PCA2", "PCA3")], col = pca_ref$Color, type = "p", radius = 1, size = 5)
# #rglwidget(scene3d(scene))
# 
# 
# plotly_df <- pca_ref
# plotly_df$Country.code..ISO.3166. <- factor(plotly_df$Country.code..ISO.3166.)
# 
# plotly_df_struct <- pca_ref_struct
# plotly_df_struct$K7 <- factor(plotly_df_struct$K7)
# 
# 
# fig.1 <- plot_ly(
#   data = plotly_df,
#   x = ~PCA1,
#   y = ~PCA2,
#   z = ~PCA3,
#   color = ~Country.code..ISO.3166.,
#   colors = ~color_palette,
#   type = 'scatter3d',
#   mode = 'markers',
#   marker = list(size = 5)
# )
# 
# # Ajouter des titres et une légende
# fig.1 <- fig.1 %>%
#   layout(
#     scene = list(
#       xaxis = list(title = 'PCA1'),
#       yaxis = list(title = 'PCA2'),
#       zaxis = list(title = 'PCA3')
#     ),
#     title = '3D PCA Plot with Color by Country',
#     legend = list(title = 'Country Code')
#   )
# 
# # Afficher le graphique interactif
# fig.1
# 
# 
# 
# 
# 
# fig.2 <- plot_ly(
#   data = plotly_df_struct,
#   x = ~PCA1,
#   y = ~PCA2,
#   z = ~PCA3,
#   color = ~K7,
#   colors = ~color_palette_struct,
#   type = 'scatter3d',
#   mode = 'markers',
#   marker = list(size = 5)
# )
# 
# # Ajouter des titres et une légende
# fig.2 <- fig.2 %>%
#   layout(
#     scene = list(
#       xaxis = list(title = paste0('PC1 (5.65% explained)')),
#       yaxis = list(title = 'PC2 (2.9% explained)'),
#       zaxis = list(title = 'PC3 (2.06% explained)')
#     ),
#     title = '3D PCA Plot with Color by pop struct',
#     legend = list(title = 'Country Code')
#   )
# 
# # Afficher le graphique interactif
# fig.2
# png_file <- paste0("all_syl_wo_ITA/faststructure_figure/PCA_K=7_EC_SYL.png")
# orca(fig.2, filename = png_file, width = 6, height = 6, units = "in", dpi = 300)
# 
# 
# pca_ref_struct$K7_numeric <- as.numeric(factor(pca_ref_struct$K7, levels = unique(pca_ref_struct$K7)))
# 
# png_file <- "all_syl_wo_ITA/faststructure_figure/PCA_K7_EC_SYL.png"
# png(png_file, width = 1600, height = 1600, res = 300)
# pdf_file <- "all_syl_wo_ITA/faststructure_figure/PCA_K7_EC_SYL.pdf"
# pdf(pdf_file, width = 6, height = 6) 
# 
# scatterplot3d(pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3, color = color_palette_struct[pca_ref_struct$K7_numeric], pch = 19, 
#               main = "3D PCA",
#               xlab = paste('PC1 (5.65% explained)'),
#               ylab = paste('PC2 (2.9% explained)'),
#               zlab = paste('PC3 (2.06% explained)'),
#               angle = 70)
# legend("topright", legend = unique(pca_ref_struct$K7),
#        col = color_palette_struct,
#        pch = 19)
# 
# dev.off()