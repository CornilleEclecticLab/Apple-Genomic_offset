rm(list = ls())

library(rgl)
library(plotly)
library(scatterplot3d)

PROJECT="tips_dom_admixed_rm_V3"
IN_FOLD=paste0("./out/", PROJECT, "/s06.plink_distance/")
OUT_FOLD=paste0("./out/", PROJECT, "/s06_5.PCA/")
METHOD <- "FS"
OUTPUT_NAME_3D="/PCA_3D_landscape_geno_K5_"
OUTPUT_NAME_2D="/PCA_2D_landscape_geno_K5_"
dir.create(OUT_FOLD,recursive = T)

pca = read.delim(paste0(IN_FOLD,"GraffiTE.merged.genotypes_",PROJECT,".eigenvec"), header = F, sep = " ")[,-1]
colnames(pca) <- c("Seq_ID", paste0("PCA",1:20))

ref = read.delim("./in/tips_dom_admixed_rm_info_sample.txt", header = T, sep = "\t", dec = ",")
ref_struct = read.delim("./in/Global_cluster_pure_wild.txt"); dim(ref_struct)
pca_ref_struct <- merge(pca, ref_struct, by="Seq_ID")

desired_levels <- c(
  "P.Lor",
  "P.Rom",
  "P.West_Eur",
  "P.Aus",
  "P.Ita",
  "wild-wild admixed"
)
color_palette_struct <- c(
  "#0072b2ff",
  "#E69F00",
  "#56b4e9ff",
  "#D55E00",
  "#009e73ff",
  "grey"   # wild-wild admixed
)
shape_palette <- c(21, 22, 23, 24, 25, 21)

## PCA 3D

# Reorder factor levels
cluster <- "Global_Cluster"
pca_ref_struct[[cluster]] <- factor(pca_ref_struct[[cluster]], levels = desired_levels)

# Numeric index for color mapping
pca_ref_struct$Cluster_num <- as.numeric(pca_ref_struct[[cluster]])

# Output file paths
png_file <- paste0(OUT_FOLD, OUTPUT_NAME_3D, METHOD,"_", PROJECT, ".png")
pdf_file <- paste0(OUT_FOLD, OUTPUT_NAME_3D, METHOD, "_", PROJECT, ".pdf")


png(png_file, width = 1600, height = 1600, res = 300)
scatterplot3d(
  pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3,
  color = color_palette_struct[pca_ref_struct$Cluster_num],
  pch = 19,
  main = "",
  xlab = 'PC1 (7.89% explained)',
  ylab = 'PC2 (5.99% explained)',
  zlab = 'PC3 (5.41% explained)',
  angle = 40
)
legend("topright",
       legend = desired_levels,
       col = color_palette_struct,
       cex = 0.5,
       pch = 19,
       bg = "white")
dev.off()

# PDF output
pdf(pdf_file, width = 6, height = 6)
scatterplot3d(
  pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3,
  color = color_palette_struct[pca_ref_struct$Cluster_num],
  pch = 19,
  main = "3D PCA",
  xlab = 'PC1 (7.89% explained)',
  ylab = 'PC2 (5.99% explained)',
  zlab = 'PC3 (5.41% explained)',
  angle = 40
)
legend("topright",
       legend = desired_levels,
       col = color_palette_struct,
       cex = 0.6,
       pch = 19,
       bg = "white")
dev.off()


## PCA 2D

# Ensure the group column is a factor with desired levels
pca_ref_struct <- merge(pca, ref_struct, by="Seq_ID")
pca_ref_struct[[cluster]] <- factor(pca_ref_struct[[cluster]], levels = desired_levels)

# Name color and shape vectors
names(color_palette_struct) <- desired_levels
names(shape_palette) <- desired_levels

# Plot
p <- ggplot(pca_ref_struct, aes(
  x = PCA1,
  y = PCA2,
  color = .data[[cluster]],
  shape = .data[[cluster]],
  fill = .data[[cluster]]
)) +
  geom_point(size = 3, stroke = 0.5) +
  scale_color_manual(values = rep("black", length(desired_levels))) +
  scale_fill_manual(values = color_palette_struct) +
  scale_shape_manual(values = shape_palette) +
  labs(
    title = paste0(OUTPUT_NAME_2D, METHOD),
    x = "PC1 (7.89% explained)",
    y = "PC2 (5.99% explained)",
    color = "Group",
    shape = "Group",
    fill = "Group"
  ) +
  theme_bw() +
  theme(legend.position = "right")

# Save as SVG
ggsave(paste0(OUT_FOLD, OUTPUT_NAME_2D, METHOD,"_", PROJECT, ".svg"), plot = p, width = 6, height = 6)

# Save as PNG (300 dpi for publication-quality)
ggsave(paste0(OUT_FOLD, OUTPUT_NAME_2D, METHOD,"_", PROJECT, ".png"), plot = p, width = 6, height = 6, dpi = 300)

