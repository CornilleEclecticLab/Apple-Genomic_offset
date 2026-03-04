rm(list = ls())
source(file = "bin/s02.load_function.R")

dir.create(OUTPUT_FOLDER,recursive = T)
pca = read.delim(PCA_file, header = TRUE, sep = "\t")
colnames(pca) <- c("Pop", "Seq_ID", paste0("PCA",1:(length(colnames(pca))-2)))
# colnames(pca) <- c("Seq_ID", paste0("PCA",1:(length(colnames(pca))-1)))
dim(pca)
ref_struct = read.csv(ref_struct_file); dim(ref_struct)
pca_ref_struct <- merge(pca, ref_struct, by="Seq_ID")

## PCA 3D

# Reorder factor levels
cluster <- "cluster"
pca_ref_struct[[cluster]] <- factor(pca_ref_struct[[cluster]], levels = desired_levels)

# Numeric index for color mapping
pca_ref_struct$Cluster_num <- as.numeric(pca_ref_struct[[cluster]])

# Output file paths
png_file <- paste0(OUTPUT_FOLDER, OUTPUT_NAME_3D, METHOD,"_", PROJECT, ".png")
pdf_file <- paste0(OUTPUT_FOLDER, OUTPUT_NAME_3D, METHOD, "_", PROJECT, ".pdf")


png(png_file, width = 1600, height = 1600, res = 300)
scatterplot3d(
  pca_ref_struct$PCA1, pca_ref_struct$PCA2, pca_ref_struct$PCA3,
  color = color_palette_struct[pca_ref_struct$Cluster_num],
  pch = 19,
  main = "",
  xlab = 'PC1 (9.80% explained)',
  ylab = 'PC2 (8.21% explained)',
  zlab = 'PC3 (4.75% explained)',
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
  xlab = 'PC1 (9.80% explained)',
  ylab = 'PC2 (8.21% explained)',
  zlab = 'PC3 (4.75% explained)',
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
    x = "PC1",
    y = "PC2",
    color = "Group",
    shape = "Group",
    fill = "Group"
  ) +
  theme_bw() +
  theme(legend.position = "right")

# Save as SVG
ggsave(paste0(OUTPUT_FOLDER, OUTPUT_NAME_2D, METHOD,"_", PROJECT, ".svg"), plot = p, width = 6, height = 6)

# Save as PNG (300 dpi for publication-quality)
ggsave(paste0(OUTPUT_FOLDER, OUTPUT_NAME_2D, METHOD,"_", PROJECT, ".png"), plot = p, width = 6, height = 6, dpi = 300)

