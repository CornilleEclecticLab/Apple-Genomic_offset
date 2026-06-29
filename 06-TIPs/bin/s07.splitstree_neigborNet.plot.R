rm(list = ls())

library(vcfR)
library(adegenet)
#install.packages("poppr")
library(poppr)
library(phangorn)
library(RColorBrewer)


PROJECT <- "tips_dom_admixed_rm_V3"
OUT_FOLD <- paste0("./out/",PROJECT,"/s07.splitstree/")
VCF <- paste0("GraffiTE.merged.genotypes_",PROJECT,".vcf.gz")

dir.create(OUT_FOLD, recursive = TRUE)


filename <- basename(sub("\\.vcf.gz$", "", VCF))

vcf <- read.vcfR(
  paste0("./out/",PROJECT,"/s01.filter_vcf/",VCF)
)
genlight_obj <- vcfR2genlight(vcf)
dist <- bitwise.dist(genlight_obj)
dist_matrix <- as.matrix(dist)
n_samples <- nrow(dist_matrix)

# Create the Nexus-style header
first_line <- paste0("\t", n_samples, paste(rep("\t", n_samples - 1), collapse = ""))
outfile <- paste0(OUT_FOLD,filename,".bitwise.dist")

# Write the first line
file_conn <- file(outfile, "w")
writeLines(first_line, file_conn)
close(file_conn)

# Now write the matrix with proper format
write.table(dist_matrix,
            file = outfile,
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = TRUE,  # Needed for Nex output
            append = TRUE)



# open dist by splitstree 4
# save it with nex file


net <- read.nexus.networx(
  paste0(OUT_FOLD, filename,".bitwise.dist.nex")
)

group_ind <- read.table(
  paste0(OUT_FOLD, "Group_genetic.ind.Global_cluster_pure_wild.tsv"),
  header = TRUE, stringsAsFactors = FALSE, sep="\t"
)
group_shape <- read.table(
   paste0(OUT_FOLD, "Group_genetic.shape.Global_cluster_pure_wild.tsv"),
  header = TRUE, stringsAsFactors = FALSE, sep="\t"
)

group_color <- read.table(
   paste0(OUT_FOLD, "Group_genetic.color.Global_cluster_pure_wild.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)


####

# Map your shapes (1–4) to filled versions (21–25)
#group_shape$shape <- c(21, 22, 23, 24, 25)  # circle, square, diamond, triangle

# Merge color and shape info again
group_info <- merge(group_shape, group_color, by = "Group")

tip_info <- merge(group_ind, group_info, by = "Group")
rownames(tip_info) <- tip_info$Ind

tip_pch <- tip_info[net$tip.label, "shape"]
tip_col <- tip_info[net$tip.label, "Color_HEX"]

# Add one NA at start for phangorn indexing
tip_pch_add <- c(NA, tip_pch)
tip_col_add <- c(NA, tip_col)

png(paste0(OUT_FOLD, filename, ".bitwise.dist.neighborNet.color.png"),
    width = 8, height = 8, units = "in", res = 300)

plot(net,
     edge.width = 0.6,
     cex = 0.5,
     show.tip.label = FALSE,
     angle = 0
)

# Use filled points
tiplabels(
  pch = as.numeric(tip_pch_add),
  cex = 0.8,
  bg = tip_col_add,
  col = "black"
)

legend("bottomleft",
       legend = group_info$Group,
       pch = group_info$shape,
       pt.bg = group_info$Color_HEX,
       col = "black",
       ncol = 2,
       bty = "n",
       title = "Group"
)

dev.off()

svg(paste0(OUT_FOLD, filename, ".bitwise.dist.neighborNet.color.svg"),
    width = 8, height = 8
)

plot(net,
     edge.width = 0.6,
     cex = 0.5,
     show.tip.label = FALSE,
     angle = 0
)

# Use filled points
tiplabels(
  pch = as.numeric(tip_pch_add),
  cex = 0.8,
  bg = tip_col_add,
  col = "black"
)

legend("bottomleft",
       legend = group_info$Group,
       pch = group_info$shape,
       pt.bg = group_info$Color_HEX,
       col = "black",
       ncol = 2,
       bty = "n",
       title = "Group"
)

dev.off()
