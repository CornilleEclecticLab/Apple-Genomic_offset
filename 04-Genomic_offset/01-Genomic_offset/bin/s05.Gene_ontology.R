rm(list = ls())
# conda activate r-go-env

# ==================================================
# 0. Load libraries
# ==================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(Biostrings)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(ggplot2)
  library(GOSemSim)
  library(enrichplot)
  library(ggpubr)
  library(stringr)
  library(ggupset)
  library(IRanges)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
})

# install.packages("./in/org.MdomesticaGDDH13.eg.db", repos=NULL, type="sources")
# install.packages('GOplot')

library(org.Mdomestica.eg.db)
library(GOplot)

# ==================================================
# 1. Read and merge SNP position files
# ==================================================
snp_dir <- "./out/dom_admixed_rm/s02.lfmm2/file/"

snp_files <- list.files(
  path      = snp_dir,
  pattern   = "chr_chr.*_SNPs_associated_with_all_bio_VCF_",
  recursive = TRUE,
  full.names = TRUE
)

stopifnot(length(snp_files) > 0)

# Helper to read and tag a file
read_file <- function(file, method) {
  chr <- str_extract(file, "chr_chr[0-9]+") %>%
    sub("^chr_", "", .) %>%
    sub("^chr", "Chr", .)

  read_csv(file, show_col_types = FALSE) %>%
    transmute(
      SNPs       = SNPs,
      chromosome = chr,
      !!method   := TRUE
    )
}

# Separate files
bonf_files <- snp_files[str_detect(snp_files, "Bonferroni")]
fdr_files  <- snp_files[str_detect(snp_files, "FDR")]

# Read all files
bonf_df <- map_dfr(bonf_files, read_file, method = "Bonferroni")
fdr_df  <- map_dfr(fdr_files,  read_file, method = "FDR")

# Merge Bonferroni + FDR
all_snps_pos <- full_join(
  bonf_df,
  fdr_df,
  by = c("SNPs", "chromosome")
) %>%
  mutate(
    Bonferroni = ifelse(is.na(Bonferroni), FALSE, Bonferroni),
    FDR        = ifelse(is.na(FDR), FALSE, FDR)
  ) %>%
  arrange(chromosome, SNPs)

# ==================================================
# 2. Read and parse GFF3 annotation
# ==================================================
gff <- read.table(
  "./in/gene_models_20170612.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE,
  fill = TRUE
)

colnames(gff) <- c(
  "chr", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

gff <- gff %>%
  mutate(
    chr = sub("^chr", "Chr", chr),
    gene_id = case_when(
      type == "gene" ~ sub(".*ID=gene:([^;]+).*", "\\1", attributes),
      TRUE           ~ sub(".*Parent=gene:([^;]+).*", "\\1", attributes)
    )
  )

# ==================================================
# 3. Build genomic ranges for overlap
# ==================================================
gff_ranges <- gff %>%
  dplyr::select(chr, start, end, type, gene_id) %>%
  group_by(chr) %>%
  summarise(
    ranges  = list(IRanges(start, end)),
    type    = list(type),
    gene_id = list(gene_id),
    .groups = "drop"
  )

# ==================================================
# 4. Overlap SNPs with GFF features
# ==================================================
find_snps_in_gff <- function(snps_df, gff_chr) {

  snp_ir <- IRanges(start = snps_df$SNPs, end = snps_df$SNPs)
  hits   <- findOverlaps(snp_ir, gff_chr$ranges[[1]])

  if (!length(hits)) return(NULL)

  tibble(
    SNPs          = snps_df$SNPs[queryHits(hits)],
    chromosome    = snps_df$chromosome[queryHits(hits)],
    pvalues       = snps_df$pvalues[queryHits(hits)],
    fscores       = snps_df$fscores[queryHits(hits)],
    adj.r.squared = snps_df$adj.r.squared[queryHits(hits)],
    Bonferroni    = snps_df$Bonferroni[queryHits(hits)],
    FDR           = snps_df$FDR[queryHits(hits)],
    gff_type      = gff_chr$type[[1]][subjectHits(hits)],
    gene_id       = gff_chr$gene_id[[1]][subjectHits(hits)]
  )
}


snps_in_gff <- bind_rows(lapply(unique(all_snps_pos$chromosome), function(chr_name) {

  chr_snps <- dplyr::filter(all_snps_pos, chromosome == chr_name)
  chr_gff  <- dplyr::filter(gff_ranges, chr == chr_name)

  if (nrow(chr_gff) == 0) return(NULL)

  find_snps_in_gff(chr_snps, chr_gff)
}))

# ==================================================
# 6. Define GO universe from FASTA
# ==================================================
fasta <- readDNAStringSet("./in/GDDH13_1-1_mrna.fasta")

universe_genes <- unique(
  sub(".*gene[:=]([^;| ]+).*", "\\1", names(fasta))
)

length(universe_genes)

# ==================================================
# 5. Extract candidate genes
# ==================================================
out_dir <- "./out/dom_admixed_rm/s05.Gene_ontology/by_chr/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

semdata <- godata(  OrgDb = org.Mdomestica.eg.db,
                   keytype = "GID",
                   ont = "BP"
)

df_go_bonf <- tibble()
df_go_fdr  <- tibble()

for (chr in unique(snps_in_gff$chromosome)) {

  message("chr: ", chr)

  # ======================
  # Bonferroni
  # ======================
  gene_list_bonf <- snps_in_gff %>%
    filter(
      gff_type == "gene",
      chromosome == chr,
      Bonferroni == TRUE
    ) %>%
    pull(gene_id) %>%
    unique()

  if (length(gene_list_bonf) > 0) {

    go_bonf <- enrichGO(
      gene          = gene_list_bonf,
      universe      = universe_genes,
      OrgDb         = org.Mdomestica.eg.db,
      keyType       = "GID",
      ont           = "BP",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
    #go_bonf_2 <- simplify(go_bonf, semData = semdata)

    if (!is.null(go_bonf) && nrow(as.data.frame(go_bonf)) > 0) {

      ggsave(
        filename = file.path(out_dir, paste0("GO_BP_dotplot_Bonf_", chr, ".png")),
        plot     = dotplot(go_bonf, showCategory = 15),
        width    = 8, height = 6, dpi = 300
      )

      df_go_bonf <- bind_rows(
        df_go_bonf,
        as.data.frame(go_bonf) %>%
          mutate(chromosome = chr, correction = "Bonferroni")
      )
    }
  }

  if (nrow(as.data.frame(go_bonf)) >= 3) {

    ego_bonf_sim <- pairwise_termsim(
      go_bonf,
      method  = "Wang",
      semData = semdata
    )

    # Enrichment map
    ggsave(
      filename = file.path(out_dir, paste0("GO_BP_emap_Bonf_", chr, ".png")),
      plot     = emapplot(ego_bonf_sim, showCategory = 30),
      width    = 14,
      height   = 11,
      dpi      = 300
    )

    # Tree plot
    n_terms <- nrow(as.data.frame(ego_bonf_sim))
    n_clusters <- min(5, n_terms)

    ggsave(
      filename = file.path(out_dir, paste0("GO_BP_tree_Bonf_", chr, ".png")),
      plot     = treeplot(
        ego_bonf_sim,
        cluster.params = list(n = n_clusters)
      ),
      width    = 14,
      height   = 11,
      dpi      = 300
    )

  }

  # ======================
  # FDR
  # ======================
  gene_list_fdr <- snps_in_gff %>%
    filter(
      gff_type == "gene",
      chromosome == chr,
      FDR == TRUE
    ) %>%
    pull(gene_id) %>%
    unique()

  if (length(gene_list_fdr) > 0) {

    go_fdr <- enrichGO(
      gene          = gene_list_fdr,
      universe      = universe_genes,
      OrgDb         = org.Mdomestica.eg.db,
      keyType       = "GID",
      ont           = "BP",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
    #go_fdr_2 <- simplify(go_fdr, semData = semdata)

    if (!is.null(go_fdr) && nrow(as.data.frame(go_fdr)) > 0) {

      ggsave(
        filename = file.path(out_dir, paste0("GO_BP_dotplot_FDR_", chr, ".png")),
        plot     = dotplot(go_fdr, showCategory = 15),
        width    = 8, height = 6, dpi = 300
      )

      df_go_fdr <- bind_rows(
        df_go_fdr,
        as.data.frame(go_fdr) %>%
          mutate(chromosome = chr, correction = "FDR")
      )
    }
  }

  if (nrow(as.data.frame(go_fdr)) >= 3) {

    ego_fdr_sim <- pairwise_termsim(
      go_fdr,
      method  = "Wang",
      semData = semdata
    )

    # Enrichment map
    ggsave(
      filename = file.path(out_dir, paste0("GO_BP_emap_FDR_", chr, ".png")),
      plot     = emapplot(ego_fdr_sim, showCategory = 30),
      width    = 14,
      height   = 11,
      dpi      = 300
    )

    # Tree plot
    n_terms <- nrow(as.data.frame(ego_fdr_sim))
    n_clusters <- min(5, n_terms)

    ggsave(
      filename = file.path(out_dir, paste0("GO_BP_tree_Bonf_", chr, ".png")),
      plot     = treeplot(
        ego_fdr_sim,
        cluster.params = list(n = n_clusters)
      ),
      width    = 14,
      height   = 11,
      dpi      = 300
    )
  }

}

dir.create("./out/dom_admixed_rm/s05.Gene_ontology/file/", recursive=T)
write.csv(df_go_bonf,
          "./out/dom_admixed_rm/s05.Gene_ontology/file/Bonferroni_by_chr_Geno_Ontology.csv",
          row.names = F)
write.csv(df_go_fdr,
          "./out/dom_admixed_rm/s05.Gene_ontology/file/FDR_by_chr_Geno_Ontology.csv",
          row.names = F)


# ========================
### full genome correction
# ========================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(Biostrings)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(ggplot2)
  library(GOSemSim)
  library(enrichplot)
  library(ggpubr)
  library(stringr)
  library(ggupset)
  library(IRanges)
  library(readr)
})

library(org.Mdomestica.eg.db)

# =========================
# Input files
# =========================
bonf_file <- "./out/dom_admixed_rm/s02.lfmm2/file/All_chr_SNPs_associated_with_all_bio_VCF_position_Bonferroni.csv"
fdr_file  <- "./out/dom_admixed_rm/s02.lfmm2/file/All_chr_SNPs_associated_with_all_bio_VCF_position_FDR.csv"

bonf_df <- read_csv(bonf_file, show_col_types = FALSE) %>%
  mutate(Bonferroni = TRUE)

fdr_df <- read_csv(fdr_file, show_col_types = FALSE) %>%
  mutate(FDR = TRUE)

final_df <- full_join(
  bonf_df[, c("SNPs", "chromosome", "pvalues", "fscores", "adj.r.squared", "Bonferroni")],
  fdr_df[, c("SNPs", "chromosome", "FDR")],
  by = c("SNPs", "chromosome")
) %>%
  mutate(
    Bonferroni = ifelse(is.na(Bonferroni), FALSE, Bonferroni),
    FDR        = ifelse(is.na(FDR), FALSE, FDR),
    chromosome = sub("^chr", "Chr", chromosome)
  )

# =========================
# GFF parsing
# =========================
gff <- read.table(
  "./in/gene_models_20170612.gff3",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE,
  fill = TRUE
)

colnames(gff) <- c(
  "chr", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

gff <- gff %>%
  mutate(
    chr = sub("^chr", "Chr", chr),
    gene_id = case_when(
      type == "gene" ~ sub(".*ID=gene:([^;]+).*", "\\1", attributes),
      TRUE           ~ sub(".*Parent=gene:([^;]+).*", "\\1", attributes)
    )
  )

gff_ranges <- gff %>%
  dplyr::select(chr, start, end, type, gene_id) %>%
  group_by(chr) %>%
  summarise(
    ranges  = list(IRanges(start, end)),
    type    = list(type),
    gene_id = list(gene_id),
    .groups = "drop"
  )

find_snps_in_gff <- function(snps_df, gff_chr) {

  snp_ir <- IRanges(start = snps_df$SNPs, end = snps_df$SNPs)
  hits   <- findOverlaps(snp_ir, gff_chr$ranges[[1]])

  if (!length(hits)) return(NULL)

  tibble(
    SNPs          = snps_df$SNPs[queryHits(hits)],
    chromosome    = snps_df$chromosome[queryHits(hits)],
    pvalues       = snps_df$pvalues[queryHits(hits)],
    fscores       = snps_df$fscores[queryHits(hits)],
    adj.r.squared = snps_df$adj.r.squared[queryHits(hits)],
    Bonferroni    = snps_df$Bonferroni[queryHits(hits)],
    FDR           = snps_df$FDR[queryHits(hits)],
    gff_type      = gff_chr$type[[1]][subjectHits(hits)],
    gene_id       = gff_chr$gene_id[[1]][subjectHits(hits)]
  )
}

snps_in_gff <- bind_rows(lapply(unique(final_df$chromosome), function(chr_name) {

  chr_snps <- filter(final_df, chromosome == chr_name)
  chr_gff  <- filter(gff_ranges, chr == chr_name)

  if (nrow(chr_gff) == 0) return(NULL)

  find_snps_in_gff(chr_snps, chr_gff)
}))

# =========================
# Gene lists
# =========================
gene_list_bonf <- snps_in_gff %>%
  filter(gff_type == "gene", Bonferroni) %>%
  pull(gene_id) %>%
  unique()

gene_list_fdr <- snps_in_gff %>%
  filter(gff_type == "gene", FDR) %>%
  pull(gene_id) %>%
  unique()

# =========================
# Universe genes
# =========================
fasta <- readDNAStringSet("./in/GDDH13_1-1_mrna.fasta")

universe_genes <- unique(
  sub(".*gene[:=]([^;| ]+).*", "\\1", names(fasta))
)

# =========================
# GO enrichment
# =========================
go_bp_bonf <- enrichGO(
  gene          = gene_list_bonf,
  universe      = universe_genes,
  OrgDb         = org.Mdomestica.eg.db,
  keyType       = "GID",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_bp_fdr <- enrichGO(
  gene          = gene_list_fdr,
  universe      = universe_genes,
  OrgDb         = org.Mdomestica.eg.db,
  keyType       = "GID",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# =========================
# Semantic data
# =========================
semdata <- godata(
  OrgDb  = org.Mdomestica.eg.db,
  keytype = "GID",
  ont     = "BP"
)

safe_pairwise_termsim <- function(ego, semdata) {

  if (is.null(ego)) return(NULL)
  df <- as.data.frame(ego)
  if (nrow(df) < 3) return(NULL)

  out <- tryCatch(
    pairwise_termsim(ego, method = "Wang", semData = semdata),
    error = function(e) NULL
  )

  if (is.null(out)) return(NULL)
  if (is.null(out@termsim)) return(NULL)
  if (is.null(dimnames(out@termsim))) return(NULL)

  out
}

# =========================
# Output
# =========================
out_dir <- "./out/dom_admixed_rm/s05.Gene_ontology/full_genome/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Bonferroni ----
if (!is.null(go_bp_bonf) && nrow(as.data.frame(go_bp_bonf)) > 0) {

  ggsave(
    file.path(out_dir, "GO_BP_Bonferroni_dotplot.png"),
    dotplot(go_bp_bonf, showCategory = 15),
    width = 8, height = 6, dpi = 300
  )

  ego_bonf_sim <- safe_pairwise_termsim(go_bp_bonf, semdata)

  if (!is.null(ego_bonf_sim)) {

    ggsave(
      file.path(out_dir, "GO_BP_Bonferroni_emap.png"),
      emapplot(ego_bonf_sim, showCategory = 30),
      width = 14, height = 11, dpi = 300
    )

    n_terms <- nrow(as.data.frame(ego_bonf_sim))
    n_clusters <- min(5, n_terms)

    ggsave(
      file.path(out_dir, "GO_BP_Bonferroni_tree.png"),
      treeplot(
        ego_bonf_sim,
        cluster.params = list(n = n_clusters)
      ),
      width = 14, height = 11, dpi = 300
    )
  }
}

# ---- FDR ----
if (!is.null(go_bp_fdr) && nrow(as.data.frame(go_bp_fdr)) > 0) {

  suppressWarnings(
    ggsave(
      file.path(out_dir, "GO_BP_FDR_dotplot.png"),
      dotplot(go_bp_fdr, showCategory = 15),
      width = 8, height = 6, dpi = 300
    )
  )

  ego_fdr_sim <- safe_pairwise_termsim(go_bp_fdr, semdata)

  if (!is.null(ego_fdr_sim)) {

    suppressWarnings(
      ggsave(
        file.path(out_dir, "GO_BP_FDR_emap.png"),
        emapplot(ego_fdr_sim, showCategory = 30),
        width = 14, height = 11, dpi = 300
      )
    )

    n_terms <- nrow(as.data.frame(ego_fdr_sim))
    n_clusters <- min(5, n_terms)

    suppressWarnings(
      ggsave(
        file.path(out_dir, "GO_BP_FDR_tree.png"),
        treeplot(
          ego_fdr_sim,
          cluster.params = list(n = n_clusters)
        ),
        width = 14, height = 11, dpi = 300
      )
    )
  }
}

