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

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)
# install.packages("./in/org.MdomesticaGDDH13.eg.db", repos=NULL, type="sources")
# install.packages('GOplot')

library(org.Mdomestica.eg.db)
library(GOplot)


## objectif: Look at SNPs (FDR) position in the genome and count how many are in :
# - intergenitic region
# - genes
# - 2 kb before a gene

# --- 1. Settings and Data Loading ---
PROJECT <- "dom_admixed_rm"
snp_path <- paste0("./out/", PROJECT, "/s02.lfmm2_V3/file/All_SNPs_Filtered_Results.csv")
gff_path <- "./in/gene_models_20170612.gff3"

data_snp <- read.csv(snp_path) %>%
  filter(FDR_signif == TRUE) %>%
  mutate(chromosome = sub("^chr", "Chr", chromosome))

gff_raw <- read.table(
  gff_path, sep = "\t", header = FALSE, 
  comment.char = "#", quote = "", stringsAsFactors = FALSE, fill = TRUE
)

colnames(gff_raw) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# --- 2. Prepare Annotation Ranges ---
# Standardize GFF and extract Gene IDs
gff_base <- gff_raw %>%
  mutate(
    chr = sub("^chr", "Chr", chr),
    gene_id = case_when(
      type == "gene" ~ sub(".*ID=gene:([^;]+).*", "\\1", attributes),
      TRUE           ~ sub(".*Parent=gene:([^;]+).*", "\\1", attributes)
    )
  )

# Create the 2kb Upstream regions based on strand
# --- 1. Préparation des zones de distance (Buffers) ---

# Zone 2kb (on réutilise l'existant mais on l'élargit potentiellement aux deux côtés si besoin)
# Ici, on va créer des zones "Flanking" (amont et aval) pour être complet
flanking_2kb <- gff_base %>%
  filter(type == "gene") %>%
  mutate(
    start_2 = pmax(1, start - 2000),
    end_2   = end + 2000,
    type    = "Intergenic_0_2kb"
  ) %>%
  dplyr::select(chr, start = start_2, end = end_2, type, gene_id)

# Zone 5kb
flanking_5kb <- gff_base %>%
  filter(type == "gene") %>%
  mutate(
    start_5 = pmax(1, start - 5000),
    end_5   = end + 5000,
    type    = "Intergenic_2_5kb"
  ) %>%
  dplyr::select(chr, start = start_5, end = end_5, type, gene_id)

# --- 2. Mise à jour de la hiérarchie ---
# On combine les régions internes et les nouvelles zones de distance
all_reference_regions <- bind_rows(
  gff_base %>% dplyr::select(chr, start, end, type, gene_id),
  flanking_2kb,
  flanking_5kb
)

# --- 3. Optimized Overlap Function ---
# We use a simple function that takes a SNP dataframe and a Reference dataframe
perform_overlap <- function(snps, refs) {
  # Group by chromosome to speed up IRanges
  all_hits <- lapply(intersect(unique(snps$chromosome), unique(refs$chr)), function(c) {
    s_sub <- snps %>% filter(chromosome == c)
    r_sub <- refs %>% filter(chr == c)
    
    if(nrow(s_sub) == 0 | nrow(r_sub) == 0) return(NULL)
    
    s_ir <- IRanges(start = s_sub$SNPs, end = s_sub$SNPs)
    r_ir <- IRanges(start = r_sub$start, end = r_sub$end)
    
    hits <- findOverlaps(s_ir, r_ir)
    
    if(length(hits) == 0) return(NULL)
    
    # Merge SNP data with the matched reference type
    cbind(s_sub[queryHits(hits), ], 
          ref_type = r_sub$type[subjectHits(hits)], 
          gene_id  = r_sub$gene_id[subjectHits(hits)])
  })
  bind_rows(all_hits)
}

# Run the overlap
raw_overlaps <- perform_overlap(data_snp, all_reference_regions)

# --- 4. Priority Ranking and Final Classification ---
# Nouvelle priorité : 
# On veut d'abord savoir si c'est DANS le gène (CDS, UTR, etc.)
# PUIS si c'est dans les 2kb, PUIS dans les 5kb.
priority_map <- c(
  "CDS"              = 1, 
  "five_prime_UTR"   = 2, 
  "three_prime_UTR"  = 3, 
  "exon"             = 4, 
  "ncRNA"            = 5,
  "mRNA"             = 6, 
  "gene"             = 7, 
  "Intergenic_0_2kb" = 8,
  "Intergenic_2_5kb" = 9
)

# --- 3. Exécution de l'Overlap et Classification ---
raw_overlaps <- perform_overlap(data_snp, all_reference_regions)

annotated_snps <- raw_overlaps %>%
  mutate(rank = priority_map[ref_type]) %>%
  group_by(chromosome, SNPs) %>%
  slice_min(rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(final_type = case_when(
    # Regroupement "Intra-gene"
    ref_type %in% c("CDS", "five_prime_UTR", "three_prime_UTR", "exon", "ncRNA", "mRNA", "gene") ~ "Intra-gene",
    # Conservation des zones de distance
    ref_type == "Intergenic_0_2kb" ~ "Intergenic < 2kb",
    ref_type == "Intergenic_2_5kb" ~ "Intergenic 2-5kb",
    TRUE ~ ref_type
  ))

# --- 4. Statistiques finales avec Intergenic > 5kb ---
total_sig_snps <- nrow(data_snp)
hit_snps_count <- n_distinct(paste(annotated_snps$chromosome, annotated_snps$SNPs))

summary_counts <- annotated_snps %>%
  count(final_type, name = "count") %>%
  add_row(final_type = "Intergenic > 5kb", count = total_sig_snps - hit_snps_count) %>%
  mutate(percentage = (count / total_sig_snps) * 100)

print(summary_counts)

order_levels <- c("Intra-gene", "Intergenic < 2kb", "Intergenic 2-5kb", "Intergenic > 5kb")
summary_counts$final_type <- factor(summary_counts$final_type, levels = order_levels)

p_final <- ggplot(summary_counts, aes(x = final_type, y = count, fill = final_type)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.7) +
  # Ajout du nombre de SNPs au-dessus de chaque barre
  geom_text(aes(label = comma(count)), vjust = -0.5, fontface = "bold", size = 4.5) +
  scale_fill_manual(values = c("#264653", "#2a9d8f", "#e9c46a", "#f4a261")) +
  theme_minimal() +
  labs(
    title = "Significant SNPs: Proximity to Annotated Genes",
    subtitle = paste0("Total SNPs analyzed: ", comma(total_sig_snps)),
    x = "Genomic Localization",
    y = "Number of SNPs"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 11, color = "black")
  )

ggsave("SNP_Distribution.png", p_final, width = 10, height = 6, dpi = 300)


library(ggplot2)
library(dplyr)
library(scales)

# --- 1. Préparation des données pour le graphique ---
# On s'assure que les chromosomes sont dans le bon ordre (Chr01, Chr02...)
# et non dans l'ordre alphabétique (Chr01, Chr10, Chr11...).
# On suppose que vos chromosomes s'appellent "Chr01", "Chr02", etc.
level_order <- sprintf("Chr%02d", 1:max(as.numeric(sub("Chr", "", unique(annotated_snps$chromosome)))))

plot_map_data <- annotated_snps %>%
  mutate(chromosome = factor(chromosome, levels = level_order)) %>%
  filter(ref_type=="Intergenic_0_2kb")%>%
  # Pour geom_linerange, on définit une petite hauteur fixe (ymin/ymax)
  # pour que chaque SNP ressemble à une barre verticale sur la ligne du Chr.
  mutate(ymin = 0, ymax = 1)

# --- 2. Création de la figure (Carte Chromosomique) ---
karyoplot <- ggplot(plot_map_data, aes(x = SNPs, ymin = ymin, ymax = ymax, color = final_type)) +
  # geom_linerange trace une ligne verticale entre ymin et ymax à la position x (SNPs)
  geom_linerange(size = 0.5, alpha = 0.8) +
  # Facet_grid crée une ligne séparée par chromosome.
  # scales = "free_x" permet à chaque Chr d'avoir sa propre échelle de longueur.
  # space = "free_x" ajuste la largeur de la facette à la longueur réelle du Chr.
  facet_grid(chromosome ~ ., scales = "free_x", space = "free_x", switch = "y") +
  # Palette de couleurs Set1 (distincte)
  scale_color_brewer(palette = "Set1", name = "Localization") +
  # Formatage de l'axe X en Mégabases (Mb) pour plus de lisibilité
  scale_x_continuous(labels = function(x) paste0(round(x/1e6, 1), " Mb")) +
  theme_minimal() +
  labs(
    title = "Genome-Wide Distribution of Significant SNPs",
    subtitle = paste0("Total prioritized SNPs: ", comma(nrow(plot_map_data))),
    x = "Physical Position along Chromosome",
    y = ""
  ) +
  theme(
    # Supprime les labels et grilles de l'axe Y devenus inutiles
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # Améliore l'affichage des labels de chromosomes sur le côté
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

# --- 3. Affichage et Sauvegarde ---
print(karyoplot)

# Sauvegarde en haute résolution (recommandé car il y a beaucoup de données)
ggsave("SNP_Genome_Map_Intergenic_0_2kb.png", karyoplot, width = 12, height = 8, dpi = 300)




total_snps <- sum(summary_counts$count)

# Calcul des pourcentages et préparation pour le graphique
plot_data <- summary_counts %>%
  mutate(percentage = (count / total_snps) * 100) %>%
  arrange(desc(percentage)) %>%
  # Création d'une seule barre (x = "")
  mutate(group = "SNPs")

# 2. Création de la figure
p <- ggplot(plot_data, aes(x = group, y = percentage, fill = reorder(final_type, -percentage))) +
  # Barre empilée
  geom_bar(stat = "identity", width = 0.6, color = "white") +
  # Ajout des labels de pourcentage (uniquement si > 2% pour la lisibilité)
  geom_text(aes(label = ifelse(percentage > 2, paste0(round(percentage, 1), "%"), "")),
            position = position_stack(vjust = 0.5), 
            color = "white", fontface = "bold", size = 4) +
  # Dégradé de couleur professionnel
  scale_fill_viridis_d(option = "D", name = "Categories") +
  # Passage en horizontal
  coord_flip() +
  # Thème épuré
  theme_minimal() +
  labs(
    title = paste0("Genomic Distribution of Significant SNPs (Total: ", comma(total_snps), ")"),
    subtitle = "Percentage of SNPs by functional annotation",
    x = "",
    y = "Percentage (%)"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 9)
  )

# 3. Affichage et sauvegarde
print(p)

# Optionnel : Sauvegarder en haute résolution
ggsave("SNP_Distribution_Gradient.png", p, width = 10, height = 4, dpi = 300)


genes_from_mrna <- annotated_snps %>%
  # On cible les mRNA, mais vous pourriez aussi inclure CDS/UTR/Intron
  filter(final_type == "mRNA" | ref_type == "mRNA") %>%
  # Nettoyage du gene_id au cas où il reste des préfixes "gene:"
  mutate(gene_id = sub("gene:", "", gene_id)) %>%
  dplyr::select(chromosome, SNPs, gene_id, pvalues, FDR_signif)

# 2. Créer une liste unique des gènes candidats
# Utile pour les analyses d'enrichissement GO ou KEGG
unique_candidate_genes <- genes_from_mrna %>%
  group_by(gene_id) %>%
  summarise(
    n_snps = n(),
    best_pvalue = min(pvalues),
    chromosome = dplyr::first(chromosome)
  ) %>%
  arrange(best_pvalue)

# Affichage des premiers gènes trouvés
print(head(unique_candidate_genes))


# library(ggplot2)
# library(dplyr)
# library(scales)

# chr_lengths <- gff_base %>%
#   group_by(chr) %>%
#   summarise(length = max(end), .groups = "drop")

# # Join with annotated SNPs
# plot_data <- annotated_snps %>%
#   left_join(chr_lengths, by = c("chromosome" = "chr")) %>%
#   mutate(chromosome = factor(chromosome, levels = sprintf("Chr%02d", 1:max(as.numeric(sub("Chr", "", unique(chromosome)))))))

# # Define the 3 distinct groups for the 3 rings
# plot_data <- plot_data %>%
#   mutate(track = case_when(
#     final_type == "Intra-gene" ~ "Track 1: Intra-gene",
#     final_type %in% c("Intergenic < 2kb", "Intergenic 2-5kb") ~ "Track 2: Intergenic < 5kb",
#     final_type == "Intergenic > 5kb" ~ "Track 3: Intergenic > 5kb"
#   ))

# # 2. Create the Circular Plot
# # We use geom_point with high transparency (alpha) to show density
# circle_dist_plot <- ggplot(plot_data, aes(x = SNPs, y = track, color = final_type)) +
#   # We use small tiles or points to represent SNP density
#   geom_jitter(size = 0.1, alpha = 0.4, width = 0) +
#   # Facet by chromosome to create the circular segments
#   facet_wrap(~chromosome, ncol = 6) + # Initially rectangular to check
#   # Convert to Circular
#   coord_polar() +
#   scale_color_manual(values = c(
#     "Intra-gene" = "#264653",
#     "Intergenic < 2kb" = "#2a9d8f",
#     "Intergenic 2-5kb" = "#e9c46a",
#     "Intergenic > 5kb" = "#e76f51"
#   )) +
#   theme_minimal() +
#   labs(
#     title = "Genome-Wide Triple Ring SNP Distribution",
#     subtitle = "Outer: Intra-gene | Middle: < 5kb | Inner: > 5kb",
#     color = "Category"
#   ) +
#   theme(
#     axis.title = element_blank(),
#     axis.text.y = element_blank(),
#     axis.text.x = element_text(size = 7),
#     panel.grid.major.y = element_line(color = "grey90"),
#     legend.position = "bottom",
#     strip.text = element_text(face = "bold")
#   )

# ggsave("SNP_circle_plot.png", circle_dist_plot, width = 10, height = 10, dpi = 300)

# --- 1. Préparation des Coordonnées Génomiques Globales ---
# Pour mettre tous les chromosomes dans un seul cercle, on doit les placer bout à bout.

# Calculer la longueur de chaque chromosome
chr_info <- gff_base %>%
  group_by(chr) %>%
  summarise(length = max(end), .groups = "drop") %>%
  # S'assurer de l'ordre Chr01, Chr02...
  mutate(chr_num = as.numeric(sub("Chr", "", chr))) %>%
  arrange(chr_num) %>%
  mutate(
    # Calculer le décalage (offset) pour chaque chromosome
    tot = cumsum(as.numeric(length)) - length,
    center = tot + length/2
  )

# Transformer les positions des SNPs en coordonnées globales (Genome-wide)
# plot_data_global <- annotated_snps %>%
#   left_join(chr_info, by = c("chromosome" = "chr")) %>%
#   mutate(global_pos = SNPs + tot) %>%
#   # Définir les 3 couches (Rings)
#   mutate(track_id = case_when(
#     final_type == "Intra-gene" ~ 3,         # Cercle extérieur
#     final_type %in% c("Intergenic < 2kb", "Intergenic 2-5kb") ~ 2, # Milieu
#     final_type == "Intergenic > 5kb" ~ 1    # Cercle intérieur
#   ))

plot_data_global <- annotated_snps %>%
  left_join(chr_info, by = c("chromosome" = "chr")) %>%
  mutate(global_pos = SNPs + tot) %>%
  # Définir les 3 couches (Rings)
  mutate(track_id = case_when(
    final_type == "Intra-gene" ~ 1,         # Cercle extérieur
    # final_type %in% c("Intergenic < 2kb", "Intergenic 2-5kb") ~ 2, # Milieu
    # final_type == "Intergenic > 5kb" ~ 1    # Cercle intérieur
  ))

# --- 2. Création du Circos Plot avec ggplot2 ---
library(ggplot2)

circos_plot <- ggplot() +
  # 1. Dessiner les lignes de séparation des chromosomes
  geom_vline(data = chr_info, aes(xintercept = tot), color = "grey80", linetype = "dashed") +
  
  # 2. Tracer les SNPs sur 3 niveaux (tracks)
  # On utilise geom_segment pour faire des "traits" radiaux
  geom_segment(data = plot_data_global, 
               aes(x = global_pos, xend = global_pos, 
                   y = track_id - 0.4, yend = track_id + 0.4, 
                   color = final_type),
               alpha = 0.5, size = 0.2) +
  
  # 3. Ajouter les noms des chromosomes à l'extérieur
  geom_text(data = chr_info, aes(x = center, y = 3.8, label = chr), 
            size = 3, fontface = "bold") +
  
  # 4. Paramètres de l'axe et transformation circulaire
  scale_y_continuous(limits = c(0, 4)) + 
  scale_color_manual(values = c(
    "Intra-gene" = "#E63946"       # Rouge
    # "Intergenic < 2kb" = "#457B9D", # Bleu
    # "Intergenic 2-5kb" = "#A8DADC",# Bleu clair
    # "Intergenic > 5kb" = "#1D3557"  # Bleu marine
  ), name = "Annotation") +
  coord_polar(theta = "x") +
  
  # 5. Esthétique
  theme_void() +
  labs(
    title = "Global Circos Plot of Significant SNPs",
    subtitle = "Outer Ring: Intra-gene | Middle: < 5kb | Inner: > 5kb"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom"
  )

# Affichage
print(circos_plot)

# Sauvegarde recommandée en PDF ou PNG haute résolution pour garder les détails
ggsave("Global_Circos_SNP_Map_intregene.png", circos_plot, width = 10, height = 10, dpi = 600)
## same thinks but per bio variables


# ==================================================
# 1. Read and merge SNP position files
# ==================================================
# snp_dir <- "./out/dom_admixed_rm/s02.lfmm2/file/"

# snp_files <- list.files(
#   path      = snp_dir,
#   pattern   = "chr_chr.*_SNPs_associated_with_all_bio_VCF_",
#   recursive = TRUE,
#   full.names = TRUE
# )

# stopifnot(length(snp_files) > 0)

# # Helper to read and tag a file
# read_file <- function(file, method) {
#   chr <- str_extract(file, "chr_chr[0-9]+") %>%
#     sub("^chr_", "", .) %>%
#     sub("^chr", "Chr", .)

#   read_csv(file, show_col_types = FALSE) %>%
#     transmute(
#       SNPs       = SNPs,
#       chromosome = chr,
#       !!method   := TRUE
#     )
# }

# # Separate files
# bonf_files <- snp_files[str_detect(snp_files, "Bonferroni")]
# fdr_files  <- snp_files[str_detect(snp_files, "FDR")]

# # Read all files
# bonf_df <- map_dfr(bonf_files, read_file, method = "Bonferroni")
# fdr_df  <- map_dfr(fdr_files,  read_file, method = "FDR")

# # Merge Bonferroni + FDR
# all_snps_pos <- full_join(
#   bonf_df,
#   fdr_df,
#   by = c("SNPs", "chromosome")
# ) %>%
#   mutate(
#     Bonferroni = ifelse(is.na(Bonferroni), FALSE, Bonferroni),
#     FDR        = ifelse(is.na(FDR), FALSE, FDR)
#   ) %>%
#   arrange(chromosome, SNPs)

# ==================================================
# 2. Read and parse GFF3 annotation
# ==================================================
# gff <- read.table(
#   "./in/gene_models_20170612.gff3",
#   sep = "\t",
#   header = FALSE,
#   comment.char = "#",
#   quote = "",
#   stringsAsFactors = FALSE,
#   fill = TRUE
# )

# colnames(gff) <- c(
#   "chr", "source", "type", "start", "end",
#   "score", "strand", "phase", "attributes"
# )

# gff <- gff %>%
#   mutate(
#     chr = sub("^chr", "Chr", chr),
#     gene_id = case_when(
#       type == "gene" ~ sub(".*ID=gene:([^;]+).*", "\\1", attributes),
#       TRUE           ~ sub(".*Parent=gene:([^;]+).*", "\\1", attributes)
#     )
#   )

# ==================================================
# 3. Build genomic ranges for overlap
# ==================================================
# gff_ranges <- gff %>%
#   dplyr::select(chr, start, end, type, gene_id) %>%
#   group_by(chr) %>%
#   summarise(
#     ranges  = list(IRanges(start, end)),
#     type    = list(type),
#     gene_id = list(gene_id),
#     .groups = "drop"
#   )

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







#####

rm(list = ls())

PROJECT <- "dom_admixed_rm"
snp_path <- paste0("./out/", PROJECT, "/s02.lfmm2_V3/file/All_SNPs_Filtered_Results.csv")
gff_path <- "./in/gene_models_20170612.gff3"

data_snp <- read.csv(snp_path) %>%
  filter(FDR_signif == TRUE) %>%
  mutate(chromosome = sub("^chr", "Chr", chromosome))

gff_raw <- read.table(
  gff_path, sep = "\t", header = FALSE, 
  comment.char = "#", quote = "", stringsAsFactors = FALSE, fill = TRUE
)

colnames(gff_raw) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# --- 2. Prepare Annotation Ranges ---
# Standardize GFF and extract Gene IDs
gff_base <- gff_raw %>%
  mutate(
    chr = sub("^chr", "Chr", chr),
    gene_id = case_when(
      type == "gene" ~ sub(".*ID=gene:([^;]+).*", "\\1", attributes),
      TRUE           ~ sub(".*Parent=gene:([^;]+).*", "\\1", attributes)
    )
  )

gene_bounds <- gff_base %>%
  filter(type == "gene") %>% filter(chr != "Chr00") %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(
    # La fin du gène précédent
    prev_end = lag(end),
    # L'espace commence à la fin du gène précédent + 1
    inter_start = prev_end + 1,
    # L'espace finit au début du gène actuel - 1
    inter_end = start - 1
  ) %>%
  # Garder seulement les lignes où il y a un espace (inter_start < inter_end)
  filter(!is.na(inter_start) & inter_start < inter_end) %>%
  mutate(
    gap_size = inter_end - inter_start + 1,
    # Étiquetage selon vos critères
    intergenic_label = case_when(
      gap_size <= 2000 ~ "Intergenic_0_2kb",
      gap_size <= 5000 ~ "Intergenic_2_5kb",
      TRUE             ~ "Intergenic_above_5kb"
    )
  ) %>%
  ungroup()

intergenic_ranges <- gene_bounds %>%
  dplyr::select(chr, start = inter_start, end = inter_end, gap_size, intergenic_label)

all_reference_regions <- bind_rows(
    gff_base %>% filter(chr != "Chr00") %>% 
      dplyr::select(chr, start, end, type, gene_id), 
    intergenic_ranges %>% 
      mutate(type = intergenic_label) # On utilise le label (0-2kb, etc.) comme 'type'
  ) %>%
  # 2. Groupement par chromosome
  group_by(chr) %>%
  # 3. Tri par position de début (start)
  arrange(start, .by_group = TRUE) %>%
  ungroup()


perform_overlap <- function(snps, refs) {
  # Group by chromosome to speed up IRanges
  all_hits <- lapply(intersect(unique(snps$chromosome), unique(refs$chr)), function(c) {
    s_sub <- snps %>% filter(chromosome == c)
    r_sub <- refs %>% filter(chr == c)
    
    if(nrow(s_sub) == 0 | nrow(r_sub) == 0) return(NULL)
    
    s_ir <- IRanges(start = s_sub$SNPs, end = s_sub$SNPs)
    r_ir <- IRanges(start = r_sub$start, end = r_sub$end)
    
    hits <- findOverlaps(s_ir, r_ir)
    
    if(length(hits) == 0) return(NULL)
    
    # Merge SNP data with the matched reference type
    cbind(s_sub[queryHits(hits), ], 
          ref_type = r_sub$type[subjectHits(hits)], 
          gene_id  = r_sub$gene_id[subjectHits(hits)])
  })
  bind_rows(all_hits)
}

# Run the overlap
raw_overlaps <- perform_overlap(data_snp, all_reference_regions)

priority_map <- c(
  "CDS"              = 1, 
  "five_prime_UTR"   = 2, 
  "three_prime_UTR"  = 3, 
  "exon"             = 4, 
  "ncRNA"            = 5,
  "mRNA"             = 6, 
  "gene"             = 7, 
  "Intergenic_0_2kb" = 8,
  "Intergenic_2_5kb" = 9,
  "Intergenic_above_5kb" = 10
)

annotated_snps <- raw_overlaps %>%
  mutate(rank = priority_map[ref_type]) %>%
  group_by(chromosome, SNPs) %>%
  slice_min(rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(final_type = case_when(
    # Regroupement "Intra-gene"
    ref_type %in% c("CDS", "five_prime_UTR", "three_prime_UTR", "exon", "ncRNA", "mRNA", "gene") ~ "Intra-gene",
    # Conservation des zones de distance
    ref_type == "Intergenic_0_2kb" ~ "Intergenic < 2kb",
    ref_type == "Intergenic_2_5kb" ~ "Intergenic 2-5kb",
    ref_type == "Intergenic_above_5kb" ~ "Intergenic > 5kb",
    TRUE ~ ref_type
  ))

## ongoing








## per bio

bio_files <- list.files("./out/dom_admixed_rm/s02.lfmm2_V3/file/", pattern = "bio.*.csv", full.name=TRUE)

gff_path <- "./in/gene_models_20170612.gff3"

gff_raw <- read.table(
  gff_path, sep = "\t", header = FALSE, 
  comment.char = "#", quote = "", stringsAsFactors = FALSE, fill = TRUE
)

colnames(gff_raw) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# --- 2. Prepare Annotation Ranges ---
# Standardize GFF and extract Gene IDs
gff_base <- gff_raw %>%
  mutate(
    chr = sub("^chr", "Chr", chr),
    gene_id = case_when(
      type == "gene" ~ sub(".*ID=gene:([^;]+).*", "\\1", attributes),
      TRUE           ~ sub(".*Parent=gene:([^;]+).*", "\\1", attributes)
    )
  )

gene_bounds <- gff_base %>%
  filter(type == "gene") %>% filter(chr != "Chr00") %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(
    # La fin du gène précédent
    prev_end = lag(end),
    # L'espace commence à la fin du gène précédent + 1
    inter_start = prev_end + 1,
    # L'espace finit au début du gène actuel - 1
    inter_end = start - 1
  ) %>%
  # Garder seulement les lignes où il y a un espace (inter_start < inter_end)
  filter(!is.na(inter_start) & inter_start < inter_end) %>%
  mutate(
    gap_size = inter_end - inter_start + 1,
    # Étiquetage selon vos critères
    intergenic_label = case_when(
      gap_size <= 2000 ~ "Intergenic_0_2kb",
      gap_size <= 5000 ~ "Intergenic_2_5kb",
      TRUE             ~ "Intergenic_above_5kb"
    )
  ) %>%
  ungroup()

intergenic_ranges <- gene_bounds %>%
  dplyr::select(chr, start = inter_start, end = inter_end, gap_size, intergenic_label)

all_reference_regions <- bind_rows(
    gff_base %>% filter(chr != "Chr00") %>% 
      dplyr::select(chr, start, end, type, gene_id), 
    intergenic_ranges %>% 
      mutate(type = intergenic_label) # On utilise le label (0-2kb, etc.) comme 'type'
  ) %>%
  # 2. Groupement par chromosome
  group_by(chr) %>%
  # 3. Tri par position de début (start)
  arrange(start, .by_group = TRUE) %>%
  ungroup()


perform_overlap <- function(snps, refs) {
  # Group by chromosome to speed up IRanges
  all_hits <- lapply(intersect(unique(snps$chromosome), unique(refs$chr)), function(c) {
    s_sub <- snps %>% filter(chromosome == c)
    r_sub <- refs %>% filter(chr == c)
    
    if(nrow(s_sub) == 0 | nrow(r_sub) == 0) return(NULL)
    
    s_ir <- IRanges(start = s_sub$SNPs, end = s_sub$SNPs)
    r_ir <- IRanges(start = r_sub$start, end = r_sub$end)
    
    hits <- findOverlaps(s_ir, r_ir)
    
    if(length(hits) == 0) return(NULL)
    
    # Merge SNP data with the matched reference type
    cbind(s_sub[queryHits(hits), ], 
          ref_type = r_sub$type[subjectHits(hits)], 
          gene_id  = r_sub$gene_id[subjectHits(hits)])
  })
  bind_rows(all_hits)
}

priority_map <- c(
  "CDS"              = 1, 
  "five_prime_UTR"   = 2, 
  "three_prime_UTR"  = 3, 
  "exon"             = 4, 
  "ncRNA"            = 5,
  "mRNA"             = 6, 
  "gene"             = 7, 
  "Intergenic_0_2kb" = 8,
  "Intergenic_2_5kb" = 9,
  "Intergenic_above_5kb" = 10
)

all_bio_results <- list()

for (file in bio_files) {
  # Extraire le nom court (ex: bio1) du nom de fichier
  bio_name <- str_extract(basename(file), "bio[0-9]+")
  message(">>> Processing: ", bio_name)
  
  data_snp <- read.csv(file) %>%
    filter(FDR_signif == TRUE) %>% ## faire avec Bonf
    mutate(chromosome = sub("^chr", "Chr", chromosome))
  
  if(nrow(data_snp) == 0) {
    message("No significant SNPs for ", bio_name)
    next
  }

  raw_overlaps <- perform_overlap(data_snp, all_reference_regions)

  annotated_snps <- raw_overlaps %>%
    mutate(rank = priority_map[ref_type]) %>%
    group_by(chromosome, SNPs) %>%
    slice_min(rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      final_type = case_when(
        ref_type %in% c("CDS", "five_prime_UTR", "three_prime_UTR", "exon", "ncRNA", "mRNA", "gene") ~ "Intra-gene",
        ref_type == "Intergenic_0_2kb" ~ "Intergenic < 2kb",
        ref_type == "Intergenic_2_5kb" ~ "Intergenic 2-5kb",
        ref_type == "Intergenic_above_5kb" ~ "Intergenic > 5kb",
        TRUE ~ ref_type
      ),
      bio_source = bio_name # Ajout de l'identifiant Bio
    )
  
  all_bio_results[[bio_name]] <- annotated_snps
}

# Combiner tout en un seul grand dataframe
combined_results <- bind_rows(all_bio_results)

combined_results <- combined_results %>%
  mutate(gene_id_clean = str_extract(gene_id, "MD\\d{2}G\\d+"))


# load SNP for all_bio in GEA
data_snp_all_bio <- read.csv("./out/dom_admixed_rm/s02.lfmm2_V3/file/All_SNPs_Filtered_Results.csv") %>%
    filter(FDR_signif == TRUE) %>% ## faire avec Bonf
    mutate(chromosome = sub("^chr", "Chr", chromosome))

full_bio_overlapp <- perform_overlap(data_snp_all_bio, all_reference_regions)

annotated_snps_all_bio <- full_bio_overlapp %>%
    mutate(rank = priority_map[ref_type]) %>%
    group_by(chromosome, SNPs) %>%
    slice_min(rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      final_type = case_when(
        ref_type %in% c("CDS", "five_prime_UTR", "three_prime_UTR", "exon", "ncRNA", "mRNA", "gene") ~ "Intra-gene",
        ref_type == "Intergenic_0_2kb" ~ "Intergenic < 2kb",
        ref_type == "Intergenic_2_5kb" ~ "Intergenic 2-5kb",
        ref_type == "Intergenic_above_5kb" ~ "Intergenic > 5kb",
        TRUE ~ ref_type
      ),
      bio_source = "Full_geno" # Ajout de l'identifiant Bio
    )

annotated_snps_all_bio <- annotated_snps_all_bio %>%
  mutate(gene_id_clean = str_extract(gene_id, "MD\\d{2}G\\d+"))


# Load required libraries
library(UpSetR)
library(dplyr)
library(tidyr)

# 1. Define the p-value variables and threshold
p_bio_cols <- paste0("p_bio", sprintf("%02d", 1:19))

sig_threshold <- 0.05 # Adjust this threshold as needed (e.g., 0.01 or 0.001)

# 2. Create a list where each element contains unique genes significant for that variable
# We exclude NA gene IDs (intergenic SNPs)
gene_list <- lapply(p_bio_cols, function(col) {
  combined_results %>%
    filter(ref_type == "CDS") %>%
    filter(!!sym(col) < sig_threshold) %>%
    filter(!is.na(gene_id_clean)) %>%
    pull(gene_id_clean) %>%
    unique()
})



# Name the list elements after the p_bio variables
names(gene_list) <- p_bio_cols

gene_list$Full_bio <- unique(annotated_snps_all_bio[annotated_snps_all_bio$ref_type == "CDS", ]$gene_id_clean)

## Blast

table_flowering <- read.table("./in/flowering_GDDH_florid.txt", sep="\t", header=T)

library(ggVennDiagram)
library(ggplot2)

# 1. Create a subset list with only the 4 groups you want to compare
# Mapping the names from your summary to your request
selected_sets <- list(
  bio19 = gene_list$p_bio19,
  bio08 = gene_list$p_bio08,
  bio12 = gene_list$p_bio12,
  Full_geno = gene_list$Full_bio # Adjust this if the name is different in your environment
)

# 2. Generate the Venn Diagram
p <- ggVennDiagram(selected_sets) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Venn Diagram: bio19, 08, 12, and Full_geno",
       caption = "Data source: gene_list") +
  theme_void() + 
  # Désactive le découpage du texte qui dépasse
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "none", 
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10, "pt")
  )

ggsave(
  filename = "venn_diagram_bio19_8_12_Full_geno.png", 
  plot = p, 
  width = 8, 
  height = 6, 
  units = "in", 
  dpi = 300
)

## Blast
table_flowering <- read.table("./in/flowering_GDDH_florid.txt", sep="\t", header=T)
# Extract the intersection data
venn_data <- process_region_data(Venn(selected_sets))
# Ensure the flowering gene list is a vector for easy comparison
flowering_genes <- as.character(table_flowering$GD_gene)

# Create a clean results list
comparison_results <- venn_data %>%
  rowwise() %>%
  mutate(
    # Count how many genes are in this specific intersection
    total_genes = length(item),
    # Find which genes in this intersection are in the flowering table
    matching_flowering = list(intersect(item, flowering_genes)),
    # Count the matches
    count_flowering = length(matching_flowering)
  ) %>%
  ungroup()

# View the summary of intersections and flowering gene hits
print(comparison_results[, c("name", "total_genes", "count_flowering")])

all_genes_with_groups <- comparison_results %>%
  dplyr::select(name, item) %>%
  unnest(item) %>%
  rename(GD_gene = item, Venn_Group = name)

# 2. Join with your flowering table to get the 'flowering_gene' column
# We use left_join to keep only the genes that are in your Venn sets
final_gene_summary <- all_genes_with_groups %>%
  inner_join(table_flowering, by = "GD_gene") %>%
  # Select only the columns you requested
  dplyr::select(GD_gene, Venn_Group, flowering_gene, TAIR)

# 3. Print the first few rows to check
print(head(final_gene_summary, 20))

print(final_gene_summary[final_gene_summary$Venn_Group %in% c("Full_geno", "bio19/Full_geno", "bio08/Full_geno", "bio12/Full_geno", "bio19/bio08/Full_geno",
"bio19/bio12/Full_geno", "bio08/bio12/Full_geno", "bio19/bio08/bio12/Full_geno"),], n=178)

print(combined_results %>% 
  filter(gene_id_clean %in% final_gene_summary[final_gene_summary$Venn_Group %in% c(
    "Full_geno", "bio19/Full_geno", "bio08/Full_geno", "bio12/Full_geno", 
    "bio19/bio08/Full_geno", "bio19/bio12/Full_geno", "bio08/bio12/Full_geno", 
    "bio19/bio08/bio12/Full_geno"),]$GD_gene) %>% 
  dplyr::select(SNPs, p_bio08, p_bio12, p_bio19, chromosome, ref_type) %>%
  filter(!is.na(p_bio08) | !is.na(p_bio12) | !is.na(p_bio19)), n=1741)

fasta <- readDNAStringSet("./in/GDDH13_1-1_mrna.fasta")

universe_genes <- unique(
  sub(".*gene[:=]([^;| ]+).*", "\\1", names(fasta))
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

semdata <- godata(
  OrgDb  = org.Mdomestica.eg.db,
  keytype = "GID",
  ont     = "BP"
)

# Initialize the plot storage lists
dotplot_list <- list()
emapplot_list <- list()
treeplot_list <- list()

for (group in unique(final_gene_summary$Venn_Group)){
#for (group in venn_data$name){
  print(group)
  #current_genes <- venn_data[venn_data$name == group,]$item[[1]]
  current_genes <- final_gene_summary[final_gene_summary$Venn_Group == group,]$GD_gene
  
  message(paste("--- Processing Intersection:", group, "| Genes:", length(current_genes), "---"))
  
  # Run Enrichment
  go_bp <- enrichGO(
    gene          = current_genes, 
    universe      = universe_genes,
    OrgDb         = org.Mdomestica.eg.db, 
    keyType       = "GID", 
    ont           = "BP", 
    pvalueCutoff  = 0.05, 
    qvalueCutoff  = 0.05, 
    readable      = TRUE,
    pAdjustMethod = "fdr"
  ); print(go_bp)
  
  if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
    
    # Extract gene count for title to verify it matches UpSet
    n_genes <- length(current_genes)
    plot_title <- paste(group, "(n =", n_genes, ")")
    
    # A. DOTPLOT (Remove legend immediately)
    dotplot_list[[group]] <- dotplot(go_bp, showCategory = 10) + 
                             labs(title = plot_title) +
                             theme(legend.position = "none") 
    
    # B. EMAPPLOT & TREEPLOT (Need > 2 terms for similarity calculation)
    if (nrow(as.data.frame(go_bp)) > 2) {
      
      # Assuming you have safe_pairwise_termsim defined, otherwise use pairwise_termsim(go_bp)
      go_bp_sim <- pairwise_termsim(go_bp) 
      
      emapplot_list[[group]] <- emapplot(go_bp_sim, showCategory = 20) + 
                                labs(title = plot_title) +
                                theme(legend.position = "none")
      
      n_terms <- nrow(as.data.frame(go_bp_sim))
      if (n_terms >= 2) {
        n_clusters <- min(5, n_terms - 1)
        treeplot_list[[group]] <- suppressWarnings(
          treeplot(go_bp_sim, cluster.params = list(n = max(1, n_clusters)))
        ) + labs(title = plot_title) +
          theme(legend.position = "none") 
      }
    }
  }
}





# Convert to data frame if it isn't already
df <- as.data.frame(fromList(gene_list))

# Total genes in each p_bio group
total_per_group <- colSums(df[, grep("p_bio", names(df))])

# Count overlaps with Full_bio
overlaps_full_bio <- sapply(df[, grep("p_bio", names(df))], function(x) sum(x == 1 & df$Full_bio == 1))

# Count genes in p_bio but NOT in Full_bio
unique_to_group <- sapply(df[, grep("p_bio", names(df))], function(x) sum(x == 1 & df$Full_bio == 0))

library(ggplot2)

# 1. Get the list of variable names
bio_names <- names(overlaps_full_bio)

# 2. Calculate the total size to order the bars from largest to smallest
total_counts <- overlaps_full_bio + unique_to_group
ordered_bios <- names(sort(total_counts, decreasing = FALSE)) # Ascending so largest is on top when flipped

# 3. Create the summary dataframe for ggplot
summary_df <- data.frame(
  Bio = factor(rep(bio_names, 2), levels = ordered_bios), # Order by total size
  Status = factor(
    rep(c("In Full_bio", "Not in Full_bio"), each = length(bio_names)),
    levels = c("Not in Full_bio", "In Full_bio") # Stack "In Full_bio" on top/right
  ),
  Count = c(overlaps_full_bio, unique_to_group)
)

# 4. Generate the stacked bar chart
png("All_Bios_vs_Full_bio_Histogram.png", width = 1800, height = 2000, res = 200)

ggplot(summary_df, aes(x = Bio, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  coord_flip() + # Flips the chart horizontally to make the 19 names readable
  # Add exact numbers inside the bars
  geom_text(aes(label = ifelse(Count > 0, Count, "")), # Only labels > 0 to avoid clutter
            position = position_stack(vjust = 0.5), 
            size = 4, color = "black") +
  # Match your existing color scheme
  scale_fill_manual(values = c("Not in Full_bio" = "darkorange", "In Full_bio" = "steelblue")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Gene Counts per Bioclimatic Variable",
    subtitle = "Separated by overlap with Full_bio reference list",
    x = "Bioclimatic Variable",
    y = "Number of Genes",
    fill = "Reference Status"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.text.y = element_text(face = "bold") # Make the bio names bold
  )

dev.off()

library(UpSetR)

# 1. Get all bio variables and sort them alphabetically so they go 01 to 19
bio_cols <- sort(grep("p_bio", names(df), value = TRUE))

# 2. Build the exact sequence of intersections you want to show
intersections_to_show <- list()

for(bio in bio_cols) {
  # Add "Bio Only" intersection
  intersections_to_show[[length(intersections_to_show) + 1]] <- bio
  
  # Add "Bio + Full_bio" intersection
  intersections_to_show[[length(intersections_to_show) + 1]] <- c(bio, "Full_bio")
}

# 3. Define the sets to include in the plot (Horizontal rows)
# Putting Full_bio first ensures it stays at the top of the matrix
target_sets <- c("Full_bio", rev(bio_cols)) 

# 4. Generate the ordered UpSet plot
png("Strict_Order_Bios_vs_Full_bio_UpSet.png", width = 4500, height = 2000, res = 200)

upset(df, 
      sets = target_sets, 
      keep.order = TRUE,             # Keeps the horizontal rows exactly as defined above
      intersections = intersections_to_show, # Forces the exact left-to-right vertical bar sequence!
      
      # Since we alternate 1) Bio Only, 2) Bio+Full_bio, we can alternate the colors!
      main.bar.color = rep(c("darkorange", "steelblue"), length(bio_cols)), 
      
      sets.bar.color = "gray30",
      matrix.color = "black",
      text.scale = c(2, 1.5, 1.5, 1.5, 1.5, 1.3),
      point.size = 3, 
      line.size = 1)

dev.off()

# 1. Identify your p_bio group columns
p_bio_cols <- names(df)[grep("p_bio", names(df))]

# 2. CALCULATE SIZES AND SORT: 
# We get the sum of each column, sort them descending, and keep the names.
# We keep "Full_bio" out of the sort or put it at the top manually.
group_sizes <- colSums(df[, p_bio_cols])
ordered_p_bio <- names(sort(group_sizes, decreasing = TRUE))

# Combine with Full_bio at the top (or bottom, depending on your preference)
ordered_sets <- c("Full_bio", ordered_p_bio)

# 3. Construct the intersections list based on the NEW order
# This ensures the vertical bars follow the same horizontal order
intersections_to_show <- list()
for(col in ordered_p_bio) {
  intersections_to_show[[length(intersections_to_show) + 1]] <- c(col, "Full_bio")
  intersections_to_show[[length(intersections_to_show) + 1]] <- col
}

n_only=0;n_full_bio=0
for (i in intersections_to_show) { if (length(i) >1) { n_full_bio <- n_full_bio+1 } else { n_only <- n_only+1 }}

# 4. Generate the plot
png("Strict_Order_Bios_vs_Full_bio_UpSet.png", width = 4500, height = 2000, res = 200)

upset(df, 
      sets = ordered_sets,           # <--- This orders the horizontal 'Set Size' bars
      intersections = intersections_to_show, 
      keep.order = TRUE,             # <--- This maintains the order of your list
      sets.bar.color = "gray30",
      main.bar.color = c(rep("steelblue", n_only), rep("darkorange",n_full_bio)),
      text.scale = c(2, 1.5, 1.5, 1.5, 1.5, 1.2),
      matrix.color = "black",
      point.size = 3, 
      line.size = 1)

dev.off()

# Get all unique genes across all lists to match the rows of the 'df'
# ---------------------------------------------------------
# 1. Build the binary matrix WITH gene IDs as row names
# ---------------------------------------------------------
all_unique_genes <- unique(unlist(gene_list))

# Create an empty dataframe of 0s
df <- data.frame(matrix(0, nrow = length(all_unique_genes), ncol = length(gene_list)))
rownames(df) <- all_unique_genes # Set row names to actual genes
colnames(df) <- names(gene_list)

# Fill with 1s where the gene is present
for (set_name in names(gene_list)) {
  df[gene_list[[set_name]], set_name] <- 1
}

# ---------------------------------------------------------
# 2. Extract EXCLUSIVE genes matching the UpSet bars
# ---------------------------------------------------------
target_gene_lists <- list()

# These are all the sets included in the plot
all_sets <- ordered_sets 

for(col in ordered_p_bio) {
  
  # Identify sets that should NOT contain the gene for these specific bars
  other_sets <- setdiff(all_sets, c(col, "Full_bio"))
  
  # --- 1. THE OVERLAP: (1 in col, 1 in Full_bio, 0 in ALL other sets) ---
  is_overlap <- df[[col]] == 1 & 
                df$Full_bio == 1 & 
                rowSums(df[, other_sets, drop = FALSE]) == 0
                
  overlap_genes <- rownames(df)[is_overlap]
  
  if(length(overlap_genes) > 20) {
    target_gene_lists[[paste0(col, "_Overlap")]] <- overlap_genes
  }
  
  # --- 2. THE UNIQUE: (1 in col, 0 in Full_bio, 0 in ALL other sets) ---
  is_unique <- df[[col]] == 1 & 
               df$Full_bio == 0 & 
               rowSums(df[, other_sets, drop = FALSE]) == 0
               
  unique_genes <- rownames(df)[is_unique]
  
  if(length(unique_genes) > 20) {
    target_gene_lists[[paste0(col, "_Unique")]] <- unique_genes
  }
}

# Make sure your output directory exists
if(!dir.exists("GO_Results_Intersections")) dir.create("GO_Results_Intersections")

fasta <- readDNAStringSet("./in/GDDH13_1-1_mrna.fasta")

universe_genes <- unique(
  sub(".*gene[:=]([^;| ]+).*", "\\1", names(fasta))
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

semdata <- godata(
  OrgDb  = org.Mdomestica.eg.db,
  keytype = "GID",
  ont     = "BP"
)

# Initialize the plot storage lists
dotplot_list <- list()
emapplot_list <- list()
treeplot_list <- list()

# --- 1. RUN GO ENRICHMENT FOR EACH INTERSECTION ---
for (label in names(target_gene_lists)) {
  current_genes <- target_gene_lists[[label]]
  
  message(paste("--- Processing Intersection:", label, "| Genes:", length(current_genes), "---"))
  
  # Run Enrichment
  go_bp <- enrichGO(
    gene          = current_genes, 
    universe      = universe_genes,
    OrgDb         = org.Mdomestica.eg.db, 
    keyType       = "GID", 
    ont           = "BP", 
    pvalueCutoff  = 0.05, 
    qvalueCutoff  = 0.05, 
    readable      = TRUE,
    pAdjustMethod = "fdr"
  )#; print(go_bp)
  
  if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
    
    # Extract gene count for title to verify it matches UpSet
    n_genes <- length(current_genes)
    plot_title <- paste(label, "(n =", n_genes, ")")
    
    # A. DOTPLOT (Remove legend immediately)
    dotplot_list[[label]] <- dotplot(go_bp, showCategory = 10) + 
                             labs(title = plot_title) +
                             theme(legend.position = "none") 
    
    # B. EMAPPLOT & TREEPLOT (Need > 2 terms for similarity calculation)
    if (nrow(as.data.frame(go_bp)) > 2) {
      
      # Assuming you have safe_pairwise_termsim defined, otherwise use pairwise_termsim(go_bp)
      go_bp_sim <- pairwise_termsim(go_bp) 
      
      emapplot_list[[label]] <- emapplot(go_bp_sim, showCategory = 20) + 
                                labs(title = plot_title) +
                                theme(legend.position = "none")
      
      n_terms <- nrow(as.data.frame(go_bp_sim))
      if (n_terms >= 2) {
        n_clusters <- min(5, n_terms - 1)
        treeplot_list[[label]] <- suppressWarnings(
          treeplot(go_bp_sim, cluster.params = list(n = max(1, n_clusters)))
        ) + labs(title = plot_title) +
          theme(legend.position = "none") 
      }
    }
  }
}

# --- 2. THE COLLECTOR LOGIC ---

# We use ncol = 3 for the grids since there are many intersections

# --- A. DOTPLOTS ---
if (length(dotplot_list) > 0) {
  message("Generating Master Dotplot Grid...")
  
  # Calculate global max count for size consistency
  all_counts <- unlist(lapply(dotplot_list, function(p) p$data$Count))
  max_c <- if(length(all_counts) > 0) max(all_counts) else 50
  
  master_legend_plot <- dotplot_list[[1]] + 
    theme(legend.position = "left") + 
    scale_color_continuous(low = "red", high = "blue", limits = c(0, 0.05), name = "p.adjust") +
    scale_size_continuous(limits = c(1, max_c), name = "Count")
  
  shared_legend <- cowplot::get_legend(master_legend_plot)
  
  combined_grid <- wrap_plots(dotplot_list, ncol = 3) + 
                   plot_annotation(title = "GO BP Dotplots: UpSet Intersections")
  
  final_dot <- wrap_plots(shared_legend, combined_grid, widths = c(1, 15))
  
  n_rows <- ceiling(length(dotplot_list) / 3)
  ggsave(file.path("GO_Results_Intersections", "All_Intersections_Dotplots.png"), 
         final_dot, width = 24, height = 5 * n_rows, dpi = 300, limitsize = FALSE)
}

# --- B. EMAPPLOTS ---
if (length(emapplot_list) > 0) {
  message("Generating Master Emapplot Grid...")
  
  master_emap_legend <- emapplot_list[[1]] + 
    theme(legend.position = "left") + 
    scale_color_continuous(low = "red", high = "blue", limits = c(0, 0.05), name = "p.adjust")
  
  shared_emap_legend <- cowplot::get_legend(master_emap_legend)
  
  combined_emap_grid <- wrap_plots(emapplot_list, ncol = 3) + 
                        plot_annotation(title = "Enrichment Maps: UpSet Intersections")
  
  final_emap <- wrap_plots(shared_emap_legend, combined_emap_grid, widths = c(1, 15))
  
  n_rows_emap <- ceiling(length(emapplot_list) / 3)
  ggsave(file.path("GO_Results_Intersections", "All_Intersections_Emaps.png"), 
         final_emap, width = 24, height = 7 * n_rows_emap, dpi = 300, limitsize = FALSE)
}

# --- C. TREEPLOTS ---
if (length(treeplot_list) > 0) {
  message("Generating Master Treeplot Grid...")
  
  master_tree_legend <- treeplot_list[[1]] + 
    theme(legend.position = "left") + 
    scale_color_continuous(low = "red", high = "blue", limits = c(0, 0.05), name = "p.adjust")
  
  shared_tree_legend <- cowplot::get_legend(master_tree_legend)
  
  combined_tree_grid <- wrap_plots(treeplot_list, ncol = 2) + # Using 2 cols for trees as they need width
                        plot_annotation(title = "Treeplots: UpSet Intersections")
  
  final_tree <- wrap_plots(shared_tree_legend, combined_tree_grid, widths = c(1, 12))
  
  n_rows_tree <- ceiling(length(treeplot_list) / 2)
  ggsave(file.path("GO_Results_Intersections", "All_Intersections_Trees.png"), 
         final_tree, width = 30, height = 8 * n_rows_tree, dpi = 300, limitsize = FALSE)
}

message("All GO plots generated and saved successfully!")

# Group by axes and take overlapping gene across Temp and Prec

# Define the column names for each category
# sprintf("%02d", 1:11) ensures we get "01", "02" instead of "1", "2"
temp_vars <- paste0("p_bio", sprintf("%02d", 1:11)) 
precip_vars <- paste0("p_bio", 12:19)               

# Create the new columns: if the row sum > 0, the gene is in the group
df$Temperature <- ifelse(rowSums(df[, temp_vars, drop = FALSE]) > 0, 1, 0)
df$Precipitation <- ifelse(rowSums(df[, precip_vars, drop = FALSE]) > 0, 1, 0)

sets_sizes <- colSums(df[, c("Full_bio", "Temperature", "Precipitation")])
ordered_sets_sizes <- names(sort(sets_sizes, decreasing = TRUE))

png("Temp_vs_Precip_vs_Full_bio.png", width = 2000, height = 1200, res = 200)

upset(df, 
      sets = ordered_sets_sizes, 
      keep.order = TRUE,
      main.bar.color = "steelblue",
      sets.bar.color = "gray30",
      matrix.color = "black",
      text.scale = 1.5,
      order.by = "freq") # Orders the bars by size

dev.off()

# Initialize a new list for these specific macro-intersections
macro_gene_lists <- list()

# 1. Temperature ONLY (Not in Precip, Not in Full_bio)
macro_gene_lists[["Temp_Only"]] <- rownames(df)[df$Temperature == 1 & df$Precipitation == 0 & df$Full_bio == 0]

# 2. Precipitation ONLY (Not in Temp, Not in Full_bio)
macro_gene_lists[["Precip_Only"]] <- rownames(df)[df$Temperature == 0 & df$Precipitation == 1 & df$Full_bio == 0]

# 3. Temp & Precip Overlap (But NOT in Full_bio)
macro_gene_lists[["Temp_Precip_Overlap"]] <- rownames(df)[df$Temperature == 1 & df$Precipitation == 1 & df$Full_bio == 0]

# 4. Temp & Full_bio Overlap (Not in Precip)
macro_gene_lists[["Temp_Full_Overlap"]] <- rownames(df)[df$Temperature == 1 & df$Precipitation == 0 & df$Full_bio == 1]

# 5. Precip & Full_bio Overlap (Not in Temp)
macro_gene_lists[["Precip_Full_Overlap"]] <- rownames(df)[df$Temperature == 0 & df$Precipitation == 1 & df$Full_bio == 1]

# 6. ALL THREE: Temp + Precip + Full_bio
macro_gene_lists[["All_Three_Overlap"]] <- rownames(df)[df$Temperature == 1 & df$Precipitation == 1 & df$Full_bio == 1]

# Filter out lists with 20 or fewer genes to save processing time
macro_gene_lists <- macro_gene_lists[sapply(macro_gene_lists, length) > 20]


library(ggplot2)

# 1. Create a summary dataframe using the lengths of your specific lists
summary_macro_df <- data.frame(
  Condition = factor(
    c("Temp spe", "Temp spe", "Prec spe", "Prec spe"),
    levels = c("Temp spe", "Prec spe") # Ensures Temp is on the left, Prec on the right
  ),
  Status = factor(
    c("Not in Full_bio", "In Full_bio", "Not in Full_bio", "In Full_bio"),
    levels = c("Not in Full_bio", "In Full_bio") # Stack order
  ),
  Count = c(
    length(macro_gene_lists[["Temp_Only"]]),           # Temp spe, Not in Full_bio
    length(macro_gene_lists[["Temp_Full_Overlap"]]),   # Temp spe, In Full_bio
    length(macro_gene_lists[["Precip_Only"]]),         # Prec spe, Not in Full_bio
    length(macro_gene_lists[["Precip_Full_Overlap"]])  # Prec spe, In Full_bio
  )
)

# 2. Generate the stacked bar chart
png("Temp_Prec_Specific_vs_Full_bio.png", width = 1600, height = 1500, res = 200)

ggplot(summary_macro_df, aes(x = Condition, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  # Add exact numbers inside the bars
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), 
            size = 6, color = "black") +
  # Colors matching the previous plots
  scale_fill_manual(values = c("Not in Full_bio" = "darkorange", "In Full_bio" = "steelblue")) +
  theme_minimal(base_size = 15) +
  labs(
    title = "Temperature vs. Precipitation Specific Genes",
    subtitle = "Separated by presence in Full_bio reference list",
    x = "Environmental Category",
    y = "Number of Genes",
    fill = "Reference Status"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.text.x = element_text(face = "bold", size = 14)
  )

dev.off()

library(patchwork)
library(cowplot)

# --- SETUP DIRECTORY ---
if(!dir.exists("GO_Results_Macro")) dir.create("GO_Results_Macro")

# Initialize storage lists for plots
dotplot_list_macro <- list()
emapplot_list_macro <- list()
treeplot_list_macro <- list()

# --- 1. RUN GO ENRICHMENT FOR MACRO INTERSECTIONS ---
for (label in names(macro_gene_lists)) {
  current_genes <- macro_gene_lists[[label]]
  
  message(paste("--- Processing Macro Group:", label, "| Genes:", length(current_genes), "---"))
  
  # Run Enrichment
  go_bp <- enrichGO(
    gene          = current_genes, 
    universe      = universe_genes,
    OrgDb         = org.Mdomestica.eg.db, 
    keyType       = "GID", 
    ont           = "BP", 
    pvalueCutoff  = 0.05, 
    qvalueCutoff  = 0.05, 
    readable      = TRUE,
    pAdjustMethod = "fdr"
  )
  
  if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
    
    # Extract exact gene count for the title
    n_genes <- length(current_genes)
    plot_title <- paste(label, "\n(n =", n_genes, ")")
    
    # A. DOTPLOT (Remove legend immediately)
    dotplot_list_macro[[label]] <- dotplot(go_bp, showCategory = 10) + 
                                   labs(title = plot_title) +
                                   theme(legend.position = "none",
                                         plot.title = element_text(size = 12, face = "bold"))
    
    # B. EMAPPLOT & TREEPLOT (Need > 2 terms to calculate similarity)
    if (nrow(as.data.frame(go_bp)) > 2) {
      
      go_bp_sim <- pairwise_termsim(go_bp) 
      
      emapplot_list_macro[[label]] <- emapplot(go_bp_sim, showCategory = 20) + 
                                      labs(title = plot_title) +
                                      theme(legend.position = "none",
                                            plot.title = element_text(size = 12, face = "bold"))
      
      n_terms <- nrow(as.data.frame(go_bp_sim))
      if (n_terms >= 2) {
        n_clusters <- min(5, n_terms - 1)
        treeplot_list_macro[[label]] <- suppressWarnings(
          treeplot(go_bp_sim, cluster.params = list(n = max(1, n_clusters)))
        ) + labs(title = plot_title) +
            theme(legend.position = "none",
                  plot.title = element_text(size = 12, face = "bold")) 
      }
    }
  } else {
    message(paste("    No significant GO terms found for", label))
  }
}

# --- 2. COLLECT AND PLOT ---
# Since there is a maximum of 6 macro lists, ncol = 3 is perfect (max 2 rows).

# A. DOTPLOTS
if (length(dotplot_list_macro) > 0) {
  message("Generating Macro Dotplot Grid...")
  
  all_counts <- unlist(lapply(dotplot_list_macro, function(p) p$data$Count))
  max_c <- if(length(all_counts) > 0) max(all_counts) else 50
  
  master_legend <- dotplot_list_macro[[1]] + 
    theme(legend.position = "left") + 
    scale_color_continuous(low = "red", high = "blue", limits = c(0, 0.05), name = "p.adjust") +
    scale_size_continuous(limits = c(1, max_c), name = "Count")
  
  shared_legend <- get_legend(master_legend)
  
  combined_grid <- wrap_plots(dotplot_list_macro, ncol = 3) + 
                   plot_annotation(title = "GO BP Dotplots: Macro Intersections (Temp vs Precip)")
  
  final_dot <- wrap_plots(shared_legend, combined_grid, widths = c(1, 15))
  
  n_rows <- ceiling(length(dotplot_list_macro) / 3)
  ggsave(file.path("GO_Results_Macro", "Macro_Dotplots.png"), 
         final_dot, width = 22, height = 5 * n_rows, dpi = 300)
}

# B. EMAPPLOTS
if (length(emapplot_list_macro) > 0) {
  message("Generating Macro Emapplot Grid...")
  
  master_emap_legend <- emapplot_list_macro[[1]] + 
    theme(legend.position = "left") + 
    scale_color_continuous(low = "red", high = "blue", limits = c(0, 0.05), name = "p.adjust")
  
  shared_emap_legend <- get_legend(master_emap_legend)
  
  combined_emap_grid <- wrap_plots(emapplot_list_macro, ncol = 3) + 
                        plot_annotation(title = "Enrichment Maps: Macro Intersections")
  
  final_emap <- wrap_plots(shared_emap_legend, combined_emap_grid, widths = c(1, 15))
  
  n_rows_emap <- ceiling(length(emapplot_list_macro) / 3)
  ggsave(file.path("GO_Results_Macro", "Macro_Emaps.png"), 
         final_emap, width = 22, height = 7 * n_rows_emap, dpi = 300)
}

# C. TREEPLOTS
if (length(treeplot_list_macro) > 0) {
  message("Generating Macro Treeplot Grid...")
  
  master_tree_legend <- treeplot_list_macro[[1]] + 
    theme(legend.position = "left") + 
    scale_color_continuous(low = "red", high = "blue", limits = c(0, 0.05), name = "p.adjust")
  
  shared_tree_legend <- get_legend(master_tree_legend)
  
  combined_tree_grid <- wrap_plots(treeplot_list_macro, ncol = 2) + # Trees need more horizontal space
                        plot_annotation(title = "Treeplots: Macro Intersections")
  
  final_tree <- wrap_plots(shared_tree_legend, combined_tree_grid, widths = c(1, 12))
  
  n_rows_tree <- ceiling(length(treeplot_list_macro) / 2)
  ggsave(file.path("GO_Results_Macro", "Macro_Trees.png"), 
         final_tree, width = 26, height = 8 * n_rows_tree, dpi = 300)
}

message("--- Macro GO plots saved to 'GO_Results_Macro' ---")


# library(UpSetR)

# # 1. Define your variable groupings (if not already in your environment)
# temp_vars <- paste0("p_bio", sprintf("%02d", 1:11)) 
# precip_vars <- paste0("p_bio", 12:19)               

# # 2. Create the new "Minimum 3" columns
# # rowSums calculates how many of the specific bios a gene is found in.
# # We then use ifelse to assign a 1 if that sum is >= 3, and a 0 otherwise.
# df$Temp_min3 <- ifelse(rowSums(df[, temp_vars, drop = FALSE]) >= 3, 1, 0)
# df$Precip_min3 <- ifelse(rowSums(df[, precip_vars, drop = FALSE]) >= 3, 1, 0)

# # 3. Define the sets we want to plot
# target_sets <- c("Full_bio", "Temp_min3", "Precip_min3")

# # 4. Order the sets by size (Horizontal bars)
# set_sizes <- colSums(df[, target_sets])
# ordered_target_sets <- names(sort(set_sizes, decreasing = TRUE))

# # 5. Generate and save the UpSet plot
# png("Temp_Precip_min3_vs_Full_bio.png", width = 2500, height = 1500, res = 200)

# upset(df, 
#       sets = ordered_target_sets, 
#       keep.order = TRUE,             # Keeps the order we just calculated
#       main.bar.color = "steelblue",
#       sets.bar.color = "gray30",
#       matrix.color = "black",
#       text.scale = c(2, 1.5, 1.5, 1.5, 1.5, 1.3), # Make labels slightly larger
#       order.by = "freq")             # Orders vertical intersection bars from largest to smallest

# dev.off()

library(UpSetR)

# 1. Define your individual variables
temp_vars <- paste0("p_bio", sprintf("%02d", 1:11)) 
precip_vars <- paste0("p_bio", 12:19)
#all_bio_vars <- c("Full_bio", temp_vars, precip_vars) # Include Full_bio if you want it
all_bio_vars <- c(temp_vars, precip_vars)

# 2. Filter the dataframe
# We only keep rows (genes) that have >= 3 in Temp OR >= 3 in Precip
temp_df <-  df[rowSums(df[, temp_vars, drop = FALSE]) >= 3, ]; print(nrow(temp_df))
prec_df <- df[rowSums(df[, precip_vars, drop = FALSE]) >= 3, ]; print(nrow(prec_df))
both_df <- df[rowSums(df[, temp_vars, drop = FALSE]) >= 3  & 
                rowSums(df[, precip_vars, drop = FALSE]) >= 3, ]; print(nrow(both_df))

library(ggplot2)

# 1. Calculate the counts for each dataframe
summary_df <- data.frame(
  Condition = factor(
    rep(c("Temp (>= 3 bios)", "Precip (>= 3 bios)", "Both"), each = 2),
    levels = c("Temp (>= 3 bios)", "Precip (>= 3 bios)", "Both") # Keeps them in this order
  ),
  Status = factor(
    rep(c("In Full_bio", "Not in Full_bio"), 3),
    levels = c("Not in Full_bio", "In Full_bio") # 'Not in Full_bio' at the bottom of the stack
  ),
  Count = c(
    sum(temp_df$Full_bio == 1), sum(temp_df$Full_bio == 0),
    sum(prec_df$Full_bio == 1), sum(prec_df$Full_bio == 0),
    sum(both_df$Full_bio == 1), sum(both_df$Full_bio == 0)
  )
)

# 2. Generate the stacked bar chart
png("Climate_Conditions_vs_Full_bio.png", width = 1800, height = 1500, res = 200)

ggplot(summary_df, aes(x = Condition, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  # Add exact numbers inside the bars for easy reading
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 5, color = "black", fontface = "bold") +
  # Use colors matching your UpSet plots
  scale_fill_manual(values = c("Not in Full_bio" = "darkorange", "In Full_bio" = "steelblue")) +
  theme_minimal(base_size = 15) +
  labs(
    title = "Genes Found in Multiple Climate Variables",
    subtitle = "Separated by presence in Full_bio reference list",
    x = "Environmental Condition",
    y = "Total Number of Genes",
    fill = "Reference Status"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.text.x = element_text(face = "bold")
  )

dev.off()

library(UpSetR)

# ---------------------------------------------------------
# 1. UpSet Plot: Temperature Condition (>= 3 Temp bios)
# ---------------------------------------------------------
# Focuses on Full_bio and the temperature variables
png("UpSet_Temp_Condition.png", width = 3500, height = 2000, res = 200)

upset(temp_df, 
      sets = c("Full_bio", temp_vars), 
      keep.order = FALSE,           # Sorts horizontal set bars by size
      nintersects = 100,             # Shows only the top 30 intersections (adjust as needed)
      order.by = "freq",            # Sorts vertical intersection bars by size
      mainbar.y.label = "Intersection Size", 
      sets.x.label = "Genes (Temp >= 3)",
      matrix.color = "black",
      main.bar.color = "black",
      sets.bar.color = "gray30",
      text.scale = c(1.8, 1.3, 1.3, 1.3, 1.5, 1.2))

dev.off()

# ---------------------------------------------------------
# 2. UpSet Plot: Precipitation Condition (>= 3 Precip bios)
# ---------------------------------------------------------
# Focuses on Full_bio and the precipitation variables
png("UpSet_Precip_Condition.png", width = 3500, height = 2000, res = 200)

upset(prec_df, 
      sets = c("Full_bio", precip_vars), 
      keep.order = FALSE, 
      nintersects = 100,
      order.by = "freq",
      mainbar.y.label = "Intersection Size", 
      sets.x.label = "Genes (Precip >= 3)",
      matrix.color = "black",
      main.bar.color = "black",
      sets.bar.color = "gray30",
      text.scale = c(1.8, 1.3, 1.3, 1.3, 1.5, 1.2))

dev.off()

# ---------------------------------------------------------
# 3. UpSet Plot: BOTH Conditions (>= 3 Temp AND >= 3 Precip)
# ---------------------------------------------------------
# Focuses on Full_bio and ALL bio variables
png("UpSet_Both_Conditions.png", width = 4500, height = 2500, res = 200)

upset(both_df, 
      sets = c("Full_bio", all_bio_vars), 
      keep.order = FALSE, 
      nintersects = 100,             # Increased since there are more variables
      order.by = "freq",
      mainbar.y.label = "Intersection Size", 
      sets.x.label = "Genes (Temp & Precip >= 3)",
      matrix.color = "black",
      main.bar.color = "black",
      sets.bar.color = "gray30",
      text.scale = c(2, 1.5, 1.5, 1.5, 1.5, 1.3))

dev.off()


robust_df <- df[rowSums(df[, temp_vars, drop = FALSE]) >= 3 | 
                rowSums(df[, precip_vars, drop = FALSE]) >= 3, ]

message(paste("Number of genes meeting the criteria:", nrow(robust_df)))

# 3. Calculate sizes and sort the sets for this specific subset of genes
set_sizes <- colSums(robust_df[, all_bio_vars])

# Optional: filter out variables that have 0 genes in this specific subset to keep the plot clean
active_vars <- names(set_sizes[set_sizes > 0])
ordered_active_sets <- names(sort(set_sizes[active_vars], decreasing = TRUE))

# 4. Generate the UpSet plot showing the individual bioclimatic variables
png("Robust_Genes_Individual_Bios_Overlap.png", width = 3500, height = 1500, res = 200)

upset(robust_df, 
      sets = ordered_active_sets, 
      keep.order = TRUE,
      nintersects = 60,              # Show the top 40 overlapping patterns
      order.by = "freq",             # Order vertical bars by size
      main.bar.color = "steelblue",
      sets.bar.color = "gray30",
      matrix.color = "black",
      text.scale = 1.5,
      point.size = 3, 
      line.size = 1)

dev.off()










library(ComplexUpset)

# Convert list to data frame
df = as.data.frame(fromList(gene_list))
set_names = names(gene_list)

# Define the logic column for coloring
df$Category <- ifelse(df$Full_bio == 1, "Overlap with Full_bio",
               ifelse(rowSums(df[, set_names]) > 1, "Other Overlaps", "No Overlap"))

# 3. Generate the UpSet plot
# 'nsets' is set to 19 to show all variables
# 1. Open the PDF device
png("Upset_Plot_Genes.png", width = 2000, height = 1000, res = 150)

# 2. Generate the plot
upset(fromList(gene_list), nsets = 19, order.by = "freq", nintersects=60)

# 3. Close the device
dev.off()



# 1. Load libraries
library(ComplexUpset)
library(ggplot2)

# 2. Prepare the Data Frame
# Convert your list to a data frame and identify the sets
df <- as.data.frame(fromList(gene_list))
set_names <- names(gene_list)

# 3. Create the Category Logic
# We calculate this row by row for the data frame
df$Category <- ifelse(
    df$Full_bio == 1, 
    "Overlap with Full_bio",
    ifelse(rowSums(df[, set_names]) > 1, "Other Overlaps", "No Overlap")
)

# 4. Open the PNG device
png("Upset_Plot_Genes_Complex.png", width = 2500, height = 1200, res = 150)

# 5. Generate the Plot
# In ComplexUpset, we use the 'mapping' inside intersection_size to color bars
upset(
    df, 
    set_names,
    n_intersections = 60,
    width_ratio = 0.1,
    base_annotations = list(
        'Intersection size' = intersection_size(
            counts = TRUE,
            mapping = aes(fill = Category) # This looks at the Category column we made
        )
    ),
    # This moves the dots/matrix to the bottom and allows styling
    annotations = list() 
) + 
scale_fill_manual(values = c(
    "Overlap with Full_bio" = "#E41A1C", # Red
    "Other Overlaps" = "#377EB8",       # Blue
    "No Overlap" = "#999999"             # Gray
)) +
theme(legend.position = "bottom")

# 6. Close the device
dev.off()



data_matrix <- fromList(gene_list)



# 2. Get all existing combinations in the data
# This finds every intersection that actually has at least one gene
all_intersections <- unique(apply(data_matrix, 1, function(x) names(x)[x == 1]))

# 1. Define the name of your target variable
target_var <- "Full_bio" 

# 2. Create a function to assign priority scores
get_priority <- function(x) {
  degree <- length(x)
  has_target <- target_var %in% x
  
  if (degree == 1) {
    return(1) # Specific/Singletons first
  } else if (has_target) {
    return(2) # Overlaps involving Full_bio second
  } else {
    return(3) # All other overlaps last
  }
}

# 3. Apply the priority and then sub-sort by degree/size
# We use a data frame to handle multiple sorting levels easily
sort_df <- data.frame(
  index = 1:length(all_intersections),
  priority = sapply(all_intersections, get_priority),
  degree = sapply(all_intersections, length)
)

# Sort by Priority (1, 2, 3) then by Degree (smallest to largest)
final_order_indices <- sort_df[order(sort_df$priority, sort_df$degree), "index"]
ordered_intersections <- all_intersections[final_order_indices]

# 4. Generate the plot
png("Upset_Ordered_By_Specificity.png", width = 5000, height = 1000, res = 150)

upset(data_matrix, 
      intersections = ordered_intersections, 
      sets = p_bio_cols,           # Keep side bars in order p_bio01 to p_bio19
      keep.order = TRUE,           # CRITICAL: Forces the plot to use our sorted list
      sets.bar.color = "gray50",
      main.bar.color = "steelblue",
      text.scale = 1.2,
      mb.ratio = c(0.6, 0.4))

dev.off()


# 1. Identify which specific intersection pattern each gene belongs to
# We create a "key" for each gene (e.g., "p_bio01:p_bio02")
intersection_keys <- apply(data_matrix, 1, function(row) {
  paste(sort(names(data_matrix)[row == 1]), collapse = ":")
})
counts_table <- table(intersection_keys)

filtered_intersections <- ordered_intersections[sapply(ordered_intersections, function(x) {
  key <- paste(sort(x), collapse = ":")
  if (key %in% names(counts_table)) {
    return(counts_table[key] >= 20)
  } else {
    return(FALSE)
  }
})]
# 4. Generate the Plot
png("Upset_Filtered_20plus.png", width = 2000, height = 1000, res = 150)

upset(data_matrix, 
      intersections = filtered_intersections, 
      sets = p_bio_cols,           # Keep side bars in order
      keep.order = TRUE,           # Keeps specifics first, then combinations
      sets.bar.color = "gray50",
      main.bar.color = "steelblue",
      text.scale = 1.2,
      mb.ratio = c(0.6, 0.4))

dev.off()

# install.packages("ComplexUpset")
library(ComplexUpset)
library(ggplot2)

# 1. Use your existing logic to get the sets (column names)
all_sets <- colnames(data_matrix)

# 2. Generate the plot using upset() from ComplexUpset
# Note: ComplexUpset uses 'intersections' to filter/sort.
p <- upset(
    data_matrix,
    all_sets,
    intersections = ordered_intersections, # Your sorted list from before
    sort_intersections = FALSE,            # THIS IS KEY: Tells it not to re-sort
    width_ratio = 0.1,
    base_annotations = list(
        'Intersection size' = intersection_size(
            counts = TRUE,
            mapping = aes(fill = "steelblue")
        ) + scale_fill_identity()
    ),
    set_descriptors = list(
        'Set size' = set_size()
    )
)

# 3. Save the plot
ggsave("Upset_Ordered_By_Specificity_Fixed.png", plot = p, width = 16, height = 8)



# 1. Load libraries
library(ComplexUpset)
library(ggplot2)

# 2. Prepare the Data Frame
df <- as.data.frame(fromList(gene_list))
set_names <- names(gene_list)

# 3. Create the Category Logic
df$Category <- ifelse(
    df$Full_bio == 1, 
    "Overlap with Full_bio",
    ifelse(rowSums(df[, set_names]) > 1, "Other Overlaps", "No Overlap")
)

# 4. Open the PNG device
png("Upset_Plot_Genes_Final.png", width = 2500, height = 1500, res = 150)

# 5. Generate the Plot
upset(
    df, 
    set_names,
    n_intersections = 60,
    width_ratio = 0.1,
    
    # Color the BARS
    base_annotations = list(
        'Intersection size' = intersection_size(
            counts = TRUE,
            mapping = aes(fill = Category)
        )
    ),
    
    # Color the POINTS (The Matrix)
    # We map the 'color' aesthetic of the matrix to our Category column
    matrix = intersection_matrix(
        geom = geom_point(
            size = 3,
            mapping = aes(color = Category)
        )
    )
) + 
# Apply the colors to both 'fill' (bars) and 'color' (points)
scale_fill_manual(values = c(
    "Overlap with Full_bio" = "#E41A1C", 
    "Other Overlaps" = "#377EB8", 
    "No Overlap" = "#999999"
)) +
scale_color_manual(values = c(
    "Overlap with Full_bio" = "#E41A1C", 
    "Other Overlaps" = "#377EB8", 
    "No Overlap" = "#999999"
)) +
theme(legend.position = "bottom")

# 6. Close the device
dev.off()




# Calcul des proportions
prop_data <- combined_results %>%
  group_by(bio_source, final_type) %>%
  tally() %>%
  group_by(bio_source) %>%
  mutate(percentage = n / sum(n) * 100,
         count = n)

plot_prop <- ggplot(prop_data, aes(x = bio_source, y = percentage, fill = final_type)) +
  geom_bar(stat = "identity", color = "white") +
  scale_fill_viridis_d(option = "D", name = "Genomic Category") +
  theme_minimal() +
  labs(title = "SNP Distribution Proportion by Bio Variable",
       x = "Environmental Variable", y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Bio_Variables_Comparison_prop.png", plot_prop, width = 15, height = 20, dpi=300)

# 1. Calculer le total par bio_source pour les étiquettes
label_data <- prop_data %>%
  group_by(bio_source) %>%
  summarise(total_snps = sum(n), .groups = "drop")

# 2. Création du graphique
plot_count <- ggplot(prop_data, aes(x = bio_source, y = n, fill = final_type)) +
  # Barre empilée
  geom_bar(stat = "identity", color = "white") +
  # Ajout du texte au-dessus de la pile complète
  geom_text(data = label_data, 
            aes(x = bio_source, y = total_snps, label = comma(total_snps)),
            inherit.aes = FALSE, # On ignore le fill pour le texte global
            vjust = -0.5, 
            fontface = "bold", 
            size = 2) +
  scale_fill_viridis_d(option = "D", name = "Genomic Category") +
  theme_minimal() +
  labs(title = "Number of SNPs by Bio Variable",
       subtitle = "Total counts displayed above each bar",
       x = "Environmental Variable", 
       y = "Number of Significant SNPs") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  )

# 3. Sauvegarde
ggsave("Bio_Variables_Comparison_count.png", plot_count, width = 12, height = 8, dpi = 300)


## get the number of genes

bio_summary <- combined_results %>%
  filter(ref_type == "CDS", !is.na(gene_id_clean)) %>%
  # On pivote toutes les colonnes p_bio vers le bas
  pivot_longer(
    cols = starts_with("p_bio"), 
    names_to = "bio_variable", 
    values_to = "p_val"
  ) %>%
  # On garde uniquement les lignes où la p-value existe
  filter(!is.na(p_val)) %>%
  # On groupe par bio_variable et on compte les gènes uniques
  group_by(bio_variable) %>%
  summarise(n_unique_genes = n_distinct(gene_id_clean),
  gene_list = list(unique(gene_id_clean)))

print(bio_summary)


fasta <- readDNAStringSet("./in/GDDH13_1-1_mrna.fasta")

universe_genes <- unique(
  sub(".*gene[:=]([^;| ]+).*", "\\1", names(fasta))
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

semdata <- godata(
  OrgDb  = org.Mdomestica.eg.db,
  keytype = "GID",
  ont     = "BP"
)

bio_var_int <- c("p_bio08", "p_bio19")

library(patchwork)
library(enrichplot)

# Main Loop per Bioclimatic Variable
for (bio in bio_var_int) {
  dotplot_list <- list()
  emapplot_list <- list()
  treeplot_list <- list()
  
  message(paste("--- Processing:", bio, "---"))
  
  for (chr in chrs_to_process) {
    message(paste("--- Processing:", chr, "---"))
    
    current_genes <- bio_summary %>% 
      filter(bio_variable == bio, chromosome == chr) %>% 
      pull(gene_list) %>% unlist()
    
    if (length(current_genes) == 0) next
    
    go_bp <- enrichGO(
      gene = current_genes, universe = universe_genes,
      OrgDb = org.Mdomestica.eg.db, keyType = "GID",
      ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
    )
    
    if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
      
      # --- STRATEGY: Remove legends immediately ---
      dotplot_list[[chr]] <- dotplot(go_bp, showCategory = 10) + 
                             labs(title = chr) +
                             theme(legend.position = "none") # Kill the legend here
      
      if (nrow(as.data.frame(go_bp)) > 2) {
        go_bp_sim <- safe_pairwise_termsim(go_bp, semdata)
        
        emapplot_list[[chr]] <- emapplot(go_bp_sim, showCategory = 20) + 
                                labs(title = chr) +
                                theme(legend.position = "none") # Kill the legend here
        
        n_terms <- nrow(as.data.frame(go_bp_sim))
        if (n_terms >= 2) {
          n_clusters <- min(5, n_terms - 1)
          treeplot_list[[chr]] <- suppressWarnings(
            treeplot(go_bp_sim, cluster.params = list(n = max(1, n_clusters)))
          ) + labs(title = chr) +
            theme(legend.position = "none") # Kill the legend here
        }
      }
    }
  }
  
  # --- 4. THE COLLECTOR LOGIC ---
  
  if (length(dotplot_list) > 0) {
    # 1. Calculate global max count for size consistency
    all_counts <- unlist(lapply(dotplot_list, function(p) p$data$Count))
    max_c <- if(length(all_counts) > 0) max(all_counts) else 50
    
    # 2. Create the Master Legend Plot
    # We ensure this plot has the legend enabled and the scales are fixed
    master_legend_plot <- dotplot_list[[1]] + 
      theme(legend.position = "left") + 
      scale_color_continuous(low = "red", high = "blue", limits = c(0, 0.05), name = "p.adjust") +
      scale_size_continuous(limits = c(1, max_c), name = "Count")
    
    # 3. Extract the legend
    # Note: Requires the cowplot package
    shared_legend <- cowplot::get_legend(master_legend_plot)
    
    # 4. Create the Grid (2 columns per row)
    # Ensure all plots in the list have their legends hidden
    dotplot_list_no_legend <- lapply(dotplot_list, function(p) p + theme(legend.position = "none"))
    
    combined_grid <- wrap_plots(dotplot_list_no_legend, ncol = 2) + 
                     plot_annotation(title = paste("Dotplots for", bio))
    
    # 5. Assemble: Legend on the LEFT, Grid on the RIGHT
    # widths = c(1, 10) gives the legend 1 part space and the grid 10 parts
    final_dot <- wrap_plots(shared_legend, combined_grid, widths = c(1, 10))
    
    # 6. Save with dynamic height
    n_rows <- ceiling(length(dotplot_list) / 2)
    ggsave(file.path("GO_Results_Grouped", paste0(bio, "_dotplots_grouped.png")), 
           final_dot, width = 18, height = 5 * n_rows, dpi = 300)
  }
  
  # --- B. EMAPPLOTS (Network Maps) ---
  if (length(emapplot_list) > 0) {
    # 1. Create Master Legend for Emap (p.adjust only)
    master_emap_legend <- emapplot_list[[1]] + 
      theme(legend.position = "left") + 
      scale_color_continuous(low = "red", high = "blue", limits = c(0, 0.05), name = "p.adjust")
    
    shared_emap_legend <- cowplot::get_legend(master_emap_legend)
    
    # 2. Prepare Grid (2 columns, legends removed)
    emap_list_no_legend <- lapply(emapplot_list, function(p) p + theme(legend.position = "none"))
    
    combined_emap_grid <- wrap_plots(emap_list_no_legend, ncol = 2) + 
                         plot_annotation(title = paste("Enrichment Maps for", bio))
    
    # 3. Assemble: Legend Left
    final_emap <- wrap_plots(shared_emap_legend, combined_emap_grid, widths = c(1, 10))
    
    # 4. Save (Emaps need more height, approx 7 inches per row)
    n_rows_emap <- ceiling(length(emapplot_list) / 2)
    ggsave(file.path("GO_Results_Grouped", paste0(bio, "_emaps_grouped.png")), 
           final_emap, width = 20, height = 7 * n_rows_emap, dpi = 300)
  }

  # --- C. TREEPLOTS (Clustered Dendrograms) ---
  if (length(treeplot_list) > 0) {
    # 1. Create Master Legend for Tree (p.adjust only)
    master_tree_legend <- treeplot_list[[1]] + 
      theme(legend.position = "left") + 
      scale_color_continuous(low = "red", high = "blue", limits = c(0, 0.05), name = "p.adjust")
    
    shared_tree_legend <- cowplot::get_legend(master_tree_legend)
    
    # 2. Prepare Grid (2 columns, legends removed)
    tree_list_no_legend <- lapply(treeplot_list, function(p) p + theme(legend.position = "none"))
    
    combined_tree_grid <- wrap_plots(tree_list_no_legend, ncol = 2) + 
                         plot_annotation(title = paste("Treeplots for", bio))
    
    # 3. Assemble: Legend Left
    final_tree <- wrap_plots(shared_tree_legend, combined_tree_grid, widths = c(1, 8))
    
    # 4. Save (Trees are dense, approx 8 inches per row)
    n_rows_tree <- ceiling(length(treeplot_list) / 2)
    ggsave(file.path("GO_Results_Grouped", paste0(bio, "_trees_grouped.png")), 
           final_tree, width = 35, height = 8 * n_rows_tree, dpi = 300)
  }
}


# 1. Install and load the library
if (!require("vcfR")) install.packages("vcfR")
library(vcfR)

# 2. Load your VCF file
# Replace 'your_file.vcf' with your actual filename
vcf <- read.vcfR("/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/06-TIPs/out/tips_dom_admixed_rm_V3/s01.filter_vcf/GraffiTE.merged.genotypes_tips_dom_admixed_rm_V3.vcf.gz")

# 3. Extract the fixed columns (CHROM, POS, ID, REF, ALT, etc.)
# This creates a matrix of the first 8 columns
vcf_fix <- as.data.frame(getFIX(vcf))

# 4. Select only CHROM and POS (and ID if you need it)
positions <- vcf_fix[, c("CHROM", "POS", "ID")] 

# 5. View the result
head(positions)

all_TIPS_regions <- positions %>%
  # 1. Format CHROM to match GFF (e.g., chr1 -> Chr01, chr10 -> Chr10)
  mutate(CHROM = str_replace(CHROM, "chr", "Chr")) %>%
  mutate(CHROM = case_when(
    # If it ends in a single digit (e.g., Chr1), add a 0 (Chr01)
    str_detect(CHROM, "^Chr\\d$") ~ str_replace(CHROM, "Chr", "Chr0"),
    # If it is Chr0, make it Chr00
    CHROM == "Chr0" ~ "Chr00",
    # Otherwise keep as is (e.g., Chr10 stays Chr10)
    TRUE ~ CHROM
  )) %>%
  # 2. Filter out the specific chromosome (now named Chr00)
  filter(CHROM != "Chr00") %>%
  # 3. Sort by Chromosome and Position
  group_by(CHROM) %>%
  arrange(POS, .by_group = TRUE) %>%
  ungroup()

# Check the new names
unique(all_TIPS_regions$CHROM)

# 1. First, ensure POS is numeric (if you haven't already)
all_TIPS_regions <- all_TIPS_regions %>%
  mutate(POS = as.integer(POS))

# 2. Join the entire combined_results with all_TIPS_regions
# This finds matches for ALL intergenic types at once
all_matches <- combined_results %>%
  inner_join(all_TIPS_regions, 
             by = c("chromosome" = "CHROM", "SNPs" = "POS"))

# 3. Create the summary table
final_summary_table <- all_matches %>%
  # Keep only the intergenic types you are interested in
  filter(final_type %in% c("Intergenic < 2kb", "Intergenic 2-5kb", "Intergenic > 5kb")) %>%
  # Count occurrences by chromosome and type
  group_by(chromosome, final_type) %>%
  tally() %>%
  # Pivot to put the intergenic spaces in columns
  pivot_wider(names_from = final_type, values_from = n, values_fill = 0) %>%
  # Optional: Sort by chromosome name
  arrange(chromosome)

# 4. View the final result
print(final_summary_table)

# Select columns starting with p_bio and count non-NA entries
p_bio_counts <- all_matches %>%
  summarise(across(starts_with("p_bio"), ~ sum(!is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "p_bio_variable", values_to = "match_count")

print(p_bio_counts)

library(ggplot2)

# Prepare data for plotting (long format)
plot_data <- all_matches %>%
  filter(final_type %in% c("Intergenic < 2kb", "Intergenic 2-5kb", "Intergenic > 5kb"))

# Create the plot
plot_TIPs <- ggplot(plot_data, aes(x = chromosome, fill = final_type)) +
  geom_bar(position = "stack") + # Use position = "dodge" for side-by-side bars
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Distribution of TIPs by Chromosome and Intergenic Distance",
       x = "Chromosome",
       y = "Count of TIP Matches",
       fill = "Distance Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path("TIPs_match_intergenic_region.png"), 
           plot_TIPs, width = 10, height = 8, dpi = 300)


all_matches <- all_matches %>%
  mutate(Type = str_extract(ID, "INS|DEL"))

#1. Prepare the data
type_plot_data <- all_matches %>%
  # Keep p_bio columns and the new Type column
  dplyr::select(starts_with("p_bio"), Type) %>%
  # Make it long format
  pivot_longer(cols = starts_with("p_bio"), names_to = "Variable", values_to = "p_val") %>%
  # Only keep rows where a p-value exists
  filter(!is.na(p_val))

# 2. Plot
plot_TIPs_type <- ggplot(type_plot_data, aes(x = Variable, fill = Type)) +
  geom_bar(position = "fill") + # "fill" shows the percentage (0 to 100%)
  theme_minimal() +
  labs(title = "Ratio of Insertions vs. Deletions per p_bio Variable",
       y = "Proportion",
       x = "Environmental Variable") +
  scale_fill_manual(values = c("INS" = "#69b3a2", "DEL" = "#404080")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path("TIPs_type_match_intergenic_region.png"), 
           plot_TIPs_type, width = 10, height = 8, dpi = 300)


library(ggplot2)
library(dplyr)

# 1. Prepare data
combined_plot_data <- all_matches %>%
  filter(final_type %in% c("Intergenic < 2kb", "Intergenic 2-5kb", "Intergenic > 5kb")) %>%
  # Ensure Type and final_type are clean
  filter(!is.na(Type), !is.na(final_type))

# 2. Create the mixed plot
# We use 'chromosome' as X, 'fill' for distance, and 'group' or 'alpha' won't work alone, 
# so we use facet_wrap or position_dodge2
plot_combined <- ggplot(combined_plot_data, aes(x = Type, fill = final_type)) +
  geom_bar(position = "stack") + 
  # This splits the view by chromosome so you can see INS vs DEL side-by-side per Chr
  facet_wrap(~chromosome, scales = "free_x") + 
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(title = "TIP Count per Chromosome: Insertion vs Deletion by Distance",
       x = "Variant Type",
       y = "Number of TIPs",
       fill = "Intergenic Category") +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 0))

# 3. Save
ggsave("TIPs_Combined_Analysis.png", plot_combined, width = 14, height = 10, dpi = 300)


library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Reshape data to long format for p_bio variables
bio_mixed_data <- all_matches %>%
  # Filter for the intergenic categories we want
  filter(final_type %in% c("Intergenic < 2kb", "Intergenic 2-5kb", "Intergenic > 5kb")) %>%
  # Select only relevant columns
  dplyr::select(ID, Type, final_type, starts_with("p_bio")) %>%
  # Pivot p_bio columns into a single column
  pivot_longer(cols = starts_with("p_bio"), 
               names_to = "Bio_Variable", 
               values_to = "p_val") %>%
  # Keep only the rows where that specific TIP matched that Bio Variable
  filter(!is.na(p_val))

# 2. Create the plot
# This uses 'Bio_Variable' on the X-axis
# It uses 'Type' (INS/DEL) to create side-by-side bars (facet)
# It uses 'final_type' to stack colors within the bars
plot_bio_mixed <- ggplot(bio_mixed_data, aes(x = Type, fill = final_type)) +
  geom_bar(position = "stack") + 
  facet_wrap(~Bio_Variable, scales = "free_y") + 
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(title = "TIP Matches per Environmental Variable",
       subtitle = "Side-by-side Insertion/Deletion stacked by Intergenic Distance",
       x = "Variant Type",
       y = "Number of Matches",
       fill = "Distance Category") +
  theme(strip.background = element_rect(fill = "gray95"),
        axis.text.x = element_text(angle = 0))

# 3. Save the plot
ggsave("TIPs_Bio_Mixed_Analysis.png", plot_bio_mixed, width = 12, height = 10, dpi = 300)




# 1. Prepare your GFF (filter for 'gene' only)
# Assuming your GFF is loaded as 'gff_data'

  # Extract Gene ID from the attribute column


genes_only <- gff_base %>%
  filter(type == "gene") %>%
  # Extract Gene ID from the attribute column
  mutate(gene_id = str_extract(attributes, "MD\\d{2}G\\d+"))

# 2. Filter TIPs for p_bio19
pbio19_tips <- all_matches %>%
  filter(!is.na(p_bio19))

# 3. Find the nearest gene
# We will use a 'cross join' or a 'rolling join'. 
# For a simple list, we can use this logic:
nearest_genes <- pbio19_tips %>%
  rowwise() %>%
  mutate(closest_gene = genes_only %>%
           filter(chromosome == chromosome) %>%
           mutate(dist = abs(start - SNPs)) %>%
           slice_min(dist, n = 1) %>%
           pull(gene_id),
         distance_to_gene = genes_only %>%
           filter(chromosome == chromosome) %>%
           mutate(dist = abs(start - SNPs)) %>%
           slice_min(dist, n = 1) %>%
           pull(dist),
         # Determine if it's upstream or downstream
         # (Simple version: assumes + strand; for - strand, logic reverses)
         position_relative = ifelse(SNPs < genes_only$start[genes_only$gene_id == closest_gene], 
                                   "Upstream", "Downstream")) %>%
  ungroup()

# View the 52 Insertions and 24 Deletions for p_bio19
head(nearest_genes)

# 1. Extraction des gènes ayant un SNP significatif dans le CDS pour p_bio19
# On utilise !is.na(p_bio19) si c'est une colonne de p-valeurs
genes_only_test <- combined_results %>%
  filter(!is.na(p_bio19) & ref_type == "CDS") %>%
  pull(gene_id_clean) %>%
  unique()

# 2. Calcul de l'intersection
# On compare directement les deux vecteurs
genes_en_commun <- intersect(nearest_genes$closest_gene, genes_only_test)

# 3. Résultat
length(genes_en_commun)


