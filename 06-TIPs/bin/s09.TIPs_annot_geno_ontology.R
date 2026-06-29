# ==============================================================================
# 0. ENVIRONMENT PREPARATION & DEPENDENCIES
# ==============================================================================
rm(list = ls())

# Core spatial, genomic, and data manipulation libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(readxl)
  library(janitor)
  library(vcfR)
  library(IRanges)
  library(GenomicRanges)
  library(ggplot2)
  library(UpSetR)
  library(grid)
  library(rBLAST)
  library(Biostrings)
  library(org.Mdomestica.eg.db)
  library(GOplot)
  library(clusterProfiler)
  library(GOSemSim)
})

# Define Project and Global Directory paths
PROJECT  <- "tips_dom_admixed_rm_V3"
OUT_FOLD <- file.path("./out", PROJECT, "s09.tips_annot_gene_ontology")
if (!dir.exists(OUT_FOLD)) dir.create(OUT_FOLD, recursive = TRUE)

# File inputs
VCF_PATH  <- file.path("./out", PROJECT, "s01.filter_vcf", paste0("GraffiTE.merged.genotypes_", PROJECT, ".vcf.gz"))
GFF_PATH  <- "./in/gene_models_20170612.gff3"
BIO_FILES <- list.files(file.path("./out", PROJECT, "s04.lfmm2_for_tips/file/"), pattern = "bio.*\\.csv$", full.names = TRUE)
ALL_BIO_PATH <- file.path("./out", PROJECT, "s04.lfmm2_for_tips/file/All_TIPs_Filtered_Results.csv")

GWAS_FOLDERS <- c(
  "./in/GraffiTE_data_Fast3VmrMLM",
  "./in/GraffiTE_data_Fast3VmrMLM-Hap",
  "./in/GraffiTE_data_Fast3VmrMLM-mQTL"
)

# Global Configuration parameters
SIG_THRESHOLD <- 0.05
PRIORITY_MAP  <- c(
  "CDS" = 1, "five_prime_UTR" = 2, "three_prime_UTR" = 3, "exon" = 4, 
  "ncRNA" = 5, "mRNA" = 6, "gene" = 7, "Intergenic_0_2kb" = 8, 
  "Intergenic_2_5kb" = 9, "Intergenic_above_5kb" = 10
)

# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

# Detect Overlaps between Variants and Genomic Annotations Using IRanges
perform_overlap <- function(tips, refs) {
  common_chrs <- intersect(unique(tips$CHROM), unique(refs$chr))
  
  all_hits <- lapply(common_chrs, function(c) {
    s_sub <- tips %>% filter(CHROM == c)
    r_sub <- refs %>% filter(chr == c)
    
    if (nrow(s_sub) == 0 || nrow(r_sub) == 0) return(NULL)
    
    s_ir <- IRanges(start = as.numeric(s_sub$POS), end = as.numeric(s_sub$POS))
    r_ir <- IRanges(start = as.numeric(r_sub$start), end = as.numeric(r_sub$end))
    
    hits <- findOverlaps(s_ir, r_ir)
    if (length(hits) == 0) return(NULL)
    
    cbind(
      s_sub[queryHits(hits), ], 
      ref_type   = r_sub$type[subjectHits(hits)], 
      gene_id    = r_sub$gene_id[subjectHits(hits)],
      intergenic = r_sub$intergenic[subjectHits(hits)]
    )
  })
  bind_rows(all_hits)
}

#' Intersect Environmentally Significant Markers Against Structural GWAS Tables
find_marker_overlaps <- function(data_list, bio_list) {
  table_names <- names(data_list)
  bio_names   <- names(bio_list)
  
  results_grid <- expand.grid(Table = table_names, Bio = bio_names, stringsAsFactors = FALSE)
  
  overlap_results <- lapply(1:nrow(results_grid), function(i) {
    t_name <- results_grid$Table[i]
    b_name <- results_grid$Bio[i]
    
    current_table  <- data_list[[t_name]]
    target_markers <- bio_list[[b_name]]
    
    if (nrow(current_table) == 0 || length(target_markers) == 0) return(NULL)
    
    if (grepl("_Hap$", t_name)) {
      matched_rows <- current_table %>% 
        filter(marker_l %in% target_markers | marker_r %in% target_markers)
    } else {
      matched_rows <- current_table %>% 
        filter(marker %in% target_markers)
    }
    
    if (nrow(matched_rows) > 0) {
      return(matched_rows %>% 
               mutate(matched_bio_variable = b_name, source_gwas_table = t_name) %>% 
               dplyr::select(source_gwas_table, matched_bio_variable, everything()))
    }
    return(NULL)
  })
  
  names(overlap_results) <- paste(results_grid$Table, results_grid$Bio, sep = "__vs__")
  compact(overlap_results)
}

# ==============================================================================
# 2. GFF PARSING AND GENOMIC ENVIRONMENT CONSTRUCTION
# ==============================================================================
message(">>> Constructing GFF Genomic Intergenic Spaces...")

gff_raw <- read.table(GFF_PATH, sep = "\t", header = FALSE, comment.char = "#", quote = "", stringsAsFactors = FALSE, fill = TRUE)
colnames(gff_raw) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

gff_base <- gff_raw %>%
  filter(chr != "Chr00") %>%
  mutate(
    chr = sub("^chr", "Chr", chr),
    gene_id = case_when(
      type == "gene" ~ sub(".*ID=gene:([^;]+).*", "\\1", attributes),
      TRUE           ~ sub(".*Parent=gene:([^;]+).*", "\\1", attributes)
    )
  )

genes_only <- gff_base %>% filter(type == "gene") %>% arrange(chr, start)

# Generate space bounds matching your chromosome structures
first_intergenic <- genes_only %>%
  group_by(chr) %>%
  dplyr::slice(1) %>%
  filter(start > 1) %>%
  mutate(
    source = "manually", type = "intergenic_space",
    end = start - 1, start = 1, score = NA, strand = NA, phase = NA, attributes = NA,
    gene_id = paste("before", gene_id, sep = "_"),
    space_size = (end - start) + 1,
    intergenic = case_when(
      space_size <= 2000 ~ "intergenic_2kb",
      space_size > 2000 & space_size <= 5000 ~ "intergenic_2kb_5kb",
      space_size >= 5001 ~ "intergenic_above_5kb",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(chr, source, type, start, end, score, strand, phase, attributes, gene_id, intergenic) %>%
  ungroup()

between_intergenic <- genes_only %>%
  group_by(chr) %>%
  mutate(int_start = end + 1, int_end = lead(start) - 1, next_gene = lead(gene_id)) %>%
  filter(!is.na(int_end) & int_end >= int_start) %>% 
  mutate(
    source = "manually", type = "intergenic_space", start = int_start, end = int_end,
    score = NA, strand = NA, phase = NA, attributes = NA,
    gene_id = paste(gene_id, next_gene, sep = "_"),
    space_size = (end - start) + 1,
    intergenic = case_when(
      space_size <= 2000 ~ "intergenic_2kb",
      space_size > 2000 & space_size <= 5000 ~ "intergenic_2kb_5kb",
      space_size >= 5001 ~ "intergenic_above_5kb",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(chr, source, type, start, end, score, strand, phase, attributes, gene_id, intergenic) %>%
  ungroup()

gff_base$intergenic <- NA_character_
final_gff <- bind_rows(gff_base, first_intergenic, between_intergenic) %>% arrange(chr, start, desc(type == "gene"))

# ==============================================================================
# 3. VCF PROCESSING & MULTI-BIOCLIMATIC OVERLAPS
# ==============================================================================
message(">>> Extracting VCF Fix Metrics...")
vcf_info <- read.vcfR(VCF_PATH, verbose = FALSE)
variants_df <- as.data.frame(getFIX(vcf_info))

all_bio_results <- list()

for (file in BIO_FILES) {
  bio_name <- str_extract(basename(file), "bio[0-9]+")
  message("Processing environment matrix: ", bio_name)
  
  data_tips <- read.csv(file) %>% filter(FDR_signif == TRUE)
  if (nrow(data_tips) == 0) next
  
  matched_positions <- variants_df %>%
    filter(ID %in% data_tips$TIPs) %>%
    dplyr::select(CHROM, POS, ID) %>%
    mutate(CHROM = sub("^chr", "Chr", CHROM))
  
  data_tips_pos <- merge(matched_positions, data_tips, by.x = "ID", by.y = "TIPs")
  raw_overlaps  <- perform_overlap(data_tips_pos, final_gff)
  
  annotated_snps <- raw_overlaps %>%
    mutate(rank = PRIORITY_MAP[ref_type]) %>%
    group_by(CHROM, ID) %>%
    slice_min(rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      final_type = case_when(
        ref_type %in% c("CDS", "five_prime_UTR", "three_prime_UTR", "exon", "ncRNA", "mRNA", "gene") ~ "Intra-gene",
        intergenic == "intergenic_2kb" ~ "Intergenic < 2kb",
        intergenic == "intergenic_2kb_5kb" ~ "Intergenic 2-5kb",
        intergenic == "intergenic_above_5kb" ~ "Intergenic > 5kb",
        TRUE ~ ref_type
      ),
      bio_source = bio_name
    )
  all_bio_results[[bio_name]] <- annotated_snps
}
combined_results <- bind_rows(all_bio_results)

# ==============================================================================
# 4. FULL GENOME POSITIONING REFERENCE DATA
# ==============================================================================
message(">>> Compiling full genomic baseline tracks...")
data_tips_all_bio <- read.csv(ALL_BIO_PATH) %>% filter(FDR_signif == TRUE)

matched_positions_all_bio <- variants_df %>%
  filter(ID %in% data_tips_all_bio$TIPs) %>%
  dplyr::select(CHROM, POS, ID) %>% 
  mutate(
    CHROM = sub("^chr", "Chr", CHROM),
    CHROM = paste0("Chr", str_pad(sub("^Chr", "", CHROM), width = 2, pad = "0"))
  ) %>%
  filter(CHROM != "Chr00")

data_tips_all_bio_pos <- merge(matched_positions_all_bio, data_tips_all_bio, by.x = "ID", by.y = "TIPs")
full_bio_overlapp     <- perform_overlap(data_tips_all_bio_pos, final_gff)

annotated_tips_all_bio <- full_bio_overlapp %>%
  mutate(rank = PRIORITY_MAP[ref_type]) %>%
  group_by(CHROM, ID) %>%
  slice_min(rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    final_type = case_when(
      ref_type %in% c("CDS", "five_prime_UTR", "three_prime_UTR", "exon", "ncRNA", "mRNA", "gene") ~ "Intra-gene",
      intergenic == "intergenic_2kb" ~ "Intergenic < 2kb",
      intergenic == "intergenic_2kb_5kb" ~ "Intergenic 2-5kb",
      intergenic == "intergenic_above_5kb" ~ "Intergenic > 5kb",
      TRUE ~ ref_type
    ),
    bio_source = "Full_geno"
  )

tips_all_bio_prep <- annotated_tips_all_bio %>%
  dplyr::select(
    ID, 
    CHROM, 
    POS, 
    p_Full_bio = pvalues,           # Renaming pvalues to p_Full_bio as requested
    qvalue = qvalue,
    Bonf_signif = Bonf_signif,
    FDR_signif = FDR_signif,
    ref_type = ref_type,
    gene_id = gene_id,
    intergenic = intergenic,
    rank = rank,
    final_type = final_type,
    bio_source = bio_source
  )

merged_bio_results <- combined_results %>%
  dplyr::full_join(
    tips_all_bio_prep, 
    by = c("ID", "CHROM", "POS", "qvalue", "Bonf_signif", "FDR_signif", "ref_type", 
           "gene_id", "intergenic", "rank", "final_type","bio_source")
  )

merged_bio_results %>%
  dplyr::select(ID, CHROM, POS, p_Full_bio, contains("p_bio")) %>%
  head()

annotated_summary_all_bio <- merged_bio_results %>%
  # 1. Pivot both individual p_bio columns AND p_Full_bio into a single long column
  tidyr::pivot_longer(
    cols = c(starts_with("p_bio"), "p_Full_bio"),
    names_to = "bioclim_variable",
    values_to = "p_value"
  ) %>%
  # 2. Keep only valid p-values (removes empty relationships)
  dplyr::filter(!is.na(p_value)) %>%
  # 3. Clean up the variable names (e.g., convert "p_Full_bio" to "Full_bio")
  dplyr::mutate(
    bioclim_variable = stringr::str_remove(bioclim_variable, "^p_")
  ) %>%
  # 4. Group by your structural TIP variants
  dplyr::group_by(ID, CHROM, POS, ref_type, intergenic) %>%
  # 5. Collapse all active variables (including Full_bio) into a single row
  dplyr::summarise(
    related_biovariables = paste(sort(unique(bioclim_variable)), collapse = ", "),
    .groups = "drop"
  )

# Preview the clean summaries
head(annotated_summary_all_bio)

# ==============================================================================
# 5. STRATIFIED EXPORT AND UPSET PLOT GENERATION
# ==============================================================================

write.csv(annotated_summary_all_bio, paste0(OUT_FOLD, "/GEA_TIPs_summary.csv"), row.names=F)

p_bio_cols <- grep("p_bio", colnames(combined_results), value = TRUE)

TIPs_list <- lapply(p_bio_cols, function(col) {
  combined_results %>%
    filter(!is.na(!!sym(col)) & as.numeric(!!sym(col)) < SIG_THRESHOLD) %>%
    pull(ID)
}) %>% setNames(p_bio_cols)

TIPs_list$Full_bio <- annotated_tips_all_bio$ID

df_upset       <- as.data.frame(fromList(TIPs_list))
active_bio_cols <- names(df_upset)[grep("p_bio", names(df_upset))]
ordered_p_bio  <- names(sort(colSums(df_upset[, active_bio_cols]), decreasing = TRUE))
ordered_sets   <- c("Full_bio", ordered_p_bio)

intersections_to_show <- list("Full_bio")
for (col in ordered_p_bio) {
  intersections_to_show[[length(intersections_to_show) + 1]] <- c(col, "Full_bio")
  intersections_to_show[[length(intersections_to_show) + 1]] <- col
}

n_only     <- sum(lengths(intersections_to_show) == 1)
n_full_bio <- sum(lengths(intersections_to_show) > 1)

# Export Charts
png(file.path(OUT_FOLD, "Strict_Order_Bios_vs_Full_bio_UpSet.png"), width = 4500, height = 2000, res = 200)
upset(df_upset, sets = ordered_sets, intersections = intersections_to_show, keep.order = TRUE, 
      sets.bar.color = "gray30", main.bar.color = c(rep("steelblue", n_only), rep("darkorange", n_full_bio)),
      text.scale = c(2, 1.5, 1.5, 1.5, 1.5, 2.5), matrix.color = "black", point.size = 3, line.size = 1)
dev.off()

svg(file.path(OUT_FOLD, "Strict_Order_Bios_vs_Full_bio_UpSet.svg"), width = 20, height = 11)
upset(df_upset, sets = ordered_sets, intersections = intersections_to_show, keep.order = TRUE, 
      sets.bar.color = "gray30", main.bar.color = c(rep("steelblue", n_only), rep("darkorange", n_full_bio)),
      text.scale = c(2, 1.5, 1.5, 1.5, 1.5, 2.5), matrix.color = "black", point.size = 3, line.size = 1)
dev.off()

# ==============================================================================
# GENE IN 2kb + Gene Ontology
# ==============================================================================

intergenic_2kb <- combined_results %>%
  filter(final_type == "Intergenic < 2kb") %>%
  mutate(
    # Extracts the 2nd element after splitting by "_" row-by-row
    gene_after = str_split_i(gene_id, "_", 2)
  )

gene2test <- intergenic_2kb$gene_after

## Flowering annotation from Xilong
table_flowering <- read.table("./in/flowering_GDDH_florid.txt", sep="\t", header=T)
# Ensure the flowering gene list is a vector for easy comparison
flowering_genes <- as.character(table_flowering$GD_gene)

# Create a clean results list
comparison_results <- data.frame(gene_id = unique(gene2test), stringsAsFactors = FALSE) %>%
  mutate(
    # Check if the gene exists in the flowering array (Returns TRUE/FALSE)
    is_flowering = gene_id %in% flowering_genes
  )

# Extract a vector containing only the matched items
matched_flowering <- intersect(gene2test, flowering_genes)

cat("Total Unique Target Genes tested:", nrow(comparison_results), "\n")
cat("Number of Matching Flowering Genes:", length(matched_flowering), "\n")
### NO RESULT

## Blast against the ... database
fasta <- readDNAStringSet("./in/GDDH13_1-1_mrna.fasta")
makeblastdb(file = "./in/ath.Malus_sylvestris.blast.allgene.cds.fa", dbtype = "nucl")

bl <- blast(db = "./in/ath.Malus_sylvestris.blast.allgene.cds.fa", type = "blastn")

group_sequences <- fasta[names(fasta) %in% gene2test]

result <- NULL
  
if (length(group_sequences) > 0) {
  cat("Running BLAST for", length(group_sequences), "sequences...\n")
  
  # Perform BLAST
  blast_results <- predict(bl, group_sequences, 
                           BLAST_args = "-evalue 1e-5 -perc_identity 80 -max_target_seqs 1")
  
  if (nrow(blast_results) > 0) {
    result <- blast_results
    result$Status <- "Hit Found"
  } 
}

# Export the single master file
write.table(result, paste0(OUT_FOLD, "/TIPs_2kb_after_gene_BLAST.tsv"), row.names = FALSE, sep="\t")

cat("\nProcessing complete. Master file saved to:", OUT_FOLD)

## Gene ontology
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

current_genes <- gene2test

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
## NO RESULTS

# ==============================================================================
# 6. STRUCTURAL GWAS OVERLAP COMPARISON & EXPORT
# ==============================================================================
message(">>> Aggregating Structural GWAS Datasets...")

grouped_data_list <- lapply(GWAS_FOLDERS, function(folder_path) {
  excel_files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(excel_files) == 0) return(NULL)
  
  lapply(excel_files, function(file) {
    read_excel(file) %>% 
      clean_names() %>% 
      mutate(trait_id = as.character(trait_id)) %>% 
      mutate(source_file = basename(file)) 
  }) %>% bind_rows()
}) %>% setNames(basename(GWAS_FOLDERS) %>% gsub('-', '_', .))

grouped_data_list <- compact(grouped_data_list)
all_overlaps_list <- find_marker_overlaps(grouped_data_list, TIPs_list)
graffite_tables   <- basename(GWAS_FOLDERS) %>% gsub('-', '_', .)

grouped_by_table <- lapply(graffite_tables, function(t_name) {
  matching_elements <- all_overlaps_list[names(all_overlaps_list) %>% grepl(paste0("^", t_name, "__vs__"), .)]
  if (length(matching_elements) == 0) return(NULL)
  
  bind_rows(matching_elements) %>%
    dplyr::select(-source_gwas_table) %>%
    rename(matched_bioclim = matched_bio_variable) %>%
    dplyr::select(matched_bioclim, everything())
}) %>% 
  setNames(graffite_tables) %>% 
  compact()

# Save final clean outputs
all_chromosomes <- lapply(grouped_by_table, function(df) {
  if (nrow(df) == 0) return(NULL)
  df %>%
    dplyr::mutate(
      chrom_clean = stringr::str_replace(chromosome, "(?i)chr", ""),
      chrom_clean = paste0("Chr", sprintf("%02d", as.numeric(chrom_clean)))
    ) %>%
    dplyr::pull(chrom_clean)
}) %>% unlist() %>% unique() %>% sort()

# Create a fixed named color palette using your preferred "turbo" color spectrum
chromosome_palette <- scales::viridis_pal(option = "turbo")(length(all_chromosomes))
names(chromosome_palette) <- all_chromosomes

lapply(names(grouped_by_table), function(table_name) {
  file_path <- file.path(OUT_FOLD, paste0(table_name, "_overlapp_GEA_GWAS.csv"))
  write.csv(grouped_by_table[[table_name]], file = file_path, row.names = FALSE)
  
  df_plot <- grouped_by_table[[table_name]] %>%
      dplyr::mutate(
        chromosome = stringr::str_replace(chromosome, "(?i)chr", ""), 
        chromosome = paste0("Chr", sprintf("%02d", as.numeric(chromosome))) 
      ) %>%
      dplyr::group_by(matched_bioclim, trait_name, chromosome) %>%
      dplyr::summarise(match_count = dplyr::n(), .groups = "drop")
    
  if(nrow(df_plot) == 0) return(NULL)
  
  p <- ggplot(df_plot, aes(x = matched_bioclim, y = match_count, fill = chromosome)) +
    geom_bar(stat = "identity", position = "stack", width = 0.75) +
    facet_wrap(~ trait_name, scales = "free_y", ncol = 2) +
    
    # CRITICAL CHANGE: Use the pre-determined fixed palette for homogeneous colors
    scale_fill_manual(values = chromosome_palette) +
    
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = paste("Distribution of Genomic Matches:", stringr::str_replace_all(table_name, "_", " ")),
      subtitle = "Counts stratified by bioclimatic variable, trait context, and genomic chromosome",
      x = "Bioclimatic Variable (matched_bioclim)", y = "Number of Matches Found", fill = "Chromosome"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = "gray80"),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  calculated_height <- max(6, ceiling(length(unique(df_plot$trait_name)) / 2) * 3.5)
  file_base <- file.path(OUT_FOLD, paste0(table_name, "_matches_distribution"))
  ggplot2::ggsave(paste0(file_base, ".png"), plot = p, width = 14, height = calculated_height, dpi = 300)
  ggplot2::ggsave(paste0(file_base, ".svg"), plot = p, width = 14, height = calculated_height)
})

cat("Done! All cleaned results successfully saved inside output directory:", OUT_FOLD, "\n")


library(dplyr)
library(stringr)

# --- 1. Clean and Standardize Data Formats ---
df_gwas <- grouped_by_table$GraffiTE_data_Fast3VmrMLM %>%
  dplyr::mutate(
    chr_clean = paste0("Chr", sprintf("%02d", as.numeric(stringr::str_replace(chromosome, "(?i)chr", "")))),
    pos_num   = as.numeric(position_bp)
  )

df_mqtl <- grouped_by_table$GraffiTE_data_Fast3VmrMLM_mQTL %>%
  dplyr::mutate(
    chr_clean = paste0("Chr", sprintf("%02d", as.numeric(stringr::str_replace(chromosome, "(?i)chr", "")))),
    pos_num   = as.numeric(position_bp)
  )

df_hap <- grouped_by_table$GraffiTE_data_Fast3VmrMLM_Hap %>%
  dplyr::mutate(
    chr_clean = paste0("Chr", sprintf("%02d", as.numeric(stringr::str_replace(chromosome, "(?i)chr", "")))),
    pos_l     = as.numeric(position_l_bp),
    pos_r     = as.numeric(position_r_bp)
  )

# --- 2. Step One: Intersect GWAS and mQTL (Standard Coordinates) ---
gwas_mqtl_intersect <- dplyr::inner_join(
  df_gwas, df_mqtl, 
  by = c("chr_clean", "pos_num", "trait_name", "matched_bioclim"),
  suffix = c("_gwas", "_mqtl")
)

# --- 3. Step Two: Cross-Reference with Haplotype Interval Data ---
# Checks if the overlapping variant falls between marker_l and marker_r positions, 
# or explicitly shares a marker boundary name.
three_way_overlap <- gwas_mqtl_intersect %>%
  dplyr::inner_join(df_hap, by = c("chr_clean", "trait_name", "matched_bioclim")) %>%
  dplyr::filter(
    (pos_num >= pos_l & pos_num <= pos_r) | 
    (marker_gwas == marker_l) | 
    (marker_gwas == marker_r)
  ) %>%
  dplyr::select(
    matched_bioclim, trait_name, chr_clean, 
    position_gwas = pos_num, marker_gwas, 
    hap_start_marker = marker_l, hap_end_marker = marker_r,
    position_l_bp, position_r_bp
  ) %>%
  dplyr::distinct()

# --- 4. Display Validation Results ---
print("=== Confirmed 3-Way Overlaps (Method + Position + Trait + Bio Variable) ===")
if(nrow(three_way_overlap) > 0) {
  print(three_way_overlap)
} else {
  print("No explicit 3-way matches found in the active rows. Check for missing rows/subsets.")
}