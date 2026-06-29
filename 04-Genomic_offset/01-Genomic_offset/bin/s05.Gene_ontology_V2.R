# ==============================================================================
# 0. Global Setup & Library Imports
# ==============================================================================
rm(list = ls())

suppressPackageStartupMessages({
  # Core Data Manipulation
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(readxl)
  library(janitor)
  
  # Genomic & Annotation Packages
  library(GenomicRanges)
  library(Biostrings)
  library(IRanges)
  library(org.Mdomestica.eg.db)
  library(AnnotationDbi)
  
  # Visualization & Enrichment
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(ggpubr)
  library(ggupset)
  library(UpSetR)
  library(grid)
  library(cowplot)
  library(ggVennDiagram)
  library(VennDiagram) # Required for process_region_data/Venn pipeline
  library(clusterProfiler)
  library(GOSemSim)
  library(enrichplot)
  library(GOplot)
  
  # Sequence Alignment
  library(rBLAST)
})

# --- Project Paths ---
PROJECT <- "dom_admixed_rm"
out_fold <- paste0("./out/", PROJECT, "/s05.Gene_Ontology_V2/")
dir.create(out_fold, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. GFF Reference Preparation
# ==============================================================================
gff_path <- "./in/gene_models_20170612.gff3"
gff_raw <- read.table(
  gff_path, sep = "\t", header = FALSE, 
  comment.char = "#", quote = "", stringsAsFactors = FALSE, fill = TRUE
)
colnames(gff_raw) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Standardize GFF and extract Gene IDs
gff_base <- gff_raw %>%
  dplyr::filter(chr != "Chr00") %>%
  dplyr::mutate(
    chr = stringr::str_replace(chr, "^chr", "Chr"),
    gene_id = dplyr::case_when(
      type == "gene" ~ stringr::str_replace(attributes, ".*ID=gene:([^;]+).*", "\\1"),
      TRUE           ~ stringr::str_replace(attributes, ".*Parent=gene:([^;]+).*", "\\1")
    )
  )

genes_only <- gff_base %>% 
  dplyr::filter(type == "gene") %>% 
  dplyr::arrange(chr, start)

# Generate space bounds matching chromosome structures
first_intergenic <- genes_only %>%
  dplyr::group_by(chr) %>%
  dplyr::slice(1) %>%
  dplyr::filter(start > 1) %>%
  dplyr::mutate(
    source = "manually", type = "intergenic_space",
    end = start - 1, start = 1, score = NA, strand = NA, phase = NA, attributes = NA,
    gene_id = paste("before", gene_id, sep = "_"),
    space_size = (end - start) + 1,
    intergenic = dplyr::case_when(
      space_size <= 2000 ~ "intergenic_2kb",
      space_size > 2000 & space_size <= 5000 ~ "intergenic_2kb_5kb",
      space_size >= 5001 ~ "intergenic_above_5kb",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(chr, source, type, start, end, score, strand, phase, attributes, gene_id, intergenic) %>%
  dplyr::ungroup()

between_intergenic <- genes_only %>%
  dplyr::group_by(chr) %>%
  dplyr::mutate(int_start = end + 1, int_end = dplyr::lead(start) - 1, next_gene = dplyr::lead(gene_id)) %>%
  dplyr::filter(!is.na(int_end) & int_end >= int_start) %>% 
  dplyr::mutate(
    source = "manually", type = "intergenic_space", start = int_start, end = int_end,
    score = NA, strand = NA, phase = NA, attributes = NA,
    gene_id = paste(gene_id, next_gene, sep = "_"),
    space_size = (end - start) + 1,
    intergenic = dplyr::case_when(
      space_size <= 2000 ~ "intergenic_2kb",
      space_size > 2000 & space_size <= 5000 ~ "intergenic_2kb_5kb",
      space_size >= 5001 ~ "intergenic_above_5kb",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(chr, source, type, start, end, score, strand, phase, attributes, gene_id, intergenic) %>%
  dplyr::ungroup()

gff_base$intergenic <- NA_character_
final_gff <- dplyr::bind_rows(gff_base, first_intergenic, between_intergenic) %>% 
  dplyr::arrange(chr, start, desc(type == "gene"))

# ==============================================================================
# 2. Overlap Function and Mapping Configuration
# ==============================================================================
perform_overlap <- function(snps, refs) {
  common_chrs <- intersect(unique(snps$chromosome), unique(refs$chr))
  
  all_hits <- lapply(common_chrs, function(c) {
    s_sub <- snps %>% dplyr::filter(chromosome == c)
    r_sub <- refs %>% dplyr::filter(chr == c)
    
    if(nrow(s_sub) == 0 || nrow(r_sub) == 0) return(NULL)
    
    s_ir <- IRanges::IRanges(start = s_sub$SNPs, end = s_sub$SNPs)
    r_ir <- IRanges::IRanges(start = r_sub$start, end = r_sub$end)
    
    hits <- findOverlaps(s_ir, r_ir)
    if(length(hits) == 0) return(NULL)
    
    cbind(s_sub[queryHits(hits), ], 
          ref_type = r_sub$type[subjectHits(hits)], 
          gene_id  = r_sub$gene_id[subjectHits(hits)],
          intergenic = r_sub$intergenic[subjectHits(hits)],
          ref_start  = r_sub$start[subjectHits(hits)],
          ref_end    = r_sub$end[subjectHits(hits)]
    )
  })
  dplyr::bind_rows(all_hits)
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

# ==============================================================================
# 3. Individual Bio Files Processing & Distribution Plot
# ==============================================================================
bio_files <- list.files("./out/dom_admixed_rm/s02.lfmm2_V3/file/", pattern = "bio.*.csv", full.names = TRUE)
all_bio_results <- list()

for (file in bio_files) {
  bio_name <- stringr::str_extract(basename(file), "bio[0-9]+")
  message(">>> Processing: ", bio_name)
  
  data_snp <- read.csv(file) %>%
    dplyr::filter(FDR_signif == TRUE) %>% 
    dplyr::mutate(chromosome = stringr::str_replace(chromosome, "^chr", "Chr"))
  
  if(nrow(data_snp) == 0) {
    message("No significant SNPs for ", bio_name)
    next
  }

  raw_overlaps <- perform_overlap(data_snp, final_gff)

  annotated_snps <- raw_overlaps %>%
    dplyr::mutate(rank = priority_map[ref_type]) %>%
    dplyr::group_by(chromosome, SNPs) %>%
    dplyr::slice_min(rank, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      final_type = dplyr::case_when(
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

combined_results <- dplyr::bind_rows(all_bio_results) %>%
  dplyr::mutate(gene_id_clean = stringr::str_extract(gene_id, "MD\\d{2}G\\d+"))

# --- Generate Stacked Bar Chart ---
df_long <- combined_results %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("p_bio"),
    names_to = "p_bio_variable",
    values_to = "p_value"
  ) %>%
  dplyr::filter(!is.na(p_value))

df_counts <- df_long %>%
  dplyr::group_by(p_bio_variable, ref_type) %>%
  dplyr::summarise(snp_count = dplyr::n_distinct(SNPs), .groups = "drop")

df_totals <- df_counts %>%
  dplyr::group_by(p_bio_variable) %>%
  dplyr::summarise(total = sum(snp_count)) %>%
  dplyr::arrange(desc(total))

df_counts <- df_counts %>%
  dplyr::mutate(p_bio_variable = factor(p_bio_variable, levels = df_totals$p_bio_variable))

p_snps <- ggplot(df_counts, aes(x = p_bio_variable, y = snp_count, fill = ref_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_y_continuous(labels = scales::comma) + 
  scale_fill_viridis_d(option = "plasma") + 
  labs(
    title = "Distribution of Significant SNPs across Bioclimatic Variables",
    subtitle = "Bars sorted by total unique SNP count and stratified by feature type",
    x = "Bioclimatic Variable",
    y = "Number of Unique SNPs",
    fill = "Reference Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggplot2::ggsave(file.path(out_fold, "SNPs_count_per_p_bio_ref_type.png"), plot = p_snps, width = 13, height = 8, dpi = 300)
ggplot2::ggsave(file.path(out_fold, "SNPs_count_per_p_bio_ref_type.svg"), plot = p_snps, width = 13, height = 8)

# ==============================================================================
# 4. Genome-wide (Full) Bio Processing & UpSet Computations
# ==============================================================================
data_snp_all_bio <- read.csv("./out/dom_admixed_rm/s02.lfmm2_V3/file/All_SNPs_Filtered_Results.csv") %>%
    dplyr::filter(FDR_signif == TRUE) %>% 
    dplyr::mutate(chromosome = stringr::str_replace(chromosome, "^chr", "Chr"))

full_bio_overlapp <- perform_overlap(data_snp_all_bio, final_gff)

annotated_snps_all_bio <- full_bio_overlapp %>%
    dplyr::mutate(rank = priority_map[ref_type]) %>%
    dplyr::group_by(chromosome, SNPs) %>%
    dplyr::slice_min(rank, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      final_type = dplyr::case_when(
        ref_type %in% c("CDS", "five_prime_UTR", "three_prime_UTR", "exon", "ncRNA", "mRNA", "gene") ~ "Intra-gene",
        intergenic == "intergenic_2kb" ~ "Intergenic < 2kb",
        intergenic == "intergenic_2kb_5kb" ~ "Intergenic 2-5kb",
        intergenic == "intergenic_above_5kb" ~ "Intergenic > 5kb",
        TRUE ~ ref_type
      ),
      bio_source = "Full_geno"
    ) %>%
    dplyr::mutate(gene_id_clean = stringr::str_extract(gene_id, "MD\\d{2}G\\d+"))

snps_all_bio_prep <- annotated_snps_all_bio %>%
  dplyr::select(
    SNPs, 
    chromosome, 
    p_Full_bio = pvalues,           # Renaming pvalues to p_Full_bio as requested
    qvalue,
    Bonf_signif,
    FDR_signif,
    ref_type,
    gene_id,
    intergenic,
    ref_start,
    ref_end,
    rank,
    final_type,
    bio_source,
    gene_id_clean
  )

merged_bio_results <- combined_results %>%
  dplyr::full_join(
    snps_all_bio_prep, 
    by = c("SNPs", "chromosome", "qvalue", "Bonf_signif", "FDR_signif", "ref_type", 
           "gene_id", "intergenic", "ref_start", "ref_end", "gene_id_clean", "rank", "final_type","bio_source")
  )

merged_bio_results %>%
  dplyr::select(SNPs, chromosome, p_Full_bio, contains("p_bio")) %>%
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
  dplyr::group_by(SNPs, chromosome, ref_type, intergenic, gene_id, gene_id_clean, ref_start, ref_end) %>%
  # 5. Collapse all active variables (including Full_bio) into a single row
  dplyr::summarise(
    related_biovariables = paste(sort(unique(bioclim_variable)), collapse = ", "),
    .groups = "drop"
  )

# Preview the clean summaries
head(annotated_summary_all_bio)

# GEA_SNPS <- bind_rows(annotated_summary, annotated_summary_full_bio)
write.csv(annotated_summary_all_bio, paste0(out_fold, "/GEA_SNPs_summary.csv"), row.names=F)

# --- List UpSet Construction ---
p_bio_cols <- paste0("p_bio", sprintf("%02d", 1:19))
sig_threshold <- 0.05 

gene_list <- lapply(p_bio_cols, function(col) {
  combined_results %>%
    dplyr::filter(ref_type == "CDS") %>%
    dplyr::filter(!!rlang::sym(col) < sig_threshold) %>%
    dplyr::filter(!is.na(gene_id_clean)) %>%
    dplyr::pull(gene_id_clean) %>%
    unique()
}) %>% setNames(p_bio_cols)

SNPs_list <- lapply(p_bio_cols, function(col) {
  combined_results %>%
    dplyr::filter(!!rlang::sym(col) < sig_threshold) %>%
    dplyr::pull(SNPs) %>%
    unique()
}) %>% setNames(p_bio_cols)

gene_list$Full_bio <- unique(annotated_snps_all_bio[annotated_snps_all_bio$ref_type == "CDS", ]$gene_id_clean)
SNPs_list$Full_bio <- annotated_snps_all_bio$SNPs

df_upset <- as.data.frame(UpSetR::fromList(gene_list))
p_bio_cols_present <- names(df_upset)[grep("p_bio", names(df_upset))]

group_sizes <- colSums(df_upset[, p_bio_cols_present])
ordered_p_bio <- names(sort(group_sizes, decreasing = TRUE))
ordered_sets <- c("Full_bio", ordered_p_bio)

intersections_to_show <- list("Full_bio")
for(col in ordered_p_bio) {
  intersections_to_show[[length(intersections_to_show) + 1]] <- c(col, "Full_bio")
  intersections_to_show[[length(intersections_to_show) + 1]] <- col
}

n_only <- sum(sapply(intersections_to_show, length) == 1)
n_full_bio <- sum(sapply(intersections_to_show, length) > 1)

my_upset_grob <- grid::grid.grabExpr(
  print(
    UpSetR::upset(df_upset, 
          sets = ordered_sets,           
          intersections = intersections_to_show, 
          keep.order = TRUE,             
          sets.bar.color = "gray30",
          main.bar.color = c(rep("steelblue", n_only), rep("darkorange", n_full_bio)),
          text.scale = c(2, 1.5, 1.5, 1.5, 1.5, 2.5),
          matrix.color = "black",
          point.size = 3, 
          line.size = 1)
  )
)

ggplot2::ggsave(file.path(out_fold, "Strict_Order_Bios_vs_Full_bio_UpSet.png"), plot = my_upset_grob, width = 22.5, height = 10, dpi = 200)
ggplot2::ggsave(file.path(out_fold, "Strict_Order_Bios_vs_Full_bio_UpSet.svg"), plot = my_upset_grob, width = 20, height = 11)

# ==============================================================================
# 5. GWAS Structural Comparison Cross-checks
# ==============================================================================
message(">>> Aggregating Structural GWAS Datasets...")
GWAS_FOLDERS <- c("./in/SNP_data_Fast3VmrMLM", "./in/SNP_data_Fast3VmrMLM-Hap")

grouped_data_list <- lapply(GWAS_FOLDERS, function(folder_path) {
  excel_files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(excel_files) == 0) return(NULL)
  
  lapply(excel_files, function(file) {
    readxl::read_excel(file) %>% 
      janitor::clean_names() %>% 
      dplyr::mutate(trait_id = as.character(trait_id), source_file = basename(file)) 
  }) %>% dplyr::bind_rows()
}) %>% setNames(basename(GWAS_FOLDERS) %>% stringr::str_replace_all('-', '_'))

grouped_data_list <- purrr::compact(grouped_data_list)

filtered_data_list <- lapply(names(grouped_data_list), function(table_name) {
  df <- grouped_data_list[[table_name]]
  if (nrow(df) == 0) return(df)
  
  condition <- stringr::str_detect(df$trait_name, "(?i)^(AGR_Seasonal|RGR_Seasonal)$")
  if (stringr::str_detect(table_name, "_Hap$")) {
    condition <- condition | stringr::str_detect(df$trait_name, "(?i)^RGR[ _]annual[ _]ww$")
  }
  df %>% dplyr::filter(!condition)
}) %>% setNames(names(grouped_data_list)) %>% purrr::compact()

# --- Marker Overlap Function ---
find_marker_overlaps <- function(data_list, bio_list) {
  results_grid <- expand.grid(Table = names(data_list), Bio = names(bio_list), stringsAsFactors = FALSE)
  
  overlap_results <- lapply(1:nrow(results_grid), function(i) {
    t_name <- results_grid$Table[i]
    b_name <- results_grid$Bio[i]
    current_table  <- data_list[[t_name]]
    target_markers <- bio_list[[b_name]]
    
    if (nrow(current_table) == 0 || length(target_markers) == 0) return(NULL)
    
    if (stringr::str_detect(t_name, "_Hap$")) {
      matched_rows <- current_table %>% 
        dplyr::mutate(pos_l_num = as.numeric(position_l_bp), pos_r_num = as.numeric(position_r_bp)) %>%
        dplyr::filter(rowSums(outer(pos_l_num, target_markers, "<=") & outer(pos_r_num, target_markers, ">=")) > 0) %>%
        dplyr::select(-pos_l_num, -pos_r_num)
    } else {
      matched_rows <- current_table %>% 
        dplyr::filter(as.numeric(position_bp) %in% target_markers)
    }
    
    if (nrow(matched_rows) > 0) {
      return(matched_rows %>% 
               dplyr::mutate(matched_bio_variable = b_name, source_gwas_table = t_name) %>% 
               dplyr::select(source_gwas_table, matched_bio_variable, dplyr::everything()))
    }
    return(NULL)
  })
  names(overlap_results) <- paste(results_grid$Table, results_grid$Bio, sep = "__vs__")
  purrr::compact(overlap_results)
}

all_overlaps_list <- find_marker_overlaps(filtered_data_list, SNPs_list)
SNP_tables <- basename(GWAS_FOLDERS) %>% stringr::str_replace_all('-', '_')

grouped_by_table <- lapply(SNP_tables, function(t_name) {
  matching_elements <- all_overlaps_list[stringr::str_detect(names(all_overlaps_list), paste0("^", t_name, "__vs__"))]
  if (length(matching_elements) == 0) return(NULL)
  
  dplyr::bind_rows(matching_elements) %>%
    dplyr::select(-source_gwas_table) %>%
    dplyr::rename(matched_bioclim = matched_bio_variable) %>%
    dplyr::select(matched_bioclim, dplyr::everything())
}) %>% setNames(SNP_tables) %>% purrr::compact()

# Export Tables & Faceted Bar Charts
lapply(names(grouped_by_table), function(table_name) {
  write.csv(grouped_by_table[[table_name]], file = file.path(out_fold, paste0(table_name, "_overlapp_GEA_GWAS.csv")), row.names = FALSE)
  
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
    scale_fill_viridis_d(option = "turbo") +
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
  file_base <- file.path(out_fold, paste0(table_name, "_matches_distribution"))
  ggplot2::ggsave(paste0(file_base, ".png"), plot = p, width = 14, height = calculated_height, dpi = 300)
  ggplot2::ggsave(paste0(file_base, ".svg"), plot = p, width = 14, height = calculated_height)
})

# ==============================================================================
# 6. Venn Subsetting and Functional Annotation Integration
# ==============================================================================
table_flowering <- read.table("./in/flowering_GDDH_florid.txt", sep="\t", header=TRUE)

selected_sets <- list(
  bio19 = gene_list$p_bio19, bio08 = gene_list$p_bio08,
  bio12 = gene_list$p_bio12, Full_geno = gene_list$Full_bio
)

p_venn <- ggVennDiagram::ggVennDiagram(selected_sets) +
  ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Venn Diagram: bio19, 08, 12, and Full_geno", caption = "Data source: gene_list") +
  ggplot2::theme_void() + 
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(legend.position = "none", plot.background = element_rect(fill = "white", color = NA))

ggplot2::ggsave(file.path(out_fold, "venn_diagram_bio19_8_12_Full_geno.png"), plot = p_venn, width = 8, height = 6, dpi = 300)
ggplot2::ggsave(file.path(out_fold, "venn_diagram_bio19_8_12_Full_geno.svg"), plot = p_venn, width = 8, height = 6)

# Extracting Venn Intersection Metrics
venn_regions <- VennDiagram::process_region_data(VennDiagram::Venn(selected_sets))
flowering_genes <- as.character(table_flowering$GD_gene)

comparison_results <- venn_regions %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    total_genes = length(item),
    matching_flowering = list(intersect(item, flowering_genes)),
    count_flowering = length(matching_flowering)
  ) %>%
  dplyr::ungroup()

all_genes_with_groups <- comparison_results %>%
  dplyr::select(name, item) %>%
  tidyr::unnest(item) %>%
  dplyr::rename(GD_gene = item, Venn_Group = name)

final_gene_summary <- all_genes_with_groups %>%
  dplyr::inner_join(table_flowering, by = "GD_gene") %>%
  dplyr::select(GD_gene, Venn_Group, flowering_gene, TAIR)

blast_gene_summary_clean <- final_gene_summary %>%
  dplyr::filter(Venn_Group %in% c("Full_geno", "bio19/Full_geno", "bio08/Full_geno", "bio12/Full_geno", 
                                  "bio19/bio08/Full_geno", "bio19/bio12/Full_geno", "bio08/bio12/Full_geno", 
                                  "bio19/bio08/bio12/Full_geno")) %>%
  tidyr::separate_rows(flowering_gene, sep = ", ")

# --- Categorization Mapping ---
functional_categories <- tibble::tribble(
  ~Category, ~Genes,
  "Photoperiod & light perception", "BBX19, CIB1, CIB4, CKB4, COL5, COL9, COP1, CRY1, CRY2, DET1, FUS2, GI, LATE, LRB1, LRB2, POB1, PFT1, MED25, PHYA, PHYB, PHYC, PHYE, RFI2, RUP2, EFO2, SPA1, SPA3, SPA4, SRR1, STO, TIC, PIF4, PIF5, PIL6, PHL",
  "Circadian clock", "CDF2, CDF3, CDF4, CDF5, ELF3, ELF4, ELF6, ELF7, ELF8, VIP6, FKF1, ADO3, LHY, LNK2, LUX, PCL1, LWD1, ATAN11, LWD2, PRR5, PRR7, TOC1, ZTL",
  "Autonomous pathway", "FCA, FLD, RSI1, FLK, FPA, FVE, MSI4, FY, LD",
  "Vernalization pathway", "VIL1, VRN5, VIL2, VEL1, VRN1, VRN2",
  "Gibberellin (GA) pathway", "FPF1, GA1, GA2, ATKS1, GA2ox1, GA2ox2, GA2ox6, GA2ox8, GA3OX1, GA3OX2, GA20ox1, GA20ox2, GA20OX3, GAI, RGA2, GID1B, GID1C, RGA1, RGL1, RGL2",
  "Aging pathway", "SPL3, SPL4, SPL5, SPL9, SPL15, TOE1, RAP2.7",
  "Floral integrators & meristem identity", "AGL15, AGL16, AGL18, AGL19, AGL24, AGL6, AP1, AP2, ATC, BRC1, TCP18, FD, flowering_gene, FT, FUL, AGL8, LFY, SOC1, AGL20, SVP, TFL1",
  "Epigenetic, chromatin & histones", "AFR1, AFR2, ARP4, ARP6, SUF3, ESD1, ASH2R, TRO, ATPRMT10, ATX1, SDG27, ATX2, ATXR7, SDG25, BAF60, CHC1, BRM, CHR2, CLF, SDG1, EBS, EFS, SDG8, EMF1, EMF2, CYR1, ESD7, TIL1, ABO4, FIE1, FIS3, HDA5, HDA6, RTS1, HDA9, HTA9, HTA11, HUA2, HUB1, HUB2, INO80, JMJ14, PKDM7B, JMJ4, JMJ15, MEE27, JMJ30, JMJ32, JMJD5, LDL1, LDL2, MBD9, MRG1, MRG2, MSI1, PIE1, PRMT4A, REF6, JMJ12, SDG7, ASSH3, SDG26, ASHH1, SEF, SWC6, SKB1, PRMT5, SPT16, SSRP1, HMG, NFD, SUF4, SUVR5, TFL2, LHP1, UBC1, UBC2, ULT1, VIM1, ORTH2, VIP3, SKI8, VIP4, VIP5, WDR5A, YAF9A, TAF14B, GAS41, LIF2",
  "RNA processing & export", "ABH1, CBP80, AtNDX, CBP20, CSTF64, CSTF77, DCL1, DCL3, DCL4, ESD4, FIO1, FLX, HLP1, MOS1, PHP, CDC73, PRP8, SUS2, EMB33, RRP6L2, SUS4",
  "Mediator complex & transcription", "MED12, CCT, CRP, MED13, GCT, MAB2, MED16, SFR6, MED18, NF-YA1, HAP2B, NF-YB1, HAP3A, NF-YB3, HAP3C, NF-YC3, NF-YC4, NF-YC9, HAP5C, TPL",
  "Regulation", "AHL22, AS1, ATH1, BFT, FBH1, FBH3, AKS1, BHLH122, FBH4, AKS3, GIS5, EMB2780, HB16, ATHB16, HULK2, ICE1, SCRM, MYB30, MYR2, SNZ, TEM1, EDF1, TEM2, RAV2, VOZ1",
  "Protein degradation (ubiquitination)", "CUL3A, CUL4, ESD6, HOS1, OTS1, ULP1D, PUB13, SIZ1, UBP12, UBP13, UBP26",
  "Transport & mobile signals", "APL, WDY, FE, FTIP1",
  "Metabolism & stress signaling", "ADG1, APS1, AKIN10, SNRK1.1, CAD2, RML1, GSH1, CPK6, EFM, FES1, FRL1, SUF8, GLK1, GPRI1, HAM1, PGI1, PGM1, RCD1, SRO1, TPS1, UGT87A2"
)

category_lookup <- functional_categories %>%
  tidyr::separate_rows(Genes, sep = "; |, | ") %>%
  dplyr::mutate(Genes = trimws(Genes)) %>%
  dplyr::distinct()

final_summary <- blast_gene_summary_clean %>%
  dplyr::left_join(category_lookup, by = c("flowering_gene" = "Genes")) %>%
  dplyr::select(GD_gene, Venn_Group, flowering_gene, Category, TAIR) %>%
  dplyr::filter(Venn_Group != "Full_geno")

write.table(final_summary, file.path(out_fold, "Xilong_flowering_annotation_final_summary.txt"), sep = "\t", row.names = FALSE)

# Functional Composition Bar Chart
group_totals <- final_summary %>% dplyr::group_by(Venn_Group) %>% dplyr::summarise(total = dplyr::n())

p_func <- ggplot(final_summary, aes(y = Venn_Group, fill = Category)) +
  geom_bar(stat = "count", position = "stack", alpha = 0.9) +
  geom_text(data = group_totals, aes(x = total, y = Venn_Group, label = paste0("n = ", total)),
            inherit.aes = FALSE, hjust = -0.1, size = 2, fontface = "bold") +
  theme_minimal() +
  scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = "Functional Composition & Gene Counts per Venn Group", x = "Number of Genes", y = "Venn Group", fill = "Functional Category") +
  theme(legend.position = "right", legend.text = element_text(size = 8), panel.grid.minor = element_blank())

ggplot2::ggsave(file.path(out_fold, "Functional_category_distribution.png"), p_func, width = 10, height = 5, dpi = 300)
ggplot2::ggsave(file.path(out_fold, "Functional_category_distribution.svg"), p_func, width = 10, height = 5)

# ==============================================================================
# 7. Functional Enrichment Analysis Profiling (GO ClusterProfiler)
# ==============================================================================
fasta <- Biostrings::readDNAStringSet("./in/GDDH13_1-1_mrna.fasta")
universe_genes <- unique(stringr::str_replace(names(fasta), ".*gene[:=]([^;| ]+).*", "\\1"))

semdata <- GOSemSim::godata(OrgDb = org.Mdomestica.eg.db, keytype = "GID", ont = "BP")

dotplot_list <- list()
emapplot_list <- list()
treeplot_list <- list()
out_cat <- "_Blast"

for (group in unique(final_gene_summary$Venn_Group)) {
  current_genes <- final_gene_summary %>% dplyr::filter(Venn_Group == group) %>% dplyr::pull(GD_gene)
  message(paste("--- Processing Intersection:", group, "| Genes:", length(current_genes), "---"))
  
  go_bp <- clusterProfiler::enrichGO(
    gene = current_genes, universe = universe_genes, OrgDb = org.Mdomestica.eg.db, 
    keyType = "GID", ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
    readable = TRUE, pAdjustMethod = "fdr"
  )
  
  if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
    plot_title <- paste(group, "(n =", length(current_genes), ")")
    dotplot_list[[group]] <- enrichplot::dotplot(go_bp, showCategory = 10) + labs(title = plot_title) + theme(legend.position = "none")
    
    if (nrow(as.data.frame(go_bp)) > 2) {
      go_bp_sim <- enrichplot::pairwise_termsim(go_bp)
      emapplot_list[[group]] <- enrichplot::emapplot(go_bp_sim, showCategory = 20) + labs(title = plot_title) + theme(legend.position = "none")
      
      n_terms <- nrow(as.data.frame(go_bp_sim))
      if (n_terms >= 2) {
        n_clusters <- min(5, n_terms - 1)
        treeplot_list[[group]] <- suppressWarnings(enrichplot::treeplot(go_bp_sim, cluster.params = list(n = max(1, n_clusters)))) + 
          labs(title = plot_title) + theme(legend.position = "none")
      }
    }
  }
}

# --- Grid Assembly & Multiplot Panel Saving ---
if (length(dotplot_list) > 0) {
  all_counts <- unlist(lapply(dotplot_list, function(p) p$data$Count))
  max_c <- if(length(all_counts) > 0) max(all_counts) else 50
  master_legend_plot <- dotplot_list[[1]] + theme(legend.position = "left") + 
    scale_color_continuous(low = "red", high = "blue", limits = c(0, 0.05), name = "p.adjust") +
    scale_size_continuous(limits = c(1, max_c), name = "Count")
  shared_legend <- cowplot::get_legend(master_legend_plot)
  combined_grid <- patchwork::wrap_plots(dotplot_list, ncol = 3) + patchwork::plot_annotation(title = "GO BP Dotplots: UpSet Intersections")
  final_dot <- patchwork::wrap_plots(shared_legend, combined_grid, widths = c(1, 15))
  ggplot2::ggsave(file.path(out_fold, paste0("All_Intersections_Dotplots", out_cat, ".png")), final_dot, width = 24, height = 5 * ceiling(length(dotplot_list) / 3), dpi = 300, limitsize = FALSE)
}

if (length(emapplot_list) > 0) {
  master_emap_legend <- emapplot_list[[1]] + theme(legend.position = "left") + scale_color_continuous(low = "red", high = "blue", limits = c(0, 0.05), name = "p.adjust")
  shared_emap_legend <- cowplot::get_legend(master_emap_legend)
  combined_emap_grid <- patchwork::wrap_plots(emapplot_list, ncol = 3) + patchwork::plot_annotation(title = "Enrichment Maps: UpSet Intersections")
  final_emap <- patchwork::wrap_plots(shared_emap_legend, combined_emap_grid, widths = c(1, 15))
  ggplot2::ggsave(file.path(out_fold, paste0("All_Intersections_Emaps", out_cat, ".png")), final_emap, width = 24, height = 7 * ceiling(length(emapplot_list) / 3), dpi = 300, limitsize = FALSE)
}

if (length(treeplot_list) > 0) {
  master_tree_legend <- treeplot_list[[1]] + theme(legend.position = "left") + scale_color_continuous(low = "red", high = "blue", limits = c(0, 0.05), name = "p.adjust")
  shared_tree_legend <- cowplot::get_legend(master_tree_legend)
  combined_tree_grid <- patchwork::wrap_plots(treeplot_list, ncol = 2) + patchwork::plot_annotation(title = "Treeplots: UpSet Intersections")
  final_tree <- patchwork::wrap_plots(shared_tree_legend, combined_tree_grid, widths = c(1, 12))
  ggplot2::ggsave(file.path(out_fold, paste0("All_Intersections_Trees", out_cat, ".png")), final_tree, width = 30, height = 8 * ceiling(length(treeplot_list) / 2), dpi = 300, limitsize = FALSE)
}

# ==============================================================================
# 8. Local Blast Target Run Validation
# ==============================================================================
rBLAST::makeblastdb(file = "./in/ath.Malus_sylvestris.blast.allgene.cds.fa", dbtype = "nucl")
bl <- rBLAST::blast(db = "./in/ath.Malus_sylvestris.blast.allgene.cds.fa", type = "blastn")
all_blast_list <- list()

for (group in venn_regions$name) {
  cat("\nProcessing Group:", group, "\n")
  current_gene_ids <- final_gene_summary %>% dplyr::filter(Venn_Group == group) %>% dplyr::pull(GD_gene)
  group_sequences <- fasta[names(fasta) %in% current_gene_ids]
  current_result <- NULL
  
  if (length(group_sequences) > 0) {
    blast_results <- stats::predict(bl, group_sequences, BLAST_args = "-evalue 1e-5 -perc_identity 80 -max_target_seqs 1")
    if (nrow(blast_results) > 0) {
      current_result <- blast_results
      current_result$Venn_Group <- group
      current_result$Status <- "Hit Found"
    } 
  }
  
  if (is.null(current_result)) {
    current_result <- data.frame(Venn_Group = group, Status = "No Result / No sequences found", stringsAsFactors = FALSE)
  }
  all_blast_list[[group]] <- current_result
}

final_combined_results <- dplyr::bind_rows(all_blast_list)
write.table(final_combined_results, paste0(out_fold, "All_Groups_BLAST_Combined.tsv"), row.names = FALSE, sep = "\t")
cat("\nProcessing complete. Master file saved to:", out_fold)