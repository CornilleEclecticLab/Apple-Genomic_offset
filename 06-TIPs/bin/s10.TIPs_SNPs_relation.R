rm(list=ls())

library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

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

## GEA
GEA_SNPS <- read.csv("/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/04-Genomic_offset/01-Genomic_offset/out/dom_admixed_rm/s05.Gene_Ontology_V2/GEA_SNPs_summary.csv")
GEA_TIPS <- read.csv("./out/tips_dom_admixed_rm_V3/s09.tips_annot_gene_ontology/GEA_TIPs_summary.csv")

# ==============================================================================
# 0. Data Preparation & Standardization
# ==============================================================================
# Standardizing Chromosome names to 'Chr01', 'Chr02', etc.
clean_chrom <- function(chr_col) {
  nums <- stringr::str_extract(chr_col, "\\d+")
  paste0("Chr", sprintf("%02d", as.numeric(nums)))
}

snps_clean <- GEA_SNPS %>%
  dplyr::mutate(
    CHROM = clean_chrom(chromosome),
    POS   = as.numeric(SNPs)
  )

tips_clean <- GEA_TIPS %>%
  dplyr::mutate(
    CHROM = clean_chrom(CHROM),
    POS   = as.numeric(POS)
  )
# ==============================================================================
# 1. Standardize and Pre-process Input Datasets
# ==============================================================================
# Ensure genomic positions are numeric to prevent character mathematical bugs
tips_numeric <- tips_clean %>%
  dplyr::mutate(POS_TIP = as.numeric(POS))

snps_numeric <- snps_clean %>%
  dplyr::mutate(
    POS_SNP   = as.numeric(POS),
    ref_start = as.numeric(ref_start),
    ref_end   = as.numeric(ref_end)
  )

# ==============================================================================
# 2. Compute Spatial & Contextual Relationships
# ==============================================================================
WINDOW_SIZE <- 10000

spatial_matches <- tips_numeric %>%
  # Merge datasets on matching chromosomes
  dplyr::inner_join(
    snps_numeric %>% 
      dplyr::select(
        CHROM, 
        POS_SNP, 
        ref_type_SNP = ref_type, 
        gene_id, 
        gene_id_clean, 
        ref_start, 
        ref_end, 
        related_biovariables_SNPs = related_biovariables
      ) %>% 
      dplyr::distinct(),
    by = "CHROM",
    relationship = "many-to-many"
  ) %>%
  # Flag proximal constraints
  dplyr::mutate(
    Inside_Gene = (POS_TIP >= ref_start & POS_TIP <= ref_end),
    Near_Gene   = (POS_TIP >= (ref_start - WINDOW_SIZE) & POS_TIP <= (ref_end + WINDOW_SIZE))
  ) %>%
  dplyr::filter(Near_Gene == TRUE) %>%
  # Formulate physical tracking and orientation details
  dplyr::mutate(
    Distance_To_Start = POS_TIP - ref_start,
    Distance_To_End   = POS_TIP - ref_end,
    
    Orientation = dplyr::case_when(
      Inside_Gene           ~ "Inside_Gene",
      Distance_To_Start < 0 ~ "Upstream",
      TRUE                  ~ "Downstream"
    ),
    
    Flank_Dist_bp = dplyr::case_when(
      Inside_Gene ~ 0,
      Orientation == "Upstream"   ~ abs(Distance_To_Start),
      Orientation == "Downstream" ~ abs(Distance_To_End)
    ),
    
    Distance_Class = dplyr::case_when(
      Inside_Gene ~ "Inside_Gene",
      Flank_Dist_bp < 2000 ~ "less_than_2kb",
      Flank_Dist_bp >= 2000 & Flank_Dist_bp <= 5000 ~ "between_2_5kb",
      Flank_Dist_bp > 5000 ~ "above_5kb"
    ),
    
    Spatial_Context = dplyr::case_when(
      Inside_Gene ~ "Inside Reference Window",
      TRUE        ~ paste0(Distance_Class, "_", Orientation)
    )
  )

# ==============================================================================
# 3. Fine-Grained Multidimensional Filtering
# ==============================================================================
filtered_matches <- spatial_matches %>%
  dplyr::filter(
    # A. Target TIP Context: Must be intergenic space flanking within 2kb of the gene
    ref_type == "intergenic_space",
    Spatial_Context %in% c("less_than_2kb_Upstream", "less_than_2kb_Downstream"),
    
    # B. Target SNP Context: Single Nucleotide variant must reside inside the gene transcript layers
    ref_type_SNP %in% c("CDS", "exon", "five_prime_UTR", "three_prime_UTR", "mRNA", "ncRNA"),
    
    # C. Environmental Co-adaptation: Vectorized verification of shared biovariables
    purrr::map2_lgl(related_biovariables, related_biovariables_SNPs, function(tip, snp) {
      tip_split <- stringr::str_split(tip, ",\\s*")[[1]]
      snp_split <- stringr::str_split(snp, ",\\s*")[[1]]
      length(intersect(tip_split, snp_split)) > 0
    })
  )

# Inspect the standardized results
head(filtered_matches)
OUT_FOLD <- "./out/tips_dom_admixed_rm_V3/s10.TIPs_SNPs_relation/"
dir.create(OUT_FOLD, recursive=T)
write.csv(filtered_matches, paste0(OUT_FOLD, "GEA_SNPs_TIPs_relation.csv"), row.names=F)

## and what now ??
## Caractérisation fine des Éléments Transposables ?
## Analyses d'enrichissement fonctionnel
## Corrélation avec les Variables Bioclimatiques -- répartition

## Caractérisation des TIPs en amont ou en aval
upset_matrix_tips <- filtered_matches %>%
  # Select only the unique combinations of TIP ID and its biovariables
  dplyr::select(ID, related_biovariables) %>%
  dplyr::distinct() %>%
  # Split variables by comma (keeping "Full_bio" completely intact as a text block)
  tidyr::separate_rows(related_biovariables, sep = ",\\s*") %>%
  # Mark presence and drop duplicates
  dplyr::mutate(Present = 1) %>%
  dplyr::distinct() %>%
  # Reshape into wide format (0 or 1 matrix)
  tidyr::pivot_wider(names_from = related_biovariables, values_from = Present, values_fill = 0) %>%
  as.data.frame()


upset_matrix_snps <- filtered_matches %>%
  # Select only the unique combinations of TIP ID and its biovariables
  dplyr::select(POS_SNP, related_biovariables_SNPs) %>%
  dplyr::distinct() %>%
  # Split variables by comma (keeping "Full_bio" completely intact as a text block)
  tidyr::separate_rows(related_biovariables_SNPs, sep = ",\\s*") %>%
  # Mark presence and drop duplicates
  dplyr::mutate(Present = 1) %>%
  dplyr::distinct() %>%
  # Reshape into wide format (0 or 1 matrix)
  tidyr::pivot_wider(names_from = related_biovariables_SNPs, values_from = Present, values_fill = 0) %>%
  as.data.frame()

upset_matrix_gene <- filtered_matches %>%
  # Select only the unique combinations of TIP ID and its biovariables
  dplyr::select(gene_id_clean, related_biovariables_SNPs) %>%
  dplyr::distinct() %>%
  # Split variables by comma (keeping "Full_bio" completely intact as a text block)
  tidyr::separate_rows(related_biovariables_SNPs, sep = ",\\s*") %>%
  # Mark presence and drop duplicates
  dplyr::mutate(Present = 1) %>%
  dplyr::distinct() %>%
  # Reshape into wide format (0 or 1 matrix)
  tidyr::pivot_wider(names_from = related_biovariables_SNPs, values_from = Present, values_fill = 0) %>%
  as.data.frame()

# Isolate columns to plot (every column except the TIP ID)
bioclim_columns_tips <- setdiff(colnames(upset_matrix_tips), "ID")
# Isolate columns to plot (every column except the POS_SNP)
bioclim_columns_snps <- setdiff(colnames(upset_matrix_snps), "POS_SNP")
# Isolate columns to plot (every column except the POS_SNP)
bioclim_columns_gene <- setdiff(colnames(upset_matrix_gene), "gene_id_clean")

# ==============================================================================
# 2. Generate and Save the UpSet Plot
# ==============================================================================
png(filename = paste0(OUT_FOLD,"TIP_ID_Biovariable_UpSet.png"), width = 2800, height = 2000, res = 300)

UpSetR::upset(
  upset_matrix_tips, 
  sets = bioclim_columns_tips,
  order.by = "freq", 
  decreasing = TRUE,
  main.bar.color = "#2c3e50", 
  sets.bar.color = "#e74c3c", 
  text.scale = 1.3
)

dev.off()

png(filename = paste0(OUT_FOLD,"SNP_ID_Biovariable_UpSet.png"), width = 2800, height = 2000, res = 300)

UpSetR::upset(
  upset_matrix_snps, 
  sets = bioclim_columns_snps,
  order.by = "freq", 
  decreasing = TRUE,
  main.bar.color = "#2c3e50", 
  sets.bar.color = "#e74c3c", 
  text.scale = 1.3
)

dev.off()

png(filename = paste0(OUT_FOLD,"Gene_Biovariable_UpSet.png"), width = 2800, height = 2000, res = 300)

UpSetR::upset(
  upset_matrix_gene, 
  sets = bioclim_columns_gene,
  order.by = "freq", 
  decreasing = TRUE,
  main.bar.color = "#2c3e50", 
  sets.bar.color = "#e74c3c", 
  text.scale = 1.3
)

dev.off()

library(UpSetR)
library(cowplot)
library(ggplotify)

# ==============================================================================
# 1. Generate and Convert UpSet Plots to ggplot Objects
# ==============================================================================
g_tips <- ggplotify::as.ggplot(
  UpSetR::upset(upset_matrix_tips, sets = bioclim_columns_tips, order.by = "freq", 
                decreasing = TRUE, main.bar.color = "#2c3e50", sets.bar.color = "#e74c3c", text.scale = 2.2)
)

g_snps <- ggplotify::as.ggplot(
  UpSetR::upset(upset_matrix_snps, sets = bioclim_columns_snps, order.by = "freq", 
                decreasing = TRUE, main.bar.color = "#2c3e50", sets.bar.color = "#e74c3c", text.scale = 2.2)
)

g_gene <- ggplotify::as.ggplot(
  UpSetR::upset(upset_matrix_gene, sets = bioclim_columns_gene, order.by = "freq", 
                decreasing = TRUE, main.bar.color = "#2c3e50", sets.bar.color = "#e74c3c", text.scale = 2.2)
)

# ==============================================================================
# 2. Combine and Save Side-by-Side
# ==============================================================================
png(filename = paste0(OUT_FOLD, "Combined_Biovariable_UpSet_Panels.png"), 
    width = 10000, height = 4000, res = 300)

cowplot::plot_grid(
  g_tips, g_snps, g_gene, 
  ncol = 3, 
  labels = c("A) TIPs", "B) SNPs", "C) Genes"),
  label_size = 18,
  label_fontface = "bold"
)

dev.off()

## GWAS SNPs and TIPs, do we found the gwas results in our list of TIPs/SNPs in potential gene relationship

# Define your input file paths
TIPs_GEA_GWAS_FILES <- c(
  './out/tips_dom_admixed_rm_V3/s09.tips_annot_gene_ontology/GraffiTE_data_Fast3VmrMLM_Hap_overlapp_GEA_GWAS.csv',
  './out/tips_dom_admixed_rm_V3/s09.tips_annot_gene_ontology/GraffiTE_data_Fast3VmrMLM_mQTL_overlapp_GEA_GWAS.csv',
  './out/tips_dom_admixed_rm_V3/s09.tips_annot_gene_ontology/GraffiTE_data_Fast3VmrMLM_overlapp_GEA_GWAS.csv'
)

SNPs_GEA_GWAS_FILES <- c(
  '../04-Genomic_offset/01-Genomic_offset/out/dom_admixed_rm/s05.Gene_Ontology_V2/SNP_data_Fast3VmrMLM_Hap_overlapp_GEA_GWAS.csv',
  '../04-Genomic_offset/01-Genomic_offset/out/dom_admixed_rm/s05.Gene_Ontology_V2/SNP_data_Fast3VmrMLM_overlapp_GEA_GWAS.csv'
)

# ==============================================================================
# Master Core Parsing Function
# ==============================================================================
extract_gea_gwas_metrics <- function(file_path) {
  # Guard check: Ensure file physically exists before reading
  if (!file.exists(file_path)) {
    warning(paste("File missing, skipping target scope:", file_path))
    return(NULL)
  }
  
  # Read incoming CSV safely forcing positions to character formats
  df <- readr::read_csv(
    file_path, 
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )
  
  if (nrow(df) == 0) return(NULL)
  
  # Determine if this file uses the Haplotype (Hap) layout architecture
  is_hap_file <- "marker_l" %in% colnames(df)
  
  if (is_hap_file) {
    # --------------------------------------------------------------------------
    # Case A: Process Haplotype blocks (Split Left & Right entries into long rows)
    # --------------------------------------------------------------------------
    # Extract Left Markers
    left_entries <- df %>%
      dplyr::select(
        biovariable  = matched_bioclim,
        trait_id,
        trait_name,
        marker       = marker_l,
        chromosome,
        position_bp  = position_l_bp,
        source_file
      ) %>%
      dplyr::mutate(hap_side = "Left")
    
    # Extract Right Markers
    right_entries <- df %>%
      dplyr::select(
        biovariable  = matched_bioclim,
        trait_id,
        trait_name,
        marker       = marker_r,
        chromosome,
        position_bp  = position_r_bp,
        source_file
      ) %>%
      dplyr::mutate(hap_side = "Right")
    
    # Combine Left and Right blocks together
    parsed_df <- dplyr::bind_rows(left_entries, right_entries)
    
  } else {
    # --------------------------------------------------------------------------
    # Case B: Process Standard Files (Single explicit Marker layouts)
    # --------------------------------------------------------------------------
    parsed_df <- df %>%
      dplyr::select(
        biovariable = matched_bioclim,
        trait_id,
        trait_name,
        marker,
        chromosome,
        position_bp,
        source_file
      ) %>%
      dplyr::mutate(hap_side = "Single_Marker")
  }
  
  # Track origin profile type metadata
  parsed_df <- parsed_df %>%
    dplyr::mutate(
      input_file_origin = basename(file_path),
      is_hap_format     = is_hap_file
    )
  
  return(parsed_df)
}

# ==============================================================================
# Execute Extractions Across All Input Targets
# ==============================================================================
# Parse all TIP variants
all_extracted_tips <- TIPs_GEA_GWAS_FILES %>%
  purrr::map(extract_gea_gwas_metrics) %>%
  purrr::keep(~ !is.null(.x)) %>%
  dplyr::bind_rows()

# Parse all SNP variants
all_extracted_snps <- SNPs_GEA_GWAS_FILES %>%
  purrr::map(extract_gea_gwas_metrics) %>%
  purrr::keep(~ !is.null(.x)) %>%
  dplyr::bind_rows()

# Clean up final variables (removing redundant row duplicates if left/right are identical)
all_extracted_tips <- dplyr::distinct(all_extracted_tips)
all_extracted_snps <- dplyr::distinct(all_extracted_snps)

# ==============================================================================
# Audit Output Previews
# ==============================================================================
print("=== Previewing Extracted TIPs Metadata ===")
print(all_extracted_tips)

print("=== Previewing Extracted SNPs Metadata ===")
print(all_extracted_snps)

# ==============================================================================
# Overlapp TIPs SNPs relation
# ==============================================================================

intersect(filtered_matches$ID, all_extracted_tips$marker) ## n=2
## [1] "MdoR9A.svim_asm.DEL.14570" "MdoRA.svim_asm.INS.2031" 
intersect(filtered_matches$POS_SNP,as.numeric(all_extracted_snps$position_bp)) ## n=0

### GO & Blast for gene

table_flowering <- read.table("./in/flowering_GDDH_florid.txt", sep="\t", header=TRUE)
flowering_genes <- as.character(table_flowering$GD_gene)

matching_gene <- intersect(flowering_genes, unique(filtered_matches$gene_id_clean)) #n=4
# [1] "MD02G1261100" "MD07G1180900" "MD08G1197300" "MD12G1197200"

flowering_annot <- table_flowering %>%
  filter(GD_gene %in% matching_gene)

filtered_matches_annot <- filtered_matches %>%
filter(gene_id_clean %in% matching_gene)

merge_annot_filtered <- merge(filtered_matches_annot, flowering_annot, by.x='gene_id_clean', by.y='GD_gene')

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

final_summary <- merge_annot_filtered %>%
  dplyr::left_join(category_lookup, by = c("flowering_gene" = "Genes")) %>%
  dplyr::select(gene_id_clean, TIP = ID, 
                pos_tip = POS, Distance_To_Start, Distance_To_End, Orientation, Spatial_Context,
                bio_tips = related_biovariables, pos_snp = POS_SNP, ref_type_SNP, bio_snp = related_biovariables_SNPs,
                gene_id_clean, TAIR, TAIR2, flowering_gene, Category
                )

write.table(final_summary, file.path(OUT_FOLD, "Flowering_gene_annnot_GEA_SNP_TIPs_relation_summary.txt"), sep = "\t", row.names = FALSE)

# # Functional Composition Bar Chart
# group_totals <- final_summary %>% dplyr::group_by(Venn_Group) %>% dplyr::summarise(total = dplyr::n())

# p_func <- ggplot(final_summary, aes(y = Venn_Group, fill = Category)) +
#   geom_bar(stat = "count", position = "stack", alpha = 0.9) +
#   geom_text(data = group_totals, aes(x = total, y = Venn_Group, label = paste0("n = ", total)),
#             inherit.aes = FALSE, hjust = -0.1, size = 2, fontface = "bold") +
#   theme_minimal() +
#   scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
#   scale_fill_viridis_d(option = "turbo") +
#   labs(title = "Functional Composition & Gene Counts per Venn Group", x = "Number of Genes", y = "Venn Group", fill = "Functional Category") +
#   theme(legend.position = "right", legend.text = element_text(size = 8), panel.grid.minor = element_blank())

# ggplot2::ggsave(file.path(out_fold, "Functional_category_distribution.png"), p_func, width = 10, height = 5, dpi = 300)
# ggplot2::ggsave(file.path(out_fold, "Functional_category_distribution.svg"), p_func, width = 10, height = 5)

# ==============================================================================
# 7. Functional Enrichment Analysis Profiling (GO ClusterProfiler)
# ==============================================================================
fasta <- Biostrings::readDNAStringSet("./in/GDDH13_1-1_mrna.fasta")
universe_genes <- unique(stringr::str_replace(names(fasta), ".*gene[:=]([^;| ]+).*", "\\1"))

semdata <- GOSemSim::godata(annoDb = org.Mdomestica.eg.db, keytype = "GID", ont = "BP")

dotplot_list <- list()
emapplot_list <- list()
treeplot_list <- list()
out_cat <- "_Blast"

# for (group in unique(final_gene_summary$Venn_Group)) {

current_genes <- unique(filtered_matches$gene_id_clean)
  
go_bp <- clusterProfiler::enrichGO(
  gene = current_genes, universe = universe_genes, OrgDb = org.Mdomestica.eg.db, 
  keyType = "GID", ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
  readable = TRUE, pAdjustMethod = "fdr"
) # No enrichment

# ==============================================================================
# 8. Local Blast Target Run Validation
# ==============================================================================
rBLAST::makeblastdb(file = "./in/ath.Malus_sylvestris.blast.allgene.cds.fa", dbtype = "nucl")
bl <- rBLAST::blast(db = "./in/ath.Malus_sylvestris.blast.allgene.cds.fa", type = "blastn")
all_blast_list <- list()

group_sequences <- fasta[names(fasta) %in% current_genes]
current_result <- NULL
  
if (length(group_sequences) > 0) {
  blast_results <- stats::predict(bl, group_sequences, BLAST_args = "-evalue 1e-5 -perc_identity 80 -max_target_seqs 1")
  if (nrow(blast_results) > 0) {
    current_result <- blast_results
    #current_result$Venn_Group <- group
    current_result$Status <- "Hit Found"
  } 
}
  
if (is.null(current_result)) {
  current_result <- data.frame(Status = "No Result / No sequences found", stringsAsFactors = FALSE)
}
# all_blast_list[[group]] <- current_result

filtered_matches_blast <- filtered_matches %>%
filter(gene_id_clean %in% current_result$qseqid)

merge_blast_filtered <- merge(filtered_matches_blast, current_result, by.x='gene_id_clean', by.y='qseqid')

merge_blast_filtered <- filtered_matches_blast %>%
  dplyr::left_join(current_result, by = c("gene_id_clean" = "qseqid")) %>%
  dplyr::select(gene_id_clean, TIP = ID, 
                pos_tip = POS, Distance_To_Start, Distance_To_End, Orientation, Spatial_Context,
                bio_tips = related_biovariables, pos_snp = POS_SNP, ref_type_SNP, bio_snp = related_biovariables_SNPs,
                gene_id_clean, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
                )



intersect(merge_blast_filtered$gene_id_clean, final_summary$gene_id_clean)


