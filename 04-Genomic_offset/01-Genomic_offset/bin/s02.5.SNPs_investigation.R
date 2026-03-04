library(vcfR)
library(dplyr)

# 1. SET YOUR TARGETS
# Names of individuals you want to compare (must match VCF column names exactly)
target_indivs <- c("E15", "E16", "P35B2_L1P41", "P38B4_L1P48", "P78D1_L1P1")

# Load your SNP position file
snp_info <- read.csv("./out/dom_admixed_rm/s02.lfmm2/file/All_chr_SNPs_associated_with_all_bio_VCF_position_FDR.csv", header = TRUE)

library(vcfR)
library(dplyr)

# 1. LOAD YOUR TARGETS
snp_info <-  read.csv("./out/dom_admixed_rm/s02.lfmm2/file/All_chr_SNPs_associated_with_all_bio_VCF_position_FDR.csv", header = TRUE)
target_indivs <- c("E15", "E16", "P35B2_L1P41", "P38B4_L1P48", "P78D1_L1P1")

# 2. CREATE A "REGION" FILTER
# We only want to read the parts of the 10M SNP file that match our target chromosomes
target_chroms <- unique(snp_info$chromosome)

# 3. USE A LOOP OR CONCATENATED COMMAND TO READ BY CHROMOSOME
# This prevents loading all 10M SNPs at once
results_list <- list()

for(chr in target_chroms){
    message("Processing ", chr, "...")
    
    # Subset SNP info for this chromosome
    chr_snps <- snp_info %>% filter(chromosome == chr)
    
    # Read ONLY this chromosome from the VCF (requires tabix)
    # If tabix isn't available, we filter during reading
    temp_vcf <- read.vcfR("./in/chr_concat.geno20_maf005.anno_filtered_dom_admixed_rm.vcf.gz", verbose = FALSE) 
    
    # Filter individuals immediately to save RAM
    gt <- extract.gt(temp_vcf, element = "GT")
    gt <- gt[, target_indivs, drop = FALSE]
    
    # Filter by specific positions
    vcf_pos <- getPOS(temp_vcf)
    vcf_chrom <- getCHROM(temp_vcf)
    
    keep_idx <- which(vcf_chrom == chr & vcf_pos %in% chr_snps$SNPs)
    results_list[[chr]] <- gt[keep_idx, , drop = FALSE]
}

# 4. COMBINE
final_gt_matrix <- do.call(rbind, results_list)

# 2. LOAD VCF
vcf <- read.vcfR("./in/chr_concat.geno20_maf005.anno_filtered_dom_admixed_rm.vcf.gz")

# 3. FILTER BY INDIVIDUALS FIRST
# We extract genotypes only for the columns in our target list
gt <- extract.gt(vcf)

# Check if all your target individuals exist in the VCF
missing_indivs <- setdiff(target_indivs, colnames(gt))
if(length(missing_indivs) > 0) {
  stop("These individuals are not in the VCF: ", paste(missing_indivs, collapse=", "))
}

gt_filtered_ind <- gt[, target_indivs, drop = FALSE]

# 4. FILTER BY SNP POSITION
# Create keys for the VCF and your CSV to find the overlap
vcf_keys <- paste(getCHROM(vcf), getPOS(vcf), sep = "_")
snp_info <- snp_info %>%
  mutate(match_key = paste(chromosome, SNPs, sep = "_"))

# Keep only the rows that match your SNP file
final_gt <- gt_filtered_ind[vcf_keys %in% snp_info$match_key, , drop = FALSE]

# 5. CONVERT TO NUMERIC (0 = Ref, 1 = Var, NA = Missing)
is_variant <- apply(final_gt, 2, function(x) {
  ifelse(x == "0/0" | x == "0|0", 0, 
         ifelse(grepl("[1-9]", x), 1, NA))
})

# 6. RESULTS
comparison_df <- as.data.frame(is_variant)
print(comparison_df)