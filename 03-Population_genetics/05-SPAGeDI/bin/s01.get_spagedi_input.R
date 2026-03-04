rm(list=ls())
library(vcfR)

# ================== INPUT FILES ==================
vcf_file  <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/01-Population_structure/03-membership_coefficient_investigating/out/All_data_wo_orientalis/s04.filter_vcf/chr_concat.geno20_maf005.anno.syno.LDpruned_dom_admixed_rm.vcf.gz"
meta_file <- "./in/data_tree_dom_admixed_rm_V2_grouped_site.txt"
out_file  <- "./out/spagedi.dat"

group_levels <- c(
  "Aaby","AlDom","Hvol","Jur_PdD","Lor_BasRh","Norh","NPDA",
  "Site_13_14_15_16","Site0_1_2","Skov","STO","syl_be_WE","Yv_DeSe"
)

# ================== METADATA ==================
meta <- read.table(meta_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
meta <- meta[, c("Seq_ID","Site_grouped","Grouped_Latitude","Grouped_Longitude")]
colnames(meta) <- c("ID","Group","Lat","Long")

# Fix European decimals
meta$Lat  <- as.numeric(gsub(",", ".", meta$Lat))
meta$Long <- as.numeric(gsub(",", ".", meta$Long))

# Map groups to integers
meta$GROUP_NUM <- match(meta$Group, group_levels)
stopifnot(!any(is.na(meta$GROUP_NUM)))

# ================== READ VCF ==================
vcf <- read.vcfR(vcf_file)
gt  <- extract.gt(vcf, element="GT")   # SNP x IND
ref <- vcf@fix[,"REF"]
alt <- vcf@fix[,"ALT"]

# SNP names
snp_names <- paste(vcf@fix[,"CHROM"], vcf@fix[,"POS"], sep="_")

# ================== MATCH INDIVIDUAL ORDER ==================
ind_ids <- colnames(gt)
meta <- meta[match(ind_ids, meta$ID), ]
stopifnot(all(meta$ID == ind_ids))

# ================== GENOTYPE CONVERSION ==================
allele_code <- c(A=1, C=2, G=3, T=4)

n_ind <- ncol(gt)
n_snp <- nrow(gt)
spagedi_mat <- matrix(0, nrow=n_ind, ncol=n_snp)

for(i in 1:n_snp) {
  
  if(!(ref[i] %in% names(allele_code) && alt[i] %in% names(allele_code))) {
    spagedi_mat[,i] <- 0
    next
  }
  
  for(j in 1:n_ind) {
    g <- gt[i,j]
    
    if(is.na(g) || g %in% c("./.", ".|.")) {
      spagedi_mat[j,i] <- 0
    } else {
      a <- strsplit(g, "[/|]")[[1]]
      a <- ifelse(a=="0", ref[i], alt[i])
      a <- sort(a)
      spagedi_mat[j,i] <- as.integer(paste0(allele_code[a[1]], allele_code[a[2]]))
    }
  }
}

# ================== WRITE SPAGEDI FILE ==================
con <- file(out_file, "w")

# ---- Line 1 ----
writeLines(paste(n_ind, length(unique(meta$GROUP_NUM)), -2, n_snp, 1, 2, sep="\t"), con)

# ---- Line 2: automatic distance classes ----
writeLines("-7", con)

# ---- Line 3: column names ----
header_line <- c("ID","GROUP","LAT","LONG", snp_names)
writeLines(paste(header_line, collapse="\t"), con)

# ---- Individuals ----
for(i in 1:n_ind) {
  line <- c(meta$ID[i], meta$GROUP_NUM[i], meta$Lat[i], meta$Long[i], spagedi_mat[i,])
  writeLines(paste(line, collapse="\t"), con)
}

# ---- END ----
writeLines("END", con)
close(con)

cat("SPAGeDi file written:", out_file, "\n")


# ================== WRITE ONE SPAGEDI FILE PER GROUP ==================

unique_groups <- unique(meta$Group)

for(gname in unique_groups) {
  
  message("Processing group: ", gname)
  
  # Subset individuals
  idx <- which(meta$Group == gname)
  meta_sub <- meta[idx, ]
  gt_sub   <- spagedi_mat[idx, , drop=FALSE]
  
  # Skip if <2 individuals
  if(nrow(meta_sub) < 2) {
    warning("Skipping ", gname, " (less than 2 individuals)")
    next
  }
  
  # SPAGeDi group numbering must start at 1
  meta_sub$GROUP_NUM <- 1
  
  n_ind_sub <- nrow(meta_sub)
  n_snp     <- ncol(gt_sub)
  
  out_file_sub <- paste0("./out/spagedi_", gname, ".dat")
  con <- file(out_file_sub, "w")
  
  # ---- Line 1 ----
  writeLines(paste(n_ind_sub, 1, -2, n_snp, 1, 2, sep="\t"), con)
  
  # ---- Line 2 ----
  writeLines("-7", con)
  
  # ---- Line 3 ----
  header_line <- c("ID","GROUP","LAT","LONG", snp_names)
  writeLines(paste(header_line, collapse="\t"), con)
  
  # ---- Individuals ----
  for(i in 1:n_ind_sub) {
    line <- c(meta_sub$ID[i], meta_sub$GROUP_NUM[i],
              meta_sub$Lat[i], meta_sub$Long[i],
              gt_sub[i, ])
    writeLines(paste(line, collapse="\t"), con)
  }
  
  writeLines("END", con)
  close(con)
  
  cat("Written:", out_file_sub, "\n")
}

### IN SPAGEDI ###

