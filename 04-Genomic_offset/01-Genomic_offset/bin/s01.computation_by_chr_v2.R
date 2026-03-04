#!/usr/bin/env Rscript

## ===============================
##   s01.computation_by_chr.R
## ===============================

# ---- Clean environment ----
rm(list = ls())

# ---- Arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript s01.computation_by_chr_v2.R <vcf_file> <project_name> <K_value> <lfmm_file> <position_data>")
}

vcf_file <- args[1]
PROJECT <- args[2]
K <- as.numeric(args[3])
lfmm_file <- args[4]
position.data <- args[5]

file_name <- basename(sub("\\.vcf\\.gz$", "", vcf_file))
input_folder <- "./in/"
computational_folder <- paste0("./out/", PROJECT, "/s01.Computation/")
nb_rep <- 20

cat("=====================================\n")
cat("Running computation for:", file_name, "\n")
cat("Project:", PROJECT, "\n")
cat("K value:", K, "\n")
cat("LFMM file:", lfmm_file, "\n")
cat("Position file:", position.data, "\n")
cat("=====================================\n\n")

# ---- Load dependencies ----
source("./bin/fonction.R")
#source("./bin/config.R")

# ---- Define folders ----
computational_folder <- paste0("./out/", PROJECT, "/s01.computation/", file_name, "/")

# ---- Check inputs ----
message("Checking input files...")
if (!file.exists(lfmm_file)) stop("LFMM file not found: ", lfmm_file)
if (!file.exists(vcf_file)) stop("VCF file not found: ", vcf_file)
message("...Done.")

# ---- Create output folders ----
message("Creating output folders...")
create_directory(computational_folder)
work.directory <- file.path(computational_folder, "file/")
fig_output <- file.path(computational_folder, "figure/")
create_directory(work.directory)
create_directory(fig_output)
message("...Done.")

# ---- Extract SNP names ----
message("Extracting SNP names...")
snp_names <- data.frame(extract_SNPs_names(vcf_file, file_name))
message("...Done.")

# ---- LFMM → GENO ----
message("Converting LFMM to GENO...")
lfmm2geno(lfmm_file, paste0(work.directory, file_name, ".geno"))
geno.file <- paste0(work.directory, file_name, ".geno")
geno.data <- read.geno(geno.file)
message("...Done.")

# ---- Filtering ----
message("Filtering data...")
all_equal <- apply(geno.data, 2, function(col) length(unique(col)) == 1)
one_different <- apply(geno.data, 2, check_one_different)
column_to_remove <- which(all_equal | one_different)
message("Columns to remove: ", length(column_to_remove))

if (length(column_to_remove) != 0 ){
geno.data_filt <- geno.data[, -column_to_remove]
snp_names_filt <- snp_names[-column_to_remove, ]
write.csv(snp_names_filt, file = paste0(work.directory, file_name, "_snp_names_filtered.csv"), row.names = FALSE)
} else {
message("There is no column to remove!")
geno.data_filt <- geno.data
snp_names_filt <- snp_names
write.csv(snp_names_filt, file = paste0(work.directory, file_name, "_snp_names_filtered.csv"), row.names = FALSE)
}

write.geno(geno.data_filt, paste0(work.directory, file_name, "_filt.geno"))
geno2lfmm(paste0(work.directory, file_name, "_filt.geno"), paste0(work.directory, file_name, "_filt.lfmm"))
message("...Done filtering.")

# ---- Imputation ----
message("Running imputation...")
project.missing <- snmf(
  paste0(work.directory, file_name, "_filt.lfmm"),
  K = K,
  entropy = TRUE,
  repetitions = nb_rep,
  project = "new"
)
best_run <- which.min(cross.entropy(project.missing, K = K))
## get population structure for each chromosome to see if it's conserved
barchart(
  project.missing,
  K = K,
  run = best_run,
  border = NA,
  space = 0,
  sort.by.Q = F,
  main = paste("K =", K)
) -> bp
print(bp$order)

#position.data <- "/scratch/mc10739/02-Wild_apple_genomic_offset/Apple-Genomic_offset/01-Population_structure/01-snmf/in/data_tree_dom_admixed_rm_V2.txt"
list_accession = read.delim(file = position.data, header = T, sep='\t', dec = ",")
list_accession_ordered <- list_accession %>%
mutate(Order = factor(Order_vcf, levels = bp$order)) %>%
arrange(Order_vcf)
get_Q = Q(project.missing, K = K, run = best_run)
get_Q_bind <- cbind(get_Q, list_accession_ordered)
get_Q_bind <- get_Q_bind %>% arrange(Order_plot)

# Identify V columns dynamically
v_columns <- grep("^V", names(get_Q_bind), value = TRUE)
print(length(v_columns))

# Convert to long format
data_long <- pivot_longer(get_Q_bind, cols = starts_with("V"), names_to = "Variable", values_to = "Proportion")
data_long$Seq_ID <- factor(data_long$Seq_ID, levels = get_Q_bind$Seq_ID)
df_wide <- pivot_wider(data_long, names_from = Variable, values_from = Proportion)
write.csv(df_wide, paste0(work.directory,"/K",K,"_value_sorted.csv"),row.names = FALSE)

# Create the barplot
color_palette <- c("V1" = "blue", "V2" = "red", "V3" = "green", "V4" = "purple", 
                   "V5" = "orange", "V6" = "brown", "V7" = "pink", "V8" = "cyan", 
                   "V9" = "yellow", "V10" = "gray")

plot <- ggplot(data_long, aes(x = Seq_ID, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Format y-axis as percentages
  labs(title = paste("K =", K), x = "Individual", y = "Proportion") +
  scale_fill_manual(values = color_palette) +  # Specify colors dynamically
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #legend.position = "none"
  )

png_file <- paste0(fig_output,"/SNPs_based_structure_",file_name,".png")
ggsave(filename = png_file, plot = plot, width = 8, height = 5, units = "in", dpi = 300)
## 

impute(
  project.missing,
  paste0(work.directory, file_name, "_filt.lfmm"),
  method = "mode",
  K = K,
  run = best_run
)

dat.imp <- read.lfmm(paste0(work.directory, file_name, "_filt.lfmm_imputed.lfmm"))
table(dat.imp)

lfmm2geno(
  paste0(work.directory, file_name, "_filt.lfmm_imputed.lfmm"),
  paste0(work.directory, file_name, "_filt.lfmm_imputed.geno")
)
message("...Done imputation.")

message("✅ Finished computation for: ", file_name)
