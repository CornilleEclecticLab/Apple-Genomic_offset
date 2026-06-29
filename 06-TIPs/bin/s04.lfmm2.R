#!/usr/bin/env Rscript
rm(list = ls())
source("./bin/fonction.R")
source("./bin/config.R")

# output_folder <- paste0("./out/", Project, "/s02.lfmm/")
# input_folder <- "./in/"
# lfmm <- "chr_concat.geno20_maf005.anno_filtered_dom_admixed_rm.lfmm"
# env.1 <- "Bioclim_orchard_chelsa_current.csv"
# vcf <- "chr_concat.geno20_maf005.anno_filtered.vcf.gz"
# computation_folder <- paste0("./out/", Project,"/s01.computation/file/")
# Y from s01.computation
# K=5
# Y <- read.geno(paste0(computation_folder, file_name, "_filt.lfmm_imputed.geno"))


## Input checking
message("Input file checking...")
#check_input_files() # Need to check if the file exist, not needed in the input 
message("...Done.")

## Output folder creation
message("Folders creation...")
create_directory(lfmm_folder)
work.directory <- paste0(lfmm_folder, "file/")
fig_output <- paste0(lfmm_folder, "figure/")
env_output <- paste0(lfmm_folder, 'env/')
create_directory(work.directory)
create_directory(fig_output)
create_directory(env_output)
message("...Done")

## Paths
work.directory = paste0(lfmm_folder, "file/")
fig_output = paste0(lfmm_folder, "figure/")
file_name <- basename(sub("\\.lfmm$", "", lfmm))

## Read tips and geno data
#tip_names_filt <- data.frame(extract_tips_names(vcf))
tip_names_filt <- read.csv(paste0(computational_folder,"file/", file_name, "_TIPs_names_filtered.csv"))[,1]
print(head(tip_names_filt))

if (file.exists(paste0(computational_folder,"file/", file_name, "_filt.lfmm_imputed.geno"))){

  Y <- read.geno(paste0(computational_folder,"file/", file_name, "_filt.lfmm_imputed.geno"))

} else {

  Y <- read.geno(paste0(computational_folder,"file/", file_name, "_filt.geno"))

}

message("lfmm dimension :"); print(dim(lfmm <- read.lfmm(paste0(computational_folder,"file/", file_name, "_filt.lfmm"))))


####----- Make env file -----
message("Creating env files...") # make a loop and create it if it doesn't exist

info_ref_all_SYL = read.table(ref, header = TRUE, sep = '\t')
sorted_info <- info_ref_all_SYL[order(info_ref_all_SYL$Order_vcf), ]

head(sorted_info)

# Define names
env_name <- basename(sub("\\.csv$", "", env.1))
out_file <- file.path(env_output, paste0(env_name, ".env"))

# Check if file already exists

env <- read.delim(env.1, header = TRUE, sep = ",", dec = ".")

env_ordered <- env[match(sorted_info$Seq_ID, env$Seq_ID), ]
head(env_ordered)

env.format.1 <- as.data.frame(scale(env_ordered[, grep("^bio", colnames(env))], center = FALSE, scale = FALSE))
write.env(env.format.1, out_file)

message("...Done.")


####----- LFMM2 -----
message("Running LFMM2...")

# Read environmental data
X <- read.env(paste0(env_output, env_name, ".env"))
mod.lfmm2 <- lfmm2(input = Y, env = X, K = K)

library(ggplot2)
library(qvalue)
library(gridExtra)

# ===== Plot latent factors (mod.lfmm2@U) =====
df_U <- as.data.frame(mod.lfmm2@U)
colnames(df_U) <- paste0("Factor", seq_len(ncol(df_U)))

p_U <- ggplot(df_U, aes(x = Factor1, y = Factor2)) +
  geom_point(color = "grey30", size = 1) +
  labs(x = "Factor 1", y = "Factor 2",
       title = "LFMM2 latent factors") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(fig_output, paste0(file_name, "_plot_mod.lfmm2.png")),
  plot = p_U, width = 8, height = 5, units = "in", dpi = 300
)
ggsave(
  filename = file.path(fig_output, paste0(file_name, "_plot_mod.lfmm2.svg")),
  plot = p_U, width = 8, height = 5, units = "in"
)

# ===== Manhattan plot =====
pv <- lfmm2.test(object = mod.lfmm2, input = Y, env = X, full = TRUE)

df_pv <- data.frame(
  TIPs = tip_names_filt,
  pvalues = as.numeric(pv$pvalues),
  fscores = as.numeric(pv$fscores),
  adj.r.squared = as.numeric(pv$adj.r.squared)
)

plot_histo(df_pv, "TIPs_pvalue_distribution", file.path(fig_output,"TIPs_full_test_pvalue_distribution.png"))

p_manhattan <- ggplot(df_pv, aes(x = seq_along(TIPs), y = -log10(pvalues))) +
  geom_point(color = "grey30", size = 0.6) +
  geom_hline(yintercept = -log10(0.1 / length(pv$pvalues)),
             color = "orange", linetype = "solid", linewidth = 0.8) +
  labs(
    title = "Manhattan plot of log10 p-value of each tip for all climate variables",
    x = "tip Index",
    y = expression(-log[10](pvalue))
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(fig_output, paste0(file_name, "_manhattan_plot.lfmm2.png")),
  plot = p_manhattan, width = 8, height = 5, units = "in", dpi = 300
)
ggsave(
  filename = file.path(fig_output, paste0(file_name, "_manhattan_plot.lfmm2.svg")),
  plot = p_manhattan, width = 8, height = 5, units = "in"
)

# Save table outputs
write.csv(df_pv, file = paste0(work.directory, file_name, "_GEAs_full.csv"), row.names = FALSE)

GEA_tips_bonf <- which(-log10(pv$pvalues) > -log10(0.1 / length(pv$pvalues)))
write.csv(GEA_tips_bonf, file = paste0(work.directory, file_name, "_TIPS_GEA_Bonf.csv"), row.names = FALSE)

q_obj <- qvalue::qvalue(pv$pvalues, fdr.level = 0.05)
GEA_tips_FDR <- which(q_obj$significant)
write.csv(GEA_tips_FDR, file = paste0(work.directory, file_name, "_TIPS_GEA_FDR.csv"), row.names = FALSE)

# ===== Variable-specific plots =====
pv.2 <- lfmm2.test(object = mod.lfmm2, input = Y, env = X, full = FALSE)
df_all <- data.frame()
name_plot <- colnames(env.format.1)

plot_list <- list()

for (i in seq_along(name_plot)) {
  df <- data.frame(
    tips = tip_names_filt,
    pval = as.numeric(pv.2$pvalues[i, ]),
    qval = as.numeric(qvalue(pchisq((pv.2$zscores[i, ] ^ 2) / pv.2$gif,
                                    df = 1, lower = FALSE))$qvalue),
    zsc = as.numeric(pv.2$zscores[i, ])
  )

  # ✅ keep names as is for plotting
  df_named <- df
  names(df_named)[2:4] <- c(paste0("pvalues_", name_plot[i]),
                            paste0("qvalue_", name_plot[i]),
                            paste0("zscores_", name_plot[i]))

  if (nrow(df_all) == 0) {
    df_all <- df_named
  } else {
    df_all <- merge(df_all, df_named, by = "tips")
  }

  # ✅ ggplot uses df (with pval, not renamed)
  p <- ggplot(df, aes(x = seq_along(tips), y = -log10(pval))) +
    geom_point(color = "grey30", size = 0.6) +
    geom_hline(yintercept = -log10(0.1 / length(pv.2$pvalues[i, ])),
               linetype = "dashed", color = "orange") +
    labs(
      title = name_plot[i],
      x = "tip index",
      y = expression(-log[10](pvalue))
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.text = element_text(size = 7)
    )

  plot_list[[i]] <- p
}


# Arrange plots in a grid (4x5)
combined_plot <- marrangeGrob(plot_list, nrow = 4, ncol = 5)

ggsave(
  filename = file.path(fig_output, paste0(file_name, "_GEAs_plot_all_bio.png")),
  plot = combined_plot, width = 8, height = 5, units = "in", dpi = 300
)
ggsave(
  filename = file.path(fig_output, paste0(file_name, "_GEAs_plot_all_bio.svg")),
  plot = combined_plot, width = 8, height = 5, units = "in"
)

# Save variable-specific results
write.csv(df_all, file = paste0(work.directory, file_name, "_GEAs_bio_sep.csv"), row.names = FALSE)
## maybe add the bio sep FDR and FULL
message("...Done.")
