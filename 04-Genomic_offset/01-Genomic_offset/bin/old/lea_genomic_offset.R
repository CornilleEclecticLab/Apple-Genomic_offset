rm(list = ls())
source("./bin/fonction.R")
source("./bin/config.R")

## Input checking
message("Input file checking...")
check_input_files()
message("...Done.")

## Input file reading...
info_ref_all_SYL = read.table(paste0(input_folder,ref), header = T, sep = '\t')  

## Output folder creation...
message("Folders creation...")
create_directory(output_folder)
work.directory = paste0(output_folder,"file/")
fig_output = paste0(output_folder,"figure/")
create_directory(work.directory)
create_directory(fig_output)
message("...Done")

## SNPs name extraction
snp_names <- data.frame(extract_SNPs_names(paste0(input_folder,vcf)))

## file name
#EC_syl/input/EC_syl.fullgenome.light.reformat.lfmm
file_name <- basename(sub("\\.lfmm$", "", paste0(input_folder,lfmm))) # Extract part after last '/' and remove lfmm

## lfmm to geno format
message("Converting lfmm to geno...")
lfmm2geno(paste0(input_folder,lfmm), paste0(work.directory,file_name,".geno"))
geno.file = paste0(work.directory,file_name,".geno")
geno.data = read.geno(geno.file)
message("...Done.")

####----- FILTERING -----
message("FILTERING ...")
# The filtering is going to count and remove column of the geno file with constant value
# only one unique value equal for all the genotype or two value but the second is 9 meaning 
# missing value

all_equal <- apply(geno.data, 2, function(col) {
  length(unique(col)) == 1
})

one_different <- apply(geno.data, 2, check_one_different)
#one_different_1 <- apply(geno.file, 2, check_one_different_1)
#which(one_different)#; which(one_different_1)

column_to_remove <-
  which(all_equal |
          one_different) # column where all value are constant
print(paste("There is", length(column_to_remove), "column(s) to remove"))

geno.data_filt <- geno.data[, -column_to_remove]
snp_names_filt <- snp_names[- column_to_remove, ]

#dim(geno.data)
#dim(geno.data_filt)

write.geno(geno.data_filt, paste0(work.directory,file_name,"_filt.geno"))
geno2lfmm(paste0(work.directory,file_name,"_filt.geno"), paste0(work.directory,file_name,"_filt.lfmm"))

message("...Done.")

####----- Imputation -----
message("Imputation...")
project.missing = snmf(
  paste0(work.directory,file_name,"_filt.lfmm"),
  K = K,
  entropy = TRUE,
  repetitions = nb_rep,
  project = "new"
)
best_run = which.min(cross.entropy(project.missing, K = K))
# Impute the missing genotypes
impute(
  project.missing,
  paste0(work.directory,file_name,"_filt.lfmm"),
  method = 'mode',
  K = K,
  run = best_run
)
message("...Done.")

# Proportion of correct imputation results
dat.imp = read.lfmm(paste0(work.directory,file_name,"_filt.lfmm_imputed.lfmm"))
table(dat.imp)# after the impute we normally don't have missing value (9) anymore

message("lfmm imputated file convertion to geno...")
lfmm2geno(
  paste0(work.directory,file_name,"_filt.lfmm_imputed.lfmm"),
  paste0(work.directory,file_name,"_filt.lfmm_imputed.geno")
)
message("...Done.")

####-----Make env file -----
message("env file creation...")
env <- read.delim(env.1, header = T, sep = '\t', dec = ".")
env.format <- env[c(3:length(env))]
#env.format <- env[c(4:9)]
env.format.1 <- scale(env[c(3:length(env))], center = F, scale = F)
write.env(env.format.1, paste0(input_folder,env.1,'.env'))

next_env <- read.delim(env.2, header = T, sep = ',', dec = '.')
new_env.format <- next_env[c(2:length(next_env))]
new_env.format.1 <- scale(next_env[c(2:length(next_env))], center = F, scale = F)
write.env(new_env.format.1, paste0(input_folder,env.2,'.env'))
message("...Done.")

####----- lfmm2 -----

message("Lfmm2...")
Y <- read.geno(paste0(work.directory,file_name,"_filt.lfmm_imputed.geno"))
X <- read.env(paste0(input_folder,env.1,'.env'))

mod.lfmm2 <- lfmm2(input = Y, env = X, K = K)

png(paste0(fig_output,file_name,'_plot_mod.lfmm2.png'), width = 8, height = 5, units = "in", res = 300)
plot(mod.lfmm2@U, col = "grey", pch = 19,
     xlab = "Factor 1",
     ylab = "Factor 2")
dev.off()

# look more at this plot
png(paste0(fig_output,file_name,'_manhattan_plot.lfmm2.png'), width = 8, height = 5, units = "in", res = 300)
pv <- lfmm2.test(object = mod.lfmm2,
                 input = Y,
                 env = X,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19, xlab = "SNP Index", ylab = "-log10(pvalue)", main = "Manhattan plot of log10 p-value of each SNP for all climate variables")
abline(h = -log10(0.1/length(pv$pvalues)), lty = 1, col = "orange", lwd = 2)
dev.off()

df_pv <- data.frame(
  SNPs = snp_names_filt,
  pvalues = as.numeric(pv$pvalues),
  fscores = as.numeric(pv$fscores),
  adj.r.squared = as.numeric(pv$adj.r.squared)
)

write.csv(df_pv, file = paste0(work.directory,file_name,"_GEAs_full.csv"),row.names = F)

GEA_snp_all <- which(-log10(pv$pvalues) > -log10(0.1/length(pv$pvalues))); length(GEA_snp_all)
write.csv(GEA_snp_all, file = paste0(work.directory,file_name,"_SNPs_associated_with_all_bio.csv"), row.names = F)

#Computing P-values and plotting their minus log10 values
pv.2 <- lfmm2.test(object = mod.lfmm2,
                   input = Y,
                   env = X,
                   full = FALSE)

df_all <- data.frame()

png(paste0(fig_output,file_name,'_GEAs_plot_all_bio.png'), width = 8, height = 5, units = "in", res = 300)
par(mfrow = c(4,5), mar = c(2, 2, 2, 2))
name_plot <- colnames(env.format)

for (i in 1:length(env.format)) {
  print(name_plot[i])
  
  GEA_snp_spe <-
    which(-log10(pv.2$pvalues) > -log10(0.1 / length(pv.2$pvalues)))
  write.csv(
    data.frame(SNP = GEA_snp_spe),
    file = paste0(work.directory, file_name, "_", name_plot[i], ".csv"),
    row.names = F
  )
  
  df <- data.frame(
    SNPs = snp_names_filt,
    pval = as.numeric(pv.2$pvalues[i, ]),
    qval = as.numeric(qvalue(pchisq((pv.2$zscores[i, ] ^ 2) / pv.2$gif,
                                    df = 1,
                                    lower = FALSE
    ))$qvalue),
    zsc = as.numeric(pv.2$zscores[i, ])
  )
  
  # Rename columns dynamically
  names(df)[2:4] <-
    c(
      paste0("pvalues_", name_plot[i]),
      paste0("qvalue_", name_plot[i]) ,
      paste0("zscores_", name_plot[i])
    )
  
  if (dim(df_all)[1] == 0) {
    df_all <- df
  } else {
    df_all <- merge(df_all, df, by = "SNPs")
  }
  
  print(dim(df_all))
  
    plot(
      -log10(pv.2$pvalues[i,]),
      col = "grey",
      cex = .6,
      pch = 19,
      main = paste0(name_plot[i]),
      xlab = "SNP index",
      ylab = "-log(pvalue)"
    )
    abline(
      h = -log10(0.1 / length(pv.2$pvalues[i,])),
      lty = 2,
      col = "orange"
    )
    
}

dev.off()
par(mfrow = c(1,1), mar = c(5,5,5,5))
print(dim(df_all))
write.csv(df_all, file = paste0(work.directory,file_name,"_GEAs_bio_sep.csv"), row.names = F)
message("...Done.")

####----- prediction ----

message("Prediction...")
X.pred <- read.env(paste0(input_folder,env.2,'.env'))

g.gap.scaled <- genetic.gap(input = Y, env = X , pred.env = X.pred, scale = T, K = K)

g.gap.scaled$vectors[,1:19]^2

par(mfrow = c(1,3))
barplot(g.gap.scaled$eigenvalues,
        col = "orange",
        xlab = "Axes",
        ylab = "Eigenvalues")
Delta = X[,1:19] - X.pred[,1:19]
squared.env.dist = rowSums(Delta^2)
#env.dist = rowSums(Delta)

#plot(squared.env.dist, g.gap.scaled$offset, cex = .6)
## il faudrait qu'on puisse colorier les points en fonction du pays d'origine




# ggplot(data = data,
#        aes(x = {{site_var}},
#            y = {{trunk_var}})) + geom_boxplot() + labs(x = labx, y = laby, title = "Boxplot") + theme_bw()

data_test = cbind(squared.env.dist,g.gap.scaled$offset, g.gap.scaled$distance)
colnames(data_test) = c("squared.env.dist","offset","offset_Relstad")
plot(log(data_test[,1]), data_test[,2])
plot(log(data_test[,1]), data_test[,3])
dim(data_test)

sorted_info =  info_ref_all_SYL[order(info_ref_all_SYL$ALL_SYL_wo_admixed),]


png(paste0(fig_output,file_name,"_Genomic_offset_sampling_site_vs_saclay_s01_(offset).png"), width = 10, height = 10, units = "in", res = 300)
ggplot(
  data_test,
  aes(
    x = log(data_test[, 1]),
    y = data_test[, 2] ,
    color = sorted_info$Country.code..ISO.3166.,
    label = sorted_info$Site
  )
) + geom_point(size = 5) + geom_text_repel(size = 3, max.overlaps = 40) + geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(x = "Log(Environmental distance)", y = "Offset value", color = "Country", title = "sampling site vs on Saclay s01")   +
  theme_bw() + theme(legend.text = element_text(size = 12), legend.title = element_text(size = 16), axis.title = element_text(size = 15))
dev.off()

png(paste0(fig_output,file_name,"_Genomic_offset_sampling_site_vs_saclay_s01_(offset_Relstad).png"), width = 10, height = 10, units = "in", res = 300)
ggplot(
  data_test,
  aes(
    x = log(data_test[, 1]),
    y = data_test[, 3] ,
    color = sorted_info$Country.code..ISO.3166.,
    label = sorted_info$Site
  )
) + geom_point(size = 5) + geom_text_repel(size = 3, max.overlaps = 40) + geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(x = "Log(Environmental distance)", y = "Offset value - Relstad", color = "Country", title = "sampling site vs on Saclay s01")   +
  theme_bw() + theme(legend.text = element_text(size = 12), legend.title = element_text(size = 16), axis.title = element_text(size = 15))
dev.off()


### with candidate loci
GEA_snp <- read.csv(paste0(work.directory,file_name,"_SNPs_associated_with_all_bio.csv"), header = T)
print(GEA_snp)
list_GEA_snp <- GEA_snp$x
print(list_GEA_snp)
g.gap.scaled.bis <- genetic.gap(input = Y, env = X , pred.env = X.pred, scale = T, K = K, candidate.loci = list_GEA_snp)

par(mfrow = c(1,3))
barplot(g.gap.scaled.bis$eigenvalues,
        col = "orange",
        xlab = "Axes",
        ylab = "Eigenvalues")


data_test_bis = cbind(squared.env.dist,g.gap.scaled.bis$offset, g.gap.scaled.bis$distance)
colnames(data_test_bis) = c("squared.env.dist","offset", "offset_Relstad")
plot(log(data_test_bis[,1]), data_test_bis[,2])
plot(log(data_test_bis[,1]), data_test_bis[,3])
dim(data_test)


png(paste0(fig_output,file_name,"_Genomic_offset_sampling_site_vs_saclay_s01_GEA_SNP_(offset).png"), width = 10, height = 10, units = "in", res = 300)
ggplot(
  data_test_bis,
  aes(
    x = log(data_test_bis[, 1]),
    y = data_test_bis[, 2] ,
    color = sorted_info$Country.code..ISO.3166.,
    label = sorted_info$Site
  )
) + geom_point(size = 5) + geom_text_repel(size = 3, max.overlaps = 40) + geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(x = "Log(Environmental distance)", y = "Offset value", color = "Country", title = "all syl wo ITA on sampling site vs on Saclay s01 (GEA SNP)")   +
  theme_bw() + theme(legend.text = element_text(size = 12), legend.title = element_text(size = 16), axis.title = element_text(size = 15))
dev.off()

png(paste0(fig_output,file_name,"_Genomic_offset_sampling_site_vs_saclay_s01_GEA_SNP_(offset_Relstad).png"), width = 10, height = 10, units = "in", res = 300)
ggplot(
  data_test_bis,
  aes(
    x = log(data_test_bis[, 1]),
    y = data_test_bis[, 3] ,
    color = sorted_info$Country.code..ISO.3166.,
    label = sorted_info$Site
  )
) + geom_point(size = 5) + geom_text_repel(size = 3, max.overlaps = 40) + geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  labs(x = "Log(Environmental distance)", y = "Offset value - Relstad", color = "Country", title = "all syl wo ITA on sampling site vs on Saclay s01 (GEA SNP)")   +
  theme_bw() + theme(legend.text = element_text(size = 12), legend.title = element_text(size = 16), axis.title = element_text(size = 15))
dev.off()
message("...Done.")