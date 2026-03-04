## config file

# Parameters
Project <- "dom_admixed_rm_PCA" # all
K=5 # all
nb_rep=20 # computation

# folder
input_folder <- "./in/"
computational_folder <- paste0("./out/", Project, "/s01.Computation/")
lfmm_folder <- paste0("./out/", Project, "/s02.lfmm2/")
go_folder <- paste0("./out/", Project, "/s03.GO/")

# file
lfmm <- "../00-File_preparation/out/dom_admixed_rm/s01.lfmm_file/chr_chr01.geno20_maf005.anno_filtered_dom_admixed_rm.lfmm" #lfmm file, built from a vcf and turn into a .012 file (only a change in the extension is needed)
# full pathway need to be available
env.1 <- "../../02-Climate_data/out/dom_admixed_rm_PCA/s01.PCA_and_correlation/PCA/PCA_coord_wild_tree.csv" #local environment
# full pathway need to be available
env.2 <- "./in/env_Saclay_Chelsa"
# full pathway need to be available
vcf <- "./in/chr_concat.geno20_maf005.anno_filtered_dom_admixed_rm.vcf.gz" #vcf from which the lfmm was built, use in the GEA for SNP id
# full pathway need to be available
ref <- "./in/data_tree_dom_admixed_rm.txt"
# in file, but full pathway need to be available






## work folder
# Project <- "dom_admixed_rm"
# input_folder <- "./in/" #input folder, in which you need to have all your input file
# output_folder <- paste0("./out/",Project,"/") #output folder, not necessary need to be built (two folder will be create in it :
# 		 #file for file output
# 		 #figure for figure output 

# ## input file (only the name)
# lfmm <- "" #lfmm file, built from a vcf and turn into a .012 file (only a change in the extension is needed)
# env.1 <- "" #local environment
# env.2 <- "" #new environment
# vcf <- "" #vcf from which the lfmm was built, use in the GEA for SNP id
# ref <- "" #data need for legend

# ## parameters
# K <- 5 #number of ancestral population (STRUCTURE)
# nb_rep <- 20 #number of repetition


# ##------- independant script

# ## full pipeline
# lfmm <- "chr_concat.geno20_maf005.anno_filtered_dom_admixed_rm.lfmm" #lfmm file, built from a vcf and turn into a .012 file (only a change in the extension is needed)
# env.1 <- "Bioclim_orchard_chelsa_current.csv" #local environment
# env.2 <- "env_Saclay_Chelsa"
# vcf <- "chr_concat.geno20_maf005.anno_filtered.vcf.gz" #vcf from which the lfmm was built, use in the GEA for SNP id
# ref <- "dom_admixed_rm_info_sample.txt" #data need for legend

# ##Lfmm2
# ## input file (only the name)
# ##Offset prediction
# #input folder is the same
# # Y.offset_file = "EC_syl.fullgenome.light.reformat_filt.lfmm_imputed.geno"
# # X.offset_file = "env_EC_syl"
# # #ref is the same
# # X.pred.folder = "new_env/"
# # GEA_snps_file = "EC_syl.fullgenome.light.reformat_SNPs_associated_with_all_bio.csv"