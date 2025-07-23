library(vegan)
library(dplyr)
library(FactoMineR)
library(corrplot)
library(factoextra)
library(ggplot2)

Project="dom_admixed_rm"
out_folder=paste0("./out/",Project,"/s03.IBE/")

dir.create(out_folder,recursive = T)

## genomic distance
gen <- read.csv(paste0("./out/",Project,"/s01.Metrics/matrix_fst.csv"), sep = ",", dec = ".", header = T)
gen_dist <- as.dist(gen)
gen_matrix <- as.matrix(gen_dist)

## geographic distance
geo <- read.csv(paste0("./out/",Project,"/s01.Metrics/matrix_dist.csv"), sep = ",", dec = ",", header = T)
geo_dist <- as.dist(geo)
geo_matrix <- as.matrix(geo_dist)


# Applying a Principal Coordinates of Neighbourhood Matrix (PCNM)
# following Borcard & Legendre 2002):

geo_pcnm <- pcnm(geo_dist)

min <- geo_pcnm$threshold #truncation distance
nb_ev_geo <- length(geo_pcnm$values) #number of eigenvalues I have 54 eigenvalue

#Selecting only the significantly positive axis
select_geo <- which(geo_pcnm$values >= 0) #I will select only positive eigenvalue
geo_pcnm_pos <-  as.data.frame(geo_pcnm$vectors)[, select_geo] #29 on 54 with dataframe transformation


env <- read.csv(paste0('../../02-Climate_data/out/',Project,'/Bioclim_sample_dom_admixed_rm_chelsa_current.csv'), sep = '\t', dec=".", header = T)
env <- env %>%
  rename_with(~ paste0(., "_current"), .cols = starts_with("Bio"))
rownames(env) <- env$Group
env <- env[-c(1)]

LGM <- read.csv(paste0('../../02-Climate_data/out/',Project,'/Bioclim_sample_dom_admixed_rm_chelsa_LGM.csv'), sep = "\t", dec = ".", header = T)
LGM <- LGM %>%
  rename_with(~ paste0(., "_LGM"), .cols = starts_with("Bio"))
rownames(LGM) <- LGM$Group
LGM <- LGM[-c(1)]

env_pcnm <- merge(env, geo_pcnm_pos, by="row.names", all=TRUE, sort=FALSE)
rownames(env_pcnm) <- env_pcnm$Row.names
env_pcnm <- env_pcnm[-c(1)]

env_pcnm_LGM <- merge(env_pcnm, LGM, by="row.names", all = TRUE, sort=FALSE)

env_pcnm_LGM <- env_pcnm_LGM[-c(1)]

write.csv(env_pcnm_LGM, paste0(out_folder,"env_pcnm_LGM.csv"), row.names = F)

corrplot(round(cor(env_pcnm_LGM[1:19]), digits=2))
corrplot(round(cor(env_pcnm_LGM[20:23]), digits = 2))
corrplot(round(cor(env_pcnm_LGM[24:42]), digits = 2))

PCA(env_pcnm_LGM)


#----dbRDA avec capscale,

rdaNullModel <- capscale(gen_dist ~ 1, data=env_pcnm_LGM)

scale_env_pcnm_LGM <- as.data.frame(scale(env_pcnm_LGM))# to scale each variable
rdaFullModel <- capscale(gen_dist ~ scale(Bio1_current) + scale(Bio2_current) + scale(Bio3_current) +
                           scale(Bio4_current) + scale(Bio5_current) + scale(Bio6_current) + scale(Bio7_current) +
                           scale(Bio8_current) + scale(Bio9_current) + scale(Bio10_current) + scale(Bio11_current) + 
                           scale(Bio12_current) + scale(Bio13_current) + scale(Bio14_current) + scale(Bio15_current) +
                           scale(Bio16_current) + scale(Bio17_current) + scale(Bio18_current) + scale(Bio19_current) +
                           scale(PCNM1) + scale(PCNM2) + scale(PCNM3) + scale(PCNM4) + 
                           scale(Bio1_LGM) + scale(Bio2_LGM) + scale(Bio3_LGM)+  scale(Bio4_LGM)+ scale(Bio5_LGM)+ 
                           scale(Bio6_LGM)+ scale(Bio7_LGM)+ scale(Bio8_LGM)+ scale(Bio9_LGM)+ scale(Bio10_LGM)+ 
                           scale(Bio11_LGM)+ scale(Bio12_LGM) + scale(Bio13_LGM) + scale(Bio14_LGM) + scale(Bio15_LGM)+ 
                           scale(Bio16_LGM) + scale(Bio17_LGM) + scale(Bio18_LGM) + scale(Bio19_LGM), 
                         data=scale_env_pcnm_LGM)#; anova(rdaFullModel)
# round(extractAIC(rdaFullModel), digits = 3)
# vif.cca(rdaFullModel)

rda1 <- ordistep(rdaNullModel, scope=formula(rdaFullModel),
                 perm.max = 10000, trace=TRUE)
vif.cca(rda1)
plot(rda1)
anova(rda1)
anova(rda1, by="mar")
#varpart(gen_dist, ~ scale(PCNM1), scale_env_pcnm_LGM)

anova.cca(rda1, by="mar")
extractAIC(rda1)

#Runing best model
rda2 <- capscale(gen_dist ~ Bio4_current, 
                 data=scale_env_pcnm_LGM)
extractAIC(rda2)
anova(rda2)
anova(rda2, by="mar")

result <- anova(rda2, by = "mar")
# Save the output to a text file
capture.output(result, file = paste0(out_folder,"IBE_best_model_anova_result.txt"))
