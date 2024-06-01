
library(readxl)
library(gtsummary)
library(ggplot2)


dat=read.csv('lung_genes_sbna.csv')

head(dat)



ggplot(dat, aes(x=mutant, y=network_score))+geom_jitter()+geom_boxplot(color="red")+facet_grid(~cancer_role)

#Edit cancer roles manually to remove mixed categorizations
dat$cancer_role_editted=dat$cancer_role
dat$cancer_role_editted[which(dat$gene=="TP53")]="TSG"
dat$cancer_role_editted[which(dat$gene=="NOTCH1")]="oncogene"
dat$cancer_role_editted[which(dat$gene=="ERBB4")]="oncogene"
dat$cancer_role_editted[which(dat$cancer_role=="oncogene, fusion")]="oncogene"
table(dat$cancer_role_editted,dat$cancer_role)
ggplot(dat, aes(x=mutant, y=network_score))+geom_jitter()+geom_boxplot(color="red")+facet_grid(~cancer_role_editted)
ggplot(dat, aes(x=mutant, y=SecondOrderIntermodularDegree_AVERAGE))+geom_jitter()+geom_boxplot(color="red")+facet_grid(~cancer_role_editted)
ggplot(dat, aes(x=mutant, y=NodeEdgeBetweennessSTRIDE_sidechain_MAX))+geom_jitter()+geom_boxplot(color="red")+facet_grid(~cancer_role_editted)

wilcox.test(dat$network_score[which(dat$cancer_role_editted=="TSG" & dat$mutant=="N")],dat$network_score[which(dat$cancer_role_editted=="TSG" & dat$mutant=="Y")])
wilcox.test(dat$network_score[which(dat$cancer_role_editted=="oncogene" & dat$mutant=="N")],dat$network_score[which(dat$cancer_role_editted=="oncogene" & dat$mutant=="Y")])
wilcox.test(dat$network_score[which(dat$cancer_role_editted=="fusion" & dat$mutant=="N")],dat$network_score[which(dat$cancer_role_editted=="fusion" & dat$mutant=="Y")])




table(dat$cancer_role)
unique(dat$gene[which(dat$cancer_role=="oncogene")])
unique(dat$gene[which(
  dat$cancer_role=="fusion")])
unique(dat$gene[which(dat$cancer_role=="TSG")])
unique(dat$gene[which(dat$cancer_role=="oncogene, fusion")])
[1] "BRAF"  "FGFR2" "ERBB2" "HGF"   "BIRC6" "NTRK2" "ROS1"  "HIP1"  "ALK" 
unique(dat$gene[which(dat$cancer_role=="oncogene, TSG")])
[1] "ERBB4"
unique(dat$gene[which(dat$cancer_role=="oncogene, TSG, fusion")])
[1] "NOTCH1" "TP53"  

