# Plotting number of significantly changed reactions and overrepresented pathways for both datasets.

cat("\014") #clear the console (CTRL+L)
rm(list=ls()) #remove all objects(variables)

# Libraries
library(tidyverse)
library(ggplot2)
library(readxl)

Table2_barplot <- read_excel(".\\Table2_barplot_AD_ROSMAP_INIT.xlsx") # Data for bar plot. It was prepared for iMAT and INIT based results.

data1<- Table2_barplot
data1$dataset<- as.factor(data1$dataset)
ad<- filter(data1, dataset=="ROSMAP")
ad<-mutate(ad,log=log2(number))
luad<-filter(data1, dataset=="LUAD")
luad<-mutate(luad,log=log2(number))

# ROSMAP-iMAT
ggplot(ad, aes(y=log, x=factor(normalization, level=c('DESeq2','DESeq2_cov','edgeR','edgeR_cov', 'GeTMM','GeTMM_cov', 'FPKM','FPKM_cov','TPM','TPM_cov')), fill= factor(condition,levels=c("pathway","reaction")))) + 
  geom_bar(position="stack",stat="identity") +   #position="dodge"
  scale_fill_manual(values=c('#bdd7e7', '#6baed6'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="ROSMAP cohort", fill = "")
ggsave(file="ROSMAP_iMAT_bar.png", width=7, height=6, dpi=1200)

# TCGA-iMAT
ggplot(luad, aes(y=log, x=factor(normalization, level=c('DESeq2','DESeq2_cov','edgeR','edgeR_cov', 'GeTMM','GeTMM_cov', 'FPKM','FPKM_cov','TPM','TPM_cov')), fill= factor(condition,levels=c("pathway","reaction")))) + 
  geom_bar(position="stack",stat="identity") +   #position="dodge"
  scale_fill_manual(values=c('#bdd7e7', '#6baed6'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="TCGA cohort", fill = "")
ggsave(file="TCGA_iMAT_bar.png", width=7, height=6, dpi=1200)

# ROSMAP-INIT
ggplot(ad, aes(y=log, x=factor(normalization, level=c('DESeq2_cov','edgeR_cov', 'GeTMM_cov', 'FPKM_cov','TPM_cov')), fill= factor(condition,levels=c("pathway","reaction")))) + 
  geom_bar(position="stack",stat="identity") +   #position="dodge"
  scale_fill_manual(values=c('#bdd7e7', '#6baed6'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="ROSMAP cohort", fill = "")
ggsave(file="ROSMAP_INIT_bar.png", width=7, height=6, dpi=1200)
