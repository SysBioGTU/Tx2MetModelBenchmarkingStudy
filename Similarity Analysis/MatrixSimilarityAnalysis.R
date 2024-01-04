# You can use this code to generate similarity heatmaps in the manuscript (in Figure 2 and Figure 3).
setwd("./")
library(openxlsx)
library(readxl)
library(ggheatmap)
library(ggplot2)
library(pheatmap)

# Jaccard Similarity
jaccard_similarity <- function(matrix1, matrix2) {
  intersection <- sum(matrix1 & matrix2)
  union <- sum(matrix1 | matrix2)
  jaccard_similarity <- intersection / union
  return(jaccard_similarity)
}

jaccard_index_rxn <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

## LUAD matrix similarity ##
sheet_names <- getSheetNames("LUAD_binary.xlsx")

coefs<-matrix(nrow=5, ncol=5)
for (i in 1:5){
  met1=read_xlsx("LUAD_binary.xlsx", sheet=i,col_names = TRUE) 
  for (j in 1:5){
    met2=read_xlsx("LUAD_binary.xlsx", sheet=j, col_names = TRUE)
    coefs[i,j] <-jaccard_similarity(met1, met2)
  }}
colnames(coefs)<-c("DeSeq2", "edgeR", "FPKM", "GeTMM", "TPM") 
rownames(coefs)<-c("DeSeq2", "edgeR", "FPKM", "GeTMM", "TPM")  

write.csv(data.frame(coefs), "LUAD_coefs_binarymatrix_personalized.csv")
pheatmap(coefs, color = hcl.colors(7, "YlOrBr", rev=T), angle_col = 45)
rm(coefs)

## LUAD reaction similarity ##

sheet_names<-getSheetNames("TCGA_LUAD_Summarised.xlsx")
coefs<-matrix(nrow=5, ncol=5)
for (i in 1:5){
  met1=read_xlsx("TCGA_LUAD_Summarised.xlsx", sheet=i)
  for (j in 1:5){
    met2=read_xlsx("TCGA_LUAD_Summarised.xlsx", sheet=j)
    coefs[i,j] <-jaccard_index_rxn(met1$Reaction, met2$Reaction)
  }}
colnames(coefs)<-c("DeSeq2", "edgeR", "FPKM", "GeTMM", "TPM") 
rownames(coefs)<-c("DeSeq2", "edgeR", "FPKM", "GeTMM", "TPM")  

pheatmap(coefs, color = hcl.colors(7, "YlOrBr", rev=T), angle_col = 45)

## LUAD pathway similarity ##
sheet_names<-getSheetNames("LUAD_overrepresentation_hypergeo.xlsx") 

coefs<-matrix(nrow=5, ncol=5)
for (i in 1:5){
  met1=read_xlsx("LUAD_overrepresentation_hypergeo.xlsx", sheet=i)
  for (j in 1:5){
    met2=read_xlsx("LUAD_overrepresentation_hypergeo.xlsx", sheet=j)
    coefs[i,j] <-jaccard_index_rxn(met1$Pathway, met2$Pathway)
  }
}
colnames(coefs)<-c("DeSeq2", "edgeR", "FPKM", "GeTMM", "TPM") 
rownames(coefs)<-c("DeSeq2", "edgeR", "FPKM", "GeTMM", "TPM")  

pheatmap(coefs, color = hcl.colors(7, "YlOrBr", rev=T), angle_col = 45)

## ROSMAP matrix similarity##
setwd("./Benchmark_normalization")
load("./ROSMAP_binary_data.RData")
ROSMAP.AD.Binary_01<- readRDS("ROSMAP_Noise_01.RDS")
ROSMAP.AD.Binary_005<- readRDS("ROSMAP_Noise_005.RDS")

data<-list(original=ROSMAP.AD.Binary,noise01=ROSMAP.AD.Binary_01,noise005=ROSMAP.AD.Binary_005)
coefs<-matrix(nrow=length(ROSMAP.AD.Binary), ncol=length(ROSMAP.AD.Binary))
for (i in 1:length(ROSMAP.AD.Binary)) {
  for (j in 1:length(ROSMAP.AD.Binary)) {
    mat1 <- ROSMAP.AD.Binary[[i]][,-1]; mat1 <- matrix(data = sapply(mat1,as.numeric), ncol = ncol(mat1), nrow = nrow(mat1))
    mat2 <- ROSMAP.AD.Binary[[j]][,-1]; mat2 <- matrix(data = sapply(mat2,as.numeric), ncol = ncol(mat2), nrow = nrow(mat2))
    intersection <- sum(mat1 & mat2)
    union <- sum(mat1 | mat2)
    coefs[i,j] <- intersection / union
  }
}
colnames(coefs)<-c("DeSeq2", "Deseq_cov","FPKM","FPKM_cov","GeTMM","GeTMM_cov","TPM","TPM_cov", "edgeR","edgeR_cov") 
rownames(coefs)<-c("DeSeq2", "Deseq_cov","FPKM","FPKM_cov","GeTMM","GeTMM_cov","TPM","TPM_cov", "edgeR","edgeR_cov") 

pheatmap(coefs, color = hcl.colors(7, "YlOrBr", rev=T), angle_col = 45)


## ROSMAP rxn similarity##
sheet_names<-getSheetNames("ROSMAP_Summarised.xlsx")
coefs<-matrix(nrow=10, ncol=10)
for (i in 1:10){
  met1=read_xlsx("ROSMAP_Summarised.xlsx", sheet=i,col_names = TRUE) 
  for (j in 1:10){
    met2=read_xlsx("ROSMAP_Summarised.xlsx", sheet=j, col_names = TRUE)
    coefs[i,j] <-jaccard_index_rxn(met1$Reaction, met2$Reaction)
  }
}
colnames(coefs)<-c("DeSeq2", "Deseq_cov","FPKM","FPKM_cov","GeTMM","GeTMM_cov","TPM","TPM_cov", "edgeR","edgeR_cov") 
rownames(coefs)<-c("DeSeq2", "Deseq_cov","FPKM","FPKM_cov","GeTMM","GeTMM_cov","TPM","TPM_cov", "edgeR","edgeR_cov") 

pheatmap(coefs, color = hcl.colors(7, "YlOrBr", rev=T), angle_col = 45)


## ROSMAP pathway similarity##
sheet_names<-getSheetNames("ROSMAP_overrepresentation_hypergeo_new.xlsx")
coefs<-matrix(nrow=10, ncol=10)
for (i in 1:10){
  met1=read_xlsx("ROSMAP_overrepresentation_hypergeo_new.xlsx", sheet=i,col_names = TRUE) 
  for (j in 1:10){
    met2=read_xlsx("ROSMAP_overrepresentation_hypergeo_new.xlsx", sheet=j, col_names = TRUE)
    coefs[i,j] <-jaccard_index_rxn(met1$Pathway, met2$Pathway)
  }
}
colnames(coefs)<-c("DeSeq2", "Deseq_cov","FPKM","FPKM_cov","GeTMM","GeTMM_cov","TPM","TPM_cov", "edgeR","edgeR_cov") 
rownames(coefs)<-c("DeSeq2", "Deseq_cov","FPKM","FPKM_cov","GeTMM","GeTMM_cov","TPM","TPM_cov", "edgeR","edgeR_cov") 

pheatmap(coefs, color = hcl.colors(7, "YlOrBr", rev=T), angle_col = 45)  

## ROSMAP AD-related gene similarity##
setwd("./Benchmark_normalization")

data <-as.data.frame(read_xlsx("AD_BinarizedGeneMatches_v3.xlsx", sheet = 6,col_names = TRUE))
rownames(data)<-data$Gene
data_new<- data[-c(1,3,6,9,12,15)]

sheet_names<-colnames(data_new)
coefs<-matrix(nrow=10, ncol=10)
for (i in 1:10){
  for (j in 1:10){
    coefs[i,j] <- jaccard_similarity(data_new[,i], data_new[,j])
  }
}
colnames(coefs)<-c("DeSeq2", "Deseq_cov","FPKM","FPKM_cov","GeTMM","GeTMM_cov","TPM","TPM_cov", "edgeR","edgeR_cov") 
rownames(coefs)<-c("DeSeq2", "Deseq_cov","FPKM","FPKM_cov","GeTMM","GeTMM_cov","TPM","TPM_cov", "edgeR","edgeR_cov") 

pheatmap(coefs, color = hcl.colors(7, "YlOrBr", rev=T), angle_col = 45)  

write.csv(data.frame(coefs), "./ROSMAP_coefs_binarymatrix_personalized.csv")
rm(coefs)

## LUAD cancer-related gene similarity##
setwd("./Benchmark_normalization")
data <- read_xlsx("Cancer_BinarizedGeneMatches_v2.xlsx", sheet = 1,col_names = TRUE)
data_new <- data[-c(1)]
sheet_names<-colnames(data_new)
coefs<-matrix(nrow=5, ncol=5)
for (i in 1:5){
  for (j in 1:5){
    coefs[i,j] <-jaccard_similarity(data_new[,i], data_new[,j])
  }
}
colnames(coefs)<-c("DeSeq2", "edgeR", "FPKM", "GeTMM", "TPM")   
rownames(coefs)<-c("DeSeq2", "edgeR", "FPKM", "GeTMM", "TPM")  

write.csv(data.frame(coefs), "ROSMAP_coefs_binarymatrix_personalized.csv")
pheatmap(coefs, color = hcl.colors(7, "YlOrBr", rev=T), angle_col = 45)
rm(coefs)