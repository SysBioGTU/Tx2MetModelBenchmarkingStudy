

rm(list = ls())
gc()
library(dplyr); library(tidyverse); library(tibble); library(tidyr);
library(edgeR); library(DESeq2); library(limma); 
load("./rosmapnewaalldata.RData")
protein_coding_rosmap <- openxlsx::read.xlsx("./protein_coding_rosmap.xlsx")
rosmapCount <- r1[,-c(2:4)]; rm(r1)
colnames(rosmapCount)[1] <- 'GeneID'
rosmapCount <- aggregate(.~GeneID,rosmapCount,max) 

dim(rosmapCount)
rosmapCount <- rosmapCount %>% 
  dplyr::select(-as.character(c('380_120503','367_120502','500_120515'))) |>
  remove_rownames() |> 
  column_to_rownames(var = 'GeneID')

dim(rosmapCount)

# Filter based on CPM  

ROSMAPNormalized <- list()
rawcounts<- matrix(data = sapply(rosmapCount,as.numeric), ncol = ncol(rosmapCount),nrow = nrow(rosmapCount))
colnames(rawcounts) <- colnames(rosmapCount)
rownames(rawcounts) <- rownames(rosmapCount)
#Calculate CPM values and apply cpm treshold
myCPM <- cpm(rawcounts,log=F)
thresh <- myCPM > 0.1 ## library size of ROSMAP ~50M. 5 count/ 50 = 0.1 cpm treshold
keep <- rowSums(thresh) >= 165 #Number of samples of the condition with the least number of samples  

summary(keep) 

rawcounts <-rawcounts[keep,] 
max(rawcounts); min(rawcounts)


# Access gene length information
exon_lengths_kb1 <- openxlsx::read.xlsx("./ExonLengths_GRCh38.v37.xlsx")
ex_length <- exon_lengths_kb1$exon_length[match(rownames(rawcounts),exon_lengths_kb1$X1)]
rawcounts.Len <- cbind(ex_length,rawcounts)|> as.data.frame()
colnames(rawcounts.Len)[1] <- 'Length'

#remove genes which have not got gene length information
rawcounts.Len <- rawcounts.Len[!is.na(rawcounts.Len$Length),]
dim(rawcounts.Len)

# metadata
metadatarosmap <- openxlsx::read.xlsx("./metadatarosmap.xlsx")
metadatarosmap <- metadatarosmap[which(metadatarosmap$specimenID %in% colnames(rawcounts.Len)),]

# Specify sample groups
group <- as.factor(metadatarosmap$condition)

AD <- metadatarosmap$specimenID[metadatarosmap$condition=='AD']
Controls <- metadatarosmap$specimenID[metadatarosmap$condition=='control']
# calculate RPK
df <- rawcounts.Len[,2:ncol(rawcounts.Len)]/rawcounts.Len[,1]

#GeTMM
getmm.norm <- DGEList(counts=df,group=group)
getmm.norm <- calcNormFactors(getmm.norm)
getmm.norm <- cpm(getmm.norm)
max(getmm.norm); min(getmm.norm)
getmm.norm <- getmm.norm |> as.data.frame() |> rownames_to_column(var = 'GeneID')

getmm.norm <- getmm.norm|> 
  dplyr::filter(GeneID %in% protein_coding_rosmap$ensembl_gene_id)
ROSMAPNormalized[['GeTMM_normalized']] <- getmm.norm
write.csv(getmm.norm,"./GeTMM_normalized_ROSMAP.csv")

df <- list()
df[['Control']] <- getmm.norm|> select(c(GeneID,Controls))
df[['AD']] <- getmm.norm|> select(c(GeneID,AD))
openxlsx::write.xlsx(df,'./GeTMM_normalized_ROSMAP.xlsx')

# EdgeR normalization
#Use raw counts instead of rpk values for edgeR normalization
df <- rawcounts |> 
  as.data.frame()

edger.norm <- DGEList(counts=df,group=group)
edger.norm <- calcNormFactors(edger.norm)
edger.norm <- cpm(edger.norm)
max(edger.norm); min(edger.norm)
edger.norm <- edger.norm |> 
  as.data.frame() |> 
  rownames_to_column(var = 'GeneID')

edger.norm <- edger.norm |> 
  dplyr::filter(GeneID %in% protein_coding_rosmap$ensembl_gene_id)
ROSMAPNormalized[['edgeR_normalized']] <- edger.norm
write.csv(edger.norm,"./edgeR_normalized_ROSMAP.csv")

df <- list()
df[['Control']] <- edger.norm|> select(c(GeneID,Controls))
df[['AD']] <- edger.norm|> select(c(GeneID,AD))
openxlsx::write.xlsx(df,'./edgeR_normalized_ROSMAP.xlsx')

# Deseq Normalization

df <- rawcounts |> as.data.frame()
coldata <- data.frame(row.names=colnames(df), group)

#construct DEseq dataset object
deseq2.norm <- DESeqDataSetFromMatrix(countData = df, colData=coldata, design=~group)
#Obtain normalized counts
deseq2.norm <- estimateSizeFactors(deseq2.norm)
deseq2.norm <- counts(deseq2.norm, normalized=TRUE) 
max(deseq2.norm); min(deseq2.norm)

deseq2.norm <- deseq2.norm |> as.data.frame() |> rownames_to_column(var = 'GeneID')

deseq2.norm <- deseq2.norm[deseq2.norm$GeneID %in% protein_coding_rosmap$ensembl_gene_id,]
ROSMAPNormalized[['DESEq2_normalized']] <- norm.counts.deseq
write.csv(norm.counts.deseq,"./DESEq2_normalized_ROSMAP.csv")

df <- list()
df[['Control']] <- deseq2.norm|> select(c(GeneID,Controls))
df[['AD']] <- deseq2.norm|> select(c(GeneID,AD))
openxlsx::write.xlsx(df,'./DESEq2_normalized_ROSMAP.xlsx')

# FPKM normalization
geneLen.ROSMAP <- exon_lengths_kb1$exon_length[match(rownames(rawcounts),exon_lengths_kb1$X1)]

df <- (rawcounts.Len |> select(-Length))+1
denCount <- rawcounts.Len$Length*apply(df,1,sum)
fpkm.norm <- (df/denCount)*10^6 
max(fpkm.norm); min(fpkm.norm)
fpkm.norm <- fpkm.norm|> 
  as.data.frame() |> 
  rownames_to_column(var = 'GeneID')

fpkm.norm <- fpkm.norm[fpkm.norm$GeneID %in% protein_coding_rosmap$ensembl_gene_id,]
ROSMAPNormalized[['FPKM_normalized']] <- fpkm.norm
write.csv(fpkm.norm,"./FPKM_normalized_ROSMAP.csv")

df <- list()
df[['Control']] <- fpkm.norm |> select(c(GeneID,Controls))
df[['AD']] <- fpkm.norm |> select(c(GeneID,AD))
openxlsx::write.xlsx(df,'./FPKM_normalized_ROSMAP.xlsx')

# TPM normalization
df <- rawcounts.Len |> select(-Length); df <- df+1

A.Matrix <- df/as.numeric(geneLen.ROSMAP)
tpm.norm <- (A.Matrix/apply(A.Matrix,1,sum))*10^6
max(tpm.norm); min(counts2TPM1)
tpm.norm <- tpm.norm |> 
  as.data.frame() |> 
  rownames_to_column(var = 'GeneID')

tpm.norm <- tpm.norm |> 
  dplyr::filter(GeneID%in% protein_coding_rosmap$ensembl_gene_id)
ROSMAPNormalized[['TPM_normalized']] <- tpm.norm
write.csv(tpm.norm,"./TPM_normalized_ROSMAP.csv")
summary(tpm.norm)

df <- list()
df[['Control']] <- tpm.norm |> select(c(GeneID,Controls))
df[['AD']] <- tpm.norm |> select(c(GeneID,AD))
openxlsx::write.xlsx(df,'./TPM_normalized_ROSMAP.xlsx')

saveRDS(ROSMAPNormalized,'./ROSMAPNormalized.RDS')

# Covariate adjustment - sex, age, pmi
gc()
covAdjustment <- function(exprData, metaData) {
  metaData <- metadatarosmap |> dplyr::filter(specimenID %in% colnames(exprData))
  condition.num <-as.numeric(as.factor(metaData$condition))
  sex.num <- as.numeric(as.factor(metaData$msex))
  age <- as.data.frame(metaData$age_death)
  age[age=="90+"] <- "91" 
  age <- as.numeric(age$`metaData$age_death`)
  pmi <- as.numeric(metaData$pmi)
  coef <- lm(t(exprData) ~ sex.num + age + condition.num + pmi)$coefficients
  logdata2 <- exprData
  for(i in 1:ncol(logdata2)){
    logdata2[,i] <- logdata2[,i]- (coef[2,]*sex.num[i]+coef[3,]*age[i]+coef[4,]*pmi[i])
  }
  return(logdata2)
}

transReportCovariateAdjusted <- lapply(seq_along(ROSMAPNormalized), function(x){
  df <- ROSMAPNormalized[[x]] |> 
    remove_rownames() |> 
    column_to_rownames(var ='GeneID')
  df <- covAdjustment(df, metaData = metadatarosmap)
  df2 <- list()
  df2[['Control']] <- df |> select(Controls) |> rownames_to_column(var = 'GeneID')
  df2[['AD']] <- df |> select(AD) |> rownames_to_column(var = 'GeneID')
  openxlsx::write.xlsx(df2,paste0('./',names(ROSMAPNormalized)[x],'_logsuz.xlsx'))
})
library(plyr)

# Summary statistics
df <- join_all(list(
  Control = openxlsx::read.xlsx("./DESEq2_normalized_loglu.xlsx", sheet = 1),
  AD = openxlsx::read.xlsx("./DESEq2_normalized_loglu.xlsx", sheet = 2)
),
              by = 'GeneID') |>
  as.data.frame()|>
  remove_rownames()|>
  column_to_rownames(var = 'GeneID')|>
  as.matrix()|>
  t()
boxplot(df)
saveRDS(list(logSuz=logSuz,logLu=logLu), file = './ROSMAP_CovAdj.RDS')
gc()
