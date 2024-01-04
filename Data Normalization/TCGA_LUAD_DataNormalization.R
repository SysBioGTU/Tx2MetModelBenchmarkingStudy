

rm(list = ls())
gc()
library(readr); library(openxlsx); library(tidyverse)
protein_coding_rosmap <- openxlsx::read.xlsx("H:/benchmarking_study/ROSMAP/Repeat 26.05.2023/protein_coding_rosmap.xlsx")
files <- list.files(path = 'H:/benchmarking_study/LUAD/Gene Expression Data/', pattern = '*.csv')
fileLoc <- 'H:/benchmarking_study/LUAD/Gene Expression Data/'

dataFiles <- lapply(files, function(x){
  df <- readr::read_delim(paste0(fileLoc,x),delim = ";", escape_double = FALSE, trim_ws = TRUE)
  df <- df[df$ensembl_id %in% protein_coding_rosmap$ensembl_gene_id,]
  return(df)
})
names(dataFiles) <- paste0(gsub('.csv','',files))

# LUAD metadata
library(readr)
TCGA_LUAD_Metadata <- read_csv("H:/benchmarking_study/LUAD/LUAD - Repeat/TCGA_LUAD_Metadata.csv")
TCGA_LUAD_Metadata$LUAD.data.barcode <- gsub('\\-','.',TCGA_LUAD_Metadata$LUAD.data.barcode)

Controls <- TCGA_LUAD_Metadata$LUAD.data.barcode[which(TCGA_LUAD_Metadata$LUAD.data.definition=='Solid Tissue Normal' & 
                                       TCGA_LUAD_Metadata$LUAD.data.barcode %in% colnames(dataFiles$Deseq2_normalized_LUAD))]
Cancer <- TCGA_LUAD_Metadata$LUAD.data.barcode[which(TCGA_LUAD_Metadata$LUAD.data.definition=='Primary solid Tumor' & 
                                     TCGA_LUAD_Metadata$LUAD.data.barcode %in% colnames(dataFiles$Deseq2_normalized_LUAD))]

# save LUAD data to Excel file
for(i in 1:length(dataFiles)){
  df <- list()
  df[['Control']] <- dataFiles[[i]] |> as.data.frame() |> dplyr::select(ensembl_id,Controls)
  df[['Cancer']] <- dataFiles[[i]] |> as.data.frame() |> dplyr::select(ensembl_id,Cancer)
  openxlsx::write.xlsx(df, file = paste0(saveLoc,names(dataFiles)[i],'.xlsx'))
}
gc()