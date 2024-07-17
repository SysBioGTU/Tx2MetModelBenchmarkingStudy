

cat("\014")
rm(list=ls())

#%%%%% GENE-BASED ANALYSIS | Disease Enrichment analysis
# Libraries:
library(enrichR)
library(readxl)
library(rJava) 
library(dplyr)
library(openxlsx)


#----------------------------------------------
# Disease enrichment analysis for the gene lists involved in the regulated reactions
method <- "INIT" # Model reconstruction method: iMAT or INIT
setwd(paste0("H:/Benchmarking Study/ROSMAP/Gene_analyses/Differential_model_contents/", method))
sheet_name <- excel_sheets(path = "ROSMAP_sigGeneSymbol.xlsx")

# # To see all databases:
# listEnrichrDbs()

# Enrichment analysis step
db <- c('DisGeNET')

wb <- createWorkbook()
for(i in 1:length(sheet_name)){
  # Reading the genes in the differential reactions
  geneInfo <- read_excel("ROSMAP_sigGeneSymbol.xlsx", sheet = sheet_name[i])
  geneInfo <- as.character(geneInfo$'Gene Symbol')
  #................................................
  # Enrichr --> Enriched disease terms
  enrichedDiseases <- enrichr(geneInfo, databases = db[1])[[1]] %>% 
    dplyr::filter(.$Adjusted.P.value <0.05) %>% 
    dplyr::select(Term, Overlap,Adjusted.P.value, Genes)
  #................................................
  # Writing the enrichment results
  addWorksheet(wb, sheet_name[i])
  writeData(wb, sheet_name[i], enrichedDiseases)
}

saveWorkbook(wb, paste0("H:/Benchmarking Study/ROSMAP/Gene_analyses/Disease_enrichment/AD_Disease_Enrichment_DisGeNET_", method, ".xlsx"))


