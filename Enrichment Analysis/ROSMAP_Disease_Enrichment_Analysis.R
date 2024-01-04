

cat("\014")
rm(list=ls())

#%%%%% GENE-BASED ANALYSIS II
library(enrichR)
library(readxl)
library(rJava) 
library(dplyr)
library(openxlsx)

#----------------------------------------------
# Disease Enrichment Analysis for the gene lists involved in the regulated reactions
setwd("./ROSMAP_All")
sheet_name <- excel_sheets(path = "ROSMAP_GeneSymbol.xlsx")

# # To see all databases:
# listEnrichrDbs()
db <- c('DisGeNET', 'ClinVar_2019')

wb <- createWorkbook()
for(i in 1:length(sheet_name)){
  # Read the data
  geneInfo <- read_excel("ROSMAP_GeneSymbol.xlsx", sheet = sheet_name[i])
  geneInfo <- as.character(geneInfo$'Gene Symbol')
  #................................................
  # Enrichr --> Disease Enrichment Analysis
  enrichedDiseases <- enrichr(geneInfo, databases = db[1])[[1]] %>% # db[1]: DisGeNET  &  db[2]: ClinVar_2019
    dplyr::filter(.$Adjusted.P.value <0.05) %>% 
    dplyr::select(Term, Overlap,Adjusted.P.value, Genes)
  #................................................
  # Write the results
  addWorksheet(wb, sheet_name[i])
  writeData(wb, sheet_name[i], enrichedDiseases)
}

saveWorkbook(wb, './Result_files/Disease_Enrichment/ROSMAP_AD_Disease_Enrichment_DisGeNET.xlsx')