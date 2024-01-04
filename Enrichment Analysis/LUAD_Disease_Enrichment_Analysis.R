

cat("\014")
rm(list=ls())

#%%%%% GENE-BASED ANALYSIS II
# remove.packages("xlsx")
# install.packages("rJava"); any(grepl("rJava",installed.packages())) 
# install.packages(pkgs = "xlsx")

#install.packages(pkgs = "xlsx")
library(enrichR)
library(readxl)
library(rJava) 
library(dplyr)

# detach("package:openxlsx", unload=TRUE)
# install.packages("openxlsx")
library(openxlsx)


#----------------------------------------------
# Disease Enrichment Analysis for the gene lists involved in the regulated reactions
setwd("./LUAD_ALL")
sheet_name <- excel_sheets(path = "TCGA_LUAD_GeneSymbol.xlsx")

# # To see all databases:
# listEnrichrDbs()
db <- c('DisGeNET', 'ClinVar_2019')

wb <- createWorkbook()
for(i in 1:length(sheet_name)){
  # Read the data
  geneInfo <- read_excel("TCGA_LUAD_GeneSymbol.xlsx", sheet = sheet_name[i])
  geneInfo <- as.character(geneInfo$'Gene Symbol')
  #................................................
  # Enrichr --> Disease Enrichment Analysis
  enrichedDiseases <- enrichr(geneInfo, databases = db[2])[[1]] %>% # db[1]: DisGeNET  &  db[2]: ClinVar_2019
    dplyr::filter(.$Adjusted.P.value <0.05) %>% 
    dplyr::select(Term, Overlap,Adjusted.P.value, Genes)
  #................................................
  # Write the results
  addWorksheet(wb, sheet_name[i])
  writeData(wb, sheet_name[i], enrichedDiseases)
}

#saveWorkbook(wb, 'H:/benchmarking_study/LUAD/+++Gene_analyses_muberra/Result_files/Disease_enrichment/TCGA_LUAD_Disease_Enrichment_DisGeNET.xlsx')
saveWorkbook(wb, './TCGA_LUAD_Disease_Enrichment_clinVar2019.xlsx')



