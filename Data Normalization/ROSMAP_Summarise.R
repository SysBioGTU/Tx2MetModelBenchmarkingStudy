
setwd('./')

rm(list = ls())
gc()
library(tidyverse);library(tibble);library(openxlsx)
library(readxl);library(stringr)
Human_GEM <- read_excel("./Human-GEM.xlsx")
loc <- './ROSMAP/iMAT Reports/'
files <- list.files(paste0(loc), pattern = '*.xlsx')

# read and annotate pathways and metabolic equations
iMATResults <- lapply(files, function(x){
  df <- read.xlsx(paste0(loc,x))
  df$Pathway <- Human_GEM$SUBSYSTEM[match(df$Reaction,Human_GEM$ID)]
  df$Equation <- Human_GEM$EQUATION[match(df$Reaction,Human_GEM$ID)]
  df$GPR <- Human_GEM$`GENE ASSOCIATION`[match(df$Reaction,Human_GEM$ID)]
  return(df[df$FPValue<0.05,])
})
names(iMATResults) <- paste0(gsub('.xlsx|_MappedResult','',files))

# Save results
loc <- './ROSMAP/iMAT Summaries/'
openxlsx::write.xlsx(x = iMATResults, file = paste0(loc,'ROSMAP_Summarised.xlsx'), asTable = T)

# Count pathways
iMATResultsPathways <- lapply(iMATResults, function(x){
  df <- table(x$Pathway) |> as.data.frame()
  colnames(df) <- c('Pathway','Frequency')
  return(df[order(df$Frequency, decreasing = T),])
})
openxlsx::write.xlsx(x = iMATResultsPathways, file = paste0(loc,'ROSMAP_SummarisedPathways.xlsx'), asTable = T)
