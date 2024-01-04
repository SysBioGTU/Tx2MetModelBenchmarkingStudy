library(readxl)
library(tidyr)
library(vroom)
library(dplyr)
library(openxlsx)

## Read Human Genome Scale Metabolic Model
humangem <- openxlsx::read.xlsx("Human-GEM.xlsx") 
## The number of reactions belonging to each subsystem is calculated
groups_global<-humangem %>% group_by(SUBSYSTEM) %>% summarise(Freq=n())

########### ROSMAP ###########

## Read differentially altered pathway list. Each sheet in Excel shows a different normalization method
sheet_names <- getSheetNames("ROSMAP_SummarisedPathways_new.xlsx")

## Run overrepresentation analysis in a for loop for each normalization method seperately
rxnsize<-data.frame()
for (j in 1:length(sheet_names)){
  ## J is to select different sheets (normalization methods) 
  paths <- openxlsx::read.xlsx("ROSMAP_SummarisedPathways_new.xlsx",sheet =j )
  ## Each sheet includes two colums. First is pathway names and second is number of differentially changed reaction numbers that belong to this pathway.
  ## Find total number of differentially changed reactions
  rxnsize[1,j]<- sum(paths[,2])
  
  p_value_hypergeo<-NULL
  expected<-NULL
  sign<-NULL
  ## Calculate overrepresantation or underrrepresentation p-value for each pathway
  for (i in 1:nrow(paths)){
    ## Calculate expected value first
    expected[i]<-as.numeric(groups_global[which(as.data.frame(groups_global[,1])==as.character(paths[i,1])),2]/sum(groups_global[,2])*sum(paths[,2]))
    if( expected[i]>paths[i,2]){
      ## if the number of differentially altered reaction number is smaller than the expected value calculate underrepresentation p-value using hypergeometric test
      p_value_hypergeo[i]<-phyper(as.numeric(paths[i,2]), as.numeric(groups_global[which(as.data.frame(groups_global[,1])==as.character(paths[i,1])),2]), sum(groups_global[,2])-as.numeric(groups_global[which(as.data.frame(groups_global[,1])==as.character(paths[i,1])),2]),sum(paths[,2]), lower.tail = FALSE, log.p = FALSE)
      sign[i]<-"-"
    }
    if( expected[i]<paths[i,2]){
      ## if the number of differentially altered reaction number is greater than the expected value calculate overrepresentation p-value using hypergeometric test
      p_value_hypergeo[i]<-phyper(as.numeric(paths[i,2])-1, as.numeric(groups_global[which(as.data.frame(groups_global[,1])==as.character(paths[i,1])),2]), sum(groups_global[,2])-as.numeric(groups_global[which(as.data.frame(groups_global[,1])==as.character(paths[i,1])),2]),sum(paths[,2]), lower.tail = FALSE, log.p = FALSE)
       sign[i]<-"+"
    }}
  
  ##Create result data frame 
  hyper <- data.frame(paths,expected,sign,p_value_hypergeo)%>% 
    filter (p_value_hypergeo < 0.05)%>%  ##Filter according to the p-value<0.05 cut-off to select pathways that change significantly
    filter (sign == "+") %>% ## Select only overrepresentation
    arrange(p_value_hypergeo) ## Sort by p-value
  ## Write results to Excel
  openxlsx::write.xlsx(hyper, "ROSMAP_overrepresentation_hypergeo_new.xlsx", sheetName = sheet_names[j],append = T, row.names = F)
  
  }

############ LUAD ############

## Read differentially altered pathway list. Each sheet in Excel shows a different normalization method
sheet_names<-getSheetNames("TCGA_LUAD_SummarisedPathways.xlsx")

## Run overrepresentation analysis in a for loop for each normalization method seperately
rxnsize<-data.frame()
for (j in 1:length(sheet_names)){
  ## J is to select different sheets (normalization methods) 
  paths <- openxlsx::read.xlsx("TCGA_LUAD_SummarisedPathways.xlsx",sheet =j )
  ## Each sheet includes two colums. First is pathway names and second is number of differentially changed reaction numbers that belong to this pathway.
  ## Find total number of differentially changed reactions
  rxnsize[1,j]<- sum(paths[,2])
  
  p_value_hypergeo<-NULL
  expected<-NULL
  sign<-NULL
  ## Calculate overrepresantation or underrrepresentation p-value for each pathway
  for (i in 1:nrow(paths)){
    ## Calculate expected value first
    expected[i]<-as.numeric(groups_global[which(as.data.frame(groups_global[,1])==as.character(paths[i,1])),2]/sum(groups_global[,2])*sum(paths[,2]))
    if( expected[i]>paths[i,2]){
      ## if the number of differentially altered reaction number is smaller than the expected value calculate underrepresentation p-value using hypergeometric test
      p_value_hypergeo[i]<-phyper(as.numeric(paths[i,2]), as.numeric(groups_global[which(as.data.frame(groups_global[,1])==as.character(paths[i,1])),2]), sum(groups_global[,2])-as.numeric(groups_global[which(as.data.frame(groups_global[,1])==as.character(paths[i,1])),2]),sum(paths[,2]), lower.tail = FALSE, log.p = FALSE)
      sign[i]<-"-"
    }
    if( expected[i]<paths[i,2]){
      ## if the number of differentially altered reaction number is greater than the expected value calculate overrepresentation p-value using hypergeometric test
      p_value_hypergeo[i]<-phyper(as.numeric(paths[i,2])-1, as.numeric(groups_global[which(as.data.frame(groups_global[,1])==as.character(paths[i,1])),2]), sum(groups_global[,2])-as.numeric(groups_global[which(as.data.frame(groups_global[,1])==as.character(paths[i,1])),2]),sum(paths[,2]), lower.tail = FALSE, log.p = FALSE)
      sign[i]<-"+"
    }}
  
  ##Create result data frame 
  hyper <- data.frame(paths,expected,sign,p_value_hypergeo)%>% 
    filter (p_value_hypergeo < 0.05)%>%  ##Filter according to the p-value<0.05 cut-off to select pathways that change significantly
    filter (sign == "+") %>% ## Select only overrepresentation
    arrange(p_value_hypergeo) ## Sort by p-value
  ## Write results to Excel
  xlsx::write.xlsx(hyper, "LUAD_overrepresentation_hypergeo.xlsx", sheetName = sheet_names[j],append = T, row.names = F)
  
}