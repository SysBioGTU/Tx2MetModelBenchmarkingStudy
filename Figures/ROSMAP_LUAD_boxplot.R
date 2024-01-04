
# These codes can be used to generate boxplots for Figure 2 in the manuscript 

rm(list = ls())
gc()

library(tidyverse)
library(ggpubr)
library(ggbreak) 
library(colorBlindness)
library(hrbrthemes)
library(ggtext)
library(glue)
library(viridis)

setwd('./iMAT Results/')

filesR <- list.files(path = '.', pattern = '*.xlsx')
df <- lapply(filesR, function(x){
  sheetNames <- openxlsx::getSheetNames(paste0(x))
  df <- lapply(sheetNames, function(y){
    openxlsx::read.xlsx(xlsxFile = paste0(x), sheet = paste0(y))
  })
  names(df) <- paste0(sheetNames)
  return(df)
})
names(df) <- paste0(gsub('.xlsx','',filesR))

df2 <- lapply(df, function(x){
  do.call(rbind,lapply(x, function(y){nrow(y)}))|>
    as.data.frame()|>dplyr::mutate(Method = names(x))|>
    dplyr::rename(NumberOfReactions = 1)
})
# Binary data
# ROSMAP
gc()
filesR <- list.files(path = './ROSMAP_Binary/', pattern = '*.xlsx')
df <- lapply(filesR, function(x){
  sheetNames <- openxlsx::getSheetNames(paste0('./ROSMAP_Binary/',x))
  df1 <- openxlsx::read.xlsx(xlsxFile = paste0('./ROSMAP_Binary/',x), sheet = paste0(sheetNames[1]))|>
    dplyr::summarise(across(2:165, sum))
  colnames(df1) <- paste0('Control_',1:ncol(df1))
  df2 <- openxlsx::read.xlsx(xlsxFile = paste0('./ROSMAP_Binary/',x), sheet = paste0(sheetNames[2]))|>
    dplyr::summarise(across(2:403, sum))
  colnames(df2) <- paste0('AD_',1:ncol(df2))
  return(cbind(df1,df2))
})
names(df) <- paste0(gsub('.xlsx','',filesR))
gc()
df1 <- df[c(2,1,4,3,8,7,6,5,10,9)] 
df1 <- do.call(rbind,df1)|>
  as.data.frame()|>
  rownames_to_column(var = 'Method')|>
  reshape2::melt()|>
  dplyr::rename(GROUP = 2, NumberOfReactions = 3)|>
  dplyr::mutate(GROUP = case_when(grepl('^Control',GROUP, ignore.case=T)==T ~ 'CONTROL',
                                  grepl('^AD',GROUP, ignore.case=T)==T ~ 'AD',
                                  .default = GROUP),
                Method = gsub('logsuz','cov',Method),
                Method = gsub('_ROSMAP','',Method),
                Method = gsub('_normalized','',Method),
                Method = factor(Method, levels = unique(Method), ordered = TRUE))


ggplot(df1, aes(x=reorder(Method,GROUP), y = NumberOfReactions, fill = GROUP)) + 
  geom_boxplot()+
  xlab('')+
  ylab('Number of Reactions')+
  scale_fill_manual(values=c('#bdd7e7', '#6baed6'))+
  
  theme_ipsum_rc(base_family = "Roboto Condensed",
                 base_size = 12,
                 plot_title_family = "Roboto Condensed",
                 plot_title_size = 12,
                 plot_title_face = "bold",
                 plot_title_margin = 10,
                 subtitle_family = if (.Platform$OS.type == "windows") "Roboto Condensed" else
                   "Roboto Condensed Light",
                 subtitle_size = 13,
                 subtitle_face = "plain",
                 subtitle_margin = 15,
                 strip_text_family = "Roboto Condensed",
                 strip_text_size = 12,
                 strip_text_face = "plain",
                 caption_family = if (.Platform$OS.type == "windows") "Roboto Condensed" else
                   "Roboto Condensed Light",
                 caption_size = 9,
                 caption_face = "plain",
                 caption_margin = 10,
                 axis_text_size = 12,
                 axis_title_family = "Roboto Condensed",
                 axis_title_size = 15,
                 axis_title_face = "plain",
                 axis_title_just = "rt",
                 plot_margin = margin(30, 30, 30, 30),
                 panel_spacing = grid::unit(2, "lines"),
                 grid_col = "#cccccc",
                 grid = TRUE,
                 axis_col = "#cccccc",
                 axis = FALSE,
                 ticks = FALSE)  +
  theme(axis.title = element_text(vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(vjust = 0.5, hjust=0.5))+
  labs(title = "ROSMAP")+
  facet_grid(vars(GROUP))
ggsave(paste0('Figures/ROSMAP_ReactionBins.tiff'), 
       units="in", width=6.0, height=8.5, dpi=600, compression = 'lzw')


# LUAD
gc()
filesR <- list.files(path = './LUAD_Binary/', pattern = '*.xlsx')
df <- lapply(filesR, function(x){
  sheetNames <- openxlsx::getSheetNames(paste0('./LUAD_Binary/',x))
  df1 <- openxlsx::read.xlsx(xlsxFile = paste0('./LUAD_Binary/',x), sheet = paste0(sheetNames[1]))|>
    dplyr::summarise(across(2:59, sum))
  colnames(df1) <- paste0('Control_',1:ncol(df1))
  df2 <- openxlsx::read.xlsx(xlsxFile = paste0('./LUAD_Binary/',x), sheet = paste0(sheetNames[2]))|>
    dplyr::summarise(across(2:525, sum))
  colnames(df2) <- paste0('Cancer_',1:ncol(df2))
  return(cbind(df1,df2))
})
names(df) <- paste0(gsub('.xlsx','',filesR))
gc()
df1 <- df[c(1,2,4,3,5)]
df1 <- do.call(rbind,df1)|>
  as.data.frame()|>
  rownames_to_column(var = 'Method')|>
  reshape2::melt()|>
  dplyr::rename(GROUP = 2, NumberOfReactions = 3)|>
  dplyr::mutate(GROUP = case_when(grepl('^Control',GROUP, ignore.case=T)==T ~ 'CONTROL',
                                  grepl('^Cancer',GROUP, ignore.case=T)==T ~ 'CANCER',
                                  .default = GROUP))|>
  dplyr::mutate(GROUP = case_when(grepl('^Control',GROUP, ignore.case=T)==T ~ 'CONTROL',
                                  grepl('^Cancer',GROUP, ignore.case=T)==T ~ 'CANCER',
                                  .default = GROUP),
                Method = gsub('_normalized_LUAD','',Method),
                Method = factor(Method, levels = unique(Method), ordered = TRUE))


ggplot(df1, aes(x=reorder(Method,GROUP), y = NumberOfReactions, fill = GROUP)) + 
  geom_boxplot()+
  xlab('')+
  ylab('Number of Reactions')+
  scale_fill_manual(values=c('#bdd7e7', '#6baed6'))+
  
  theme_ipsum_rc(base_family = "Roboto Condensed",
                 base_size = 12,
                 plot_title_family = "Roboto Condensed",
                 plot_title_size = 12,
                 plot_title_face = "bold",
                 plot_title_margin = 10,
                 subtitle_family = if (.Platform$OS.type == "windows") "Roboto Condensed" else
                   "Roboto Condensed Light",
                 subtitle_size = 13,
                 subtitle_face = "plain",
                 subtitle_margin = 15,
                 strip_text_family = "Roboto Condensed",
                 strip_text_size = 12,
                 strip_text_face = "plain",
                 caption_family = if (.Platform$OS.type == "windows") "Roboto Condensed" else
                   "Roboto Condensed Light",
                 caption_size = 9,
                 caption_face = "plain",
                 caption_margin = 10,
                 axis_text_size = 12,
                 axis_title_family = "Roboto Condensed",
                 axis_title_size = 15,
                 axis_title_face = "plain",
                 axis_title_just = "rt",
                 plot_margin = margin(30, 30, 30, 30),
                 panel_spacing = grid::unit(2, "lines"),
                 grid_col = "#cccccc",
                 grid = TRUE,
                 axis_col = "#cccccc",
                 axis = FALSE,
                 ticks = FALSE)  +
  theme(axis.title = element_text(vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(vjust = 0.5, hjust=1))+
  labs(title = "LUAD")+
  facet_grid(vars(GROUP))
ggsave(paste0('Figures/LUAD_ReactionBins.tiff'), 
       units="in", width=6.0, height=8.5, dpi=600, compression = 'lzw')

