

cat("\014") #clear the console (CTRL+L)
rm(list=ls()) #remove all objects(variables)

# Libraries
library(ggplot2)
library(tidyverse) 
library(reshape2)
#...............................
#Font type
library(extrafont) #text font
loadfonts(device = "win")


# GENE-BASED ANALYSIS | Plotting Enriched Disease Terms

# ** Loading the AD-/LUAD-associated enriched terms (iMAT)
# iMAT models --> Genes in the differential reactions --> Enrichr (DisGeNET)

# X-axis: normalization method, Y-axis: disease term, size: -log10(p-value)
setwd("H:/Benchmarking Study/ROSMAP/Gene_analyses/Disease_enrichment")
disease_res <- read.csv(file = 'AD_LUAD_enrichedPaths_selected_final.csv', header = TRUE, sep = ";")

disease_Long <- gather(disease_res, key="Mutations", value="Size",
                       c('RLE',	'RLEcov',	'TMM',	'TMMcov',	'GeTMM',	'GeTMMcov',	'FPKM',	'FPKMcov',	'TPM',	'TPMcov'))

# Visualization of the enrichment results
custom_colors <- rep(c("#db8e1e"), 10)

plt <- ggplot(disease_Long,aes(x=Mutations, y=Pathways, color = Mutations)) +
  geom_point(aes(size=Size)) +  scale_y_discrete(limits = rev) +
  labs(x=NULL, y=NULL) +
  theme(text=element_text(family="Calibri", color="black", size=17.5),
        axis.text.x=element_text(family="Calibri", face = "plain", color="black", size=17.5),
        axis.text.y=element_text(family="Calibri", face = "plain", color="black", size=17.5)) + 
  facet_grid(Source ~ ., scales = "free_y")

plt + scale_color_manual(values = custom_colors)


# #------------------------------------------------------------
# ** Loading the AD-associated enriched terms (INIT)
# INIT models --> Genes in the differential reactions --> Enrichr (DisGeNET)

# X-axis: normalization method, Y-axis: disease term, size: -log10(p-value)
setwd("H:/Benchmarking Study/ROSMAP/Gene_analyses/Disease_enrichment")
disease_res <- read.csv(file = 'AD_INIT_enrichedPaths_selected_final.csv', header = TRUE, sep = ";")

disease_Long <- gather(disease_res, key="Mutations", value="Size",
                       c('RLEcov',	'TMMcov',	'GeTMMcov',	'FPKMcov',	'TPMcov'))

# Visualization of the enrichment results
custom_colors <- rep(c("#db8e1e"), 10)

plt <- ggplot(disease_Long,aes(x=Mutations, y=Pathways, color = Mutations)) +
  geom_point(aes(size=Size)) +  scale_y_discrete(limits = rev) +
  labs(x=NULL, y=NULL) + theme_bw() +
  theme(text=element_text(family="Calibri", color="black", size=17.5),
        axis.text.x=element_text(family="Calibri", face = "plain", color="black", size=17.5),
        axis.text.y=element_text(family="Calibri", face = "plain", color="black", size=17.5)) + 
  facet_grid(Source ~ ., scales = "free_y")

plt + scale_color_manual(values = custom_colors)

