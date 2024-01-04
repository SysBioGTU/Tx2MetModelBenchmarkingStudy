#You can use tihs code to generate radar charts for gene-based analysis in the manuscript (Figure 3).

#install.packages("radarchart")
library(radarchart)  

# *** Radar charts for the regulated reactions & the related genes:
setwd("./Benchmarking_figures")

# Normalized ROSMAP
color <- grDevices::col2rgb(c("#005b96", "#cb3234"))

labels = c("FPKM", "TPM", "edgeR", "DESeq2", "GeTMM")
data1 <- data.frame("Reactions" = c(1621, 3083, 397, 490, 393))
data2 <- data.frame("Genes" = c(1045, 1327, 269, 256, 353))
chartJSRadar(scores=c(data1, data2), labs=labels, polyAlpha=0.2, showLegend=F, colMatrix = color,
             width = '100%', height = '400px', labelSize = 18, addDots = TRUE, lineAlpha = 1)

#_________________________________________
# Normalized + Covariate Adjusted ROSMAP
labels = c("FPKM", "TPM", "edgeR", "DESeq2", "GeTMM")
data1 <- data.frame("Reactions" = c(703, 685, 262, 301, 343))
data2 <- data.frame("Genes" = c(371, 458, 186, 204, 200))
chartJSRadar(scores=c(data1, data2), labs=labels, polyAlpha=0.2, showLegend=F, colMatrix = color,
             width = '100%', height = '400px', labelSize = 18, addDots = TRUE, lineAlpha = 1)


#.................................................................
# LUAD
labels = c("FPKM", "TPM", "edgeR", "DESeq2", "GeTMM")
data1 <- data.frame("Reactions" = c(2555, 3115, 1338, 1326, 1498))
data2 <- data.frame("Genes" = c(1124, 1245, 807, 780, 860))
chartJSRadar(scores=c(data1, data2), labs=labels, polyAlpha=0.2, showLegend=F, colMatrix = color,
             width = '100%', height = '400px', labelSize = 18, addDots = TRUE, lineAlpha = 1)


