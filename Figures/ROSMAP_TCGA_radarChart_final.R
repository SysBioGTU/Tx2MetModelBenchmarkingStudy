
#install.packages("radarchart")
library(radarchart)  

# *** RADAR CHARTS
# Plots for the numbers of differential reactions & the related genes:
labels = c("FPKM", "TPM", "TMM", "RLE", "GeTMM")
color <- grDevices::col2rgb(c("#005b96", "#da680f"))

# Normalized AD (ROSMAP) - iMAT
data1 <- data.frame("Reactions" = c(1621, 3083, 397, 490, 393))
data2 <- data.frame("Genes" = c(1045, 1327, 269, 256, 353))
chartJSRadar(scores=c(data1, data2), labs=labels, polyAlpha=0.2, showLegend=F, colMatrix = color,
             width = '100%', height = '400px', labelSize = 18, addDots = TRUE, lineAlpha = 1)

# Normalized + Covariate-Adjusted AD (ROSMAP) - iMAT
data1 <- data.frame("Reactions" = c(703, 685, 262, 301, 343))
data2 <- data.frame("Genes" = c(371, 458, 186, 204, 200))
chartJSRadar(scores=c(data1, data2), labs=labels, polyAlpha=0.2, showLegend=F, colMatrix = color,
             width = '100%', height = '400px', labelSize = 18, addDots = TRUE, lineAlpha = 1)

# Normalized + Covariate-Adjusted AD (ROSMAP) - INIT
data1 <- data.frame("Reactions" = c(1498, 2902, 1311, 1455, 1310))
data2 <- data.frame("Genes" = c(975, 1365, 606, 729, 621))
chartJSRadar(scores=c(data1, data2), labs=labels, polyAlpha=0.2, showLegend=F, colMatrix = color,
             width = '100%', height = '400px', labelSize = 18, addDots = TRUE, lineAlpha = 1)


#.................................................................
# Normalized LUAD (TCGA) - iMAT
labels = c("FPKM", "TPM", "TMM", "RLE", "GeTMM")
data1 <- data.frame("Reactions" = c(2555, 3115, 1338, 1326, 1498))
data2 <- data.frame("Genes" = c(1124, 1245, 807, 780, 860))
chartJSRadar(scores=c(data1, data2), labs=labels, polyAlpha=0.2, showLegend=F, colMatrix = color,
             width = '100%', height = '400px', labelSize = 18, addDots = TRUE, lineAlpha = 1)

# Normalized + Covariate-Adjusted LUAD - iMAT
labels = c("FPKM", "TPM", "TMM", "RLE", "GeTMM")
data1 <- data.frame("Reactions" = c(2799, 2831, 2373, 2216, 2254))
data2 <- data.frame("Genes" = c(1124, 1101, 829, 867, 801))
chartJSRadar(scores=c(data1, data2), labs=labels, polyAlpha=0.2, showLegend=F, colMatrix = color,
             width = '100%', height = '400px', labelSize = 18, addDots = TRUE, lineAlpha = 1)


