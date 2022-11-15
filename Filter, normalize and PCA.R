if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
install.packages("devtools")
devtools::install_github("LUMC/dgeAnalysis")



library(edgeR)
library(limma)
dge <- DGEList(countfinal)
dge <- dge[rowSums(abs(dge$counts)) > 1,]

#Filter genes
edgeR <- calcNormFactors(dge, method = "TMM")
cpm_perc = 25
cpm_value = 1
counts <- cpm(edgeR, log = TRUE, prior.count = 1)
selectedFeatures <- rownames(edgeR)[apply(counts, 1, function(v)
  sum(v >= cpm_value)) >= (cpm_perc / 100) * ncol(counts)]
highExprDge <- dge[selectedFeatures, , keep.lib.sizes = FALSE]

#Normalize
normDge <- calcNormFactors(highExprDge, method = "TMM")
normDge$counts <- cpm(normDge, log = TRUE, prior.count = 1)

#PCA after filter
df = t(normDge)
pc <- prcomp(df, center=TRUE, scale=TRUE)
library(factoextra)
fviz_eig(pc, choice= c("variance", "eigenvalue"), ncp =13, addlabels=TRUE)
summary(pc)$importance

library(devtools)
library("dgeAnalysis")

nieuwe_pca <- pca_data(normDge)
