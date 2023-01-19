#Pre-processing the data, filtering out the low expressing genes and making pc's of the data


#Load required packages
library(dplyr)
library(tibble)
library(edgeR)
library(limma)
library(factoextra)

#Packages suggested by Casper, needed to create the PCs.
#install.packages("devtools")
library(devtools)
#devtools::install_github("LUMC/dgeAnalysis")
library("dgeAnalysis")


names(counts) <- counts[1,]
clinical <- subset(clinical, select = -c(X, OTHER_PATIENT_ID))
clinical <- clinical[order(clinical$PATIENT_ID),]

countno <- counts[,-1]
rownames(countno) <- counts[,1]
countfinal <- countno[-1, ]
colnames(countfinal) <- countno[1, ]
countfinal <- subset(countfinal, select=-c(gene_name,gene_type))

rownames(clinical)<- clinical[,1]
all(rownames(clinical) %in% colnames(countfinal))
countfinal <- countfinal[,rownames(clinical)]
all(rownames(clinical) == colnames(countfinal))

ncol(countfinal) == nrow(clinical)

countfinal <- mutate_all(countfinal, function(x) as.numeric(as.character(x)))


df <- tibble::rownames_to_column(countfinal, "VALUE")

#Create DGE object

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
normDge_counts <- normDge$counts
normDge_counts <- normDge[-c(0,1,2,3,4),]
normDge_counts <- normDge_counts[["counts"]]


#PCA after filter
df = t(normDge_counts)
pc <- prcomp(df, center=TRUE, scale=TRUE)

#Making a plot of the pcs with the highest variances
fviz_eig(pc, choice= c("variance", "eigenvalue"), ncp =30, addlabels=TRUE)
summary(pc)$importance

#Make pcs of the data

nieuwe_pca <- pca_data(normDge)

nieuwe_pca <- data.frame(nieuwe_pca[,-1], row.names=nieuwe_pca[,1])

#Select the amount of PCs with the highest variance 
pca_data <- nieuwe_pca[, c(0:7)]
