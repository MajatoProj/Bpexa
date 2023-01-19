# Title: Full_code of the Bpexa clustering project
# Authors: Jeremy den Ouden, Marloes van Luik, David Labee and Samuel Rebel
# Commissioner:  Hailiang (Leon) Mei
# Date: 19-01-2023

#Required packages
library("ggfortify")
library(dplyr)
library(plyr)
library(ggplot2)
library(umap)
library(magrittr)
library(factoextra)
library(dbscan)
library(fpc)
library(dplyr)
library(tibble)
library(edgeR)
library(limma)
library(devtools)
library("dgeAnalysis")
library(mclust)
library(cluster)

#Data inlezen
counts <- read.csv("Primary_Tumors_unstranded_counts.csv", header = FALSE, stringsAsFactors=FALSE)
clinical <- read.csv("clinical3.csv", comment.char="#", sep = ";")

#Make a subset of clinical
#clinical_age <- data.frame(AGE=c(90,89,88,87,86,85,84,83,82,81))
#clinical_age <- data.frame(AGE=c(80,79,78,77,76,75,74,73,72,71))
#clinical_age <- data.frame(AGE=c(70,69,68,67,66,65,64,63,62,61))
#clinical_age <- data.frame(AGE=c(60,59,58,57,56,55,54,53,52,51))
#clinical_age <- data.frame(AGE=c(50,49,48,47,46,45,44,43,42,41))
clinical_age <- data.frame(AGE=c(40,39,38,37,36,35,34,33,32,31,30,29,28,27,26))

#Select the wanted subtypes
clinical_subtype <- data.frame(SUBTYPE=c("BRCA_Basal", "BRCA_LumA", "BRCA_LumB", "BRCA_Her2", "BRCA_Normal" ))
clinical_er <- data.frame(ER_STATUS=c("Positive", "Negative" ))
clinical <- match_df(clinical, clinical_age)
clinical <- match_df(clinical, clinical_subtype)
clinical <- match_df(clinical, clinical_er)

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

#hclust
hcl <- eclust(pca_data, "hclust", k=2)
fviz_silhouette(hcl)

#cluster hcl
cluster_data <- data.frame(hcl$cluster)
clinical_cluster <- cbind(clinical, cluster_data)

clinical_cluster$hcl.cluster <- as.character(clinical_cluster$hcl.cluster)


##umap
umap_results <- umap::umap(df, intgroup=c("hcl.cluster"), returnData=TRUE)
umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("PATIENT_ID") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(clinical_cluster, by = "PATIENT_ID")


ggplot(umap_plot_df,aes(x = X1,y = X2,
                        color = hcl.cluster,                        
                        shape = ER_STATUS)) + geom_point(size=2) +   ggtitle("Age <- 40 colored on ER+/ER-")

#PAM clustering
#bepalen optimaal aantal clusters voor pam
fviz_nbclust(pca_data, FUNcluster=cluster::pam, k.max = 5)
fviz_nbclust(pca_data, FUNcluster=cluster::pam, method="gap_stat", k.max  = 7)+ theme_classic()

#clusteren pam
pam1<-eclust(pca_data, "pam", k=2)
fviz_silhouette(pam1)

#cluster pam1
cluster_data <- data.frame(pam1$cluster)
clinical_cluster <- cbind(clinical, cluster_data)

clinical_cluster$pam1.cluster <- as.character(clinical_cluster$pam1.cluster)


##umap
umap_results <- umap::umap(df, intgroup=c("pam1.cluster"), returnData=TRUE)
umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("PATIENT_ID") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(clinical_cluster, by = "PATIENT_ID")


ggplot(umap_plot_df,aes(x = X1,y = X2,
                        color = pam1.cluster,                        
                        shape = ER_STATUS)) + geom_point(size = 2) +   ggtitle("Age <- 40 colored on ER+/ER-")

#Clustering K-means and visualizing with UMAP

pca_data_f <- na.omit(pca_data)
pca_data_f[Reduce(`&`, lapply(pca_data_f, function(x) !is.na(x)  & is.finite(x))),]
fviz_nbclust(pca_data_f, FUNcluster=kmeans, k.max = 8) 
fviz_nbclust(pca_data_f, FUNcluster=kmeans, method="gap_stat", k.max = 8)+ theme_classic()

#clusteren met kmeans
km1<-eclust(pca_data, "kmeans", hc_metric="eucliden",k=2)
fviz_silhouette(km1)

cluster_data <- data.frame(km1$cluster)
clinical_cluster <- cbind(clinical, cluster_data)

clinical_cluster$km1.cluster <- as.character(clinical_cluster$km1.cluster)


##umap
umap_results <- umap::umap(df, intgroup=c("km1.cluster"), returnData=TRUE)
umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("PATIENT_ID") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(clinical_cluster, by = "PATIENT_ID")


ggplot(umap_plot_df,aes(x = X1,y = X2,
                        color = km1.cluster,                        
                        shape = ER_STATUS)) + geom_point(size = 2) +   ggtitle("Age <- 40 colored on ER+/ER-")

#diana 
dia <- eclust(pca_data, "diana", k=2)
fviz_silhouette(dia)

#cluster diana
cluster_data <- data.frame(dia$cluster)
clinical_cluster <- cbind(clinical, cluster_data)

clinical_cluster$dia.cluster <- as.character(clinical_cluster$dia.cluster)

##umap
umap_results <- umap::umap(df, intgroup=c("dia.cluster"), returnData=TRUE)
umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("PATIENT_ID") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(clinical_cluster, by = "PATIENT_ID")

ggplot(umap_plot_df,aes(x = X1,y = X2,
                        color = dia.cluster,                        
                        shape = ER_STATUS)) + geom_point(size=2) +   ggtitle("Age <- 40 colored on ER+/ER-")

### Clustering scoring metrics pca
fviz_nbclust(pca_data, kmeans, method = 'wss')
fviz_nbclust(pca_data, kmeans, method = 'silhouette')

### Checking optimal eps value
kNNdistplot(pca_data, k=2)
abline(h=60, lty=2)

### Clustering
set.seed(123)
#f <- fpc::dbscan(pca_data, eps=46, MinPts = 4)
d <- dbscan::dbscan(pca_data, eps = 65, minPts = 3)
fviz_cluster(d, pca_data, geom = "points")

### Obtain quantitavtie score silhouette
D <- daisy(pca_data)
plot(silhouette(d$cluster, D), col=1:5)


score <- silhouette(d$cluster, D)
score <- data.frame(score)
sum(score$sil_width) / length(score$sil_width)


cluster_data <- data.frame(d$cluster)
clinical_cluster <- cbind(clinical, cluster_data)

clinical_cluster$d.cluster <- as.character(clinical_cluster$d.cluster)

##umap
umap_results <- umap::umap(df, intgroup=c("d.cluster"), returnData=TRUE)
umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("PATIENT_ID") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(clinical_cluster, by = "PATIENT_ID")

ggplot(umap_plot_df,aes(x = X1,y = X2,
                        color = d.cluster,                        
                        shape = ER_STATUS)) + geom_point(size = 2) + ggtitle("Age <- 40 colored on ER+/ER-")

#Adjusted rand index

#pam
if (any(pam1$cluster == 1)) {
  # Change the value to 0
  pam1$adjrand[pam1$cluster == 1] <- 'Positive'
}
if (any(pam1$cluster == 2)) {
  # Change the value to 0
  pam1$adjrand[pam1$cluster == 2] <- 'Negative'
}
#calculate the adjusted rand index
adjusted_pam <- adjustedRandIndex(clinical$ER_STATUS, pam1$adjrand)


#kmeans
if (any(km1$cluster == 1)) {
  # Change the value to 0
  km1$adjrand[km1$cluster == 1] <- 'Positive'
}
if (any(km1$cluster == 2)) {
  # Change the value to 0
  km1$adjrand[km1$cluster == 2] <- 'Negative'
}
#calculate the adjusted rand index
adjusted_kmeans <- adjustedRandIndex(clinical$ER_STATUS, km1$adjrand)


#Diana
if (any(dia$cluster == 1)) {
  # Change the value to 0
  dia$adjrand[dia$cluster == 1] <- 'Positive'
}
if (any(dia$cluster == 2)) {
  # Change the value to 0
  dia$adjrand[dia$cluster == 2] <- 'Negative'
}
#calculate the adjusted rand index
adjusted_diana <- adjustedRandIndex(clinical$ER_STATUS, dia$adjrand)


#Dbscan is difficult to cluster, sometimes it gives more than 2 clusters. 
#If that is the case, look at umap and decide if that cluster is supposed to be ER+ or ER-.
#dbscan
if (any(d$cluster == 0)) {
  # Change the value to 0
  d$adjrand[d$cluster == 0] <- 'Positive'
}
if (any(d$cluster == 1)) {
  # Change the value to 0
  d$adjrand[d$cluster == 1] <- 'Positive'
}
if (any(d$cluster == 2)) {
  # Change the value to 0
  d$adjrand[d$cluster == 2] <- 'Negative'
}
if (any(d$cluster == 3)) {
  # Change the value to 0
  d$adjrand[d$cluster == 3] <- 'Positive'
}
#calculate the adjusted rand index
adjusted_dbscan <- adjustedRandIndex(clinical$ER_STATUS, d$adjrand)

# Make the boxplot of the Adjusted Rand Index table

scores <- read.csv("Randindex1.csv", header = TRUE,sep = ';')
scores_n <- scores[,-1]
rownames(scores_n) <- scores[,1]
cluster_scores <- scores_n[c(1:5),]
boxplot(cluster_scores, main="Rand index from clusterprograms on different subsets",
        xlab="Subsets of age", ylab= "Rand index scores", ncol =5, notch = FALSE, col="lightblue")
boxplot(t(cluster_scores), main="Rand index from clusterprograms",
        xlab="Cluster programs", ylab= "Rand index", ncol =5, notch = FALSE, col="lightblue")

# Make the boxplot of the Silhouette score table

scores <- read.csv("Silhouette_scores.csv", header = TRUE)
scores_n <- scores[,-1]
rownames(scores_n) <- scores[,1]
cluster_scores <- scores_n[c(1:5),]
boxplot(cluster_scores, main="Silhouette scores from clusterprograms on different subsets",
        xlab="Subsets of age", ylab= "Silhouette scores", ncol =5, notch = FALSE, col="lightblue")
boxplot(t(cluster_scores), main="Silhouette scores from clusterprograms",
        xlab="Cluster programs", ylab= "Silhouette scores", ncol =5, notch = FALSE, col="lightblue")




