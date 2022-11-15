library(tidyverse)  
library(cluster)    
library(factoextra)
library(FactoMineR)
pca<-prcomp(df, center=FALSE, scale.=FALSE)
results <- pca$x

#Bepalen optimaal aantal clusters

fviz_nbclust(results, FUNcluster=kmeans, k.max = 8) 
fviz_nbclust(results, FUNcluster=kmeans, method="gap_stat", k.max = 8)+ theme_classic()

#clusteren met kmeans

km1<-eclust(results, "kmeans", hc_metric="eucliden",k=2)
fviz_silhouette(km1) 

km1<-eclust(results, "kmeans", hc_metric="eucliden",k=5)
fviz_silhouette(km1) 