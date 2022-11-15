library(tidyverse)  
library(cluster)    
library(factoextra)
library(FactoMineR)

#bepalen optimaal aantal clusters voor pam
fviz_nbclust(results, FUNcluster=cluster::pam, k.max = 5)
fviz_nbclust(results, FUNcluster=cluster::pam, method="gap_stat", k.max  = 7)+ theme_classic()

#clusteren pam
pam1<-eclust(results, "pam", k=2)
fviz_silhouette(pam1)

pam1<-eclust(results, "pam", k=5)
fviz_silhouette(pam1)