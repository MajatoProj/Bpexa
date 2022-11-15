library(tidyverse)  
library(cluster)    
library(factoextra)
library(FactoMineR)

hcl <- eclust(results, "hclust", k=2)
fviz_silhouette(hcl)

hcl <- eclust(results, "hclust", k=5)
fviz_silhouette(hcl)