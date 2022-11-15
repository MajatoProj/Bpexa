library(tidyverse)  
library(cluster)    
library(factoextra)
library(FactoMineR)

dia <- eclust(results, "diana", k=2)
fviz_silhouette(dia)

dia <- eclust(results, "diana", k=5)
fviz_silhouette(dia)