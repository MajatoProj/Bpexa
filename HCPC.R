library(tidyverse)  
library(cluster)    
library(factoextra)
library(FactoMineR)

hcpc <- HCPC(pca, graph = FALSE)
fviz_cluster(hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)