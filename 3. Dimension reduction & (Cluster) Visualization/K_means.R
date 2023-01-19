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

library(umap)
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
