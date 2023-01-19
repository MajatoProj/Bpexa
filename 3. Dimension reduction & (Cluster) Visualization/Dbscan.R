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

library(umap)
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