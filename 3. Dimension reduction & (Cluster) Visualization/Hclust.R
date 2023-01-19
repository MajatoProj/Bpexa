#hclust
hcl <- eclust(pca_data, "hclust", k=2)
fviz_silhouette(hcl)

#cluster hcl
cluster_data <- data.frame(hcl$cluster)
clinical_cluster <- cbind(clinical, cluster_data)

clinical_cluster$hcl.cluster <- as.character(clinical_cluster$hcl.cluster)

library(umap) 
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