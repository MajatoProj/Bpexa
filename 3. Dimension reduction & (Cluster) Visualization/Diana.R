#diana 
dia <- eclust(pca_data, "diana", k=2)
fviz_silhouette(dia)

#cluster diana
cluster_data <- data.frame(dia$cluster)
clinical_cluster <- cbind(clinical, cluster_data)

clinical_cluster$dia.cluster <- as.character(clinical_cluster$dia.cluster)

library(umap)
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