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

library(umap)
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