if (!("umap" %in% installed.packages())) {
  # Install umap package
  BiocManager::install("umap", update = FALSE)
}


library(umap)
umap_results <- umap::umap(df, intgroup=c("SUBTYPE"), returnData=TRUE)
umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("PATIENT_ID") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(clinical, by = "PATIENT_ID")
ggplot(umap_plot_df,aes(x = X1,y = X2,
                        color = SUBTYPE,
                        shape = DSS_STATUS)) + geom_point()