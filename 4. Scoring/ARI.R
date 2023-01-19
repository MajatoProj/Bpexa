#Adjusted rand index

#pam
if (any(pam1$cluster == 1)) {
  # Change the value to 0
  pam1$adjrand[pam1$cluster == 1] <- 'Positive'
}
if (any(pam1$cluster == 2)) {
  # Change the value to 0
  pam1$adjrand[pam1$cluster == 2] <- 'Negative'
}
#calculate the adjusted rand index
adjusted_pam <- adjustedRandIndex(clinical$ER_STATUS, pam1$adjrand)


#kmeans
if (any(km1$cluster == 1)) {
  # Change the value to 0
  km1$adjrand[km1$cluster == 1] <- 'Positive'
}
if (any(km1$cluster == 2)) {
  # Change the value to 0
  km1$adjrand[km1$cluster == 2] <- 'Negative'
}
#calculate the adjusted rand index
adjusted_kmeans <- adjustedRandIndex(clinical$ER_STATUS, km1$adjrand)


#Diana
if (any(dia$cluster == 1)) {
  # Change the value to 0
  dia$adjrand[dia$cluster == 1] <- 'Positive'
}
if (any(dia$cluster == 2)) {
  # Change the value to 0
  dia$adjrand[dia$cluster == 2] <- 'Negative'
}
#calculate the adjusted rand index
adjusted_diana <- adjustedRandIndex(clinical$ER_STATUS, dia$adjrand)


#Dbscan is difficult to cluster, sometimes it gives more than 2 clusters. 
#If that is the case, look at umap and decide if that cluster is supposed to be ER+ or ER-.
#dbscan
if (any(d$cluster == 0)) {
  # Change the value to 0
  d$adjrand[d$cluster == 0] <- 'Positive'
}
if (any(d$cluster == 1)) {
  # Change the value to 0
  d$adjrand[d$cluster == 1] <- 'Positive'
}
if (any(d$cluster == 2)) {
  # Change the value to 0
  d$adjrand[d$cluster == 2] <- 'Negative'
}
if (any(d$cluster == 3)) {
  # Change the value to 0
  d$adjrand[d$cluster == 3] <- 'Positive'
}
#calculate the adjusted rand index
adjusted_dbscan <- adjustedRandIndex(clinical$ER_STATUS, d$adjrand)