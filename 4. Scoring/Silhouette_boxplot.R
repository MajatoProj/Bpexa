# Make the boxplot of the Silhouette score table

scores <- read.csv("Silhouette_scores.csv", header = TRUE)
scores_n <- scores[,-1]
rownames(scores_n) <- scores[,1]
cluster_scores <- scores_n[c(1:5),]
boxplot(cluster_scores, main="Silhouette scores from clusterprograms on different subsets",
        xlab="Subsets of age", ylab= "Silhouette scores", ncol =5, notch = FALSE, col="lightblue")
boxplot(t(cluster_scores), main="Silhouette scores from clusterprograms",
        xlab="Cluster programs", ylab= "Silhouette scores", ncol =5, notch = FALSE, col="lightblue")