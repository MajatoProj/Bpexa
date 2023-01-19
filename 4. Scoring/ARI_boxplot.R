# Make the boxplot of the Adjusted Rand Index table

scores <- read.csv("Randindex1.csv", header = TRUE,sep = ';')
scores_n <- scores[,-1]
rownames(scores_n) <- scores[,1]
cluster_scores <- scores_n[c(1:5),]
boxplot(cluster_scores, main="Rand index from clusterprograms on different subsets",
        xlab="Subsets of age", ylab= "Rand index scores", ncol =5, notch = FALSE, col="lightblue")
boxplot(t(cluster_scores), main="Rand index from clusterprograms",
        xlab="Cluster programs", ylab= "Rand index", ncol =5, notch = FALSE, col="lightblue")