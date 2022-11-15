install.packages("Rtsne")
install.packages("tidyverse")
library(tidyverse)
library(Rtsne)

tsne_results <- Rtsne(countfinal, check_duplicates = FALSE)
par(mfrow=c(1,2))
plot(tsne_results$Y, col = "blue", pch = 19, cex = 1.5)
plot(tsne_results$Y, col = "black", bg= as.factor(km1$cluster), pch = 21, cex = 1.5)


tsne_results <- Rtsne(df, check_duplicates = FALSE)
par(mfrow=c(1,2))
plot(tsne_results$Y, col = "blue", pch = 19, cex = 1.5)
plot(tsne_results$Y, col = "black", bg= as.factor(km1$cluster), pch = 21, cex = 1.5)