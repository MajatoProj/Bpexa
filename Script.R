if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("DESeq2")

Primary_Tumors_unstranded_counts <- read_csv("C:/Users/jerem/Desktop/Bpexa/Primary_Tumors_unstranded_counts.csv", skip = 4)
library(DESeq2)
