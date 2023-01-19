
### Load required packages
library(plyr)

#Data inlezen
counts <- read.csv("Primary_Tumors_unstranded_counts.csv", header = FALSE, stringsAsFactors=FALSE)
clinical <- read.csv("clinical3.csv", comment.char="#", sep = ";")

#Make a subset of clinical
clinical_age <- data.frame(AGE=c(90,89,88,87,86,85,84,83,82,81))
#clinical_age <- data.frame(AGE=c(80,79,78,77,76,75,74,73,72,71))
#clinical_age <- data.frame(AGE=c(70,69,68,67,66,65,64,63,62,61))
#clinical_age <- data.frame(AGE=c(60,59,58,57,56,55,54,53,52,51))
#clinical_age <- data.frame(AGE=c(50,49,48,47,46,45,44,43,42,41))
#clinical_age <- data.frame(AGE=c(40,39,38,37,36,35,34,33,32,31,30,29,28,27,26))

#Select the wanted subtypes
clinical_subtype <- data.frame(SUBTYPE=c("BRCA_Basal", "BRCA_LumA", "BRCA_LumB", "BRCA_Her2", "BRCA_Normal" ))
clinical_er <- data.frame(ER_STATUS=c("Positive", "Negative" ))
clinical <- match_df(clinical, clinical_age)
clinical <- match_df(clinical, clinical_subtype)
clinical <- match_df(clinical, clinical_er)

