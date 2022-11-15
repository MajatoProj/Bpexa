counts <- read.csv("~/GitHub/Bpexa/Primary_Tumors_unstranded_counts.csv", header = FALSE, stringsAsFactors=FALSE)
clinicalnew <- read.csv("~/GitHub/Bpexa/clinical3.txt", comment.char="#")

names(counts) <- counts[1,]
clinical <- subset(clinicalnew, select = -c(X, OTHER_PATIENT_ID))
clinical <- clinical[order(clinical$PATIENT_ID),]

countno <- counts[,-1]
rownames(countno) <- counts[,1]
countfinal <- countno[-1, ]
colnames(countfinal) <- countno[1, ]
countfinal <- subset(countfinal, select=-c(gene_name,gene_type))

rownames(clinical)<- clinical[,1]
all(rownames(clinical) %in% colnames(countfinal))
countfinal <- countfinal[,rownames(clinical)]
all(rownames(clinical) == colnames(countfinal))
library(DESeq2)
ncol(countfinal) == nrow(clinical)
library(dplyr)
countfinal <- mutate_all(countfinal, function(x) as.numeric(as.character(x)))

library(tibble)
df <- tibble::rownames_to_column(countfinal, "VALUE")