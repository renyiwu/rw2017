# Mes-13 RNA-seq analysis no replicates with edgeR

# Install edgeR from bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")

# Load edgeR
library(edgeR)

# Set working director to where count files stored
getwd()
setwd("C:/Users/davichen/Desktop/kidney.ro1")
getwd()

# Import featurecounts output

dt1 <- read.table("mes13_featurecounts_122017.csv", 
                  sep = "\t",
                  header = TRUE,
                  stringsAsFactors = FALSE)
# Remove columns to select only Geneid and sample columns

dt2 <- dt1[,c(1,7:13)]

# Rename samples to WR1 to WR7

colnames(dt2)[2:8] <- paste("WJ",
                            1:7,
                            sep = "")

# Focus only on LG, HG, and MITC

dt3 <- dt2[,c(1:4)]

bcv <- 0.1
y <- DGEList(counts = dt3[,2:4], group =2:4, remove.zeros = TRUE)
et <- exactTest(y, dispersion=bcv^2)


