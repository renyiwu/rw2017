# DEGseq for RNA-seq data without replicates. R Wu April 2018
# source("https://bioconductor.org/biocLite.R")
# biocLite("DEGseq")

library(DEGseq)

geneExpdt1 <- read.table("data/RNA_JB6_Oct2017/RW_all_primary.dedup_hisat2_new.csv", sep = "\t", header = T)
dt <- geneExpdt1[,c(1, 7:13)]
dt1 <- dt[!rowSums(dt[2:7] < 5),] # keep rows with more than 5 reads per sample for all samples.
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 5, groupLabel1 = "TPA",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 4, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_JB6_Oct2017/DEGseq/RW4-RW3_TPA1-Control",
       rawCount = TRUE
       )


