# DEGseq for RNA-seq data without replicates. R Wu April 2018
# source("https://bioconductor.org/biocLite.R")
# biocLite("DEGseq")

library(DEGseq)

geneExpdt1 <- read.table("data/RNA_rat_Kagen_2018/featurecounts.results.rat.csv", sep = "\t", header = T)
colnames(geneExpdt1)
dt <- geneExpdt1[,c("Geneid", "B.dedup.bam", "T.dedup.bam", "C.T.dedup.bam")]
colnames(dt) <- c("Geneid", "Control", "Taxol", "Taxol-Curcumin")
dt1 <- dt[!rowSums(dt[2:4] < 5),] # keep rows with more than 5 reads per sample for all samples.
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 3, groupLabel1 = "Taxol",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_rat_Kagen_2018/Taxol-Control",
       rawCount = TRUE
       )

DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 4, groupLabel1 = "Taxol-Curcumin",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_rat_Kagen_2018/Taxol-Curcumin-Control",
       rawCount = TRUE
        )
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 4, groupLabel1 = "Taxol-Curcumin",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 3, groupLabel2 = "Taxol",
       method = "MARS",
       outputDir = "data/RNA_rat_Kagen_2018/Taxol-Curcumin-Taxol",
       rawCount = TRUE
)
