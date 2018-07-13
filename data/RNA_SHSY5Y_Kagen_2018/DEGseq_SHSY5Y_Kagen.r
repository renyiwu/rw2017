# DEGseq for RNA-seq data without replicates. R Wu April 2018
# source("https://bioconductor.org/biocLite.R")
# biocLite("DEGseq")

library(DEGseq)

geneExpdt1 <- read.table("data/RNA_SHSY5Y_Kagen_2018/featurecounts.results.human.csv", sep = "\t", header = T)
colnames(geneExpdt1)
dt <- geneExpdt1[,c("Geneid", "X1.Blank.dedup.bam", "X4.Cur.2uM.dedup.bam", "X7.PTX.5.ng.mL.dedup.bam", "X9.Cur.PTX.5.ng.mL.dedup.bam")]


colnames(dt) <- c("Geneid", "Control", "Curcumin", "PTX", "Curcumin-PTX")

dt1 <- dt[!rowSums(dt[2:4] < 5),] # keep rows with more than 5 reads per sample for all samples.
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 3, groupLabel1 = "Curcumin",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_SHSY5Y_Kagen_2018/Curcumin-Control",
       rawCount = TRUE
       )

DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 4, groupLabel1 = "PTX",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_SHSY5Y_Kagen_2018/PTX-Control",
       rawCount = TRUE
        )
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 5, groupLabel1 = "Curcumin-PTX",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_SHSY5Y_Kagen_2018/Curcumin_PTX-Control",
       rawCount = TRUE
)

DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 5, groupLabel1 = "Curcumin-PTX",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 3, groupLabel2 = "Curcumin",
       method = "MARS",
       outputDir = "data/RNA_SHSY5Y_Kagen_2018/Curcumin_PTX-Curcumin",
       rawCount = TRUE
)

DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 5, groupLabel1 = "Curcumin-PTX",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 4, groupLabel2 = "PTX",
       method = "MARS",
       outputDir = "data/RNA_SHSY5Y_Kagen_2018/Curcumin_PTX-PTX",
       rawCount = TRUE
)