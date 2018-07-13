# DEGseq for RNA-seq data without replicates. R Wu April 2018
# source("https://bioconductor.org/biocLite.R")
# biocLite("DEGseq")

library(DEGseq)

geneExpdt1 <- read.table("data/RNA_hct116_Qian_2018/featurecounts.results.human.csv", sep = "\t", header = T)
dt <- geneExpdt1[,c("Geneid", "Con.dedup.bam", "F8.dedup.bam", "L30.dedup.bam")]
colnames(dt) <- c("Geneid", "Control", "Formula", "Luteolin")
dt1 <- dt[!rowSums(dt[2:4] < 5),] # keep rows with more than 5 reads per sample for all samples.
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 3, groupLabel1 = "Formula",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_hct116_Qian_2018/DEGseq/Formula-Control",
       rawCount = TRUE
       )

DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 4, groupLabel1 = "Luteolin",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_hct116_Qian_2018/DEGseq/Luteolin-Control",
       rawCount = TRUE
        )

