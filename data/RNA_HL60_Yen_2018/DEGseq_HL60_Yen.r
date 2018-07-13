# DEGseq for RNA-seq data without replicates. R Wu April 2018
# source("https://bioconductor.org/biocLite.R")
# biocLite("DEGseq")

library(DEGseq)

geneExpdt1 <- read.table("data/RNA_hct116_Qian_2018/featurecounts.results.human.csv", sep = "\t", header = T)
dt <- geneExpdt1[,c("Geneid", "C.dedup.bam", "RA.dedup.bam", "UA.dedup.bam", "SFN.dedup.bam", "URA.dedup.bam", "SRA.dedup.bam")]

colnames(dt) <- c("Geneid", "Control", "Radiation", "UA", "SFN", "Ra-UA", "Ra-SFN")
dt1 <- dt[!rowSums(dt[2:4] < 5),] # keep rows with more than 5 reads per sample for all samples.

DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 3, groupLabel1 = "Radiation",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_HL60_Yen_2018/DEGseq/Radiation-Control",
       rawCount = TRUE
       )
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 4, groupLabel1 = "UA",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_HL60_Yen_2018/DEGseq/UA-Control",
       rawCount = TRUE
        )
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 5, groupLabel1 = "SFN",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_HL60_Yen_2018/DEGseq/SFN-Control",
       rawCount = TRUE
        )#
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 6, groupLabel1 = "Ra-UA",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 4, groupLabel2 = "Radiation",
       method = "MARS",
       outputDir = "data/RNA_HL60_Yen_2018/DEGseq/RaUA-Radiation",
       rawCount = TRUE
)#
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 7, groupLabel1 = "Ra-SFN",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 4, groupLabel2 = "Radiation",
       method = "MARS",
       outputDir = "data/RNA_HL60_Yen_2018/DEGseq/RaSFN-Radiation",
       rawCount = TRUE
)#
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 7, groupLabel1 = "Ra-SFN",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_HL60_Yen_2018/DEGseq/RaSFN-Control",
       rawCount = TRUE
)#
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 6, groupLabel1 = "Ra-UA",
        geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "Control",
        method = "MARS",
        outputDir = "data/RNA_HL60_Yen_2018/DEGseq/RaUA-Control",
        rawCount = TRUE
        )#