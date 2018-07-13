# DEGseq for RNA-seq data without replicates. R Wu April 2018
# source("https://bioconductor.org/biocLite.R")
# biocLite("DEGseq")

library(DEGseq)

geneExpdt1 <- read.table("data/RNA_JB6_Oct2017/RW_all_primary.dedup_hisat2_new.csv", sep = "\t", header = T)
dt <- geneExpdt1[,c(1, 7:13)]
dt1 <- dt[!rowSums(dt[2:7] < 5),] # keep rows with more than 5 reads per sample for all samples.
# TPA_1, mITC, Control, TPA_0, CA, FX, CDDO; TPA_1 is the new TPS treated samples by Renyi. TPA_0 is the old by shan.
colnames(dt1)
DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 5, groupLabel1 = "TPA",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 4, groupLabel2 = "Control",
       method = "MARS",
       outputDir = "data/RNA_JB6_Oct2017/DEGseq/RW4-RW3_TPA0-Control",
       rawCount = TRUE
       )

DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 6, groupLabel1 = "Corosolic Acid",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 5, groupLabel2 = "TPA",
       method = "MARS",
       outputDir = "data/RNA_JB6_Oct2017/DEGseq/RW5-RW4_CA-TPA0",
       rawCount = TRUE
        )

DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 3, groupLabel1 = "mITC",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 5, groupLabel2 = "TPA",
       method = "MARS",
       outputDir = "data/RNA_JB6_Oct2017/DEGseq/RW2-RW4_mITC-TPA0",
       rawCount = TRUE
        )

DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 8, groupLabel1 = "CDDO",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 5, groupLabel2 = "TPA",
       method = "MARS",
       outputDir = "data/RNA_JB6_Oct2017/DEGseq/RW7-RW4_CDDO-TPA0",
       rawCount = TRUE
)


DEGexp(geneExpMatrix1 = dt1, geneCol1 = 1, expCol1 = 3, groupLabel1 = "mITC",
       geneExpMatrix2 = dt1, geneCol2 = 1, expCol2 = 2, groupLabel2 = "TPA",
       method = "MARS",
       outputDir = "data/RNA_JB6_Oct2017/DEGseq/RW2-RW1_mITC-TPA1",
       rawCount = TRUE
        )

####
# Merge RNA and Methyl data
# 1. load TPA vs control RNA data (TPA0)
rna1 <- read.table("data/RNA_JB6_Oct2017/DEGseq/RW4-RW3_TPA0-Control.csv", sep = "\t", header = T)
rna1 <- rna1[c(1,5,7,8)]
colnames(rna1) <- c("gene", "log2FC.TPA-Control", "p.TPA-Control", "q.TPA-Control")

rna2 <- read.table("data/RNA_JB6_Oct2017/DEGseq/RW2-RW4_mITC-TPA0.csv", sep = "\t", header = T)
rna2 <- rna2[c(1,5,7,8)]
colnames(rna1) <- c("gene", "log2FC.mITC-TPA", "p.mITC-TPA", "q.mITC-TPA")

rna3 <- read.table("data/RNA_JB6_Oct2017/DEGseq/RW5-RW4_CA-TPA0.csv", sep = "\t", header = T)
rna3 <- rna3[c(1,5,7,8)]
colnames(rna1) <- c("gene", "log2FC.CA-TPA", "p.CA-TPA", "q.CA-TPA")



dna <- read.table("data/Methyl_JB6_Oct2017/use.JB6_methyl_oct_combined_mm10_ca_MITC-CA_fr3s2c5_anno.csv",
                  header = T,
                  sep = "\t")


# m1 <- merge(dna, rna1)
dna1 <- merge(dna, rna1, all.x = T)
# m1r <- merge(dna, rna1, all.y = T)
dna12 <-merge(dna1, rna2, all.x = T)
