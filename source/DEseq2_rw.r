# DEseq2 practice, R Wu. April 2018

# source("https://bioconductor.org/biocLite.R")
# biocLite("rnaseqGene")

library("DESeq2")
library("data.table")

dt <- fread("data/RNA_cur_8wks/RNA_8wks_primary.picard-dedup.csv")
dt <- dt[,-(2:5)]

countdata <- dt[,-(1:2)]
rownames(countdata) <- dt$Geneid
colnames(countdata) <- paste("C", c(01,26,29,02,42,46,65,70,75,79), sep = "")
countdata <- as.matrix(countdata)

coldata <- data.frame(group = c("Control", "AOM+DSS", "AOM+DSS", "Control",
                                "AOM+DSS+Cur.", "AOM+DSS+Cur.", "DSS", "DSS",
                                "DSS+Cur.", "DSS+Cur."))
rownames(coldata) <- colnames(countdata)
#colnames(coldata) <- "group"

ddsmat <- DESeqDataSetFromMatrix(countdata, coldata, design = ~ group)
ddsmat
mcols(ddsmat)
keep <- rowSums(counts(ddsmat)) >= 5 
# or,  ddsmat$group <- relevel (ddsmat$group, ref = "Control")
ddsmat$group <- factor(ddsmat$group, levels = c("Cnotrol", "DSS", "DSS+Cur.", "AOM+DSS", "AOM+DSS+Cur."))
ddsmat1 <- ddsmat[keep,]

dds <- DESeq(ddsmat1)
res <- results(dds)
res1 <- results(dds, contrast = c("group", "DSS", "Control"))

plotMA(res1)
plotCounts(dds, gene = which.max(res1$padj), intgroup = "group")

colData(dds)
plotPCA(rlog(dds), intgroup = "group")



# hierarchical clustering

counts(ddsmat1)
dtt <- t(as.matrix(counts(ddsmat1)))
ds <- dist(dtt)
h <- hclust(ds)
plot(h, xlab = "sample", asp = 2)
# Parameters for plotting
plot(h, hang = -1, xlab = "sample", ylab = "distance", cex = 1.2, lwd = 2, col = "black")
#rect.hclust(h, 5)
dev.off()
