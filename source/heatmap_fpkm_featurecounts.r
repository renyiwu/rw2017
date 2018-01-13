source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)
library(data.table)
library("gplots")
# load matrix of read counts
df <- read.table("data/shan_rna/RW_all_primary.dedup.csv", header = T, row.names = "Geneid")
df <- df[,-(1:4)]
colnames(df) <- c("length", paste("RW", 1:7, sep = ""))
df_7 <- df[-1]
# create sample matrix
mat <- matrix(paste("S",1:7, sep = ""), nrow = 7)
colnames(mat) <- "condition"
mat
rownames(mat) <- colnames(df[-1])
mat_7 <- mat
dds <- DESeqDataSetFromMatrix(countData=df_7, colData=mat_7, design= ~ condition)
mcols(dds)$basepairs <- df$length
dds_7 <- DESeq(dds)

write.table(fpkm(dds_7), file = "data/shan_rna/RW_all.fpkm.csv", sep = "\t", quote = T, col.names = NA)

#
fpkm_7 <- fpkm(dds_7)
fpkm_71 <- fpkm_7+1 #Add pseudo 1 to all values to avoid log10(0) error.
fpkm_72 <- log10(fpkm_7+1) #log10 transformation
mat_f <- fpkm_72[!rowSums(fpkm_72 <1),] #Only keep rows which have values greater than or equal to 2 (in all groups/columns)
write.table(log10(fpkm(dds_7)+1), file = "data/shan_rna/RW_all_log10_fpkm.csv", sep = "\t", quote = T, col.names = NA)
write.table(mat_f, file = "data/shan_rna/RW_all_log10_fpkm_gt1.csv", sep = "\t", quote = T, col.names = NA)

#tmp <- read.table(file = "data/shan_rna/RW_all_log10_fpkm_gt11.csv", sep = "\t", header = T, row.names = "X")
#heatmap
lmat <- rbind(c(4,1),c(3,1),c(2,1))
lmat
dev.off()
heatmap.2(mat_f,col=redgreen(255),dendrogram = "none",Colv = FALSE,Rowv = FALSE, 
          trace = "none",
          density.info = "none", key = FALSE, labRow = NA, labCol = NA,
          lmat= lmat,  lwid = c(0,4), lhei = c(4,4,4))
heatmap.2(mat,col=bluered(255),dendrogram = "none",Colv = FALSE, 
          density.info = "none", key = FALSE, labCol = NA,
          trace = "none",
          lmat= lmat,  lwid = c(0,4), lhei = c(4,4,4))
heatmap.2(mat_f)


dds$condition <- relevel(dds$condition, ref="S4")
# filter for genes with at least one count in at least two samples:
dds <- dds[ rowSums(counts(dds) >= 1) >= 2, ]  
#Run deseq
dds_7 <- DESeq(dds)
# pairwise comparisons
res1 <- results(dds_7,contrast=c("condition", "S2", "S4"))
#write out
write.table(res1, file="data/shan_rna/S2_S4_FX-TPA.csv", sep="\t", quote=T, col.names=NA)
