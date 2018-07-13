source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
library("data.table")
library("gplots")
library("pheatmap")
# load matrix of read counts

df <- read.table("data/RNA_JB6_Oct2017/RW_all_primary.dedup_hisat2_new.csv", header = T, row.names = "Geneid")
colnames(df)
df1 <- df[,c("Length", "RW03.dedup.bam", "RW04.dedup.bam", "RW07.dedup.bam", "RW02.dedup.bam")]
colnames(df1) <- c("Length", "Control", "TPA", "TPA_CDDO", "TPA_MITC")
 
df11 <- df1[-1]
mat <- matrix(relevel(factor(colnames(df11)), ref = "Control"))
colnames(mat) <- "condition"
dds <- DESeqDataSetFromMatrix(countData=df11, colData=mat, design= ~ condition)
mcols(dds)$basepairs <- df$Length

dds <- dds[rowSums(counts(dds) >= 1) >= 3, ]

dds$condition <- relevel(dds$condition, ref = "Control")
dds$condition
# Control      TPA TPA_CDDO TPA_MITC 
# Control      TPA TPA_CDDO TPA_MITC 
# Levels: Control TPA TPA_CDDO TPA_MITC

dds <- DESeq(dds)

resultsNames(dds)
# [1] "Intercept"                     "condition_TPA_vs_Control"      "condition_TPA_CDDO_vs_Control"
# [4] "condition_TPA_MITC_vs_Control"

res <- results(dds, contrast = c("condition", "TPA", "Control"))
res[order(res$pvalue),]
res[res$pvalue < 0.05,] # 183

sum(res$pvalue < 0.05, na.rm = TRUE)
# 183

plotMA(res, ylim = c(-2,2))

# log2 moderation
resLFC2 <- lfcShrink(dds, 2) # TPA vs Control. or use resLFC <- lfcShrink(dds, 3) for TPA_CDDO vs Control.
resLFC3 <- lfcShrink(dds, 3) # TPA vs Control. or use resLFC <- lfcShrink(dds, 3) for TPA_CDDO vs Control.
 
resLFC_TC <- lfcShrink(dds, contrast = c("condition", "TPA", "Control")) #
resLFC_TCC <- lfcShrink(dds, contrast = c("condition", "TPA_CDDO", "Control")) #
resLFC_TMC <- lfcShrink(dds, contrast = c("condition", "TPA_MITC", "Control")) #


plotMA(resLFC, ylim = c(-2,2)) #more useful


plotMA(resLFC3, ylim = c(-2,2))
plotMA(resLFC_TC, ylim = c(-2,2))
plotMA(resLFC_TCC, ylim = c(-2,2))
plotMA(resLFC_TMC, ylim = c(-2,2))

tiff("data/RNA_JB6_Oct2017/kidney_R01/resLFC_TMC.tiff",
     width = 480, height = 480, units = "px", pointsize = 12,
     bg = "white",
     type = "quartz")
dev.off()

# 2.2 Plot counts
plotCounts(dds, gene = which.min(res$pvalue), intgroup = "condition", col = 1:4)

plotCounts(dds, gene = which.max(res$log2FoldChange), intgroup = "condition", col = 1:4,  xlab = "Groups")

plotCounts(dds, gene = "Nfe2l2", intgroup = "condition", col = 1:4, xlab = "Groups")


abline(v = 3, h = 2000)
dev.off()





#write.table(fpkm(dds1), file = "data/RNA_human_cells_Najoo_2018/rpkm.SUM159_RNA-seq_Suh_04122018.csv", sep = "\t", quote = T, col.names = NA)

rpkm1 <- fpkm(dds)
rpkm1 <- rpkm1[!rowSums(rpkm1 < 10),]

# test[,`:=`(mean_test = apply(.SD, 1, mean), sd_test = apply(.SD, 1, sd),by=id,.SDcols=c('A','B','C','D')]

test1 <- as.data.table(rpkm1)
colnames(test1)
test1 <- cbind(GeneID=rownames(rpkm1), test1)
test2 <- test1[,`:=`(mean_test = apply(.SD, 1, mean), sd_test = apply(.SD, 1, sd)),by=GeneID, .SDcols=c("Control", "TPA", "TPA_CDDO", "TPA_MITC")]
test2$sdm <- test2$sd_test/test2$mean_test
test3 <- test2[order(-sdm)]
test4 <- test3[1:50,]

#heatmap
df1 <- test4[,2:5]
rownames(df1) <- test4$GeneID
df2 <- log2(df1)


pheatmap(df2, cluster_rows= TRUE,
         show_rownames=T,
         cluster_cols=F,
         border_color = NA)
# save as top 50 regulated genes


select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,"condition"])
pheatmap(assay(ntd)[select,],
         cluster_rows=T,
         show_rownames=T,
         cluster_cols=FALSE,
         border_color = NA,
         scale = "row")

# save as top 50 expressed genes.


heatmap.2(as.matrix(df2))




# 
mat_f <- as.matrix(df2)
lmat <- rbind(c(4,1),c(3,1),c(2,1))
lmat
dev.off()
heatmap.2(mat_f,col=redgreen(255),dendrogram = "none",Colv = FALSE,# Rowv = FALSE, 
          trace = "none",
          density.info = "none", key = FALSE, labRow = NA, labCol = NA,
          lmat= lmat,  lwid = c(0,4), lhei = c(4,4,4))
heatmap.2(mat_f,col=redgreen(255),dendrogram = "none",Colv = FALSE, 
          density.info = "none",# key = FALSE,# labCol = NA,
          trace = "none")#,
          #lmat= lmat,  lwid = c(0,4), lhei = c(4,4,4))
heatmap.2(mat_f)

dev.off()

###PCA
rpkm1
dt1 <- t(rpkm1)
dt2 <- data.frame(dt1)
dt3 <- dt2
dt3$sample <- colnames(rpkm1)

autoplot(prcomp(dt2), data = dt3, colour = 'sample')
summary(prcomp(dt2))
# Importance of components:
#         PC1       PC2       PC3       PC4       PC5
# Standard deviation     3758.0150 2120.7684 759.54156 425.41523 1.377e-11
# Proportion of Variance    0.7288    0.2321   0.02977   0.00934 0.000e+00
# Cumulative Proportion     0.7288    0.9609   0.99066   1.00000 1.000e+00




# more colors
# colfunc <- colorRampPalette(c("black", "red"))
colfunc <- colorRampPalette(c("black", "green"))
# heatmap.2(data_matrix,col=colfunc(15),scale="row", trace="none")
heatmap.2(mat_f,col=colfunc(255),dendrogram = "none",Colv = FALSE, 
          density.info = "none",# key = FALSE,# labCol = NA,
          trace = "none")


rpkm2 <- t(rpkm1)
rpkm2_1 <- rpkm2
rpkm2_1$groups <- colnames(rpkm1)
nn <- cbind(groups = colnames(rpkm1))

autoplot(prcomp(rpkm2), data = nn, colour = "groups")



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
heatmap.2(mat_f,col=bluered(255),dendrogram = "none",Colv = FALSE, 
          density.info = "none", key = FALSE, labCol = NA,
          trace = "none",
          lmat= lmat,  lwid = c(0,4), lhei = c(4,4,4))
heatmap.2(mat_f)

rpkm2 <- log10(rpkm1)
mat_f <- rpkm2
rpkm


dds$condition <- relevel(dds$condition, ref="S4")
# filter for genes with at least one count in at least two samples:
dds <- dds[ rowSums(counts(dds) >= 1) >= 2, ]  
#Run deseq
dds_7 <- DESeq(dds)
# pairwise comparisons
res1 <- results(dds_7,contrast=c("condition", "S2", "S4"))
#write out
write.table(res1, file="data/shan_rna/S2_S4_FX-TPA.csv", sep="\t", quote=T, col.names=NA)
