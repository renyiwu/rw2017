##DESeq2 for gene express. (originally by John)
#Renyi W. 2017-7-13
library(DESeq2)
library(data.table)
#

#Or load file with this command:
tab <- read.table("data/RNA_cur_8wks/RNA_8wks_primary.dedup.csv", header = T, row.names = "Geneid")
tab <- tab[-(1:4)] #remove columns 1 to 4. Other approaches may also work.
#Assign column names
colnames (tab) <- c("Length", paste("C", c("01",26,29,"02",42,46,65,70,75,79), sep = ""))
# create sample matrix
mat <- matrix(c("Control", "AOM+DSS", "AOM+DSS", "Control", "AOM+DSS+Cur.", "AOM+DSS+Cur.", "DSS", "DSS","DSS+Cur.", "DSS+Cur.", rep(1,3), "2", rep(1,6)), nrow = 10)

#
colnames(mat) <- c("Condition", "Batch")
mat
rownames(mat) <- colnames(tab[-1])


#
tab_10 <- tab[-1]
dds <- DESeqDataSetFromMatrix(countData = tab_10, colData = mat, design = ~ Condition + Batch)
mcols(dds)$basepairs <- tab$Length #For FKPM function.
dds$Condition <- relevel(dds$Condition, ref="Control")

# filter for genes with at least one count in at least two samples:
dds1 <- dds[rowSums(counts(dds) >= 1) >= 2, ]
#Run deseq
dds_10 <- DESeq(dds1)
# pairwise comparisons
res1 <- results(dds_10, contrast=c("Condition", "DSS", "Control"))
write.table(res1, file="data/rna_8wks/DSS_Control.csv", sep="\t", quote=F, col.names=NA)
res2 <- results(dds_10,contrast=c("Condition", "AOM+DSS", "Control"))
write.table(res2, file="data/rna_8wks/AOM+DSS_Control.csv", sep="\t", quote=F, col.names=NA)
res3 <- results(dds_10,contrast=c("Condition", "DSS+Cur.", "DSS"))
write.table(res3, file="data/rna_8wks/DSS+Cur._DSS.csv", sep="\t", quote=F, col.names=NA)
res4 <- results(dds_10,contrast=c("Condition", "AOM+DSS", "DSS"))
write.table(res4, file="data/rna_8wks/AOM+DSS_DSS.csv", sep="\t", quote=F, col.names=NA)
res5 <- results(dds_10,contrast=c("Condition", "AOM+DSS+Cur.", "AOM+DSS"))
write.table(res5, file="data/rna_8wks/AOM+DSS+Cur._AOM+DSS.csv", sep="\t", quote=F, col.names=NA)

#FMPK
dds_f<- DESeq(dds)
write.table(fpkm(dds_f), file = "data/rna_8wks/rna_8wks_all.fpkm.csv", sep = "\t", quote = T, col.names = NA)
#
##

#
#####################
#PCA of 18 weeks samples.
#to find outlier?
#####################

#log transformation of dds_18
rld <- rlog(dds_10, blind=F)

plotPCA(rld, intgroup = "Condition")
#or
plotPCA(rld, returnData = T, intgroup = "Condition")
#or
plotPCA(rld, intgroup = c("Condition","Batch"))
#or
plotMA(dds_10)#, ylim=c(-4,4))
#
#export to file
tiff(filename = "data/John/PCA_18wks_400x300.tiff", 
     width = 400, height = 300, units = "px", pointsize = 12,
     compression = "none",type = "windows",antialias = "none")
plotPCA(rld,
        ntop = 500,
        intgroup = c("condition"))
dev.off()
pwd()

  #
#
# bmp(filename = "Rplot%03d.bmp",
#     width = 480, height = 480, units = "px", pointsize = 12,
#     bg = "white", res = NA, ...,
#     type = c("cairo", "Xlib", "quartz"), antialias)
# 
# jpeg(filename = "Rplot%03d.jpeg",
#      width = 480, height = 480, units = "px", pointsize = 12,
#      quality = 75,
#      bg = "white", res = NA, ...,
#      type = c("cairo", "Xlib", "quartz"), antialias)
# 
# png(filename = "Rplot%03d.png",
#     width = 480, height = 480, units = "px", pointsize = 12,
#     bg = "white",  res = NA, ...,
#     type = c("cairo", "cairo-png", "Xlib", "quartz"), antialias)
# 
# tiff(filename = "Rplot%03d.tiff",
#      width = 480, height = 480, units = "px", pointsize = 12,
#      compression = c("none", "rle", "lzw", "jpeg", "zip", "lzw+p", "zip+p"),
#      bg = "white", res = NA,  ...,
#      type = c("cairo", "Xlib", "quartz"), antialias)

