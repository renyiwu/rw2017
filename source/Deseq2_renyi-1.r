##DESeq2 for gene express. (originally by John)
#Renyi W. 2017-7-13
library(DESeq2)
library(data.table)
#
# load matrix of read counts
tab <- fread("data/shan_rna/RW_all_primary.dedup.csv")
tab <- data.frame(tab)
#tab1 <- read.table("data/shan_rna/RW_all_primary.dedup.csv", sep = "\t", row.names = "Geneid")
#tab2 <- read.csv("data/shan_rna/RW_all_primary.dedup.csv")#, header = T, sep = "\t", row.names = "Geneid")
rownames(tab) <- tab$Geneid
<<<<<<< HEAD
tab <- tab[,-(1:5)] #remove columns 1 to 5. Other approaches may also work. tab[,1:5] <- NULL
=======

#Or load file with this command:
tab <- read.table("data/rna_8wks/RNA_8wks_primary.dedup.csv", header = T, row.names = "Geneid")
tab <- tab[,-(1:5)] #remove columns 1 to 5. Other approaches may also work.
>>>>>>> 69cd6b0262febc4ccd7bde53a80ced43d0b4bfdc
#Assign column names
colnames(tab) <- c("length", paste("RW",1:7, sep = ""))
tab_7 = tab[-1]
# create sample matrix
mat <- matrix(paste("S",1:7, sep = ""), nrow = 7, ncol = 1, byrow = F) #byrow default to F, meaning top-bottom then left-right
colnames(mat) <- "condition"
mat
rownames(mat) <- colnames(tab[-1])
mat_7 <- mat
dds <- DESeqDataSetFromMatrix(countData=tab_7, colData=mat_7, design= ~ condition)
dds$condition <- relevel(dds$condition, ref="S4")
# filter for genes with at least one count in at least two samples:
dds <- dds[ rowSums(counts(dds) >= 1) >= 2, ]  # down to 17979 genes
#Run deseq
dds_7 <- DESeq(dds)
# pairwise comparisons
res1 <- results(dds_7,contrast=c("condition", "S6", "S4"))
#write out
write.table(res1, file="data/shan_rna/S6-S4_mITC-TPA_l4.csv", sep="\t", quote=T, col.names=NA)
#
#
#
#
#
#
#
#
#
# load matrix of read counts
tabj <- read.csv("data/John/comb.count", sep="\t", row.names="Symbol")

# create sample matrix
mat <- matrix(c(rep("Control_8", 2), rep("AOM_DSS_8", 2), rep("AOM_DSS_Cur_8", 2),
  rep("Control_18", 4), rep("AOM_DSS_18", 4), rep("AOM_DSS_Cur_18", 4),
  0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1), nrow=18, ncol=2)
colnames(mat) <- c("condition", "batch")
row.names(mat) <- colnames(tab[-1])
# analyze just 18 week samples AND without C36, an outlier.
tab_18 <- tab[,c('C20', 'C14', 'C15', 'C19', 'C34', 'C40', 'C33', 'C54', 'C60', 'C55', 'C59')]
# alternative:
# tab_18 <- subset(tab, select=c('C20', 'C14', 'C15', 'C19', 'C34', 'C40', 'C33', 'C36', 'C54', 'C60', 'C55', 'C59'))
#tab_18$C36 <- NULL  # remove outlier -- see below
mat_18 <- mat[c('C20', 'C14', 'C15', 'C19', 'C34', 'C40', 'C33', 'C54', 'C60', 'C55', 'C59'),]
dds_18 <- DESeqDataSetFromMatrix(countData=tab_18, colData=mat_18, design= ~ batch + condition)
dds_18$condition <- relevel(dds_18$condition, ref="Control_18")
dds_18$
# filter for genes with at least one count in at least two samples:
dds_18 <- dds_18[ rowSums(counts(dds_18) >= 1) >= 2, ]  # down to 17979 genes
dds_18 <- DESeq(dds_18)
# pairwise comparisons
res1 <- results(dds_18,contrast=c("condition", "AOM_DSS_Cur_18", "Control_18"))  # equivalent to results(dds, contrast=c("condition", "AOM-DSS-Cur", "Control"))
res2 <- results(dds_18, contrast=c("condition", "AOM_DSS_18", "Control_18"))
res3 <- results(dds_18, contrast=c("condition", "AOM_DSS_Cur_18", "AOM_DSS_18"))
#
# Optional, subset results matrices to FDR < 10%, sort by log2FC, and write to file
#resSig <- subset(res, padj < 0.1)
#resSig <- resSig[order(resSig$log2FoldChange), ]
write.table(res1, file="data/John/AOM-DSS-Cur_Control.csv", sep="\t", quote=F, col.names=NA)
write.table(res2, file="data/John/AOM-DSS_Control.csv", sep="\t", quote=F, col.names=NA)
write.table(res3, file="data/John/AOM-DSS-Cur_AOM-DSS.csv", sep="\t", quote=F, col.names=NA)
#
#
#
#get FKPM data
tab <- read.csv("data/John/comb.count", sep="\t", row.names="Symbol")
# create sample matrix
mat <- matrix(c(rep("Control_8", 2), rep("AOM_DSS_8", 2), rep("AOM_DSS_Cur_8", 2),
                rep("Control_18", 4), rep("AOM_DSS_18", 4), rep("AOM_DSS_Cur_18", 4),
                0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1), nrow=18, ncol=2)
colnames(mat) <- c("condition", "batch")
row.names(mat) <- colnames(tab[-1]) # "-1" means not including the first entry.
tab_18 <- tab[,c('C20', 'C14', 'C15', 'C19', 'C34', 'C40', 'C33', 'C36', 'C54', 'C60', 'C55', 'C59')]
mat_18 <- mat[c('C20', 'C14', 'C15', 'C19', 'C34', 'C40', 'C33', 'C36', 'C54', 'C60', 'C55', 'C59'),]

dds_18 <- DESeqDataSetFromMatrix(countData=tab_18, colData=mat_18, design= ~ batch + condition)
# to get FPKM values, need transcript lengths
mcols(dds_18)$basepairs <- tab$Length  # 'Length' column of comb.count (produced by convert4.py)
dds_18$condition <- relevel(dds_18$condition, ref="Control_18")
dds_18 <- DESeq(dds_18)
write.table(fpkm(dds_18), file="data/John/samples_18_fpkm_all.csv", sep="\t", quote=F, col.names=NA)
#
#
#
#####################
#PCA of 18 weeks samples.
#to find outlier?
#####################
#
#
tab <- read.csv("data/John/comb.count", sep="\t", row.names="Symbol")
# create sample matrix
mat <- matrix(c(rep("Control_8", 2), rep("AOM_DSS_8", 2), rep("AOM_DSS_Cur_8", 2),
                rep("Control_18", 4), rep("AOM_DSS_18", 4), rep("AOM_DSS_Cur_18", 4),
                0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1), nrow=18, ncol=2)
colnames(mat) <- c("condition", "batch")
row.names(mat) <- colnames(tab[-1])
# analyze just 18 week samples AND without C36, an outlier.
tab_18 <- tab[,c('C20', 'C14', 'C15', 'C19', 'C34', 'C40', 'C33', 'C33', 'C54', 'C60', 'C55', 'C59')]
# alternative:
# tab_18 <- subset(tab, select=c('C20', 'C14', 'C15', 'C19', 'C34', 'C40', 'C33', 'C36', 'C54', 'C60', 'C55', 'C59'))
#tab_18$C36 <- NULL  # remove outlier -- see below
mat_18 <- mat[c('C20', 'C14', 'C15', 'C19', 'C34', 'C40', 'C33', 'C36', 'C54', 'C60', 'C55', 'C59'),]
#Optional. add another column to mat_18, namely Sample.
mat_18 <- cbind(mat_18, Sample=c('C20', 'C14', 'C15', 'C19', 'C34', 'C40', 'C33', 'C36', 'C54', 'C60', 'C55', 'C59'))
#change condition description.
mat_18[,1] <- matrix(c(rep("Control",4), rep("AOM+DSS",4), rep("AOM+DSS+Cur", 4)))
dds_18 <- DESeqDataSetFromMatrix(countData=tab_18, colData=mat_18, design= ~ batch + condition)
dds_18$condition <- relevel(dds_18$condition, ref="Control")
dds_18
#Optional. filter for genes with at least one count in at least two samples:
#dds_18 <- dds_18[ rowSums(counts(dds_18) >= 1) >= 2, ]  # down to 17979 genes
#main function
dds_18 <- DESeq(dds_18)
#
#log transformation of dds_18
rld <- rlog(dds_18, blind=F)
plotPCA(rld)
#or
pcad <- plotPCA(rld) #, returnData = T)
#or
pcad2 <- plotPCA(rld, returnData = T)
#or
plotPCA(rld, intgroup = c("condition","Sample"))
#or
plotMA(dds_18)#, ylim=c(-4,4))
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


# output files:
# res -> resSig -> AOM-DSS-Cur_Control.csv
# res2 -> res2Sig -> AOM-DSS_Control.csv
# res3 -> res3Sig -> AOM-DSS-Cur_AOM-DSS.csv
# NOTE: log2FoldChanges give 1st sample w.r.t 2nd --
#   i.e. in AOM-DSS-Cur_Control.csv, log2FC < -1 => lower expression in AOM-DSS-Cur than Control

summary(res1)
# out of 17979 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 222, 1.2% 
# LFC < 0 (down)   : 381, 2.1% 
# outliers [1]     : 15, 0.083% 
# low counts [2]   : 3835, 21% 
# sum(res$padj < 0.1, na.rm=T)
# # [1] 603   # sum of 222 up and 381 down
# > summary(res2)
# LFC > 0 (up)     : 270, 1.5% 
# LFC < 0 (down)   : 309, 1.7% 
# > summary(res3)
# LFC > 0 (up)     : 421, 2.3% 
# LFC < 0 (down)   : 671, 3.7% 
# NOTE: no filtering based on log2FC; numbers are quite low though:
# > sum(resSig$log2FoldChange < -1)
# [1] 11
# > sum(resSig$log2FoldChange > 1)
# [1] 25



#/////////////////////////////////////////////////////////////////////////////////////

# how to get FPKM values

# load data (same as before):
library(DESeq2)
tab <- read.csv("comb.count", sep="\t", row.names="Symbol")
mat <- matrix(c(rep("Control_8", 2), rep("AOM_DSS_8", 2), rep("AOM_DSS_Cur_8", 2),
  rep("Control_18", 4), rep("AOM_DSS_18", 4), rep("AOM_DSS_Cur_18", 4),
  0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1), nrow=18, ncol=2)
colnames(mat) <- c("condition", "batch")
row.names(mat) <- colnames(tab[-1])

# analyze just 8 week samples
tab_8 <- tab[,c('C1', 'C2', 'C26', 'C29', 'C42', 'C46')]
mat_8 <- mat[c('C1', 'C2', 'C26', 'C29', 'C42', 'C46'),]
dds_8 <- DESeqDataSetFromMatrix(countData=tab_8, colData=mat_8, design= ~ batch + condition)
# to get FPKM values, need transcript lengths
mcols(dds_8)$basepairs <- tab$Length  # 'Length' column of comb.count (produced by convert4.py)

dds_8$condition <- relevel(dds_8$condition, ref="Control_8")
#dds_8 <- dds_8[ rowSums(counts(dds_8) >= 1) >= 2, ]  # filter values (?)
dds_8 <- DESeq(dds_8)
write.table(fpkm(dds_8), file="samples_8_fpkm-.csv", sep="\t", quote=F, col.names=NA)

# FPKM averages can be added to matrix (actually easier to do with a data.frame)
fpkm_8 <- data.frame(fpkm(dds_8))
fpkm_8$Control <- rowMeans(fpkm_8[,c('C1', 'C2')])
fpkm_8$AOM_DSS <- rowMeans(fpkm_8[,c('C26', 'C29')])
fpkm_8$AOM_DSS_Cur <- rowMeans(fpkm_8[,c('C42', 'C46')])
write.table(fpkm_8, file='samples_8_fpkm.csv', sep='\t', quote=F, col.names=NA)



#////////////////////////////////////////////////////////////////////////////////////////////


# alternative for making list of transcript lengths: https://www.biostars.org/p/83901/
library(GenomicFeatures)
txdb <- makeTxDbFromGFF('/home/john/genomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf',
  format='gtf', organism='Mus musculus')
  # n.b. "In BioC 3.1, all the makeTranscriptDbFrom*() functions were renamed makeTxDbFrom*()."
exonsList <- exonsBy(txdb, by="gene")
# two possibilities:
dds_8@rowRanges <- exonsList                                      # add to DESeqDataSet directly(?)
geneSizes <- lapply(exonsList,function(x){sum(width(reduce(x)))}) # or create lengths to add later
#  -> errors with both approaches



/////////////////////////////////////////////////////////////////////////////////////////////

// dedup with picard
$ fol="C14"; java -jar ~/tools/picard-tools-2.5.0/picard.jar MarkDuplicates \
  I=$fol/accepted_hits.bam O=$fol/accepted_hits_dedup.bam \
  M=$fol/dedup.txt REMOVE_DUPLICATES=true
// etc. for the other samples

$ /home/john/tools/cufflinks-2.2.1.Linux_x86_64/cuffdiff -p 8 \
  -o cuffdiff_dedup2 --no-update-check \
  -L Control,AOM-DSS,AOM-DSS-Cur \
  /home/john/genomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf \
  ../../RNAseq/results/C20/accepted_hits_dedup.bam,C14_ca/accepted_hits_dedup.bam,C15_ca/accepted_hits_dedup.bam,C19_ca/accepted_hits_dedup.bam \
  ../../RNAseq/results/C34/accepted_hits_dedup.bam,../../RNAseq/results/C40/accepted_hits_dedup.bam,C33_ca/accepted_hits_dedup.bam,C36_ca/accepted_hits_dedup.bam \
  ../../RNAseq/results/C54/accepted_hits_dedup.bam,../../RNAseq/results/C60/accepted_hits_dedup.bam,C55_ca/accepted_hits_dedup.bam,C59_ca/accepted_hits_dedup.bam

// dendrogram looks better w.r.t. batch effect,
//   but no appreciable improvement in MDSplot or PCAplot


// cummeRbund
$ R
> library(cummeRbund)
> cuff <- readCufflinks(gtfFile='/home/john/genomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf', genome='mm9')
> csDendro(genes(cuff), replicates=T)

> pdf("aom_dss.pdf")
> par(mar=c(8, 4, 4, 2) + 0.1)
> csDendro(genes(cuff), replicates=T)
> dev.off()
// large batch effect


//////////////////////////////////////////////////////////////////////////

// subsetting genes
$ R
> library(cummeRbund)
> cuff <- readCufflinks(
  gtfFile='/home/john/genomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf',
  genome='mm9')

# load data
> diff <- diffData(genes(cuff))

# retrieve genes called significant by cufflinks
> diffGenes <- diff[diff$significant == 'yes',]$gene_id
> diffGenes2 <- getGenes(cuff, diffGenes)  # does not count duplicates twice --
                                           #   good, since diffGenes *does* have duplicate gene names

# k-means clustering
> ic <- csCluster(diffGenes2, k=10)
> csClusterPlot(ic)

# dendrogram, MDS, PCA
> csDendro(diffGenes2, replicates=T)
> MDSplot(diffGenes2, replicates=T)
> PCAplot(diffGenes2, replicates=T)


# can try other filtering options
# e.g. both fpkm > 1:
> diffGenes3 <- diff[diff$value_1 > 1 & diff$value_2 > 1,]$gene_id
> diffGenes4 <- getGenes(cuff, diffGenes3)
# test status is OK: diff$status == 'OK'

# repFpkm(genes(cuff))$fpkm has fpkm for each replicate individually
> repFpkm <- repFpkm(genes(cuff))
> repFpkm[repFpkm$gene_id == 'Tnf',]

# can also get matrix of fpkm values for each replicate
> repFpkmMat <- repFpkmMatrix(genes(cuff))
# subset matrix: get just rows where more than 4 replicates have fpkm >= 50
> repFpkmMat[rowSums(repFpkmMat >= 50) > 4,]
# get gene names for rows where more than 5 replicates have fpkm > 1
> diffGenes5 <- row.names(repFpkmMat[rowSums(repFpkmMat > 1) > 5,])
> diffGenes6 <- getGenes(cuff, diffGenes5)


//////////////////////////////////////////////////////////////////////////

// heatmap
$ python ~/rutgers/makeHeatmap.py ../RNAseq3/fastq/cuffdiff_dedup2/gene_exp.diff zzz9 Control
$ R
> library(gplots)
> df <- read.csv("file.txt", sep="\t", row.names=1, check.names=FALSE)
> mat <- data.matrix(df)
> heatmap.2(mat, dendrogram="row", trace="none",
  col=redgreen(100), breaks=seq(-5, 5, 0.1),  # limit scale to [-5,5]
  density.info="none", key.title=NA, key.xlab="log2-fold change",
  cexCol=1, srtCol=0, adjCol=c(0.5, 0), offsetRow=0, Colv=F,
  main="Glorious heatmap", xlab="Treatment", ylab="Gene")
# consider decreasing the scale to maximize color diffs:
# col=redgreen(80), breaks=seq(-2, 2, 0.05)



//////////////////////////////////////////////////////////////////////////



// also generate read counts (not FPKM) for feeding to DESeq2 or edgeR
// in ~/RNAseq3/:
$ /home/john/tools/subread-1.5.1-source/bin/featureCounts \
  -a /home/john/genomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf \
  -p -o combraw.count \
  ../../RNAseq/results/C20/accepted_hits_dedup.bam C14_ca/accepted_hits_dedup.bam C15_ca/accepted_hits_dedup.bam C19_ca/accepted_hits_dedup.bam \
  ../../RNAseq/results/C34/accepted_hits_dedup.bam ../../RNAseq/results/C40/accepted_hits_dedup.bam C33_ca/accepted_hits_dedup.bam C36_ca/accepted_hits_dedup.bam \
  ../../RNAseq/results/C54/accepted_hits_dedup.bam ../../RNAseq/results/C60/accepted_hits_dedup.bam C55_ca/accepted_hits_dedup.bam C59_ca/accepted_hits_dedup.bam
$ python convert3.py combraw.count comb.count
// featureCounts: -M --fraction
// sorting bam by name does not affect featureCounts output

$ R
> library(edgeR)
> tab <- read.csv("comb.count", sep="\t", row.names="Symbol")
> group <- c(rep("Control", 4), rep("AOM-DSS", 4), rep("AOM-DSS-Cur", 4))
> batch <- as.factor(c(0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1))
> dge <- DGEList(counts=tab, group=group)
# dge$samples has metainfo, incl. "lib.size"
# dge$counts has actual count data
# optional arg: samples=c("Control_0", "Control_1", ...)

# filter out genes without at least 1 cpm in at least two samples
> keep <- rowSums(cpm(dge) > 1) >= 2  # cpm(dge) performs calculations: count / total * 10^6
> dge2 <- dge[keep, , keep.lib.sizes=F]  # last term means recalculate lib.sizes in dge2$samples
> nrow(dge$counts)
[1] 24391
> nrow(dge2$counts)
[1] 13097

# normalize, estimate dispersion
> dge3 <- calcNormFactors(dge2)
> dge4 <- estimateDisp(dge3)
Design matrix not provided. Switch to the classic mode.


> mat <- model.matrix(~group + 0 + batch)
> dge5 <- estimateDisp(dge3, mat)
> fit <- glmFit(dge5, mat)

> lrt.2vs1 <- glmLRT(fit, coef=2)
> lrt.3vs1 <- glmLRT(fit, coef=3)




> logCPM <- cpm(dge,log=TRUE,prior.count=5)
> batch <- as.factor(c(0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1))
> logCPM <- removeBatchEffect(logCPM, batch=batch)
> plotMDS(logCPM, col=as.numeric(batch))

removeBatchEffect():
     This function is not intended to be used prior to linear
     modelling. For linear modelling, it is better to include the batch
     factors in the linear model.

> dge4$counts[which(row.names(dge4$counts) == "Duox2"),]
 C20  C14  C15  C19  C34  C40  C33  C36  C54  C60  C55  C59 
 801  393 1280 1110 1703 2774 1778 2523  449  282  430  589 


////////////////////////////////////////////////////////////////////////

MDS plot (http://bioinf.wehi.edu.au/folders/smchd1/Smchd1.html)

#pdf("mdsplot.pdf",width=10,height=5)
par(mfrow=c(1,2))
plotMDS(vK,label=c(1,2,3,4,5,6),col=c( "darkolivegreen4", "darkolivegreen4","deepskyblue" ,"deepskyblue" ,"deepskyblue" ,"darkolivegreen4"),cex=2,main="Neural Stem Cell")
legend("topleft",legend=c("Wild Type","Smchd1-null"),text.col=c( "darkolivegreen4","deepskyblue"))
plotMDS(v,label=c(1,2,3,4,5,6,7),xlim=c(-2.5,5),cex=2,main="Lymphoma Cell",col=c("deepskyblue","deepskyblue","deepskyblue","deepskyblue","darkolivegreen4","darkolivegreen4","darkolivegreen4"))

////////////////////////////////////////////////////////////////////////

https://support.bioconductor.org/p/59195/
   library(edgeR)
   logCPM <- cpm(y,log=TRUE,prior.count=5)
   logCPM <- removeBatchEffect(logCPM, batch=batch)

////////////////////////////////////////////////////////////////////////


// can use R version of featureCounts:
> library(Rsubread)
> fc2 <- featureCounts(c("/home/john/RNAseq3/C14/accepted_hits_dedup.bam",
  "/home/john/RNAseq3/C15/accepted_hits_dedup.bam",
  "/home/john/RNAseq3/C19/accepted_hits_dedup.bam"),
  annot.ext="/home/john/genomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf", 
  isGTFAnnotationFile=T, isPairedEnd=T)
> fc2$targets <- c("C14", "C15", "C19")  # doesn't really fix sample names


// featureCounts can be given individual samples:
$ for fol in */accepted_hits_dedup.bam; do \
  /home/john/tools/subread-1.5.1-source/bin/featureCounts \
  -a /home/john/genomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf \
  -p -o ${fol%/*}raw.count $fol; \
  done
// also done in ~/RNAseq/ for the previous samples
// then, convert featureCounts output to edgeR-compatible two-column format:
$ python convert.py C14raw.count C14.count
// mv count files (including from ~/RNAseq/) to ~/RNAseq3/count/
$ R
> library(edgeR)
# load converted count files
> dge <- readDGE(
  files=c("C20.count", "C14.count", "C15.count", "C19.count",
    "C34.count", "C40.count", "C33.count", "C36.count",
    "C54.count", "C60.count", "C55.count", "C59.count"),
  group=c("Control", "Control", "Control", "Control",
    "AOM-DSS", "AOM-DSS", "AOM-DSS", "AOM-DSS",
    "AOM-DSS-Cur", "AOM-DSS-Cur", "AOM-DSS-Cur", "AOM-DSS-Cur"),
  labels=c("Control_0", "Control_1", "Control_2", "Control_3",
    "AOM-DSS_0", "AOM-DSS_1", "AOM-DSS_2", "AOM-DSS_3",
    "AOM-DSS-Cur_0", "AOM-DSS-Cur_1", "AOM-DSS-Cur_2", "AOM-DSS-Cur_3"))


////////////////////////////////////////////////////////////////////////

# trying (and failing) to write a 'distfun' for heatmap.2 that will
#   handle 'NA' values appropriately
> rdist.alt <- function(x) {
na.id <- is.na(x[1,]) | is.na(x[2,])
x1 <- x[1,][!na.id]
x2 <- x[2,][!na.id]
return(sqrt(sum((x1 - x2)^2)))
}
# still produces error with heatmap.2:
Error in if (is.na(n) || n > 65536L) stop("size cannot be NA nor exceed 65536") : 
  missing value where TRUE/FALSE needed


////////////////////////////////////////////////////////////////////////

// read counts
$ for file in fastq/*; do echo -n `basename $file`; echo -en "\t"; zcat $file | echo $((`wc -l`/4)); done
Kong_C2_S1_R1_001.fastq.gz	35147457
Kong_C2_S1_R2_001.fastq.gz	35147457
Kong_C14_S2_R1_001.fastq.gz	31870111
Kong_C14_S2_R2_001.fastq.gz	31870111
Kong_C15_S3_R1_001.fastq.gz	29916789
Kong_C15_S3_R2_001.fastq.gz	29916789
Kong_C19_S4_R1_001.fastq.gz	33451923
Kong_C19_S4_R2_001.fastq.gz	33451923
Kong_C33_S5_R1_001.fastq.gz	30099448
Kong_C33_S5_R2_001.fastq.gz	30099448
Kong_C36_S6_R1_001.fastq.gz	30766544
Kong_C36_S6_R2_001.fastq.gz	30766544
Kong_C55_S7_R1_001.fastq.gz	30392109
Kong_C55_S7_R2_001.fastq.gz	30392109
Kong_C59_S8_R1_001.fastq.gz	29130758
Kong_C59_S8_R2_001.fastq.gz	29130758

// cuffdiff without dedup
$ /home/john/tools/cufflinks-2.2.1.Linux_x86_64/cuffdiff -p 8 \
  -o cuffdiff --no-update-check \
  -L Control,AOM-DSS,AOM-DSS-Cur \
  /home/john/genomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf \
  ../RNAseq/results/C20/accepted_hits.bam,C14/accepted_hits.bam,C15/accepted_hits.bam,C19/accepted_hits.bam \
  ../RNAseq/results/C34/accepted_hits.bam,../RNAseq/results/C40/accepted_hits.bam,C33/accepted_hits.bam,C36/accepted_hits.bam \
  ../RNAseq/results/C54/accepted_hits.bam,../RNAseq/results/C60/accepted_hits.bam,C55/accepted_hits.bam,C59/accepted_hits.bam 



/////////////////////////////////////////////////////////////////////////

from Ron:
I was wondering if we could ask Stephen to compare Tophat analysis of RNAseq
data using the "official" RUCDR settings (--library-type fr-secondstrand
--inner-dist-mean 200 and -b2-sensitive) to the default settings? 

// try --library-type fr-secondstrand
$ tophat --library-type fr-secondstrand \
  --transcriptome-index /home/john/genomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes/genes \
  -o C14_test \
  /home/john/genomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome \
  C14_R1_ca.fastq.gz C14_R2_ca.fastq.gz

// slightly fewer alignments (but that doesn't mean it's wrong):
$ cat C14_test/align_summary.txt 
Left reads:
          Input     :  31859283
           Mapped   :  31048130 (97.5% of input)
            of these:  13242039 (42.7%) have multiple alignments (13560 have >20)
Right reads:
          Input     :  31859283
           Mapped   :  30679337 (96.3% of input)
            of these:  13084723 (42.6%) have multiple alignments (13574 have >20)
96.9% overall read mapping rate.

Aligned pairs:  30198008
     of these:  12880226 (42.7%) have multiple alignments
                  756330 ( 2.5%) are discordant alignments
92.4% concordant pair alignment rate.

// ... compared to previous:
$ cat C14_ca/align_summary.txt 
Left reads:
          Input     :  31859283
           Mapped   :  31051368 (97.5% of input)
            of these:  13243578 (42.7%) have multiple alignments (13557 have >20)
Right reads:
          Input     :  31859283
           Mapped   :  30683325 (96.3% of input)
            of these:  13085849 (42.6%) have multiple alignments (13571 have >20)
96.9% overall read mapping rate.

Aligned pairs:  30204363
     of these:  12880680 (42.6%) have multiple alignments
                  731501 ( 2.4%) are discordant alignments
92.5% concordant pair alignment rate.



