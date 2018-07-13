library(DESeq2)
ls()
tab <- read.csv("comb.count", sep="\t", row.names="Symbol")
head(tab)
dim(tab)
help(matrix)
mat <- matrix(c(rep("Control", 4), rep("AOM-DSS", 4), rep("AOM-DSS-Cur", 4)))
mat
mat <- matrix(c(rep("Control", 4), rep("AOM-DSS", 4), rep("AOM-DSS-Cur", 4)), nrow=12, ncol=2)
mat
mat <- matrix(c(rep("Control", 4), rep("AOM-DSS", 4), rep("AOM-DSS-Cur", 4), 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1), nrow=12, ncol=2)
mat
help(matrix)
mat <- matrix(c(rep("Control", 4), rep("AOM-DSS", 4), rep("AOM-DSS-Cur", 4), 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1), nrow=12, ncol=2, dimnames=c(c("condition", "batch"), c("Control_0", "Control_1", "Control_2", "Control_3", "AOM-DSS_0", "AOM-DSS_1", "AOM-DSS_2", "AOM-DSS_3", "AOM-DSS-Cur_0", "AOM-DSS-Cur_1", "AOM-DSS_Cur_2", "AOM-DSS-Cur_3"))
)
mat <- matrix(c(rep("Control", 4), rep("AOM-DSS", 4), rep("AOM-DSS-Cur", 4), 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1), nrow=12, ncol=2, dimnames=list(c("condition", "batch"), c("Control_0", "Control_1", "Control_2", "Control_3", "AOM-DSS_0", "AOM-DSS_1", "AOM-DSS_2", "AOM-DSS_3", "AOM-DSS-Cur_0", "AOM-DSS-Cur_1", "AOM-DSS_Cur_2", "AOM-DSS-Cur_3"))
)
help(matrix)
mat <- matrix(c(rep("Control", 4), rep("AOM-DSS", 4), rep("AOM-DSS-Cur", 4), 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1), nrow=12, ncol=2, dimnames=list(("condition", "batch"), ("Control_0", "Control_1", "Control_2", "Control_3", "AOM-DSS_0", "AOM-DSS_1", "AOM-DSS_2", "AOM-DSS_3", "AOM-DSS-Cur_0", "AOM-DSS-Cur_1", "AOM-DSS_Cur_2", "AOM-DSS-Cur_3"))
mat <- matrix(c(rep("Control", 4), rep("AOM-DSS", 4), rep("AOM-DSS-Cur", 4), 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1), nrow=12, ncol=2, dimnames=list("condition", "batch", "Control_0", "Control_1", "Control_2", "Control_3", "AOM-DSS_0", "AOM-DSS_1", "AOM-DSS_2", "AOM-DSS_3", "AOM-DSS-Cur_0", "AOM-DSS-Cur_1", "AOM-DSS_Cur_2", "AOM-DSS-Cur_3"))
help(list)
mat
colnames(mat) <- c("condition", "batch")
mat
mat <- matrix(c(rep("Control", 4), rep("AOM-DSS", 4), rep("AOM-DSS-Cur", 4), 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1), nrow=12, ncol=2)
colnames(mat) <- c("condition", "batch")
mat
ls()
paste("Control_", rep(0:3))
row.names(mat) <- c(paste("Control_", rep(0:3)), paste("AOM-DSS_", rep(0:3)), paste("AOM-DSS-Cur_", rep(0:3))
)
mat
row.names(mat) <- c(paste("Control_", rep(0:3), sep=""), paste("AOM-DSS_", rep(0:3), sep=""), paste("AOM-DSS-Cur_", rep(0:3), sep=""))
mat
mat[1,1]
taba
tab
head(tab)
row.names(mat) <- col.names(tab)
colnames(tab)
row.names(mat) <- colnames(tab)
mat
rep(0:3)
help(paste)
dds <- DESeqDataSetFromMatrix(countData=tab, colData=mat)
dds <- DESeqDataSetFromMatrix(countData=tab, colData=mat, design=~condition)
dds
dds$condition
dds$condition <- relevel(dds$condition, ref="Control")
dds$condition
dds
help(as.factor)
help(DESeqDataSetFromMatrix)
dds <- DESeqDataSetFromMatrix(countData=tab, colData=mat, design= ~ batch + condition)
dds$condition <- relevel(dds$condition, ref="Control")
dds
dds$condition
counts(dds)
head(counts(dds))
head(rowSums(counts(dds)))
head(rowSums(counts(dds))>1)
head(rowSums(counts(dds)>1))
dds2 <- dds[ rowSums(counts(dds) >= 1) >= 2, ]
dds
dds2
counts(dds2)
head(counts(dds2))
dds[! rowSums(counts(dds) >= 1) >= 2, ]
head(counts(dds[! rowSums(counts(dds) >= 1) >= 2, ]))
head(counts(dds[! rowSums(counts(dds) >= 1) >= 2, ]), n=20)
dds
dds2
dds <- DESeq(dds)
res <- results(dds)
res
dds
res$pvalue
rowSum(res$pvalue < 0.05)
rowSums(res$pvalue < 0.05)
sum(res$pvalue < 0.05)
res$pvalue < 0.05
(res$pvalue < 0.05) == TRUE
count(res$pvalue < 0.05)
sum(res$pvalue < 0.05)
sum(res$pvalue != 'NA' & res$pvalue < 0.05)
sum(res$pvalue == 'NA')
sum(is.na(res$pvalue))
sum(!is.na(res$pvalue) & res$pvalue < 0.05)
sum(!is.na(res$pvalue) & res$pvalue >= 0.05)
length(res$pvalue)
sum(!is.na(res$pvalue) && res$pvalue < 0.05)
head(res)
summary(res)
sum(res$padj < 0.1, na.rm=T)
sum(res$padj < 0.1)
sum(!is.na(res$padj) & res$padj < 0.1)
length(p$adj)
length(res$padj)
sum(is.na(res$padj))
res
resOrdered <- res[order(res$padj),]
resOrdered
order(p$adj)
order(res$padj)
head(order(res$padj))
which(order(res$padj) == 1)
res[2565]
res[2565,]
res[2564,]
res[2566,]
resOrdered
head(resOrdered, n=20)
summary(res)
summary(resOrdered)
help(results)
resultsNames(dds)
help(results)
res2 <- results(dds, contrast=c("conditionAOM.DSS.Cur", "conditionControl"))
res2 <- results(dds, contrast=c("AOM.DSS.Cur", "Control"))
res2 <- results(dds, contrast=("AOM.DSS.Cur", "Control"))
q()
library(DESeq2)
ls()
tab <- read.csv('comb.count', sep="\t", row.names="Symbol")
dim(tab)
mat
q()
library(DESeq2)
tab <- read.csv("comb.count", sep="\t", row.names="Symbol")
dim(tab)
head(tab)
tab$C36 <- NULL
head(tab)
mat <- matrix(c(rep("Control", 4), rep("AOM-DSS", 3), rep("AOM-DSS-Cur", 4),
0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1), nrow=11, ncol=2)
mat
colnames(mat) <- c("condition", "batch")
row.names(mat) <- colnames(tab)
mat
dds <- DESeqDataSetFromMatrix(countData=tab, colData=mat, design= ~ batch + condition)
dds$condition <- relevel(dds$condition, ref="Control")
dds
ls()
dds <- dds[ rowSums(counts(dds) >= 1) >= 2, ]
dds
dds <- DESeq(dds)
res <- results(dds)
summary(res)
res2 <- results(dds, contrast=c("condition", "AOM-DSS-Cur", "Control"))
summary(res2)
res2 <- results(dds, contrast=c("condition", "AOM-DSS", "Control"))
summary(res2)
ls()
resOrdered
res3 <- results(dds, contrast=c("condition", "AOM-DSS-Cur", "AOM-DSS2"))
res3 <- results(dds, contrast=c("condition", "AOM-DSS-Cur", "AOM-DSS"))
summary(res3)
res
dim(res)
lenth(res$pvalue)
length(res$pvalue)
sum(res$padj < 0.1)
sum(res$padj < 0.1, na.rm=T)
ls)
ls()
rm(dds2)
resOrdered
rm(resOrdered)
ls()
q()
ls()
mat
res
res2
res3
res
help(results)
dds
resultsNames(dds)
library(DESeq2)
resultsNames(dds)
summary(res)
dds$batch
dds$condition
dds$sizeFactor
help(results)
summary(res)
sum(res$padj < 0.1, na.rm=T)
sum(res$padj < 0.1 & res$log2FoldChange > 1, na.rm=T)
sum(res$padj < 0.1 & res$log2FoldChange < -1, na.rm=T)
resSig <- subset(res, padj < 0.1)
resSig
sum(resSig$log2FoldChange < -1)
sum(resSig$log2FoldChange > 1)
sum(resSig$log2FoldChange <= 1 & resSig$log2FoldChange >= -1)
sum(resSig$pvalue < 0.05)
resSig[order(resSig$log2FoldChange), ]
type(resSig)
class(resSig)
resSig
max(resSig$padj)
resSig[order(resSig$padj), ]
resSig <- resSig[order(resSig$log2FoldChange), ]
resSig
write.table(resSig)
help(write.table)
resSig.all
ttf.all
ttf
format(resSig, digits=6)
dim(res)
dim(resSig)
write.csv(resSig, file="zzz.csv", quote=F)
write.csv(resSig, file="zzz.csv", sep="\t", quote=F)
write.table(resSig, file="zzz.csv", sep="\t", quote=F)
help(write.table)
write.table(resSig, file="zzz.csv", sep="\t", quote=F, row.names=T)
write.table(resSig, file="zzz.csv", sep="\t", quote=F, row.names=F)
help(write.table)
write.table(resSig, file="zzz.csv", sep="\t", quote=F, col.names=NA)
plotCounts(dds, gene='Tnf')
resSig
plotMA(res)
res2
res2Sig <- subset(res2, padj < 0.1)
res2Sig
res2Sig <- res2Sig[order(res2Sig$log2FoldChange),
]
res2Sig
write.table(res2Sig, file="AOM-DSS_Control.csv", sep="\t", quote=F, col.names=NA)
res3
res3Sig <- subset(res3, padj < 0.1)
res3Sig <- res3Sig[order(res3Sig$log2FoldChange),]
res3Sig
write.table(res3Sig, file="AOM-DSS-Cur_AOM-DSS.csv", sep="\t", quote=F, col.names=NA)
ls()
q()
ls()
tab
ls()
mat
dim(tab)
head(tab)
dds_18 <- DESeqDataSetFromMatrix(countData=tab, colData=mat, design= ~ batch + condition)
library(DESeq2)
dds_18 <- DESeqDataSetFromMatrix(countData=tab, colData=mat, design= ~ batch + condition)
mat
dds_18 <- DESeqDataSetFromMatrix(countData=tab, colData=mat, design= ~ batch + condition)
dds_18$condition <- relevel(dds_18$condition, ref="Control")
dds_18 <- DESeq(dds_18)
head(fpkm(dds_18))
mcols(dds_18)$basepairs <- tab$Length
head(fpkm(dds_18))
head(dds_18)
dds_18 <- DESeqDataSetFromMatrix(countData=tab, colData=mat, design= ~ batch + condition)
mcols(dds_18)$basepairs <- tab$Length
tab2 <- read.csv('comb.count', sep='\t', row.names='Symbol')
head(tab2)
mcols(dds_18)$basepairs <- tab2$Length
head(fpkm(dds_18))
rowMeans(fpkm(dds_18['Tnf', c('C20', 'C14', 'C15', 'C19')]))
fpkm(dds_18['Tnf', c('C20', 'C14', 'C15', 'C19')])
fpkm(dds_18['Tnf',])
dds_18['Tnf',]
rowMeans(fpkm(dds_18[, c('C20', 'C14', 'C15', 'C19')]))
rowMeans(fpkm(dds_18[, c('C20', 'C14', 'C15', 'C19')]))$Tnf
fpkmControl18 <- rowMeans(fpkm(dds_18[, c('C20', 'C14', 'C15', 'C19')]))
head(fpkmControl18)
fpkmControl18[1]
class(fpkmControl18[1])
fpkmControl18[1][1]
fpkm_18 <- data.frame(fpkm(dds_18))
fpkm_18$Control <- rowMeans(fpkm(dds_18[, c('C20', 'C14', 'C15', 'C19')]))
fpkm_18$AOMDSS <- rowMeans(fpkm(dds_18[, c('C34', 'C40', 'C33')]))
fpkm_18$AOMDSSCur <- rowMeans(fpkm(dds_18[, c('C54', 'C60', 'C55', 'C59')]))
head(fpkm_18)
fpkm_18['Tnf',]
'Tnf'
numeric('Tnf')
ls()
q()
library(DESeq2)
ls()
mat
dds_18
dds_18$condition
dds_18$batch
res
help(results)
dim(res)
ls()
res2
res3
res2
sum(abs(res2$log2FoldChange) > 1 & res2$padj < 0.05)
sum(abs(res2$log2FoldChange) > 1 & res2$padj < 0.05, na.rm=T)
sum(res2$log2FoldChange > 1 & res2$padj < 0.05, na.rm=T)
sum(res2$log2FoldChange < -1 & res2$padj < 0.05, na.rm=T)
head(res2[which(res2$log2FoldChange > 1 & res2$padj < 0.05),])
dim(res2[which(res2$log2FoldChange > 1 & res2$padj < 0.05),])
rownames(res2[which(res2$log2FoldChange > 1 & res2$padj < 0.05),])
aomdss.up <- rownames(res2[which(res2$log2FoldChange > 1 & res2$padj < 0.05),])
aomdss.up
aomdss.down <- rownames(res2[which(res2$log2FoldChange < -1 & res2$padj < 0.05),])
aomdss.down
identical(rownames(res2), rownames(res))
identical(rownames(res2), rownames(res3))
write.table('geneNames.txt', rownames(res2), sep='\t')
head(rownames(res2))
write.table(rownames(res2), 'geneNames.txt', sep='\t')
write.table(rownames(res2), 'geneNames.txt', sep='\t', quote=F, row.names=F)
aomdss.up
aomdsscur.down <- rownames(res[which(res$log2FoldChange < -1 & res$padj < 0.05),])
aomdsscur.down
aomdsscur.down <- rownames(res[which(res$log2FoldChange > 1 & res$padj < 0.05),])
aomdsscur.down <- rownames(res[which(res$log2FoldChange < -1 & res$padj < 0.05),])
aomdsscur.up <- rownames(res[which(res$log2FoldChange > 1 & res$padj < 0.05),])
aomdsscur.up
aomdss.aomdsscur.up <- rownames(res[which(res3$log2FoldChange > 1 & res3$padj < 0.05),])
aomdss.aomdsscur.down <- rownames(res[which(res3$log2FoldChange < -1 & res3$padj < 0.05),])
aomdss.aomdsscur.up
aomdss.aomdsscur.down
aomdss.up
head(res2)
res2['Tnf
',]
res2['Tnf',]
res['Tnf',]
res3['Tnf',]
mat
res2
res2@metadata
dds
colnames(dds)
head(tab)
ls()
fpkm_18
head(fpkm_18)
fpkm_18['Tnf',]
res3['Tnf',]
res2['Tnf',]
res['Tnf',]
head(rownames(res2))
which(rownames(res2) == 'Tnf')
res2[7099,]
aomdss.up
genes
genes <- read.csv('../../MethylSeq/genes.txt', sep='\t', row.names=F, header=F)
genes <- read.csv('../../MethylSeq/genes.txt', sep='\t', header=F)
head(genes)
dim(genes)
intersect(aomdss.up, genes$V1)
intersect(aomdss.down, genes$V1)
add <- read.table('../../MethylSeq/add.txt', header=F)
add
add <- read.table('../../MethylSeq/add.txt', header=T)
add
intersect(add, aomdss.up)
intersect(add, aomdss.down)
adu <- read.table('../../MethylSeq/adu.txt', header=T)
adu
intersect(adu, 'Mgat3')
intersect(add$x, aomdss.down)
intersect(add$x, aomdss.up)
intersect(adu$x, aomdss.up)
intersect(adu$x, aomdss.down)
aomdss.down
adu$x
intersect(adu$x, c('FBXw7', 'Rint1', 'Lrtm2'))
intersect(adu, c('FBXw7', 'Rint1', 'Lrtm2'))
adno <- read.table('../../MethylSeq/adno.txt', header=T)
head(adno)
dim(adno)
intersect(adno, aomdss.up)
intersect(adno$x, aomdss.up)
intersect(adno$x, aomdss.down)
ls()
res2['Tnf',]
aomdss.up
res2['Arid3a',]
curu <- read.table('../../MethylSeq/curu.txt', header=T)
curu
intersect(curu$x, aomdsscur.up)
intersect(curu$x, aomdsscur.down)
curd <- read.table('../../MethylSeq/curd.txt', header=T)
curd
intersect(curd$x, aomdsscur.up)
intersect(curd$x, aomdsscur.down)
curno <- read.table('../../MethylSeq/curno.txt', header=T)
head(curno)
dim(curno)
intersect(curno$x, aomdsscur.down)
intersect(curno$x, aomdsscur.up)
intersect(aomdss.down, genes)
head(genes)
intersect(aomdss.down, genes$V1)
intersect(aomdss.up, genes$V1)
intersect(aomdsscur.down, genes$V1)
intersect(aomdsscur.up, genes$V1)
intersect(aomdss.aomdsscur.down, genes$V1)
intersect(aomdss.aomdsscur.up, genes$V1)
doubu <- read.table('../../MethylSeq/doubu.txt', header=T)
doubd <- read.table('../../MethylSeq/doubd.txt', header=T)
doubno <- read.table('../../MethylSeq/doubno.txt', header=T)
doubu
intersect(doubu$x, aomdsscur.up)
intersect(doubu$x, aomdss.aomdsscur.up)
intersect(doubu$x, aomdss.aomdsscur.down)
intersect(doubu$x, genes$V1)
intersect(doubd$x, aomdss.aomdsscur.up)
intersect(doubd$x, aomdss.aomdsscur.down)
intersect(doubd$x, genes$V1)
intersect(doubno$x, aomdss.aomdsscur.up)
intersect(doubno$x, aomdss.aomdsscur.down)
intersect(doubno$x, genes$V1)
length(genes$V1)
length(doubno$x)
intersect(doubno$x, union(aomdss.aomdsscur.down, aomdss.aomdsscur.up))
setdiff(doubno$x, union(aomdss.aomdsscur.down, aomdss.aomdsscur.up))
head(res, n=2
)
head(res3, n=2)
aomdss.aomdsscur.down2 <- rownames(res3[which(res3$log2FoldChange < 1 & res3$padj < 0.05),])
length(aomdss.aomdsscur.down)
length(aomdss.aomdsscur.down2)
head(aomdss.aomdsscur.down)
head(aomdss.aomdsscur.down2)
dim(res)
dim(res3)
aomdss.aomdsscur.down2 <- rownames(res3[which(res3$log2FoldChange < -1 & res3$padj < 0.05),])
identical(aomdss.aomdsscur.down, aomdss.aomdsscur.down2)
sum(res3[which(res3$log2FoldChange < -1 & res3$padj < 0.05, na.rm=T)
sum(res3$log2FoldChange < -1 & res3$padj < 0.05, na.rm=T)
sum(res3$log2FoldChange > 1 & res3$padj < 0.05, na.rm=T)
sum(res3$padj < 0.1, na.rm=T)
sum(res$padj < 0.1, na.rm=T)
rownames(res2[which(res2$padj < 0.1),])
rownames(res2[which(res2$padj < 0.1),])
ls()
q()
View(tab)
View(tab)
View(tab2)
View(fpkm_18)
