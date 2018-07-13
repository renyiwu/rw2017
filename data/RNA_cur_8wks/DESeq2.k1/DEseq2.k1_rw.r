# DEseq2 practice, R Wu. April 2018

# source("https://bioconductor.org/biocLite.R")
# biocLite("rnaseqGene")

# Load libraries.
library("DESeq2")
library("data.table")

# Load data from csv file.
dt <- fread("data/RNA_cur_8wks/DESeq2.k1/RNA_cur_8wks_k1_picardup_primary.csv")
dt <- dt[,-(2:5)]

# Modify colnames and reorder.
colnames(dt) <- c("Geneid","Length", paste("C", c("01",26,29,"02",42,46,65,70,75,79), sep = ""))
dt1 <- dt[, c("Geneid","Length", paste("C", c("01", "02", 26,29,42,46,65,70,75,79), sep = ""))]

# Format counts.
cts <- as.matrix(dt1[,3:12])
rownames(cts) <- dt1$Geneid

# Format meta data.
coldata <- data.frame(condition = rep(c("Control", "DSS", "DSS+Cur.", "AOM+DSS", "AOM+DSS+Cur."), each = 2),
                      row.names = colnames(cts))
# Check 
head(cts, 2)
coldata
all(colnames(cts) %in% rownames(coldata)) #TRUE
all(colnames(cts) == rownames(coldata)) #TRUE

# Construct DESeq data table.
dds <- DESeqDataSetFromMatrix(cts, coldata, ~ condition)
mcols(dds) <- DataFrame(mcols(dds), basepairs = dt1$Length) # Add more info for RPKM calculation.
# mcols(dds) <- NULL

# Pre-filtering. Method one:
keep <- rowSums(counts(dds)) >= 10 # Keep rows that have at least 10 reads total.
dds <- dds[keep,]

# Prefiltering. Method two:
# dds2 <- dds[rowSums(counts(dds) >= 1) >= 2, ]
# keep genes with at least one count in at least two samples.


# re-order factor levels. This is not necessary if one specifies the groups 
# on which the results command will perform the calculation latter.
dds$condition <- factor (dds$condition, levels = c("Control", "DSS", "DSS+Cur.", "AOM+DSS", "AOM+DSS+Cur."))
# or use relevel
# dds$condition <- relevel(dds$condition, ref = "Control") # in this case, only one factor was re-ordered)

# If some samples were removed later, use this to remove the factors also.
# dds$condition <- droplevels(dds$condition)

# Run DESeq
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds) # be default, the last one in 'condition' vs the first one.
results(dds, contrast = c("condition", "AOM+DSS+Cur.", "Control")) # specify which vd which.

# Reorder output by p, padj, or other.
res[order(res$pvalue, res$log2FoldChange),]

summary(res)

# Check how many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm = TRUE)
# 128

# Change the FDR cutoff, alpha. If not set, it is 0.1 by default.
res5 <- results(dds, alpha = 0.05)

summary(res5)
# out of 15347 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 31, 0.2% 
# LFC < 0 (down)   : 30, 0.2% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 4463, 29% 
# (mean count < 34)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

sum(res5$padj < 0.05, na.rm = TRUE)
# 129

## 2 Exploring and exporting results

# 2.1 MA plot. M is minus. A is average. (logA-logB)/(loA+logB)
plotMA(res, ylim = c(-2,2))

# log2 moderation
resLFC <- lfcShrink(dds, coef = 2)
plotMA(resLFC, ylim = c(-2,2)) #more useful

# 2.2 Plot counts
plotCounts(dds, gene = which.min(res$padj), intgroup = "condition", col = rep(1:5, each = 2))

plotCounts(dds, gene = which.max(res$log2FoldChange), intgroup = "condition", col = rep(1:5, each = 2))

plotCounts(dds, gene = "Nfe2l2", intgroup = "condition", col = rep(1:5, each = 2), xlab = "Groups")

abline(v = 3, h = 2000)
dev.off()

# Or save the output so it can be plotted by another function. 
d <- plotCounts(dds, gene = "Nfe2l2", intgroup = "condition", returnData = TRUE)

library("ggplot2")
ggplot(d, aes(x=condition, y=count)) +
        geom_point(position = position_jitter(w = 0.1, h = 0)) +
        scale_y_log10(breaks = c(25, 100, 400)) # looks terrible.

mcols(res)$description


## 3. Write results to files
write.table(as.data.frame(res), "data/RNA_cur_8wks/DESeq2.k1/k1.AOM+DSS+Cur_Control.csv",
          sep = "\t",
          quote = F,
          col.names = NA) # leave the first cell blank and save row names.

# use subset function to filter resutls
# write.table(as.data.frame(subset(res, padj < 0.1)), "data/file.csv")

## 4 Multi-factor design.


## 5. Extract transformed values.
vsd <- vst(dds, blind = F)
rld <- rlog(dds, blind = F)
head(assay(vsd), 3)
head(assay(rld), 3)
head(assay(dds), 3)

# Check effects of transformation on the variance
ntd <- normTransform(dds)

library("vsn")
meanSdPlot((assay(ntd)))
meanSdPlot((assay(vsd)))
meanSdPlot((assay(rld)))

install.packages("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:30]
df <- as.data.frame(colData(dds)[,"condition"])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE)



# Compare sample wise distance.
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



# Hierarchical clustering
h <- hclust(sampleDists)
plot(h)
plot(h, xlab = "sample", asp = 2)
# Parameters for plotting
plot(h, hang = -1, xlab = "sample", ylab = "distance", cex = 1.2, lwd = 2, col = "black")
#rect.hclust(h, 5)
dev.off()

# PCA plot
plotPCA(vsd, intgroup = "condition")

# Plot with ggplot
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) + # Add shape = condition if desired.
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        coord_fixed()


# Plot dispersion
 plotDispEsts(dds)
 
 # Access dispersion
 dispersions(dds)
 head(dispersions(dds)[order(-dispersions(dds))], 50)

 head(mcols(dds)$dispersion)

 head(coef(dds))

 
 # Calculate fpkm 
 
 fpkm <- fpkm(dds)
 
 write.table(fpkm, "data/RNA_cur_8wks/DESeq2.k1/RNA_cur_8wks_k1_fpkm.csv",
             sep = "\t",
             quote = F,
             col.names = NA)
 
 
# to get the gene by clicking on the plot.  # Keeps crashing...
# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]

#

# countdata <- as.matrix(countdata)
# coldata <- data.frame(group = c("Control", "AOM+DSS", "AOM+DSS", "Control",
                                # "AOM+DSS+Cur.", "AOM+DSS+Cur.", "DSS", "DSS",
                                # "DSS+Cur.", "DSS+Cur."))
#rep(c("a","b"), each = 3, 2) >>> aaabbbaaabbb





 