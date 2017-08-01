library(ggplot2)
dt1 <- read.table("data/shan_rna/RW_all.fpkm.csv", header = T, sep = "\t", row.names = "X")
dt2 <- dt1[!rowSums(dt1 < 1),] #keep rows that have values greater than 1 in each cell.
#
#Hirarchical clustering
#colnames(dt2) <- c("TPA_1", "mITC", "Control", "TPA_0", "CA", "FX", "CDDO")
dtt <- t(dt2)
ds <- dist(dtt)
h <- hclust(ds)
plot(h, xlab = "sample", asp = 2)
# Parameters for plotting
plot(h, hang = -1, xlab = "sample", ylab = "distance", cex = 1.2, lwd = 2, col = "black")
rect.hclust(h, 3)
dev.off()
#
#PCA analysis
#
library(ggfortify)
#ncol(dtt)
dtp1 <- data.frame(dtt)
#class(dtp1)
dtp2 <- dtp1
dtp2$Sample <- rownames(dtp2)
#or
dtp2$Sample <- c("TPA_1", "mITC", "Control", "TPA_0", "CA", "FX", "CDDO")
rownames(dtp2) <- c("TPA_1", "mITC", "Control", "TPA_0", "CA", "FX", "CDDO")
prcomp(dtp1)
autoplot(prcomp(dtp1), data = dtp2, colour = 'Sample')
autoplot(prcomp(dtp1), data = dtp2, colour = 'Samples', label = TRUE, label.size = 3)
autoplot(prcomp(dtp1), data = dtp2, colour = 'Samples', shape = FALSE, label.size = 3) #Labels from rownames(dtp2)
autoplot(prcomp(dtp1), data = dtp2, colour = 'Samples', loadings = TRUE)
