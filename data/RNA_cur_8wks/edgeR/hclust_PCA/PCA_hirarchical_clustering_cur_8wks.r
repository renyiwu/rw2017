#plot hirarchical clustering graph and PCA -- principal component analysis. Using cur-8wks data as an example.
#R Wu 1-30-2018. Tested on SOP-1163 Ubuntu 16.04

#1 Hirarchical clustering
library(ggplot2)
dt <- read.table("data/RNA_cur_8wks/edgeR/rpkm_all_cur_8wks_edgeR_1-30-2018.csv",
                  header = T, sep = "\t", row.names = "Geneid")
dt1 <- dt[-1] 
colnames(dt1) <- substring(colnames(dt1), 1 ,3)

dt2 <- dt1[!rowSums(dt1 < 1),] #keep rows that have values greater than 1 in each cell.
#

#colnames(dt2) <- c("TPA_1", "mITC", "Control", "TPA_0", "CA", "FX", "CDDO")
dtt <- t(dt2)
ds <- dist(dtt)
h <- hclust(ds)
plot(h, xlab = "sample", asp = 2)
# Parameters for plotting
plot(h, hang = -1, xlab = "sample", ylab = "distance", cex = 1.2, lwd = 2, col = "black")
#rect.hclust(h, 5)
dev.off()





#2 PCA analysis
#
library(ggfortify)

dt <- read.table("data/RNA_cur_8wks/edgeR/rpkm_all_cur_8wks_edgeR_1-30-2018.csv",
                 header = T, sep = "\t", row.names = "Geneid")
dt1 <- dt[-1] 
colnames(dt1) <- substring(colnames(dt1), 1 ,3)

dt2 <- dt1[!rowSums(dt1 < 1),] #keep rows that have values greater than 1 in each cell.
dtt <- t(dt2)
#ncol(dtt)
dtp1 <- data.frame(dtt)
#class(dtp1)
dtp2 <- dtp1
dtp2$Sample <- rownames(dtp2)
dtp2$group <- c("Con", "AOM_DSS", "AOM_DSS", "Con", "AOM_DSS_Cur", "AOM_DSS_Cur", "DSS", "DSS", "DSS_Cur", "DSS_Cur")
#or

prcomp(dtp1)
summary(prcomp(dtp1))
# > summary(prcomp(dtp1))
# Importance of components:
#         PC1       PC2       PC3      PC4       PC5       PC6       PC7       PC8       PC9
# Standard deviation     4020.796 1.473e+03 1.309e+03 812.2966 693.26968 576.17983 451.46135 378.38388 353.16331
# Proportion of Variance    0.735 9.867e-02 7.792e-02   0.0300   0.02185   0.01509   0.00927   0.00651   0.00567
# Cumulative Proportion     0.735 8.337e-01 9.116e-01   0.9416   0.96346   0.97855   0.98782   0.99433   1.00000

autoplot(prcomp(dttp1), data = dtp2, colour = 'Sample')
autoplot(prcomp(dtp1), data = dtp2, colour = 'group')
autoplot(prcomp(dtp1), data = dtp2, colour = 'group', label = TRUE, label.size = 3)
autoplot(prcomp(dtp1), data = dtp2, colour = 'group', shape = FALSE, label.size = 3) #Labels from rownames(dtp2)
autoplot(prcomp(dtp1), data = dtp2, colour = 'group', loadings = TRUE)
dev.off()
