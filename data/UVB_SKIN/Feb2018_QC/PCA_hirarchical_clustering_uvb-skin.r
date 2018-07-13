#plot hirarchical clustering graph and PCA -- principal component analysis. Using cur-8wks data as an example.
#R Wu 1-30-2018. Tested on SOP-1163 Ubuntu 16.04

#1 Hirarchical clustering
library(ggplot2)
dt <- read.table("data/uvb-skin/combined_dmr_fraction_s32_uvb-skin_renyi_02092018.csv",
                  header = T, sep = "\t", row.names = "Geneid")

dt <- read.table("data/uvb-skin/combined_dmr_fraction_s32_uvb-skin_renyi_02092018.csv",
                 header = T, sep = "\t")
dt1 <- dt[-(1:4)]

colnames(dt1) <- substring(colnames(dt1), 2 ,10)

# dt2 <- dt1[!rowSums(dt1 < 1),] #keep rows that have values greater than 1 in each cell.
#

#colnames(dt2) <- c("TPA_1", "mITC", "Control", "TPA_0", "CA", "FX", "CDDO")
dtt <- t(dt1)
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

dt <- read.table("data/uvb-skin/rpkm_all_uvb-skin_edgeR_2-9-2018.csv",
                 header = T, sep = "\t", row.names = "Geneid")

dt <- read.table("data/uvb-skin/combined_dmr_fraction_s32_uvb-skin_renyi_02092018.csv",
                 header = T, sep = "\t")

# dt1 <- dt[-c(1:4, 5:12, 13:20, 21:28, 29:36)]  #dt1 <- dt[-c(1, 10:17, 18:33)]  #1, 2:9, 10:17, 18:25, 26:33
dt1 <- dt[-c(1:4, 5:12, 13:20, 21:28)]
#
dt1 <- na.omit(dt1)

dtt <- t(dt1)
class(dtt)
dtp1 <- data.frame(dtt)
class(dtp1)

dtp2 <- dtp1

dtp2$Sample <- rownames(dtp2)

dtp2$Group <- substring(rownames(dtp2), 2, 8)
summary(prcomp(dtp1))

# autoplot(prcomp(dtp1), data = dtp2, colour = 'Sample')
autoplot(prcomp(dtp1), data = dtp2, colour = 'Group') # pca-on-dmr-2wks-g

autoplot(prcomp(dtp1), data = dtp2, colour = 'Group', label = TRUE, label.size = 3)
autoplot(prcomp(dtp1), data = dtp2, colour = 'Group', shape = FALSE, label.size = 3) #Labels from rownames(dtp2)
autoplot(prcomp(dtp1), data = dtp2, colour = 'Group', loadings = TRUE)
dev.off()


##
dt <- read.table("data/uvb-skin/rpkm_all_uvb-skin_edgeR_2-9-2018.csv",
                 header = T, sep = "\t", row.names = "Geneid")
colnames(dt)
dt1 <- dt[-1]
#dt1 <- dt[-c(1, 10:17, 18:33)]  #1, 2:9, 10:17, 18:25, 26:33
dtt <- t(dt1)
dtp1 <- data.frame(dtt)
#class(dtp1)
dtp2 <- dtp1
dtp2$Sample <- rownames(dtp2)
dtp2$Group <- substring(rownames(dtp2), 2, 8)

# > summary(prcomp(dtp1))
# Importance of components:
#         PC1       PC2       PC3      PC4       PC5       PC6       PC7       PC8       PC9
# Standard deviation     4020.796 1.473e+03 1.309e+03 812.2966 693.26968 576.17983 451.46135 378.38388 353.16331
# Proportion of Variance    0.735 9.867e-02 7.792e-02   0.0300   0.02185   0.01509   0.00927   0.00651   0.00567
# Cumulative Proportion     0.735 8.337e-01 9.116e-01   0.9416   0.96346   0.97855   0.98782   0.99433   1.00000
autoplot(prcomp(dtp1), data = dtp2, colour = 'Group', label = TRUE, label.size = 3)
autoplot(prcomp(dtp1), data = dtp2, colour = 'Sample')
autoplot(prcomp(dtp1), data = dtp2, colour = 'Group')
##
