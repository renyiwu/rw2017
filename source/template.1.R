# R scrpits originally by John. 
#setwd("C:/Users/renyi/Google Drive/WORK/Rutgers/Kong_Lab/R")
##Run following two lines to install packages if needed --remove leading hashtags to activate  ##
# install.packages("gplots")
library("gplots")
#####read table: check.names indicate that the +_will be recognized
df<-read.table(file="data/heatmap/8 dss.txt",sep = "\t", header = TRUE, check.names = FALSE)
#####heatmap needs numeric matrix. Only use the expression as matrix
mat <- as.matrix(df[c(2:3)])
### Use the gene symbols as rowname
rownames(mat) <- df$gene
mat
# heatmap.2(mat,col=redgreen(255),dendrogram = "row",Colv = FALSE,hclustfun = hclust)
#set output format
# png(filename = "data/heatmap/heatmap1.png",
#     height = 5,
#     width = 5,
#     units = 'in',
#     res = 300)
#margins
# par(mar=c(5, 4, 4, 2) + 0.1)
#without key and labels
lmat <- rbind(c(4,1),c(3,1),c(2,1))
lmat
dev.off()
heatmap.2(mat,col=redgreen(255),dendrogram = "none",Colv = FALSE,Rowv = FALSE, 
          trace = "none",
          density.info = "none", key = FALSE, labRow = NA, labCol = NA,
          lmat= lmat,  lwid = c(0,4), lhei = c(4,4,4))

##with key and labels
heatmap.2(mat,col=redgreen(255),dendrogram = "none",Colv = FALSE,Rowv = FALSE,
          trace = "none",
          density.info = "none", key.title = "expression", keysize = 2, key.xlab =NA)
# # heatmap.2(mat,col=redgreen(255),dendrogram = "both")
# # #heatmap.2(mat, dendrogram = "row", scale = "none", margins = c(1,10), trace = "none", col = redgreen(255), density.info = "none", key = TRUE, keysize = 2, key.title = FALSE, key.xlab = expression(log[2]~(fold~change)), cexCol=1, srtCol=60)
# heatmap.2(mat, dendrogram = "row", scale = "none", margins = c(1,10), trace = "none", col = redgreen(255), density.info = "none", key = TRUE, keysize = 2, key.title = FALSE, key.xlab = expression(log[2]~expression), cexCol=1, srtCol=NULL)
heatmap.2(mat, lmat= lmat,  lwid = c(1,4), lhei = c(4,4,4), cexRow = 1, cexCol = 1)
