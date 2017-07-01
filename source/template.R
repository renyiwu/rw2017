# R scrpits originally by John. 
setwd("C:/Users/renyi/Google Drive/WORK/Rutgers/Kong_Lab/R")
##Run following two lines to install packages if needed --remove leading hashtags to activate  ##
# install.packages("gplots")
library("gplots")
#####read table: check.names indicate that the +_will be recognized
df<-read.table(file="heatmap/Genes in common 8 wk DSS.txt",sep = "\t", header = TRUE, check.names = FALSE)
df<-read.table(file="heatmap/Genes in common 8 wk AOMDSS.txt",sep = "\t", header = TRUE, check.names = FALSE)
df<-read.table(file="heatmap/Genes in common 18wk AOMDSS.txt",sep = "\t", header = TRUE, check.names = FALSE)
df<-read.table(file="heatmap/pathways 8 weeks AOMDSS.txt",sep = "\t", header = TRUE, check.names = FALSE)
df<-read.table(file="heatmap/pathways 18 weeks AOMDSS.txt",sep = "\t", header = TRUE, check.names = FALSE)
#####heatmap needs numeric matrix. Only use the expression as matrix
mat <- as.matrix(df[c(2:3)])
### Use the gene symbols as rowname
rownames(mat) <- df$pathway
heatmap.2(mat,col=redgreen(255),dendrogram = "row",Colv = FALSE,hclustfun = hclust)

#without key and labels
heatmap.2(mat,col=redgreen(255),dendrogram = "none",Colv = FALSE,Rowv = FALSE, 
          trace = "none",
          density.info = "none", key = FALSE, labRow = NA, labCol = NA)
##with key and labels
heatmap.2(mat,col=redgreen(255),dendrogram = "none",Colv = FALSE,Rowv = FALSE, 
          trace = "none",
          density.info = "none", key.title = "expression", keysize = 2, key.xlab =NA)
heatmap.2(mat,col=redgreen(255),dendrogram = "both")
#heatmap.2(mat, dendrogram = "row", scale = "none", margins = c(1,10), trace = "none", col = redgreen(255), density.info = "none", key = TRUE, keysize = 2, key.title = FALSE, key.xlab = expression(log[2]~(fold~change)), cexCol=1, srtCol=60)
heatmap.2(mat, dendrogram = "row", scale = "none", margins = c(1,10), trace = "none", col = redgreen(255), density.info = "none", key = TRUE, keysize = 2, key.title = FALSE, key.xlab = expression(log[2]~expression), cexCol=1, srtCol=NULL)
heatmap.2(mat)

