library(data.table)
library("gplots")


p12 <- read.table("data/Chao/pathways_12wks.csv", sep = "\t", header = T, row.names = "Pathways")
p20 <- read.table("data/Chao/pathways_20wks_p2.csv", sep = "\t", header = T, row.names = "Pathways")
p20_1 <- read.table("data/Chao/pathways_20wks_p2.csv", sep = "\t", header = T, row.names = "Pathways")


lmat <- rbind(c(4,1),c(3,1),c(2,1))
lmat
dev.off()

heatmap.2(as.matrix(p20),col=redgreen(255),dendrogram = "none",Colv = FALSE,Rowv = FALSE, 
          trace = "none",
          density.info = "none", key = FALSE, labRow = NA, labCol = NA,
          lmat= lmat,  lwid = c(0,4), lhei = c(4,4,4))


heatmap.2(mat_f,col=redgreen(255),dendrogram = "none",Colv = FALSE,Rowv = FALSE, 
          trace = "none",
          density.info = "none", key = FALSE, labRow = NA, labCol = NA,
          lmat= lmat,  lwid = c(0,4), lhei = c(4,4,4))
heatmap.2(mat,col=bluered(255),dendrogram = "none",Colv = FALSE, 
          density.info = "none", key = FALSE, labCol = NA,
          trace = "none",
          lmat= lmat,  lwid = c(0,4), lhei = c(4,4,4))
heatmap.2(mat_f)




# install.packages("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:30]
df <- as.data.frame(colData(dds)[,"condition"])

colnames(p20)

p20 <- p20[order(-p20$wt20wk_ko20wk),]

pheatmap(p20,
         cluster_rows= F,
         show_rownames= T,
         # treeheight_row = 0,
         cluster_cols=FALSE,
         border_color = NA,
         color = redgreen(255))

pheatmap(p20_1,
         cluster_rows= F,
         show_rownames= T,
         # treeheight_row = 0,
         cluster_cols=FALSE,
         border_color = NA,
         color = redgreen(255))
