library("gplots")
#tmp <- read.table(file = "data/shan_rna/RW_all_log10_fpkm_gt11.csv", sep = "\t", header = T, row.names = "X")
#heatmap
dt <- read.table("data/uvb-skin/rpkm_all_uvb-skin_edgeR_2-9-2018.csv", sep = "\t", header = T, row.names = "Geneid")
dt <- dt[-1]

# Filtering
dt1 <- dt[rowSums(dt <= 40) >= 32,]
dt1 <- dt1[rowSums(dt1 >= 5) >= 32,]


dt2 <- log10(dt1)
mat <- as.matrix(dt2)
lmat <- rbind(c(4,1),c(3,1),c(2,1))

lmat <- rbind(c(4,3), c(2,1))
dev.off()

tiff("data/uvb-skin/heatmap_rpkm_5-40.tiff",
     width = 1500, height = 1500,
     units = "px", pointsize = 12, bg = "white")


heatmap.2(mat,col=redgreen(255),dendrogram = "none",Colv = FALSE,Rowv = FALSE, 
          trace = "none",
          density.info = "none", key = FALSE, labRow = NA, labCol = NA,
          lmat= lmat,  lwid = c(1,4), lhei = c(4,4,4))

heatmap.2(mat,col=bluered(255),dendrogram = "none",Colv = FALSE, 
          density.info = "none", key = FALSE, labCol = NA,
          trace = "none",
          lmat= lmat,  lwid = c(0,4), lhei = c(4,4,4))

###
heatmap.2(mat, lmat= lmat,  lwid = c(1,4), lhei = c(1,4),
          trace = "none",
          col = redgreen(255)# dendrogram control
         )
          # Rowv = TRUE,
          # Colv=if(symm)"Rowv" else TRUE,
          # distfun = dist,
          # hclustfun = hclust,
          # dendrogram = c("both","row","column","none"),
          # reorderfun = function(d, w) reorder(d, w),
          # symm = FALSE,
          # 
          # # data scaling
          # scale = c("none","row", "column"),
          # na.rm=TRUE,
          # 
          # # image plot
          # revC = identical(Colv, "Rowv"),
          # add.expr,
          # 
          # # mapping data to colors
          # breaks,
          # symbreaks=any(x < 0, na.rm=TRUE) || scale!="none",
          # 
          # # colors
          # col="heat.colors",
          # 
          # # block sepration
          # colsep,
          # rowsep,
          # sepcolor="white",
          # sepwidth=c(0.05,0.05),
          # 
          # # cell labeling
          # cellnote,
          # notecex=1.0,
          # notecol="cyan",
          # na.color=par("bg"),
          # 
          # # level trace
          # trace=c("column","row","both","none"),
          # tracecol="cyan",
          # hline=median(breaks),
          # vline=median(breaks),
          # linecol=tracecol,
          # 
          # # Row/Column Labeling
          # margins = c(5, 5),
          # ColSideColors,
          # RowSideColors,
          # cexRow = 0.2 + 1/log10(nr),
          # cexCol = 0.2 + 1/log10(nc),
          # labRow = NULL,
          # labCol = NULL,
          # srtRow = NULL,
          # srtCol = NULL,
          # adjRow = c(0,NA),
          # adjCol = c(NA,0),
          # offsetRow = 0.5,
          # offsetCol = 0.5,
          # colRow = NULL,
          # colCol = NULL,
          # 
          # # color key + density info
          # key = TRUE,
          # keysize = 1.5,
          # density.info=c("histogram","density","none"),
          # denscol=tracecol,
          # symkey = any(x < 0, na.rm=TRUE) || symbreaks,
          # densadj = 0.25,
          # key.title = NULL,
          # key.xlab = NULL,
          # key.ylab = NULL,
          # key.xtickfun = NULL,
          # key.ytickfun = NULL,
          # key.par=list(),
          # 
          # # plot labels
          # main = NULL,
          # xlab = NULL,
          # ylab = NULL,
          # 
          # # plot layout
          # lmat = NULL,
          # lhei = NULL,
          # lwid = NULL,
          # 
          # # extras
          # extrafun=NULL
         
