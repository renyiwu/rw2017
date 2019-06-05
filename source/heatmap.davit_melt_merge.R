# Source: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

require(data.table)
require(ggplot2)

########


meltdata.wu <- function(file1, file2, file3){
  require(data.table)
  dt1 <- read.table(file1, header = T)
  dt1 <- cbind(Geneid = row.names(dt1), dt1[c(2,5, 6)])
  row.names(dt1) <- c()
  colnames(dt1) <- c("Geneid", "LogFC", "Pvalue", "FDR")
  dt2 <- read.table(file2, header = T)
  dt2 <- dt2[c(1,3,6,7)]
  colnames(dt2) <- c("Geneid", "LogFC", "Pvalue", "FDR")
  
  # Merge the data
  dt3 <- merge(dt1, #Deseq2
               dt2, # edgeR
               by = "Geneid",
               all = T,
               suffixes = c(".Deseq2",".edgeR"))
  
  
  # Write table
  #dir.create("data/RNA_cur_8wks/edgeR-DEseq2/")
  write.table(dt3, file3,
              col.names = T, row.names = F,
              sep = "\t", quote = F)
  ##
  return(dt3)
}
###########
f1 <- "data/RNA_cur_8wks/DEseq2/AOM+DSS_DSS.csv"
f2 <- "data/RNA_cur_8wks/edgeR/AOM_DSS-DSS.csv"
f3 <- "data/RNA_cur_8wks/edgeR-DEseq2/AOM+DSS_DSS.csv"
meltdata.wu(f1, f2, f3)




################
dt1 <- read.table("data/RNA_cur_8wks/DEseq2/AOM+DSS_Control.csv", header = T)
dt1 <- cbind(Geneid = row.names(dt1), dt1[c(2,5, 6)])
row.names(dt1) <- c()
colnames(dt1) <- c("Geneid", "LogFC", "Pvalue", "FDR")
dt2 <- read.table("data/RNA_cur_8wks/edgeR/AOM_DSS-Control.csv", header = T)
dt2 <- dt2[c(1,3,6,7)]
colnames(dt2) <- c("Geneid", "LogFC", "Pvalue", "FDR")

# Merge the data
dt3 <- merge(dt1, #Deseq2
             dt2, # edgeR
             by = "Geneid",
             all = T,
             suffixes = c(".Deseq2",".edgeR"))
 

# Write table
dir.create("data/RNA_cur_8wks/edgeR-DEseq2/")
write.table(dt3, "data/RNA_cur_8wks/edgeR-DEseq2/AOM_DSS-Control.csv",
            col.names = T, row.names = F,
            sep = "\t", quote = F)
write.csv(dt3, "data/RNA_cur_8wks/edgeR-DEseq2/AOM_DSS-Control_2.csv", row.names = F)





dt3 <- as.data.table(dt3)

# Melt the data
dt4 <- melt.data.table(data = dt3,
                       id.vars = "Geneid",
                       measure.vars = c("LogFC.Deseq2",
                                        "LogFC.edgeR"))

dt4$gene <- factor(dt4$gene)
levels(dt4$variable) <- c("AOMDSS - Negative Control",
                          "AOMDSSCur - AOMDSS")

# Heatmap
# tiff(filename = "tmp/heatmap1.tiff",
#      height = 5,
#      width = 6,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")

png(filename = "tmp/heatmap1.png",
     height = 5,
     width = 5,
     units = 'in',
     res = 300)
ggplot(data = dt4) +
  geom_tile(aes(x = variable,
                y = gene,
                fill = value)) +
  scale_fill_gradient2(low = "red", 
                       high = "green", 
                       mid = "black", 
                       midpoint = 0, 
                       limit = c(-3, 3), 
                       space = "Lab", 
                       name="Pearson\nCorrelation") +
  scale_x_discrete("Treatment") + 
  scale_y_discrete("Gene") +
  ggtitle("Hitmap")  +
  theme(axis.text.x = element_text(angle = 5))
graphics.off()

# NOTE (05/25/2017):
# i. Repeat this for:
#    a. Con vs AOMDSS 8wk study2.xlsx vs. AOMDSS vs AOMDSS Cur 8 wkstudy 2.xlsxb
#    b. Con vs DSS study2.xlsx vs. DSS vs DSS Cur study2.xlsx
# 
# ii. Do Venn Diagrams on these 3 sets (i.e. 6 files)
