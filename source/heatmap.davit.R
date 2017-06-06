# Source: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

require(data.table)
require(ggplot2)

# Read data
#DATA_HOME <- "/home/administrator/Documents/gene.seq.kong/data"
DATA_HOME <- "C:/git_local/gene.seq.kong/data"
# Positive vs. negative controls
dt1 <- fread(paste(DATA_HOME,
                   "AOMDSS-control 18 wk.csv",
                   sep = "/"))
dt1 <- dt1[, 1:3]
names(dt1) <- c("gene",
                "dlog2.ctrl.aomdss",
                "padj.ctrl.aomdss")
dt1

# Positive control vs. treatment
dt2 <- fread(paste(DATA_HOME,
                   "AOMDSSCur- AOMDSS 18 wk.csv",
                   sep = "/"))
dt2 <- dt2[, 1:3]
names(dt2) <- c("gene",
                "dlog2.aomdss.cur",
                "padj.aomds.cur")
dt2

# Merge the data
dt3 <- merge(dt1,
             dt2,
             by = "gene")
dt3

# Melt the data
dt4 <- melt.data.table(data = dt3,
                       id.vars = "gene",
                       measure.vars = c("dlog2.ctrl.aomdss",
                                        "dlog2.aomdss.cur"))

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
