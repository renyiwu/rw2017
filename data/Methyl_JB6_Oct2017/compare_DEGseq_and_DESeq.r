library("DESeq2", "data.table")

dt01 <- fread("data/RNA_JB6_Oct2017/shan_rna/pairwise/S4_S3_TPA-control.csv")
dt01 <- dt01[,c(1,3,6)]
colnames(dt01) <- c("Geneid", "des.logFC", "des.p")

dt02 <- fread("data/RNA_JB6_Oct2017/DEGseq/RW4-RW3_TPA0-Control.csv")
dt02 <- dt02[, c(1,5,7)]
colnames(dt02) <-c("Geneid", "deg.logFC", "deg.p")

# Merge data tables
dt.m <- merge(dt01, dt02, by = "Geneid", # by = intersect(names(x), names(y)) or set by.x and by.y seperately.
              all = FALSE, # or set all.x and all.y seperately
              suffixs = "" # or c(".x", ".y")
              )

sum(dt.m$des.p < 0.1, na.rm = TRUE)
# 1102
nrow(dt.m[dt.m$des.p < 0.1,])
# 1102
nrow(dt.m[dt.m$des.p < 0.05,])
# 330



nrow(dt.m[dt.m$deg.p < 0.1,])
# 8200
nrow(dt.m[dt.m$deg.p < 0.05,])
# 7638


# Melt
melt(dt.m, id = 1, measure.vars = c(2,4))
