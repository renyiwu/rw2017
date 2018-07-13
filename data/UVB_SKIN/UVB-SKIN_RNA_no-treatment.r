# DESeq2 analyses on UVB-SKIN project.
# R WU 5-1-2018

# Load libraries.
library("DESeq2")
library("data.table")
library("gplots")

dt <- fread("data/UVB_SKIN/featurescounts_uvb-skin_dedup_renyi_2-9-2018.csv")
dt <-dt[, -(2:5)]
colnames(dt)
# [1] "Geneid"              "Length"              "02w_CON_0.dedup.bam" "02w_CON_1.dedup.bam"
# [5] "02w_SFN_0.dedup.bam" "02w_SFN_1.dedup.bam" "02w_UAA_0.dedup.bam" "02w_UAA_1.dedup.bam"
# [9] "02w_UVB_0.dedup.bam" "02w_UVB_1.dedup.bam" "15w_CON_0.dedup.bam" "15w_CON_1.dedup.bam"
# [13] "15w_SFN_0.dedup.bam" "15w_SFN_1.dedup.bam" "15w_UAA_0.dedup.bam" "15w_UAA_1.dedup.bam"
# [17] "15w_UVB_0.dedup.bam" "15w_UVB_1.dedup.bam" "25t_CON_0.dedup.bam" "25t_CON_1.dedup.bam"
# [21] "25t_SFN_0.dedup.bam" "25t_SFN_1.dedup.bam" "25t_UAA_0.dedup.bam" "25t_UAA_1.dedup.bam"
# [25] "25t_UVB_0.dedup.bam" "25t_UVB_1.dedup.bam" "25w_CON_0.dedup.bam" "25w_CON_1.dedup.bam"
# [29] "25w_SFN_0.dedup.bam" "25w_SFN_1.dedup.bam" "25w_UAA_0.dedup.bam" "25w_UAA_1.dedup.bam"
# [33] "25w_UVB_0.dedup.bam" "25w_UVB_1.dedup.bam"

colnames(dt)[3:34] <- substring(colnames(dt), 1, 9)[3:34]
colnames(dt)
# [1] "Geneid"    "Length"    "02w_CON_0" "02w_CON_1" "02w_SFN_0" "02w_SFN_1" "02w_UAA_0"
# [8] "02w_UAA_1" "02w_UVB_0" "02w_UVB_1" "15w_CON_0" "15w_CON_1" "15w_SFN_0" "15w_SFN_1"
# [15] "15w_UAA_0" "15w_UAA_1" "15w_UVB_0" "15w_UVB_1" "25t_CON_0" "25t_CON_1" "25t_SFN_0"
# [22] "25t_SFN_1" "25t_UAA_0" "25t_UAA_1" "25t_UVB_0" "25t_UVB_1" "25w_CON_0" "25w_CON_1"
# [29] "25w_SFN_0" "25w_SFN_1" "25w_UAA_0" "25w_UAA_1" "25w_UVB_0" "25w_UVB_1"

# dt1 <- dt[, c("Geneid", "Length",
#               "02w_CON_0", "02w_CON_1", "02w_UVB_0", "02w_UVB_1", "02w_SFN_0", "02w_SFN_1","02w_UAA_0", "02w_UAA_1",
#               "15w_CON_0", "15w_CON_1", "15w_UVB_0", "15w_UVB_1", "15w_SFN_0", "15w_SFN_1","15w_UAA_0", "15w_UAA_1",
#               "25w_CON_0", "25w_CON_1", "25w_UVB_0", "25w_UVB_1", "25w_SFN_0", "25w_SFN_1","25w_UAA_0", "25w_UAA_1",
#               "25t_CON_0", "25t_CON_1", "25t_UVB_0", "25t_UVB_1", "25t_SFN_0", "25t_SFN_1","25t_UAA_0", "25t_UAA_1")]

dt1 <- dt[, c("Geneid", "Length",
              "02w_CON_0", "02w_CON_1", "02w_UVB_0", "02w_UVB_1",
              "15w_CON_0", "15w_CON_1", "15w_UVB_0", "15w_UVB_1",
              "25w_CON_0", "25w_CON_1", "25w_UVB_0", "25w_UVB_1",
              "25t_CON_0", "25t_CON_1", "25t_UVB_0", "25t_UVB_1")]

colnames(dt1)

# Format counts
cts <- as.matrix(dt1[,3:18])
rownames(cts) <- dt1$Geneid
summary(cts)
colnames(cts)

# Format meta-data
coldata <- data.frame(timepoint = rep(c("02wks", "15wks", "25wks", "25tmr"), each = 4),
                      condition = rep(c("control", "uvb"), each = 2, 4),
                      row.names = colnames(cts))
coldata

# Construct DEseq data table.
dds <- DESeqDataSetFromMatrix(cts, coldata, ~ timepoint + condition)
dds$group <- factor(paste0(dds$timepoint, dds$condition))
mcols(dds) <- DataFrame(mcols(dds), basepairs = dt1$Length) # Add more info for RPKM calculation.
# mcols(dds) <- NULL


# FPKM
fpkm <- as.data.frame(fpkm(dds))
fwrite(fpkm, file = "data/UVB_SKIN/May/fpkm.all.csv", sep = "\t",
       row.names = T)

# Pre-filtering. Method one:
keep <- rowSums(counts(dds)) >= 10 # Keep rows that have at least 10 reads total.
dds <- dds[keep,]
# method two
dds <- dds[rowSums(counts(dds) >= 1 ) >= 16, ]


# design 1. two factors no interaction
design(dds)
# ~timepoint + condition

colData(dds)
# > as.data.frame(colData(dds))
# timepoint condition        group sizeFactor
# 02w_CON_0     02wks   control 02wkscontrol  0.9140164
# 02w_CON_1     02wks   control 02wkscontrol  0.3804442
# 02w_UVB_0     02wks       uvb     02wksuvb  0.9240064
# 02w_UVB_1     02wks       uvb     02wksuvb  1.0257816
# 15w_CON_0     15wks   control 15wkscontrol  0.8656174
# 15w_CON_1     15wks   control 15wkscontrol  1.0988142
# 15w_UVB_0     15wks       uvb     15wksuvb  1.1369088
# 15w_UVB_1     15wks       uvb     15wksuvb  1.1418339
# 25w_CON_0     25wks   control 25wkscontrol  1.0181438
# 25w_CON_1     25wks   control 25wkscontrol  1.1259548
# 25w_UVB_0     25wks       uvb     25wksuvb  1.1637136
# 25w_UVB_1     25wks       uvb     25wksuvb  0.9961790
# 25t_CON_0     25tmr   control 25tmrcontrol  1.6479497
# 25t_CON_1     25tmr   control 25tmrcontrol  0.9445733
# 25t_UVB_0     25tmr       uvb     25tmruvb  1.1827763
# 25t_UVB_1     25tmr       uvb     25tmruvb  1.2403322


dds <- DESeq(dds) 

resultsNames(dds)
# [1] "Intercept"                "timepoint_15wks_vs_02wks" "timepoint_25tmr_vs_02wks" "timepoint_25wks_vs_02wks" "condition_uvb_vs_control"

colnames(model.matrix(~timepoint + condition, coldata))
# [1] "(Intercept)"    "timepoint15wks" "timepoint25tmr" "timepoint25wks" "conditionuvb"  


# design 2. interactions

design(dds) <- ~timepoint*condition
colData(dds)

dds <- DESeq(dds) 

resultsNames(dds)
# [1] "Intercept"                   "timepoint_15wks_vs_02wks"    "timepoint_25tmr_vs_02wks"    "timepoint_25wks_vs_02wks"   
# [5] "condition_uvb_vs_control"    "timepoint15wks.conditionuvb" "timepoint25tmr.conditionuvb" "timepoint25wks.conditionuvb" 


# uvb vs control in 15 weeks group
results(dds, list(c("condition_uvb_vs_control", "timepoint15wks.conditionuvb" )))["Xkr4",]
# Xkr4 0.6745287      -2.517546  4.288322 -0.5870701 0.5571566        NA
 


# Design 3.
design(dds) <- ~ group
colData(dds)
# > as.data.frame(colData(dds))
# timepoint condition        group sizeFactor
# 02w_CON_0     02wks   control 02wkscontrol  0.9140164
# 02w_CON_1     02wks   control 02wkscontrol  0.3804442
# 02w_UVB_0     02wks       uvb     02wksuvb  0.9240064
# 02w_UVB_1     02wks       uvb     02wksuvb  1.0257816
# 15w_CON_0     15wks   control 15wkscontrol  0.8656174
# 15w_CON_1     15wks   control 15wkscontrol  1.0988142
# 15w_UVB_0     15wks       uvb     15wksuvb  1.1369088
# 15w_UVB_1     15wks       uvb     15wksuvb  1.1418339
# 25w_CON_0     25wks   control 25wkscontrol  1.0181438
# 25w_CON_1     25wks   control 25wkscontrol  1.1259548
# 25w_UVB_0     25wks       uvb     25wksuvb  1.1637136
# 25w_UVB_1     25wks       uvb     25wksuvb  0.9961790
# 25t_CON_0     25tmr   control 25tmrcontrol  1.6479497
# 25t_CON_1     25tmr   control 25tmrcontrol  0.9445733
# 25t_UVB_0     25tmr       uvb     25tmruvb  1.1827763
# 25t_UVB_1     25tmr       uvb     25tmruvb  1.2403322


dds <- DESeq(dds) 

resultsNames(dds)
# [1] "Intercept"                          "group_02wksuvb_vs_02wkscontrol"     "group_15wkscontrol_vs_02wkscontrol"
# [4] "group_15wksuvb_vs_02wkscontrol"     "group_25tmrcontrol_vs_02wkscontrol" "group_25tmruvb_vs_02wkscontrol"    
# [7] "group_25wkscontrol_vs_02wkscontrol" "group_25wksuvb_vs_02wkscontrol" 

results(dds, contrast = c("group", "15wksuvb", "15wkscontrol"))["Xkr4",]
# Xkr4 0.6745287      -2.517569  4.288378 -0.587068  0.557158        NA
# Almost the same as in the previous design.


# 

names <- sub("_vs.*","", resultsNames(dds)) #(pattern, replacement, x)
names <- sub("group_", "", names)
names[1] <- "02wkscontrol"
names
n <- length(resultsNames(dds))
m <- 1
for (m in 1:(n-1)) {
        for (i in (m+1):n) {
                file <- paste0("data/UVB_SKIN/temp1/group_", names[i],"-",names[m], ".csv")
                write.table(results(dds, contrast = c("group", names[i], names[m])),
                            file, sep = "\t",
                            quote = F,
                            col.names = NA)
                }
        m <- m+1
}


m <- 1
for (m in 1:(n-1)) {
        for (i in (m+1):n) {
                file <- paste0("data/UVB_SKIN/temp2/group_", names[m],"-",names[i], ".csv")
                write.table(results(dds, contrast = c("group", names[m], names[i])),
                            file, sep = "\t",
                            quote = F,
                            col.names = NA)
        }
        m <- m+1
}

## stop here


# 5-8-2018

vsd <- vst(dds, blind = F)
rld <- rlog(dds, blind = F)
head(assay(vsd), 3)
head(assay(rld), 3)
head(assay(dds), 3)

# Check effects of transformation on the variance
ntd <- normTransform(dds)
head(assay(ntd), 3)

plotPCA(vsd, intgroup = c("condition", "timepoint"))


library("vsn")
meanSdPlot((assay(ntd)))
meanSdPlot((assay(vsd)))
meanSdPlot((assay(rld)))
meanSdPlot((assay(dds)))

# install.packages("pheatmap")
library("pheatmap")

df <- as.data.frame(colData(dds)[,c("timepoint", "condition")])

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:1000]

pheatmap(assay(rld)[select,],
         cluster_rows= T,
         treeheight_row = 0,
         show_rownames= F,
         cluster_cols= F,
         treeheight_col = 0,
         annotation_col = df,
         #color = palette(rainbow(1024)), #rgb((255:1)/255, (1:255)/255, 0), #heat.colors(1024), #gray(0:90 / 100), #
         border_color = NA,
         scale = "row", #"row",
         legend = T)
dev.off()

# Data for epithelial cells only.
dds1 <- dds
dds1 <- dds1[, paste0(rep(c("02w", "15w", "25w"), each = 4),
                      "_",
                      rep(c("CON", "UVB"), each = 2, 3),
                      "_",
                      rep(c("0", "1"), 6))]

# Check groups.
as.data.frame(colData(dds1))

# transform (log2(count + 1))
ntd1 <- normTransform(dds1)
# heatmap
select1 <- order(rowMeans(counts(dds1,normalized=TRUE)),
                decreasing=TRUE)[1:10000]
df1 <- as.data.frame(colData(dds1)[,c("condition", "timepoint")]) # The last column will be 
# used first to draw the column annotation by "annotation_col = ".

pheatmap(assay(ntd1)[select1,],
         cluster_rows= T,
         treeheight_row = 0,
         show_rownames=F,
         cluster_cols=F, 
         annotation_col = df1,
         annotation_colors = list(condition = c(control = "green", uvb = "#FF0000"), timepoint = c( '02wks' = "#00FFFF", '15wks' = "#FFFF00", '25wks' = "#FF00FF")), #, "#66A61E", "#FF0000")),
         #color = redblue(255), #color = heat.colors(1024), # colorpanel(n, low, mid, high) 
        # color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
         border_color = NA,
         scale = "row",
         legend = T)

####
# 4 samples in time point "25t", which are whole skin from control group and tumor in UVB group,
# had a very different expression pattern from the other 12 samples.
# This difference was due to tissue (whole skin vs epithelial cells)?

# 5-9
# dds1 <- dds1[rowSums(counts(dds) >= 1) >= 10, ] # keep rows that have a value no smaller than 1 on at least 10 columns.
design(dds1)
as.data.frame(colData(dds1))
resultsNames(dds1)
# [1] "Intercept"                          "group_02wksuvb_vs_02wkscontrol"    
# [3] "group_15wkscontrol_vs_02wkscontrol" "group_15wksuvb_vs_02wkscontrol"    
# [5] "group_25tmrcontrol_vs_02wkscontrol" "group_25tmruvb_vs_02wkscontrol"    
# [7] "group_25wkscontrol_vs_02wkscontrol" "group_25wksuvb_vs_02wkscontrol"

dds1$group <- droplevels(dds1$group) # Required
dds1 <- DESeq(dds1)




results(dds1, contrast = c("group", "02wksuvb", "02wkscontrol"))
# Xkr4             0.6745287    -0.72338412 4.2703510 -0.1693969  0.8654845         NA

# 1
res <- as.data.frame(results(dds1, contrast = c("group", "02wksuvb", "02wkscontrol")))
fwrite(res,
       "data/UVB_SKIN/May/weeks02_UVB-control.csv",
       sep = "\t",
       row.names = TRUE)


# 2

fwrite(as.data.frame(results(dds1, contrast = c("group", "15wksuvb", "15wkscontrol"))),
       "data/UVB_SKIN/May/weeks15_UVB-control.csv",
       sep = "\t",
       row.names = TRUE)

# 3
fwrite(as.data.frame(results(dds1, contrast = c("group", "25wksuvb", "25wkscontrol"))),
       "data/UVB_SKIN/May/weeks25_UVB-control.csv",
       sep = "\t",
       row.names = TRUE)


# 4
fwrite(as.data.frame(results(dds1, contrast = c("group", "25tmruvb", "25tmrcontrol"))),
       "data/UVB_SKIN/May/weeks25tmr_UVB-control.csv",
       sep = "\t",
       row.names = TRUE)

# 5
fwrite(as.data.frame(results(dds1, contrast = c("group", "15wkscontrol", "02wkscontrol"))),
       "data/UVB_SKIN/May/weeks15-02_control.csv",
       sep = "\t",
       row.names = TRUE)

# 6
fwrite(as.data.frame(results(dds1, contrast = c("group", "25wkscontrol", "02wkscontrol"))),
       "data/UVB_SKIN/May/weeks25-02_control.csv",
       sep = "\t",
       row.names = TRUE)

# 7
fwrite(as.data.frame(results(dds1, contrast = c("group", "15wksuvb", "02wksuvb"))),
       "data/UVB_SKIN/May/weeks15-02_uvb.csv",
       sep = "\t",
       row.names = TRUE)

# 8
fwrite(as.data.frame(results(dds1, contrast = c("group", "25wksuvb", "02wksuvb"))),
       "data/UVB_SKIN/May/weeks25-02_uvb.csv",
       sep = "\t",
       row.names = TRUE)

#9 whole skin vs epidermal 
res9 <- results(dds1, contrast = c("group", "25tmrcontrol", "02wkscontrol"))

sum(res9$log2FoldChange > 1 | res9$log2FoldChange < -1, na.rm = TRUE)
# 9066

sum(res9$padj  < 0.05, na.rm = TRUE)
# 8535 

sum(res9$pvalue  < 0.05, na.rm = TRUE)
# 9260 

# Check how many genes are significantly changed.
sum(res$padj < 0.1, na.rm = TRUE)
# 3429

summary(results(dds1, alpha = 0.1))

# MA plot
res1 <- results(dds1, contrast = c("group", "02wksuvb", "02wkscontrol"))
plotMA(res1)
res1LFC <- lfcShrink(dds1, coef = 2)
plotMA(res1LFC)

# 
vsd1 <- vst(dds1, blind = F)

plotPCA(vsd1, intgroup = c("timepoint", "condition"))
 









colnames(model.matrix(~timepoint*condition, coldata))
# the same as above. 

# To use another design
# design(dds) <- formula(~ condition)
# and re-run des
# dds <- DEseq(dds)

# To include interactions

dds$group <- factor(paste0(dds$timepoint, dds$condition))
colData(dds)

design(dds) <- ~ group
dds
design(dds)
resultsNames(dds)


# design <- model.matrix(~ timepoint %in% condition, coldata)
# design <- model.matrix(~ condition + condition %in% timepoint, coldata)
# design <- model.matrix(~ timepoint + condition %in% timepoint, coldata)
# design <- model.matrix(~ timepoint / condition %in% timepoint, coldata)
# colnames(design)


# dds$group <- factor(paste0(dds$genotype, dds$condition))
# design(dds) <- ~ group
# dds <- DESeq(dds)
# resultsNames(dds)
# results(dds, contrast=c("group", "IB", "IA"))


