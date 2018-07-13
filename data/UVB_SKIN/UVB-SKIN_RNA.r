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

dt1 <- dt[, c("Geneid", "Length",
              "02w_CON_0", "02w_CON_1", "02w_UVB_0", "02w_UVB_1", "02w_SFN_0", "02w_SFN_1","02w_UAA_0", "02w_UAA_1",
              "15w_CON_0", "15w_CON_1", "15w_UVB_0", "15w_UVB_1", "15w_SFN_0", "15w_SFN_1","15w_UAA_0", "15w_UAA_1",
              "25w_CON_0", "25w_CON_1", "25w_UVB_0", "25w_UVB_1", "25w_SFN_0", "25w_SFN_1","25w_UAA_0", "25w_UAA_1",
              "25t_CON_0", "25t_CON_1", "25t_UVB_0", "25t_UVB_1", "25t_SFN_0", "25t_SFN_1","25t_UAA_0", "25t_UAA_1")]

# Format counts
cts <- as.matrix(dt1[,3:34])
rownames(cts) <- dt1$Geneid
summary(cts)
colnames(cts)

# Format meta-data
coldata <- data.frame(timepoint = rep(c("02wks", "15wks", "25wks", "25tmr"), each = 8),
                      condition = rep(c("control", "uvb", "uvb+sfn", "uvb+ua"), each = 2, 4),
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

# design 1. two factors no interaction
design(dds)
# ~timepoint + condition
colData(dds)
dds <- dds[rowSums(counts(dds)) >= 32, ]
dds <- DESeq(dds) 
resultsNames(dds)
colnames(model.matrix(~timepoint + condition, coldata))
# [1] "Intercept"                    "timepoint_15wks_vs_02wks"    
# [3] "timepoint_25tmr_vs_02wks"     "timepoint_25wks_vs_02wks"    
# [5] "condition_uvb_vs_control"     "condition_uvb.sfn_vs_control"
# [7] "condition_uvb.ua_vs_control"

results(dds, contrast = c("timepoint", "15wks", "02wks" ))
# Xkr4       0.5753933     -0.1746752 1.48750583 -0.11742824 9.065207e-01           NA

results(dds, contrast = c("timepoint", "25wks", "02wks" ))
# Xkr4       0.5753933     0.38579357 1.47084834  0.26229323 7.930954e-01           NA

results(dds, contrast = c("timepoint", "25tmr", "02wks" ))
# Xkr4       0.5753933      2.2434748 1.37545496   1.6310784 1.028738e-01 1.411058e-01

## 25wks vs 15 wks
results(dds, contrast = c("timepoint", "25wks", "15wks" ))
# Xkr4       0.5753933     0.56046876 1.46618226 0.38226404 0.7022655        NA


results(dds, contrast = c("condition", "uvb", "control" ))
# Xkr4       0.5753933  -2.1844103128 1.39001411 -1.57150226  0.1160660         NA

results(dds, contrast = c("condition", "uvb.sfn", "control" ))
# Xkr4       0.5753933     -1.8212242 1.37906756 -1.3206200 0.18662811        NA

results(dds, contrast = c("condition", "uvb.ua", "control" ))
# Xkr4       0.5753933   -1.475237589 1.37447176 -1.07331241 0.2831310        NA


# design 2. interactions

design(dds) <- ~timepoint*condition
colData(dds)
dds <- DESeq(dds) 
resultsNames(dds)
# [1] "Intercept"                       "timepoint_15wks_vs_02wks"       
# [3] "timepoint_25tmr_vs_02wks"        "timepoint_25wks_vs_02wks"       
# [5] "condition_uvb_vs_control"        "condition_uvb.sfn_vs_control"   
# [7] "condition_uvb.ua_vs_control"     "timepoint15wks.conditionuvb"    
# [9] "timepoint25tmr.conditionuvb"     "timepoint25wks.conditionuvb"    
# [11] "timepoint15wks.conditionuvb.sfn" "timepoint25tmr.conditionuvb.sfn"
# [13] "timepoint25wks.conditionuvb.sfn" "timepoint15wks.conditionuvb.ua" 
# [15] "timepoint25tmr.conditionuvb.ua"  "timepoint25wks.conditionuvb.ua" 

results(dds, name = resultsNames(dds)[2])["Xkr4",]
# Xkr4 0.5753933        1.56539  4.926499  0.317749 0.7506754        NA
results(dds, name = resultsNames(dds)[3])["Xkr4",]
# Xkr4 0.5753933       2.499482  4.858168 0.5144906  0.606909        NA
results(dds, name = resultsNames(dds)[4])["Xkr4",]
# Xkr4 0.5753933       2.604486  4.859705  0.535935 0.5920035        NA
results(dds, name = resultsNames(dds)[5])["Xkr4",]
# Xkr4 0.5753933     -0.7368332  5.018886 -0.1468121 0.8832803        NA
results(dds, name = resultsNames(dds)[6])["Xkr4",]
# Xkr4 0.5753933      0.9727249  5.006858 0.1942785 0.8459578        NA
results(dds, name = resultsNames(dds)[7])["Xkr4",]
# Xkr4 0.5753933     -0.6305716  5.018886 -0.1256397 0.9000171        NA
results(dds, name = resultsNames(dds)[8])["Xkr4",]
# Xkr4 0.5753933      -1.780291  7.032917 -0.2531369 0.8001624        NA
results(dds, name = resultsNames(dds)[9])["Xkr4",]
# Xkr4 0.5753933     -0.3568628  6.910348 -0.05164181 0.9588141        NA
results(dds, name = resultsNames(dds)[10])["Xkr4",]
# Xkr4 0.5753933      -2.741569  6.986296 -0.392421 0.6947472        NA
results(dds, name = resultsNames(dds)[11])["Xkr4",]
# Xkr4 0.5753933      -3.510448   7.02434 -0.4997548 0.6172477        NA
results(dds, name = resultsNames(dds)[12])["Xkr4",]
# Xkr4 0.5753933       -2.18847  6.901617 -0.3170952 0.7511714        NA
results(dds, name = resultsNames(dds)[13])["Xkr4",]
# Xkr4 0.5753933      -4.379194  6.977661 -0.627602 0.5302647        NA
results(dds, name = resultsNames(dds)[14])["Xkr4",]
# Xkr4 0.5753933       -1.83072  7.032917 -0.2603073 0.7946267        NA
results(dds, name = resultsNames(dds)[15])["Xkr4",]
# Xkr4 0.5753933       1.015121  6.861528 0.1479438 0.8823871        NA
results(dds, name = resultsNames(dds)[16])["Xkr4",]
# Xkr4 0.5753933      -2.476443  6.986296 -0.3544715 0.7229856        NA
# 

# uvb vs control in 15 weeks group
results(dds, list(c("condition_uvb_vs_control", "timepoint15wks.conditionuvb" )))["Xkr4",]
# Xkr4 0.5753933      -2.517124  4.927161 -0.510867 0.6094442        NA
write.table(results(dds, list(c("condition_uvb_vs_control", "timepoint15wks.conditionuvb" ))),
            "data/UVB_SKIN/temp/timepoint.x.condition_15wksuvb-15wkscontrol.csv", sep = "\t", quote = F,
            col.names = NA)


# Design 3.
design(dds) <- ~ group
colData(dds)
dds <- DESeq(dds) 
resultsNames(dds)
# [1] "Intercept"                          "group_02wksuvb_vs_02wkscontrol"    
# [3] "group_02wksuvb.sfn_vs_02wkscontrol" "group_02wksuvb.ua_vs_02wkscontrol" 
# [5] "group_15wkscontrol_vs_02wkscontrol" "group_15wksuvb_vs_02wkscontrol"    
# [7] "group_15wksuvb.sfn_vs_02wkscontrol" "group_15wksuvb.ua_vs_02wkscontrol" 
# [9] "group_25tmrcontrol_vs_02wkscontrol" "group_25tmruvb_vs_02wkscontrol"    
# [11] "group_25tmruvb.sfn_vs_02wkscontrol" "group_25tmruvb.ua_vs_02wkscontrol" 
# [13] "group_25wkscontrol_vs_02wkscontrol" "group_25wksuvb_vs_02wkscontrol"    
# [15] "group_25wksuvb.sfn_vs_02wkscontrol" "group_25wksuvb.ua_vs_02wkscontrol" 

results(dds, contrast = c("group", "15wksuvb", "15wkscontrol"))["Xkr4",]
# Xkr4 0.5753933      -2.517217  4.927297 -0.5108717 0.6094409        NA
# values are almost the same as in the previous design.

write.table(results(dds, contrast = c("group", "15wksuvb", "15wkscontrol")),
            "data/UVB_SKIN/temp/group_15wksuvb-15wkscontrol.csv", sep = "\t", quote = F,
            col.names = NA)
# 

names <- sub("_vs.*","", resultsNames(dds))
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

library("vsn")
meanSdPlot((assay(ntd)))
meanSdPlot((assay(vsd)))
meanSdPlot((assay(rld)))

# install.packages("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("timepoint", "condition")])
pheatmap(assay(vsd)[select,], cluster_rows=F, show_rownames=F,
         cluster_cols=F, annotation_col = df,
         #color = heat.colors(1024), #
         border_color = NA,
         scale = "row",
         legend = T)
dev.off()

# Data for Control and UVB groups.
dds1 <- dds
dds1 <- dds1[, paste0(rep(c("02w", "15w", "25w", "25t"), each = 4),
                      "_",
                      rep(c("CON", "UVB"), each = 2, 4),
                      "_",
                      rep(c("0", "1"), 8))]

# Check groups.
as.data.frame(colData(dds1))

# transform (log2(count + 1))
ntd1 <- normTransform(dds1)
# heatmap
select1 <- order(rowMeans(counts(dds1,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df1 <- as.data.frame(colData(dds1)[,c("condition", "timepoint")]) # The last column will be 
# used first to draw the column annotation by "annotation_col = ".

pheatmap(assay(ntd1)[select1,], cluster_rows=F, show_rownames=F,
         cluster_cols=F, 
         annotation_col = df1,
         annotation_colors = list(condition = c(control = "green", uvb = "#FF0000")), #, "#66A61E", "#FF0000")),
         color = redblue(255), #color = heat.colors(1024), # colorpanel(n, low, mid, high) 
        # color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
         border_color = NA,
         scale = "row",
         legend = T)



# 5-9
dds1 <- dds1[rowSums(counts(dds) >= 1) >= 10, ] # keep rows that have a value no smaller than 1 on at least 10 columns.
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


