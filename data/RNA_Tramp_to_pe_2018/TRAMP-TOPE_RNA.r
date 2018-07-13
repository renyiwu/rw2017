# DESeq2 analyses on TRAMP PEITC, tocotrienol project. TO(cotrienol) -- PE(ITC) -> TOPE
# R WU 7-9-2018

# Load libraries.
library("DESeq2")
library("data.table")
library("gplots")

dt <- fread("data/RNA_Tramp_to_pe_2018/featurecounts.results.TR-all.csv")
dt <-dt[, -(2:5)]
colnames(dt)
# [1] "Geneid"         "Length"         "TR01.dedup.bam" "TR02.dedup.bam" "TR04.dedup.bam"
# [6] "TR05.dedup.bam" "TR07.dedup.bam" "TR09.dedup.bam" "TR10.dedup.bam" "TR11.dedup.bam"
# [11] "TR13.dedup.bam" "TR14.dedup.bam" "TR16.dedup.bam" "TR17.dedup.bam" "TR19.dedup.bam"
# [16] "TR20.dedup.bam" "TR21.dedup.bam" "TR23.dedup.bam"

colnames(dt)[3:18] <- substring(colnames(dt), 1, 4)[3:18]
colnames(dt)
# [1] "Geneid" "Length" "TR01"   "TR02"   "TR04"   "TR05"   "TR07"   "TR09"   "TR10"   "TR11"  
# [11] "TR13"   "TR14"   "TR16"   "TR17"   "TR19"   "TR20"   "TR21"   "TR23" 



# Format counts
cts <- as.matrix(dt[,3:18])
rownames(cts) <- dt$Geneid
summary(cts)
colnames(cts)

# Format meta-data
coldata <- data.frame(group = c(rep(c("TRAMP_24wks", "PEITC_24wks", "TOCO_24wks", "TRAMP_16wks", "PEITC_16wks", "TOCO_16wks", "WILD_24wks"), each = 2), "WILD_16wks", "TRAMP.T_24wks"),
                      timepoint = c(rep(c("24", "16"), each = 6), "24", "24", "16", "24"),
                      genetype = c(rep("tramp", 12), rep("wildtype", 3), "tramp"),
                      condition = c(rep(c("control", "peitc", "toco"), each = 2, 2), rep("control", 4)),
                      sampletype = c(rep("Ventral", 15), "tumor"),
                      row.names = colnames(cts))
coldata

# 		group timepoint genetype condition sampletype
# TR01   TRAMP_24wks        24    tramp   control    Ventral
# TR02   TRAMP_24wks        24    tramp   control    Ventral
# TR04   PEITC_24wks        24    tramp     peitc    Ventral
# TR05   PEITC_24wks        24    tramp     peitc    Ventral
# TR07    TOCO_24wks        24    tramp      toco    Ventral
# TR09    TOCO_24wks        24    tramp      toco    Ventral
# TR10   TRAMP_16wks        16    tramp   control    Ventral
# TR11   TRAMP_16wks        16    tramp   control    Ventral
# TR13   PEITC_16wks        16    tramp     peitc    Ventral
# TR14   PEITC_16wks        16    tramp     peitc    Ventral
# TR16    TOCO_16wks        16    tramp      toco    Ventral
# TR17    TOCO_16wks        16    tramp      toco    Ventral
# TR19    WILD_24wks        24 wildtype   control    Ventral
# TR20    WILD_24wks        24 wildtype   control    Ventral
# TR21    WILD_16wks        16 wildtype   control    Ventral
# TR23 TRAMP.T_24wks        24    tramp   control      tumor

# Construct DEseq data table.
dds <- DESeqDataSetFromMatrix(cts, coldata, ~ group)
#dds$group <- factor(paste0(dds$timepoint, dds$condition))
mcols(dds) <- DataFrame(mcols(dds), basepairs = dt$Length) # Add more info for RPKM calculation.
# mcols(dds) <- NULL


# FPKM
fpkm <- as.data.frame(fpkm(dds))
fwrite(fpkm, file = "data/RNA_Tramp_to_pe_2018/TRAMP_fpkm.all_7-9-2018.csv", sep = "\t",
       row.names = T)

# Pre-filtering. 
# Method one:
# keep <- rowSums(counts(dds)) >= 10 # Keep rows that have at least 10 reads total.
# dds <- dds[keep,]

# Prefiltering. Method two:
dds <- dds[rowSums(counts(dds) >= 5) >= 16, ]

dds
# keep genes with at least one count in at least two samples.



# Design 3.
design(dds)
# ~group
colData(dds)

dds <- DESeq(dds) 

resultsNames(dds)
# [1] "Intercept"                          "group_PEITC_24wks_vs_PEITC_16wks"   "group_TOCO_16wks_vs_PEITC_16wks"   
# [4] "group_TOCO_24wks_vs_PEITC_16wks"    "group_TRAMP_16wks_vs_PEITC_16wks"   "group_TRAMP_24wks_vs_PEITC_16wks"  
# [7] "group_TRAMP.T_24wks_vs_PEITC_16wks" "group_WILD_16wks_vs_PEITC_16wks"    "group_WILD_24wks_vs_PEITC_16wks"   

results(dds, contrast = c("group", "PEITC_24wks", "PEITC_16wks"))

# write.table(results(dds, contrast = c("group", "15wksuvb", "15wkscontrol")),
            # "data/UVB_SKIN/temp/group_15wksuvb-15wkscontrol.csv", sep = "\t", quote = F,
            # col.names = NA)
#

# names <- sub("_vs.*","", resultsNames(dds))
# names <- sub("group_", "", names)
# names[1] <- "02wkscontrol"
# names
# (names)

names <- as.character(unique(colData(dds)[, "group"]))

n <- length(resultsNames(dds))
m <- 1
for (m in 1:(n-1)) {
        for (i in (m+1):n) {
                file <- paste0("data/RNA_Tramp_to_pe_2018/temp1/group_", names[i],"-", names[m], ".csv")
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
                file <- paste0("data/RNA_Tramp_to_pe_2018/temp2/group_", names[m],"-",names[i], ".csv")
                write.table(results(dds, contrast = c("group", names[m], names[i])),
                            file, sep = "\t",
                            quote = F,
                            col.names = NA)
        }
        m <- m+1
}
# dir.create("data/RNA_Tramp_to_pe_2018/123")
## stop here

### 5-11-2018 exploring with pathways


library(readxl)
p1 <- read_excel("~/github/rw2017/data/RNA_Tramp_to_pe_2018/new/IPA/Cannonical Pathways - group TRAMP_24wks vs WILD_24wks log2FC-0.8_q-0.1_2860-genes.xls", 
                 skip = 1)
p1[1,5]

p2 <- read_excel("data/RNA_Tramp_to_pe_2018/new/IPA/Cannonical Pathways - group PEITC_24wks vs TRAMP_24wks p-0.1_2485-genes.xls",
                 skip = 1)



colnames(p1)[2:5] <- paste0(colnames(p1)[2:5], ".tramp-vs-wt")
colnames(p2)[2:5] <- paste0(colnames(p2)[2:5], ".peitc-vs-tramp")

pc <- merge(p1, p2)
pc1 <- na.omit(pc)

# Keep rows with opposite activation z scores.
keep <- (pc1$`z-score.tramp-vs-wt` * pc1$`z-score.peitc-vs-tramp`) < 0
pc2 <- pc1[keep,]
colnames(pc2)
# Re-order
pc2 <- pc2[, c("Ingenuity Canonical Pathways", "-log(p-value).tramp-vs-wt", "Ratio.tramp-vs-wt", "z-score.tramp-vs-wt",
               "-log(p-value).peitc-vs-tramp", "Ratio.peitc-vs-tramp", "z-score.peitc-vs-tramp",
               "Molecules.tramp-vs-wt", "Molecules.peitc-vs-tramp")]

# Filter by p value.
pc3 <- pc2[pc2$`-log(p-value).tramp-vs-wt` >= 2 & pc2$`-log(p-value).peitc-vs-tramp` >= 2,]
pc4 <- pc2[pc2$`-log(p-value).tramp-vs-wt` >= 2.15 & pc2$`-log(p-value).peitc-vs-tramp` >= 2.1,]
# Save to disk
fwrite(pc3, "data/RNA_Tramp_to_pe_2018/new/IPA/contrast_pathways_tramp-wt.peitc-tramp_significant.csv", sep = "\t")


###
plotCounts(dds, gene = "Cdc7", intgroup = "group")

## extract pathway names and moleculars with loop and rbind. credit to Davit

library("data.table")
pc3
class(pc3)
# "data.frame"
colnames(pc3)
tmp <- strsplit(pc3$`Molecules.peitc-vs-tramp`, split = ",")
class(tmp)
# "list"
length(tmp)
# 26
out <- list()
for (i in 1:length(tmp)){
        out[[i]] <- data.table(Pathway = rep(pc3$`Ingenuity Canonical Pathways`[i],
                                             length(tmp[[i]])),
                               gene = tmp[[i]])
}
out
tt1 <- do.call("rbind", out)


########## END rbind ####



strsplit(pc3[1, 8], ",")
tolower(strsplit(pc3[1, 8], ","))
chartr("A-Z", "a-z", strsplit(pc3[1, 8], ","))



pc3$`Ingenuity Canonical Pathways`
pc3[3,1] <- "Cell Cycle G1-S Checkpoint Regulation"
pc3[9,1] <- "LPS-IL-1 Mediated Inhibition of RXR Function"


for (j in 1:nrow(pc3)){ # nrow(pc3)
        
i <- strsplit(pc3[j, 8], ",")
p <-pc3[j, 1]
dir.create(paste0("data/RNA_Tramp_to_pe_2018/new/IPA/dds1/", p))
i2 <- i[[1]]
for (i in i2) { #strsplit(pc3[1, 8], ",")
        print(i)
        g <- (paste0(substring(i, 1,1), tolower(substring(i, 2, nchar(i)))))
        print(g)
        if (g %in% rownames(dds)){
                if (!file.exists(paste0("data/RNA_Tramp_to_pe_2018/new/IPA/dds1/", p, "/", g, ".tif"))){
                tiff(file =paste0("data/RNA_Tramp_to_pe_2018/new/IPA/dds1/", p, "/", g, ".tif"), width = 6, height = 3, units = "in", res = 96)
                 plotCounts(dds1, gene = g, intgroup = "group", col = rep(1:3, each = 2),
                            cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5, xlab = "Groups")
                  dev.off()
                }
        }
}
}


plotCounts(dds1, gene = "Fos", intgroup = "group",
           cex.lab=1.5, cex.axis=1, cex.main=2, cex.sub=3, cex = 1, las =3)


### subset DDS. 5-12-18 

colData(dds)
dds1 <- dds[, paste0("TR", c("01","02","04","05",19,20))]
dds1$group <- droplevels(dds1$group)
design(dds1)
dds1 <- DESeq(dds1)


vsd1 <- vst(dds1, blind = F)
rld1 <- rlog(dds1, blind = F)
head(assay(vsd1), 3)
head(assay(rld1), 3)
head(assay(dds1), 3)





# Check effects of transformation on the variance
ntd1 <- normTransform(dds1)

library("vsn")
meanSdPlot((assay(ntd1)))
meanSdPlot((assay(vsd1)))
meanSdPlot((assay(rld1)))

# install.packages("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(dds1,normalized=TRUE)),
                decreasing=TRUE)[1:500]
df <- as.data.frame(colData(dds1)[,"group", drop = F])
pheatmap(assay(vsd1)[select,], cluster_rows=T, show_rownames=F,
         cluster_cols=F, annotation_col = df,
         #color = heat.colors(1024), #
         border_color = NA,
         scale = "row",
         legend = T)
dev.off()

pheatmap(assay(vsd1)[select,])
#### End subset DDS



plotCounts(dds, gene = g, intgroup = "group")


length(strsplit(pc3[1, 8][1], ","))
###
colnames(dds)
rownames(counts(dds))
rownames(dds)
## Use apply to achieve this goal....

i <- strsplit(pc3[1, 8], ",")
i2 <- i[[1]]
g <- (paste0(substring(i, 1,1), tolower(substring(i, 2, nchar(i)))))

pa <-function(x){paste0(x, "-a")
        cat(x)}

lapply(i[[1]], pa)

lapply(i2, pa)

length(i[[1]])

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


