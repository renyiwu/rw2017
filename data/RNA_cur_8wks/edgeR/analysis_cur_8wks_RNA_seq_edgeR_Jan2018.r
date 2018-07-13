# https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library(edgeR)
dt <- read.table("data/RNA_cur_8wks/RNA_8wks_primary.dedup.csv", header = T)
dt <- dt[-(2:5)]
colnames (dt) <- c("Geneid", "Length", paste("C", c("01",26,29,"02",42,46,65,70,75,79), sep = ""))

#Define groups (and batches)
group <- factor(c("Con", "AOM_DSS", "AOM_DSS", "Con", "AOM_DSS_Cur", "AOM_DSS_Cur", "DSS", "DSS", "DSS_Cur", "DSS_Cur"))
group <- relevel(group, ref = "Con")
# > group
# [1] Con         AOM_DSS     AOM_DSS     Con         AOM_DSS_Cur AOM_DSS_Cur DSS         DSS         DSS_Cur    
# [10] DSS_Cur    
# Levels: Con AOM_DSS AOM_DSS_Cur DSS DSS_Cur
batch <- factor(c(rep("1",3), "2", rep("1", 6)))
# > batch
# [1] 1 1 1 2 1 1 1 1 1 1
# Levels: 1 2

y <- DGEList(counts = dt[3:12], group = group, genes = dt[1:2])

#Define design
#  design <- model.matrix(~group+batch) # This style works only with "coef" setting as in glmLRT(fit, coef = 2)
# Or 
design <- model.matrix(~0+group+batch)
rownames(design) <- colnames(dt[-(1:2)])
# design <- model.matrix(~0+group+batch) ## Here the leading "0+" tells not to produce an "Intercept" group.
# This method is not recommended if more than one variance is included in the model. Eg, groups, time points, batches effects etc.

# > design
# groupCon groupAOM_DSS groupAOM_DSS_Cur groupDSS groupDSS_Cur batch2
# C1.dedup.bam         1            0                0        0            0      0
# C26.dedup.bam        0            1                0        0            0      0
# C29.dedup.bam        0            1                0        0            0      0
# C2.dedup.bam         1            0                0        0            0      1
# C42.dedup.bam        0            0                1        0            0      0
# C46.dedup.bam        0            0                1        0            0      0
# C65.dedup.bam        0            0                0        1            0      0
# C70.dedup.bam        0            0                0        1            0      0
# C75.dedup.bam        0            0                0        0            1      0
# C79.dedup.bam        0            0                0        0            1      0
# attr(,"assign")
# [1] 1 1 1 1 1 2
# attr(,"contrasts")
# attr(,"contrasts")$group
# [1] "contr.treatment"
# 
# attr(,"contrasts")$batch
# [1] "contr.treatment"

# Data filtering based on counts (cpm, counts per million reads)
keep <- rowSums(cpm(y)>3) >= 2
#or
# keep <- rowSums(y$counts>5) >= 2
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts) #required
y <- calcNormFactors(y) #required

#Check sample clustering with PCA plots.
plotMDS(y, method="bcv", col=as.numeric(y$samples$group))
legend("bottomleft", as.character(unique(y$samples$group)), col=as.numeric(unique(y$samples$group)), pch=20)

# Estimate dispersion
y <- estimateDisp(y, design = design, robust = T)
y$common.dispersion
plotBCV(y)

# Model fitting
fit <- glmFit(y, design)

#USe coef to compare any other group(s) with the first (referece, or intercept) group.
#Can not be used to compare two non-reference gropoups.
#Reference group is the group with the name in lowest alphabetic order. 
#Can be reordered with "group <- relevel(group, "DSS")" etc
# 

## Method 1  ONLY WORKS with "design <- model.matrix(~group+batch)". NO leading "0+"
#
# > colnames(design)
# [1] "(Intercept)"      "groupAOM_DSS"     "groupAOM_DSS_Cur" "groupDSS"         "groupDSS_Cur"     "batch2"          
#  
#
lrt <- glmLRT(fit, coef = 2) #ONLY WORKS with "design <- model.matrix(~group+batch)". NO leading "0+"
tt02 <- topTags(lrt, n = nrow(y$counts), sort.by = "none") #
# Then save data with write.table etc
#
# > colnames(design)
# [1] "(Intercept)"      "groupAOM_DSS"     "groupAOM_DSS_Cur" "groupDSS"         "groupDSS_Cur"     "batch2"          


# ##########Method 2 #######################
## use contrast parameter to assign the TWO groups to be compared.
# Either one of the two samples below works.

# > colnames(design)
# [1] "groupCon"         "groupAOM_DSS"     "groupAOM_DSS_Cur" "groupDSS"         "groupDSS_Cur"    
# [6] "batch2"          

# 2.1
contrast <- makeContrasts(groupDSS_Cur-groupDSS, levels = design)
lrt1 <- glmLRT(fit, contrast = contrast)
tt1 <- topTags(lrt1, n = nrow(y$counts), sort.by = "none")

# or  2.2
lrt2 <- glmLRT(fit, contrast = c(0, 0, 1, -1, 0, 0))
tt2 <- topTags(lrt2, n = nrow(y$counts), sort.by = "none")

# > colnames(design)
# [1] "groupCon"         "groupAOM_DSS"     "groupAOM_DSS_Cur" "groupDSS"         "groupDSS_Cur"    
# [6] "batch2"

#Write data to disk.
write.table(tt2,
            file = "data/RNA_cur_8wks/tt2_DSSCur-DSS.csv",
            sep = "\t", col.names = T, row.names = F)

#1 AOM_DSS vs Control
lrt2 <- glmLRT(fit, contrast = c(-1, 1, 0, 0, 0, 0))
tt2 <- topTags(lrt2, n = nrow(y$counts), sort.by = "none")
write.table(tt2,
            file = "data/RNA_cur_8wks/edgeR/AOM_DSS-Control.csv",
            sep = "\t", col.names = T, row.names = F)

#2 DSS vs Control
lrt2 <- glmLRT(fit, contrast = c(-1, 0, 0, 1, 0, 0))
tt2 <- topTags(lrt2, n = nrow(y$counts), sort.by = "none")
write.table(tt2,
            file = "data/RNA_cur_8wks/edgeR/DSS-Control.csv",
            sep = "\t", col.names = T, row.names = F)


#3 AOM_DSS_Cur vs AOM_DSS
lrt2 <- glmLRT(fit, contrast = c(0, -1, 1, 0, 0, 0))
tt2 <- topTags(lrt2, n = nrow(y$counts), sort.by = "none")
write.table(tt2,
            file = "data/RNA_cur_8wks/edgeR/AOM_DSS_Cur-AOM_DSS.csv",
            sep = "\t", col.names = T, row.names = F)

#4 DSS_Cur vs DSS
lrt2 <- glmLRT(fit, contrast = c(0, 0, 0, -1, 1, 0))
tt2 <- topTags(lrt2, n = nrow(y$counts), sort.by = "none")
write.table(tt2,
            file = "data/RNA_cur_8wks/edgeR/DSS_Cur-DSS.csv",
            sep = "\t", col.names = T, row.names = F)

#5 AOM_DSS vs DSS
lrt2 <- glmLRT(fit, contrast = c(0, 1, 0, -1, 0, 0))
tt2 <- topTags(lrt2, n = nrow(y$counts), sort.by = "none")
write.table(tt2,
            file = "data/RNA_cur_8wks/edgeR/AOM_DSS-DSS.csv",
            sep = "\t", col.names = T, row.names = F)




##### End of Differentially Expression (DE) analysis ######

#######################
# Chapter 2.  Get RPKM

# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library(edgeR)
dt <- read.table("data/RNA_cur_8wks/RNA_8wks_primary.dedup.csv", header = T)#, row.names = "Geneid")
dt <- dt[-(2:5)]
#colnames(c) <- c("Length", paste("RW", 1:7, sep = ""))

group <- factor(c("Con", "AOM_DSS", "AOM_DSS", "Con", "AOM_DSS_Cur", "AOM_DSS_Cur", "DSS", "DSS", "DSS_Cur", "DSS_Cur"),
                c("Con", "DSS", "DSS_Cur", "AOM_DSS", "AOM_DSS_Cur"))
batch <- factor(c(rep("1",3), "2", rep("1", 6)))
y <- DGEList(counts = dt[3:12], group = group, genes = dt[1:2])

rpkm <- rpkm(y, y$genes$Length)
rpkm <- cbind(y$genes, rpkm)
write.table(rpkm,
            file = "data/RNA_cur_8wks/edgeR/rpkm_all_cur_8wks_edgeR_1-30-2018.csv",
            sep = "\t", col.names = T, row.names = F) 

# Tested on SOP-1163 Ubuntu 16.04
# R Wu. 1-30-2018