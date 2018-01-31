# https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library(edgeR)
dt <- read.table("data/RNA_cur_8wks/RNA_8wks_primary.dedup.csv", header = T)
dt <- dt[-(2:5)]

#Define groups (and batches)
group <- factor(c("Con", "AOM_DSS", "AOM_DSS", "Con", "AOM_DSS_Cur", "AOM_DSS_Cur", "DSS", "DSS", "DSS_Cur", "DSS_Cur"))
group <- relevel(group, ref = "Con")
#> group
# [1] Con         AOM_DSS     AOM_DSS     Con         AOM_DSS_Cur AOM_DSS_Cur DSS         DSS         DSS_Cur    
# [10] DSS_Cur    
# Levels: Con DSS DSS_Cur AOM_DSS AOM_DSS_Cur
batch <- factor(c(rep("1",3), "2", rep("1", 6)))
# > batch
# [1] 1 1 1 2 1 1 1 1 1 1
# Levels: 1 2

y <- DGEList(counts = dt[3:12], group = group, genes = dt[1:2])

#Define design

#  design <- model.matrix(~group+batch) # This style works only with "coef" setting as in glmLRT(fit, coef = 2)
# Or 
design <- model.matrix(~0+group+batch) # This style works only with "contrast" setting as in glmLRT(fit, contrast = c(0, -1, 1, 0, 0, 0))
rownames(design) <- colnames(dt[-(1:2)])
# design <- model.matrix(~0+group+batch) ## Here the leading "0+" tells not to produce an "Intercept" group.
# This method is not recommended if more than one variance is included in the model. Eg, groups, time points, batches effects etc.

 
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
### Method 1. ONLY WORKS with "design <- model.matrix(~group+batch)". NO leading "0+"
# seems
lrt <- glmLRT(fit, coef = 2)
tt02 <- topTags(lrt, n = nrow(y$counts), sort.by = "none")

lrt <- glmLRT(fit, coef = 3) #to Compare groups 2 - 5 to the first group.
tt23 <- topTags(lrt, n = nrow(y$counts), sort.by = "none")

#Save data with write.table etc

# Method 2 ONLY WORKS with "design <- model.matrix(~0+group+batch)". WITH leading "0+"
## or use contrast parameter to assign the TWO groups to be compared.
# Either one of the two samples below works.

# 2.1
contrast <- makeContrasts(groupDSS_Cur-groupDSS, levels = design)
lrt1 <- glmLRT(fit, contrast = contrast)
tt1 <- topTags(lrt1, n = nrow(y$counts), sort.by = "none")

# 2.2
lrt2 <- glmLRT(fit, contrast = c(0, -1, 1, 0, 0, 0))
tt2 <- topTags(lrt2, n = nrow(y$counts), sort.by = "none")


#Write data to disk.
write.table(tt2,
            file = "data/RNA_cur_8wks/tt2_DSSCur-DSS.csv",
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