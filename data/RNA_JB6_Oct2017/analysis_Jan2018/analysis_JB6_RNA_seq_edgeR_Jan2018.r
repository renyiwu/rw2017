#https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

# Method #1, using house keeping genes as a reference and calculate the dispersion.
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)
dt <- read.table("data/RNA_JB6_Oct2017/RW_all_primary.dedup_hisat2_new.csv", header = T)#, row.names = "Geneid")
dt <- dt[-(2:5)]
colnames(dt) <- c("Geneid", "Length", paste("RW", 1:7, sep = ""))
groups <- c("TPA1", "mITC", "Control", "TPA", "CA", "FX", "CDDO") #Need to match the order of samples in the above dataset
y <- DGEList(counts = dt[3:9], group = relevel(factor(groups),"Control"), #set "Control" as the reference group.
             genes = dt[1:2]) #data.frame(geneid = rownames(c), Lengh = c$Length))
y <- y[rowSums(cpm(y)>10) >= 2,] # Only keep genes that have more than 10 counts in at least 2 groups.


####
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)
dt <- read.table("data/RNA_JB6_Oct2017/RW_all_primary.dedup_hisat2_new.csv", header = T, row.names = "Geneid")
c <- dt[-(1:4)]
colnames(c) <- c("Length", paste("RW", 1:7, sep = ""))
head(c)
groups <- c("TPA", "mITC", "Control", "TPA", "CA", "FX", "CDDO")
groups2 <- c("TPA1", "TPA1", "Control", "Control", "C1", "CA", "CA")
y <- DGEList(counts=c[2:8], group = relevel(factor(groups),"TPA"), #set "TPA" as the reference group.
             genes = data.frame(geneid = rownames(c), Lengh = c$Length))
y$samples$group
dim(y)
y.full <- y
head(y$counts)
head(cpm(y))

apply(y$counts, 2, sum)
y <- 0
sum(y$counts)

keep <- rowSums(cpm(y)>10) >= 2
y <- y[keep,]
dim(y)
y
y$samples$lib.size <- colSums(y$counts)
y$samples

y <- calcNormFactors(y)

design <- model.matrix(~0+y$samples$group)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y, design = design)
estimateDisp(y, design = design, verbose = T)
plotBCV(y)
#or
estimateDisp(y, verbose = T)
plotBCV(y)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
tt <- topTags(lrt, n = nrow(y$counts), sort.by = "none")

head(tt$table)
y$design

plotMDS(y, method="bcv", col=as.numeric(y$samples$group))
legend("bottomleft", as.character(unique(y$samples$group)), col=as.numeric(unique(y$samples$group)), pch=20)

#No replicates. y1 <- estimateCommonDisp(y, verbose=T)
dev.off()
# # 
# design <-model.matrix(~groups)
# design2 <-model.matrix(~groups2)
# y <- calcNormFactors(y)
# y
# y <- estimateDisp(y, design2)
# y

#Excact Test
#Exact Tests for Differences between Two Groups of Negative-Binomial Counts
bcv <- 0.3

et <- exactTest(y, pair = c("TPA", "mITC" ))#, dispersion = bcv^2)
tt <- topTags(et, n = nrow(y$counts))
head(tt$table)
rownames(tt$table)
write.table(topTags(et, n = nrow(y$counts)),
            file = "data/RNA_JB6_Oct2017/analysis_Jan2018/t1NA-mITC-TPA1.csv",
            sep = "\t", col.names = NA, row.names = T) 
#or
write.table(topTags(et, n = nrow(y$counts)),
            file = "data/RNA_JB6_Oct2017/analysis_Jan2018/t1F-mITC-TPA1.csv",
            sep = "\t", col.names = T, row.names = F) 

# glmLRT Test
# Genewise Negative Binomial Generalized Linear Models
bcv <- 0.3
gt <- glmFit(y, dispersion = bcv^2)
gt <- glmLRT(gt)
write.table(topTags(gt, n = nrow(y$counts)),
            file = "data/RNA_JB6_Oct2017/analysis_Jan2018/test_bcv0.3-gt.csv",
            sep = "\t", col.names = NA, row.names = T) 



### Using housekeeping genes for normalization.
y1 <- y
y1$samples$group <- 1
y0 <- estimateDisp(y1["Actb",], trend = "none", tagwise = F)
#or the following two line:
housekeeping <- c("Actb", "Gapdh")
y0 <- estimateDisp(y1[housekeeping,], trend = "none", tagwise = F)

y$common.dispersion <- y0$common.dispersion

design <- model.matrix(~0+groups, data = y$samples)
fit <- glmFit(y, design)
mITCvsTPA1 <- makeContrasts(groupsmITC-groupsTPA1, levels = design)
lrt <- glmLRT(fit, contrast = mITCvsTPA1)

write.table(topTags(lrt, n = nrow(y$counts)),
            file = "data/RNA_JB6_Oct2017/analysis_Jan2018/test_housekeeping-lrt.csv",
            sep = "\t", col.names = NA, row.names = T) 





# GEt rpkm values

rpkm <- rpkm(y, y$genes$Lengh)
write.table(rpkm,
            file = "data/RNA_JB6_Oct2017/analysis_Jan2018/rpkm_all_edgeR_1-24-2018.csv",
            sep = "\t", col.names = NA, row.names = T) #set col.names to NA to ensure compatibility with MS-EXCEL
#or
write.csv(rpkm, file = "data/RNA_JB6_Oct2017/analysis_Jan2018/rpkm_all_edgeR_1-24-2018_2.csv")

