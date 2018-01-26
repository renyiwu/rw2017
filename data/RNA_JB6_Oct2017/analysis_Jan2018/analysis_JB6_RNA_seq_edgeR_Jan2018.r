#https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)
dt <- read.table("data/RNA_JB6_Oct2017/RW_all_primary.dedup_hisat2_new.csv", header = T, row.names = "Geneid")
c <- dt[-(1:4)]
colnames(c) <- c("Length", paste("RW", 1:7, sep = ""))
head(c)
groups <- c("TPA1", "mITC", "Control", "TPA0", "CA", "FX", "CDDO")
groups2 <- c("TPA1", "TPA1", "Control", "Control", "CA", "CA", "CA")
y <- DGEList(counts=c[2:8], group=factor(groups))
y
dim(y)
y.full <- y
head(y$counts)
head(cpm(y))

apply(y$counts, 2, sum)
y
sum(y$counts)

keep <- rowSums(cpm(y)>10) >= 2
y <- y[keep,]
dim(y)
y
y$samples$lib.size <- colSums(y$counts)
y$samples

y <- calcNormFactors(y)

plotMDS(y, method="bcv", col=as.numeric(y$samples$group))
legend("bottomleft", as.character(unique(y$samples$group)), col=as.numeric(y$samples$group), pch=20)

#No replicates. y1 <- estimateCommonDisp(y, verbose=T)
dev.off()
# # 
# design <-model.matrix(~groups)
# design2 <-model.matrix(~groups2)
# y <- calcNormFactors(y)
# y
# y <- estimateDisp(y, design2)
# y

bcv <- 0.3
et <- exactTest(y, dispersion = bcv^2)
write.table(topTags(et, n = nrow(y$counts)),
            file = "data/RNA_JB6_Oct2017/analysis_Jan2018/test_bcv0.3.csv",
            sep = "\t", col.names = NA, row.names = T) 

bcv <- 0.3
gt <- glmFit(y, dispersion = bcv^2)
gt <- glmLRT(gt)
write.table(topTags(gt, n = nrow(y$counts)),
            file = "data/RNA_JB6_Oct2017/analysis_Jan2018/test_bcv0.3-gt.csv",
            sep = "\t", col.names = NA, row.names = T) 

