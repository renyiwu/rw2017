#Calculation of LogFC of samples with no replicates using edgeR  #RPKM, FPKM see the end section.
#Renyi Wu 1-25-2018
#Ref: https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html


# Method #1, using house keeping genes as a reference for xalculation of the dispersion.
# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library(edgeR)
dt <- read.table("data/RNA_JB6_Oct2017/RW_all_primary.dedup_hisat2_new.csv", header = T)#, row.names = "Geneid")
dt <- dt[-(2:5)] #remove unneeded columns.
colnames(dt) <- c("Geneid", "Length", paste("RW", 1:7, sep = "")) #Change as needed.
groups <- c("TPA1", "mITC", "Control", "TPA", "CA", "FX", "CDDO") #Need to match the order of samples in the above dataset
y <- DGEList(counts = dt[3:9], group = relevel(factor(groups),"Control"), #set "Control" as the reference group.
             genes = dt[1:2]) #data.frame(geneid = rownames(c), Lengh = c$Length))
#y <- y[rowSums(cpm(y)>10) >= 2,] # Only keep genes with CPM (counts per million reads) values greater than 10 in at least 2 groups.Arbitrary setting.
y <- y[rowSums(y$counts>3) >= 2,]# Only keep genes that have more than 3 counts in at least 2 groups. Arbitrary setting.
y$samples$lib.size <- colSums(y$counts) #This step is required if filtering was performed with the last command.
y <- calcNormFactors(y) #Not required but suggested.
# Define design
design <- model.matrix(~0+y$samples$group) #Leading "~0" is for ... check the manual!
colnames(design) <- levels(y$samples$group)

y1 <- y
y1$samples$group <-1
design1 <- model.matrix(~0+y1$samples$group)
housekeeping <- c("Actb", "Gapdh") #ideally more genes are needed.
y1 <- estimateDisp(y1[y1$genes$Geneid == housekeeping,], design = design1, trend = "none", tagwise = F) #or y0 <- estimateDisp(y1[housekeeping,])
y$common.dispersion <- y1$common.dispersion

#1.1 Exact test
et <- exactTest(y, pair = c("TPA1", "mITC" )) # Need to be two groups at a time.
write.table(topTags(et, n = nrow(y$counts), sort.by = "none"), # may sort.by "Pvalue".
            file = "data/RNA_JB6_Oct2017/analysis_Jan2018/1-25/mITC-TPA1_housekeeping_exact-test.csv",
            sep = "\t", col.names = T, row.names = F) 
 
#1.2 glm model
mITCvsTPA1 <- makeContrasts(mITC-TPA1, levels = design)#Make contrast
gf <- glmFit(y, design = design)
gt <- glmLRT(gf, contrast = mITCvsTPA1)
write.table(topTags(gt, n = nrow(y$counts), sort.by = "none"),
            file = "data/RNA_JB6_Oct2017/analysis_Jan2018/1-25/mITC-TPA1_housekeeping_glmLRT-test.csv",
            sep = "\t", col.names = T, row.names = F) 

######################################################################
#Method #2, manually set a bcv value, then use bvc^2 as dispersion.
######################################################################
# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library(edgeR)
dt <- read.table("data/RNA_JB6_Oct2017/RW_all_primary.dedup_hisat2_new.csv", header = T)#, row.names = "Geneid")
dt <- dt[-(2:5)]
colnames(dt) <- c("Geneid", "Length", paste("RW", 1:7, sep = ""))
groups <- c("TPA1", "mITC", "Control", "TPA", "CA", "FX", "CDDO") #Need to match the order of samples in the above dataset
y <- DGEList(counts = dt[3:9], group = relevel(factor(groups),"Control"), #set "Control" as the reference group.
             genes = dt[1:2]) #data.frame(geneid = rownames(c), Lengh = c$Length))
#y <- y[rowSums(cpm(y)>10) >= 2,] # Only keep genes with CPM (counts per million reads) values greater than 10 in at least 2 groups.Arbitrary setting.
y <- y[rowSums(y$counts>3) >= 2,]# Only keep genes that have more than 3 counts in at least 2 groups. Arbitrary setting.
y$samples$lib.size <- colSums(y$counts) #This step is required if filtering was performed with the last command.
y <- calcNormFactors(y)
# Define design
design <- model.matrix(~0+y$samples$group)
colnames(design) <- levels(y$samples$group)

#2.1 Exact Test
bcv <- 0.1
et <- exactTest(y, pair = c("TPA1", "mITC" ), dispersion = bcv^2)
tt <- topTags(et, n = nrow(y$counts), sort.by = "none")
write.table(tt,
            file = "data/RNA_JB6_Oct2017/analysis_Jan2018/1-25/mITC-TPA1_bcv-0.1_exact-test.csv",
            sep = "\t", col.names = T, row.names = F) 


#2.2 glm model
# glmLRT Test
# Genewise Negative Binomial Generalized Linear Models
bcv <- 0.1 #set to any arbitrary value
mITCvsTPA1 <- makeContrasts(mITC-TPA1, levels = design) #make contrast
gf <- glmFit(y, design= design, dispersion = bcv^2)
gt <- glmLRT(gf, contrast = mITCvsTPA1)
write.table(topTags(gt, n = nrow(y$counts), sort.by = "none"),
            file = "data/RNA_JB6_Oct2017/analysis_Jan2018/1-25/mITC-TPA1_bcv-0.1_glmLRT-test.csv",
            sep = "\t", col.names = T, row.names = F) 



## Get RPKM values
# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library(edgeR)
dt <- read.table("data/RNA_JB6_Oct2017/RW_all_primary.dedup_hisat2_new.csv", header = T)#, row.names = "Geneid")
dt <- dt[-(2:5)]
colnames(dt) <- c("Geneid", "Length", paste("RW", 1:7, sep = ""))
groups <- c("TPA1", "mITC", "Control", "TPA", "CA", "FX", "CDDO") #Need to match the order of samples in the above dataset
y <- DGEList(counts = dt[3:9], group = relevel(factor(groups),"Control"), #set "Control" as the reference group.
             genes = dt[1:2]) #data.frame(geneid = rownames(c), Lengh = c$Length))

rpkm <- data.frame(rpkm(y, y$genes$Length))
rpkm <- cbind(Geneid = y$genes$Geneid, rpkm) #Use cbind to put "gene id" at the first column.
write.table(rpkm,
            file = "data/RNA_JB6_Oct2017/analysis_Jan2018/1-25/RPKM.csv",
            sep = "\t", col.names = T, row.names = F) 




#Tested on SOP-1163 Ubuntu 16.04
#Renyi Wu and others
#1-25-2018
