###Gene annotation with ChIPseeker.
#Renyi Wu 7-12-2017
#Credits to Davit, John...
#
#1, if input file has calculated methyl diff data, i.e. output of DMRfinder.r, use this script:
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")
# biocLite("ChIPseeker")
#biocLite("matrixStats")
#biocLite("rtracklayer")
#biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
#biocLite("Rcpp")


library(GenomicFeatures)
txdb <- makeTxDbFromGFF('Genes/genes.gtf', format='gtf', organism='Mus musculus')
library(ChIPseeker)
peak <- readPeakFile('data/yue/dmr_yue_2nd_dedup.csv')
peak <- readPeakFile('data/yue/dmr_yue_2nd_ori.csv')#results.csv is the output of DMRfinder.r 
peakAnno <- annotatePeak(peak, TxDb=txdb)
# add annotations to results.csv
peak2 <- read.table('data/yue/dmr_yue_2nd_dedup.csv', sep = "\t", header=T)
peak2 <- read.table('data/yue/dmr_yue_2nd_ori.csv', sep = "\t", header=T)
peak2 <- peak2[peak2$chr!="chrM",]
peak2$gene <- peakAnno@anno$geneId
peak2$feature <- peakAnno@anno$annotation
peak2$distance <- peakAnno@anno$distanceToTSS
#Save output file.
write.table(peak2, 'data/yue/dmr_yue_2nd_dedup_anno.csv', sep='\t', quote=T, row.names=F) #Use "quote=T" to avoid format issues.
write.table(peak2, 'data/yue/dmr_yue_2nd_ori_anno.csv', sep='\t', quote=T, row.names=F) #Use "quote=T" to avoid format issues.

#2, if input file has only combined counts of all groups, i.e., output of Combine_CpG.py, use this script:
com1 <- read.table("data/shan_methyl/jb6_combined.csv", header = T, sep = "\t")
com1 <- com1[com1$chr!="chrM",] #remove chromosome M, required.
com1$S1.mu <- com1$S1.X/com1$S1.N 
com1$S2.mu <- com1$S2.X/com1$S2.N
com1$S3.mu <- com1$S2.X/com1$S3.N
com1$S4.mu <- com1$S2.X/com1$S4.N
com1$S5.mu <- com1$S2.X/com1$S5.N
com1$S6.mu <- com1$S2.X/com1$S6.N
#Add more as necessary.
com1$S2_S1.diff <- com1$S2.mu - com1$S1.mu
com1$S3_S2.diff <- com1$S3.mu - com1$S2.mu
com1$S4_S2.diff <- com1$S4.mu - com1$S2.mu
com1$S5_S2.diff <- com1$S5.mu - com1$S2.mu
com1$S6_S2.diff <- com1$S6.mu - com1$S2.mu
#Add more as necessary.
#
write.table(com1, "data/combined-M.csv", sep='\t', quote=T, row.names=F) #Use "quote=T" to avoid format issues.
#
library(GenomicFeatures)
txdb <- makeTxDbFromGFF('/path/to/your/gtf/file/Mus_musculus_UCSC_mm9.2/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf') #, format='gtf', organism='Mus musculus')
library(ChIPseeker)
peak <- readPeakFile('data/combined-M.csv') #Or read from "data/combined.csv"
peakAnno <- annotatePeak(peak, TxDb=txdb)
# add annotations
peak2 <- read.table('data/combined-M.csv', sep = "\t", header=T)
peak2$gene <- peakAnno@anno$geneId
peak2$feature <- peakAnno@anno$annotation
peak2$distance <- peakAnno@anno$distanceToTSS
peak2$geneLenth <- peakAnno@anno$geneLength
peak2$geneStrand <- peakAnno@anno$geneStrand
#Save output file.
write.table(peak2, 'data/combined_anno_2.csv', sep='\t', quote=T, row.names=F) #Use "quote=T" to avoid format issues.
