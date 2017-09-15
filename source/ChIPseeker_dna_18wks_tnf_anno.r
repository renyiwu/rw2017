# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")
library(GenomicFeatures)
txdb <- makeTxDbFromGFF('Genes/mm9/genes.gtf', format='gtf', organism='Mus musculus')
library(ChIPseeker)
peak <- readPeakFile('data/combined_CpG_all_John.csv')
peakAnno <- annotatePeak(peak, TxDb=txdb)
peak2 <- read.table("data/combined_CpG_all_John.csv", header =T)
peak3 <- peak2[peak2$chr!="chrM",]
#
peak3$feature <- peakAnno@anno$annotation
peak3$distance <- peakAnno@anno$distanceToTSS
peak3$gene <- peakAnno@anno$geneId
write.table(peak3, "data/combined_CpG_all_John_anno.csv", sep='\t', quote=F, row.names=F)
#
peak4 <- peak3[peak3$gene == "Tnf",]
write.table(peak4, "data/combined_CpG_Tnf_John_anno.csv", sep='\t', quote=F, row.names=F)