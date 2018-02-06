# Automatically annotate DMR tables produced by DMRfinder.
#R Wu. 2-6-2018
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")

anno.konglab <- function(file){
  library(GenomicFeatures)
  library(ChIPseeker)
  txdb <- makeTxDbFromGFF('genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf', #This path is for running on the workstation in room 230. Change if necessary.
                          format='gtf',
                          organism='Mus musculus')
  peak <- readPeakFile(file)
  peakAnno <- annotatePeak(peak, TxDb=txdb)
  # add annotations
  peak2 <- read.table(file, header=T)
  peak2 <- peak2[peak2$chr!="chrM",] #This is required for some files...
  peak2$feature <- peakAnno@anno$annotation
  peak2$distance <- peakAnno@anno$distanceToTSS
  peak2$gene <- peakAnno@anno$geneId
  file_dir <- dirname(file)
  file_base <- basename(file)
  
  write.table(peak2, paste(file, "_anno.csv", sep = ""), sep='\t', quote=F, row.names=F)
  
  return()
}

anno.konglab("~/Documents/results.csv")



