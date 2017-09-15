#workflow for Methyl-seq, shan. use (fixed by John) bismark v 0.16.1-dev

# align with bismark
cd /dir/to/fastq/files/
perl ../bis-jphn/bismark --multicore 1 ~/Documents/mm9.2_for_bismark/ S*.fastq.gz >log.txt

#2 extract 
$ samtools view -h S1_bismark_bt2.bam | python ../../DMRfinder/extract_CpG_data.py -i -  -o S1.cov

###try this loop. worked.
for i in 2 3 4 5 6
do
samtools view -h S"$i"_bismark_bt2.bam | python ../../DMRfinder/extract_CpG_data.py -i -  -o S"$i".cov
done

##3, combine reads into regions. THis step finishes in seconds.
$ python ../../DMRfinder/combine_CpG_sites.py -o combined.csv S1.cov S2.cov S3.cov S4.cov S5.cov S6.cov

#4, run DMRfinder.r 
administrator@SOP-1482:~/Documents/Renyi/Shan_methyl_2017$ Rscript ../../DMRfinder/findDMRs.r -i combined.csv -o results.csv -n Control,TPA,Treatment S1,S2,S3 S4,S5,S6

Rscript ../../DMRfinder/findDMRs.r -i combined.csv -o results_all.csv -n Control,TPA,TPA+CA,TPA+FX,TPA+CDDO,TPA+mITC S1 S2 S3 S4 S5 S6

#fixed installation

#optional, install DSS (only for the first time) Use this link instructions to fix some problems.
# https://askubuntu.com/questions/359267/cannot-find-curl-config-in-ubuntu-13-04
sudo R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("DSS")

R
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")
library(GenomicFeatures)
txdb <- makeTxDbFromGFF('/home/administrator/Documents/mm9.2_for_bismark/Mus_musculus_UCSC_mm9.2/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf', format='gtf', organism='Mus musculus')
library(ChIPseeker)
peak <- readPeakFile('results.csv')
peakAnno <- annotatePeak(peak, TxDb=txdb)

# add annotations to comb_def.csv
peak2 <- read.table('results.csv', header=T)
peak2$feature <- peakAnno@anno$annotation
peak2$distance <- peakAnno@anno$distanceToTSS
peak2$gene <- peakAnno@anno$geneId
write.table(peak2, 'results_anno.csv', sep='\t', quote=F, row.names=F)
#or 
write.table(peak2, 'results_anno_T.csv', sep='\t', quote=T, row.names=F)
