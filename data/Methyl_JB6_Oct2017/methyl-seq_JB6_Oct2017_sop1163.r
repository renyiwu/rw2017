
peak <- readPeakFile('~/Methyl_Oct_JB6/mm10/cutadapt/JB6_methyl_oct_combined_mm10_ca_all-samples_fr3s2c1x2_CpG.csv')

peakAnno <- annotatePeak(peak, TxDb=txdb)
# add annotations to results.csv
peak2 <- read.table('~/Methyl_Oct_JB6/mm10/cutadapt/JB6_methyl_oct_combined_mm10_ca_all-samples_fr3s2c1x2_CpG.csv', sep = "\t", header=T)

peak2 <- peak2[peak2$chr!="chrM",]
peak2$gene <- peakAnno@anno$geneId
peak2$feature <- peakAnno@anno$annotation
peak2$distance <- peakAnno@anno$distanceToTSS
peak2$transcript <- peakAnno@anno$transcriptId
#Save output file.
write.table(peak2, '~/Methyl_Oct_JB6/mm10/cutadapt/JB6_methyl_oct_combined_mm10_ca_all-samples_fr3s2c1x2_CpG_anno.csv', sep='\t', quote=T, row.names=F) #Use "quote=T" to avoid format issues.

# python ~/tools/DMRfinder/combine_CpG_sites.py --helpUsage: python combine_CpG_sites.py  [options]  -o <output>  [<input>]+
#         [<input>]+    One or more files, each listing methylation counts
# for a particular sample
# -o <output>   Output file listing genomic regions and combined
# methylation counts for each sample
# Options:
#         To consider a particular CpG:
#         -r <int>    Min. number of counts at a position (def. 3)
# -s <int>    Min. number of samples with -r counts (def. 1)
# To analyze a region of CpGs:
#         -d <int>    Max. distance between CpG sites (def. 100)
# -c <int>    Min. number of CpGs in a region (def. 3)
# -x <int>    Max. length of a region (def. 500)
# To report a particular result:
#         -m <int>    Min. total counts in a region (def. 20)
# Other:
#         -f          Report methylation fraction for each sample
# -b          Memory-saving option (may take longer)
# -e <file>   File listing ordered chromosome names (comma-
#                                                            separated; used only with -b option)



# library("xlsx")
 
dt <- read.table("data/Methyl_JB6_Oct2017/JB6_methyl_oct_combined_mm10_ca_all-samples_fr3s2c5_anno.csv", header = T)
colnames(dt) <- c("chr", "start", "end", "CpG", "RW3", "RW4", "RW5", "RW6", "RW7", "RW1", "RW2", "gene", "feature", "distance", "transcript")
colnames(dt)
dt1 <- cbind(dt[1:4], dt[10:11], dt[5:9], dt[12:15])
colnames(dt1)
write.table(dt1, "data/Methyl_JB6_Oct2017/use.JB6_methyl_oct_combined_mm10_ca_all-samples_fr3s2c5_anno.csv", sep = "\t", quote = F, row.names = F)

#Use this: use.JB6_methyl_oct_combined_mm10_ca_MITC-CA_fr3s2c5_anno.xlsx

names(d1)[3:8] <- rep(c("23","24"), 3, each = 2 )


getwd()


# Merge data
dna <- read.table("data/Methyl_JB6_Oct2017/use.JB6_methyl_oct_combined_mm10_ca_all-samples_fr3s2c5_anno.csv",
                  header = T,
                  quote = "\"",
                  sep = "\t")
dna <- dna[order(dna[1], dna[2]),] # Sort on chr then start

dna <- cbind(DMR = 1:nrow(dna), dna)




dna1 <- dna[which(dna$gene %in% rownames(df1) & substring(dna$feature, 1, 3) == "Pro"),]
colnames(dna1)
unique(as.character(dna1$gene))

dmr <- dna1[, c("RW3", "RW4", "RW7", "RW2")]
colnames(dmr) <- c("Control", "TPA", "TPA_CDDO", "TPA_MITC")
rownames(dmr) <- paste0(1:nrow(dna1), "_", dna1$gene)
dmr <- na.omit(dmr)
pheatmap(dmr, cluster_rows= T,
         show_rownames=T,
         cluster_cols=F,
         border_color = NA,
         color = bluered(255)) # or color = gray.colors(255, start = 0.9, end = 0.3, gamma = 2.2, alpha = NULL)
colnames(dmr)

# load RNA data
load.rna <- function(file, study){
        rna <- read.table(file, header = T, sep = "\t")
        rna <- rna[c(1,5,7,8)]
        colnames(rna) <- c("gene", paste(c("log2FC", "pvalue", "qvalue"), study, sep = "."))
        return(rna)
}
rna1 <-load.rna("data/RNA_JB6_Oct2017/DEGseq/RW4-RW3_TPA0-Control.csv", "TPA-Control")
rna2 <- load.rna("data/RNA_JB6_Oct2017/DEGseq/RW5-RW4_CA-TPA0.csv", "CA-TPA")
rna3 <- load.rna("data/RNA_JB6_Oct2017/DEGseq/RW2-RW4_mITC-TPA0.csv", "MITC-TPA")

# merge
d1 <- merge(dna, rna1, all.x = T) # by = "gene"
d12 <- merge(d1, rna2, all.x = T)
d123 <- merge(d12, rna3, all.x = T)
d123 <- d123[order(d123[2], d123[3]),] # order again. may skip this line.
write.table(d123, "data/Methyl_JB6_Oct2017/use.JB6_methyl_oct_combined_mm10_ca_MITC-CA_fr3s2c5_anno_exp2.csv",
            row.names = F,
            quote = F,
            sep = "\t") # Very slow


library("data.table")
fwrite(d123, "data/Methyl_JB6_Oct2017/use.JB6_methyl_oct_combined_mm10_ca_MITC-CA_fr3s2c5_anno_expfo.csv",
            row.names = F,
            quote = F,
            sep = "\t") # Super fast.


