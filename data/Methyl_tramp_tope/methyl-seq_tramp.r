library("data.table")
w16 <- fread("~/Methyl-seq_gw/DMR/combined.TD.c5r5.dmr.16wks.csv_anno.csv")
w24 <- fread("~/Methyl-seq_gw/DMR/combined.TD.c5r5.dmr.24wks.csv_anno.csv")
wall <- fread("~/Methyl-seq_gw/DMR/combined.TD.all-default.csv")
head(wall)
w24.t <- wall[, -(5:34)]
w24.t <- na.omit(w24.t)
w24.t$tumor.24weeks <- w24.t$`TD23-X`/w24.t$`TD23-N`
tumor24 <- w24.t[,-(5:6)]
write.table(tumor24, "data/Methyl_tramp_tope/combined.TD.c5r5.w24.tumor.dmr.csv")


w16.nrf2 <- w16[gene == "Nfe2l2",]
w24.nrf2 <- w24[gene == "Nfe2l2",]

fwrite(w16.nrf2, "data/Methyl_tramp_tope/w16.nrf2.methylation.csv")
fwrite(w24.nrf2, "data/Methyl_tramp_tope/w24.nrf2.methylation.csv")

source("https://bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")
biocLite("ChIPseeker")
n



# 7-11-2018
library("data.table")
w24 <- fread("data/Methyl_tramp_tope/combined.TD.c5r5.dmr.24wks.csv_anno.csv")
colnames(w24)
# [1] "chr"                     "start"                   "end"                    
# [4] "CpG"                     "tramp_24.mu"             "peitc_24.mu"            
# [7] "toco_24.mu"              "wild_24.mu"              "tramp_24..peitc_24.diff"
# [10] "tramp_24..peitc_24.pval" "tramp_24..toco_24.diff"  "tramp_24..toco_24.pval" 
# [13] "tramp_24..wild_24.diff"  "tramp_24..wild_24.pval"  "peitc_24..toco_24.diff" 
# [16] "peitc_24..toco_24.pval"  "peitc_24..wild_24.diff"  "peitc_24..wild_24.pval" 
# [19] "toco_24..wild_24.diff"   "toco_24..wild_24.pval"   "feature"                
# [22] "distance"                "gene"

w24.TRAMP_WILD <- w24[, c("gene", "tramp_24..wild_24.diff", "tramp_24..wild_24.pval", "feature", "distance")]
w24.TRAMP_WILD <- na.omit(w24.TRAMP_WILD)

# Discard records on "Distal Intergenic"
w24.TRAMP_WILD.keep <- w24.TRAMP_WILD[w24.TRAMP_WILD$feature != "Distal Intergenic", ]

# Keep only records on "Promoter ..."
w24.TRAMP_WILD.promoter <- w24.TRAMP_WILD[substring(w24.TRAMP_WILD$feature, 1, 8) == "Promoter", ]

# keep only one record for each gene, based on records on "promoter". See explanation later"
w24.TRAMP_WILD.promoter.unique <- w24.TRAMP_WILD.promoter[!duplicated(w24.TRAMP_WILD.promoter, by = "gene"),]

# Save file
fwrite(w24.TRAMP_WILD.promoter.unique, "data/Methyl_tramp_tope/dmr_TRAMP_24wks-WILD_24_wks_promoter_unique.csv", sep = "\t")
# The diff values are reported as values in WILD group substracted by values in TRAMP group.

# Explanation below#
df <- w24.TRAMP_WILD.promoter
colnames(df)
# [1] "gene"                   "tramp_24..wild_24.diff" "tramp_24..wild_24.pval" "feature"                "distance" 
df.u1 <- df[!duplicated(df, by = "gene"),]
df.u2 <- unique(df, by = "gene")
df[gene == "Ncoa2",]
# gene tramp_24..wild_24.diff tramp_24..wild_24.pval          feature distance
# 1: Ncoa2             -0.1335201              0.0132053 Promoter (1-2kb)     1647
# 2: Ncoa2             -0.0968998              0.0009911 Promoter (<=1kb)     -627
df.u1[gene == "Ncoa2",]
# gene tramp_24..wild_24.diff tramp_24..wild_24.pval          feature distance
# 1: Ncoa2             -0.1335201              0.0132053 Promoter (1-2kb)     1647
df.u2[gene == "Ncoa2",]
# gene tramp_24..wild_24.diff tramp_24..wild_24.pval          feature distance
# 1: Ncoa2             -0.1335201              0.0132053 Promoter (1-2kb)     1647
 
#End explanation"


# try 2
library("data.table")
w24 <- fread("data/Methyl_tramp_tope/combined.TD.c5r5.dmr.fdr.24wks_anno.csv")
colnames(w24)
# [1] "chr"                     "start"                   "end"                     "CpG"                     "tramp_24.mu"            
# [6] "peitc_24.mu"             "toco_24.mu"              "wild_24.mu"              "tramp_24..peitc_24.diff" "tramp_24..peitc_24.pval"
# [11] "tramp_24..peitc_24.fdr"  "tramp_24..toco_24.diff"  "tramp_24..toco_24.pval"  "tramp_24..toco_24.fdr"   "tramp_24..wild_24.diff" 
# [16] "tramp_24..wild_24.pval"  "tramp_24..wild_24.fdr"   "peitc_24..toco_24.diff"  "peitc_24..toco_24.pval"  "peitc_24..toco_24.fdr"  
# [21] "peitc_24..wild_24.diff"  "peitc_24..wild_24.pval"  "peitc_24..wild_24.fdr"   "toco_24..wild_24.diff"   "toco_24..wild_24.pval"  
# [26] "toco_24..wild_24.fdr"    "feature"                 "distance"                "gene"


### 1, extract comparison TRAMP vs WILDTYPE
w24.TRAMP_WILD <- w24[, c("gene", "chr", "start", "end", "CpG","tramp_24.mu", "wild_24.mu", "tramp_24..wild_24.diff", "tramp_24..wild_24.pval", "tramp_24..wild_24.fdr", "feature", "distance")]
w24.TRAMP_WILD <- na.omit(w24.TRAMP_WILD)

# Keep only records on "Promoter ..."
w24.TRAMP_WILD.promoter <- w24.TRAMP_WILD[substring(w24.TRAMP_WILD$feature, 1, 8) == "Promoter", ]

# keep only one record for each gene, based on records on "promoter". See explanation later"
w24.TRAMP_WILD.promoter.unique <- w24.TRAMP_WILD.promoter[!duplicated(w24.TRAMP_WILD.promoter, by = "gene"),]

# Reverse the diff value. 
w24.TRAMP_WILD.promoter.unique$tramp_24..wild_24.diff <- -w24.TRAMP_WILD.promoter.unique$tramp_24..wild_24.diff

# Save file
fwrite(w24.TRAMP_WILD.promoter.unique, "data/Methyl_tramp_tope/dmr_TRAMP-24wks_vs_WILD-24wks_promoter_unique_fdr.csv", sep = "\t")

#

### 2, extract comparison PEITC vs TRAMP
w24.exract <- w24[, c("gene", "chr", "start", "end", "CpG", "tramp_24.mu", "peitc_24.mu", "tramp_24..peitc_24.diff", "tramp_24..peitc_24.pval", "tramp_24..peitc_24.fdr", "feature", "distance")]
w24.exract <- na.omit(w24.exract)

# Keep only records on "Promoter ..."
w24.exract.promoter <- w24.exract[substring(w24.exract$feature, 1, 8) == "Promoter", ]

# keep only one record for each gene, based on records on "promoter". See explanation later"
w24.exract.promoter.unique <- w24.exract.promoter[!duplicated(w24.exract.promoter, by = "gene"),]

# No need to reverse the diff value. 
 

# Save file
fwrite(w24.exract.promoter.unique, "data/Methyl_tramp_tope/dmr_PEITC-24wks_vs_TRAMP-24wks_promoter_unique_fdr.csv", sep = "\t")


### 3, extract comparison TOCO vs TRAMP
w24.exract <- w24[, c("gene", "chr", "start", "end", "CpG", "tramp_24.mu", "toco_24.mu", "tramp_24..toco_24.diff", "tramp_24..toco_24.pval", "tramp_24..toco_24.fdr", "feature", "distance")]
w24.exract <- na.omit(w24.exract)

# Keep only records on "Promoter ..."
w24.exract.promoter <- w24.exract[substring(w24.exract$feature, 1, 8) == "Promoter", ]

# keep only one record for each gene, based on records on "promoter". See explanation later"
w24.exract.promoter.unique <- w24.exract.promoter[!duplicated(w24.exract.promoter, by = "gene"),]

# No need to reverse the diff value. 


# Save file
fwrite(w24.exract.promoter.unique, "data/Methyl_tramp_tope/dmr_TOCO-24wks_vs_TRAMP-24wks_promoter_unique_fdr.csv", sep = "\t")


# 7-12
library("data.table")
w24 <- fread("data/Methyl_tramp_tope/combined.TD.c5r5.csv")
colnames(w24)
# [1] "chr"    "start"  "end"    "CpG"    "TD01-N" "TD01-X" "TD02-N" "TD02-X" "TD04-N" "TD04-X"
# [11] "TD05-N" "TD05-X" "TD07-N" "TD07-X" "TD09-N" "TD09-X" "TD10-N" "TD10-X" "TD11-N" "TD11-X"
# [21] "TD13-N" "TD13-X" "TD14-N" "TD14-X" "TD16-N" "TD16-X" "TD17-N" "TD17-X" "TD19-N" "TD19-X"
# [31] "TD20-N" "TD20-X" "TD21-N" "TD21-X" "TD23-N" "TD23-X"

w24.3 <- w24[, c("chr", "start", "end", "CpG")]
w24.3$TD01 <- w24$`TD01-X`/ w24$`TD01-N`
w24.3$TD02 <- w24$`TD02-X`/ w24$`TD02-N`
w24.3$TD04 <- w24$`TD04-X`/ w24$`TD04-N`
w24.3$TD05 <- w24$`TD05-X`/ w24$`TD05-N`
w24.3$TD19 <- w24$`TD19-X`/ w24$`TD19-N`
w24.3$TD20 <- w24$`TD20-X`/ w24$`TD20-N`

w24.3 <- na.omit(w24.3)
fwrite(w24.3, "data/Methyl_tramp_tope/combined.TD.c5r5.fraction.24wks.csv", sep = "\t")

w24a <- fread("data/Methyl_tramp_tope/combined.TD.c5r5.fraction.24wks_anno.csv")
colnames(w24a)
# [1] "chr"      "start"    "end"      "CpG"      "TD01"     "TD02"     "TD04"     "TD05"    
# [9] "TD19"     "TD20"     "feature"  "distance" "gene"

w24a.keep <- subset(w24a, gene %in% c("Hmox1", "Fos", "Jun", "Il1b"))
w24a.keep <- w24a.keep[w24a.keep$feature != "Distal Intergenic"]

w24ma <- w24a.keep[, c("TD04", "TD05", "TD01", "TD02", "TD19", "TD20")]
colnames(w24ma) <- c("PEITC_1", "PEITC_2", "TRAMP_1", "TRAMP_2", "Wildtype_1", "Wildtype_2")
rownames(w24ma) <- paste0(w24a.keep$gene, "_", c(1:10, 1:10, 1:2))

library("pheatmap")
pheatmap(w24ma,
         cluster_rows=F,# show_rownames=F,
         cluster_cols=F,
         #color = heat.colors(1024), #
        # border_color = NA,
         scale = "row",
         legend = T
         )
