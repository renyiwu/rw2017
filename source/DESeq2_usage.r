# Project: BGI db-db RNA-seq project
# Author: David Cheng & Davit Sargsyan
# Created: 7/5/2017       Last Update: 7/9/2017
#**************************************************
# # Install DESeq from Bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

# Load packages
require(DESeq2)

# Set working directory to where count files stored----
setwd("C:/Users/davichen/Desktop/DN_Data/db_db mouse data/de_analysis")
or
setwd("C:/git_local/DESeq2/de_analysis")  #(for Acer PC)
getwd() # to confirm wd

# Import featureCounts output count, tab-delimited text file into R (use Ron's script)
# d=read.table("filename",sep="\t",header=TRUE,stringsAsFactors=FALSE)

# // Meaning, read data from the filename (in quotes), using tab ("\t") as separator, into an object named "d".  
# The file has a first line containing the column headers (header=TRUE) and we want to read in any text (strings) 
# without turning them into statistical factors

# Load counts data from a tab separated values file----
dt1 <- read.table("combraw.count",
                  sep = "\t",
                  header = TRUE,
                  stringsAsFactors = FALSE)

# Rename the samples
colnames(dt1)[7:14] <- paste("DR", 
                             1:8, 
                             sep = "")

# Part I: Run a single model----
# Specify treatment groups
mat <- data.frame(condition = rep(c("16wC",
                                    "16wDB",
                                    "21wC",
                                    "21wDB"),
                                  each = 2),
                  batch = rep(1, 8))

# Prepare the data set
dds <- DESeq2::DESeqDataSetFromMatrix(countData = dt1[, 7:14],
                                      colData = mat,
                                      design = ~ condition)

# Differential expression analysis based on the Negative Binomial
dds <- DESeq2::DESeq(dds)
DESeq2::resultsNames(dds)

# Compare 2 groups at a time----
## a. 16-week DB vs. Control
out.16w <- DESeq2::results(dds, 
                           contrast = c("condition",
                                        "16wDB",
                                        "16wC"))
out.16w <- data.frame(dt1[, 1:6],
                      out.16w)
write.csv(out.16w, 
          file = "out.16w_all.csv")

## b. 21-week DB vs. Control
out.21w <- DESeq2::results(dds, 
                           contrast = c("condition",
                                        "21wDB",
                                        "21wC"))
out.21w <- data.frame(dt1[, 1:6], 
                      out.21w)
write.csv(out.21w, 
          file = "out.21w_all.csv")

# Part II: Alternatively, run the analysis 2 treatment groups at a time (Alternate)----
## A. Subset 16-week data
dt.16w <- dt1[, c(1:10)]

mat.16w <- data.frame(condition = rep(c("16wC",
                                        "16wDB"),
                                      each = 2),
                      batch = rep(1, 4))

# Prepare the data set
dds.16w <- DESeq2::DESeqDataSetFromMatrix(countData = dt.16w[, 7:10],
                                          colData = mat.16w,
                                          design = ~ condition)

# Differential expression analysis based on the Negative Binomial
dds.16w <- DESeq2::DESeq(dds.16w)

# Plot dispersion estimates
plotDispEsts(dds.16w)

# Pause, plot MA & PCA on 16wk raw counts here. Check if variability among replicates greater than variability among conditions----

plotMA(dds.16w, ylim=c(-2,2))
rld <- rlog(dds.16w, blind = F)
plotPCA(rld)

# Continue with DESeq2
DESeq2::resultsNames(dds.16w)

# Compare 16-week DB vs. Control
out.16w.2 <- DESeq2::results(dds.16w, 
                             contrast = c("condition",
                                          "16wDB",
                                          "16wC"))
# Plot Frequencies of p-values for 16 weeks
hist(out.16w.2$pvalue,
     col = "grey", border = "white", 
     xlab = "", ylab = "",
     main = "16 week Frequencies of p-values")

#optional analysis, do before running data.frame
summary(out.16w.2)

out.16w.2 <- data.frame(dt1[, 1:6],
                        out.16w.2)
write.csv(out.16w.2, 
          file = "out.16w.csv")

## B. Subset 21-week data
dt.21w <- dt1[, c(1:6, 11:14)]

mat.21w <- data.frame(condition = rep(c("21wC",
                                        "21wDB"),
                                      each = 2),
                      batch = rep(1, 4))

# Prepare the data set
dds.21w <- DESeq2::DESeqDataSetFromMatrix(countData = dt.21w[, 7:10],
                                          colData = mat.21w,
                                          design = ~ condition)

# Differential expression analysis based on the Negative Binomial
dds.21w <- DESeq2::DESeq(dds.21w)

# Plot dispersion estimates
plotDispEsts(dds.21w)

# Pause, plot MA & PCA on 21wk raw/log transformed counts here. Check if variability among replicates greater than variability among conditions----

rld <- rlog(dds.21w, blind = F)
plotPCA(rld)
plotMA(dds.21w, ylim=c(-2,2))

# Optional, shrink log2 fold changes
dds.21wLFC <- lfcShrink(dds.21w, coef = 2, out.21w.2=out.21w.2)
dds.21wLFC

# Continue with DESeq2
DESeq2::resultsNames(dds.21w)

# Compare 21-week DB vs. Control
out.21w.2 <- DESeq2::results(dds.21w, 
                             contrast = c("condition",
                                          "21wDB",
                                          "21wC"))

# Plot Frequencies of p-values for 21 weeks
hist(out.21w.2$pvalue,
     col = "grey", border = "white", 
     xlab = "", ylab = "",
     main = "21 week Frequencies of p-values")

# Optional analysis, do before running data.frame
summary(out.21w.2)


out.21w.2 <- data.frame(dt1[, 1:6],
                        out.21w.2)
write.csv(out.21w.2, 
          file = "out.21w2.csv")