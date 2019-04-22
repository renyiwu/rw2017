# list of useful packages
# R. Wu
# Feb 2019

# CRAN packages
install.packages("officer")
install.packages("magrittr") #the package for function "%>%"
install.packages("rvg") # for function "ph_with_vg"
install.packages("ggplot2")
install.packages("ggfortify")





# II. Bioconductor packages
install.packages("BiocManager")

# 1. DSS
# sudo apt-get install libcurl4-gnutls-dev
install.packages("RCurl")
# sudo apt-get install libxml2-dev
install.packages("XML")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomeInfoDb")
BiocManager::install("DSS")


# 2. GenomicFeatures
# sudo apt-get install libssl-dev
BiocManager::install("GenomicFeatures")

# 3. ChIPseeker
# sudo apt-get install libudunits2-dev
install.packages("units")
BiocManager::install("ChIPseeker")

# 4. DESeq
BiocManager::install("DESeq")

# and more
BiocManager::install("rtracklayer")
BiocManager::install("biomaRt")






BiocManager::install("")
BiocManager::install("")
BiocManager::install("")
BiocManager::install("")
BiocManager::install("")
BiocManager::install("")
BiocManager::install("")
BiocManager::install("")








# for R version below 3.5, use below for installation of Bioconductor packages
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")
