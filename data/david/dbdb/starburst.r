## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
# install also rvest httr
library("TCGAbiolinks")
biocLite("TCGAbiolinksGUI", dependencies = TRUE)
library(TCGAbiolinksGUI)
TCGAbiolinksGUI()

TCGAvisualize_starburst()
