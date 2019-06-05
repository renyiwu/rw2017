# Change the list accordingly to fit your need.
# R Wu
# June 2019


# 1. CRAN packages----
cran_packages <- c("shinydashboard",
                   "data.table",
                   "dplyr",
                   "DT",
                   "farver",
                   "fgsea",
                   "ggdendro",
                   "ggplot2",
                   "gridExtra",
                   "knitr",
                   "MASS",
                   "packrat",
                   "pheatmap",
                   "plotly",
                   "RColorBrewer",
                   "readxl",
                   "shiny",
                   "shinyFiles",
                   "shinythemes",
                   "shinyWidgets",
                   "stringr",
                   "tibble",
                   "units",
                   "VennDiagram",
                   "zip",
                   "tidyverse",
                   "shinyBS",
                   # "NewPackage", # uncomment this line and change the name in duoble quotes to add more package. Add more lines if desire
                   "BiocManager")

if (length(setdiff(cran_packages, rownames(installed.packages()))) > 0) {install.packages(setdiff(cran_packages, rownames(installed.packages()))) }
#
# # 2. Bioconductor packages----
bio_packages <- c("DESeq2",
                  #     "BiocInstaller", # not work with R version 3.6.0+
                  "DEGseq",
                  "GOSemSim",
                  "ChIPseeker",
                  "TxDb.Mmusculus.UCSC.mm10.knownGene",
                  "DSS",
                  "farver",
                  "units",
                  "vsn",
                  "fgsea",
                  "org.Mm.eg.db",
                  # "NewPackage", # uncomment this line and change the name in duoble quotes to add more package. Add more lines if desire
                  "dada2")

if (length(setdiff(bio_packages, rownames(installed.packages()))) > 0) {BiocManager::install(setdiff(bio_packages, rownames(installed.packages())), update = F)}

