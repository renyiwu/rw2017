# Install packages in batch from CRAN and bioconductor
# R WU
# April 2019


# 1. CRAN packages----
cran_packages <- c("data.table",
                   "tidyverse",
                   "stringr",
                   "tibble",
                   "units",
                   "VennDiagram",
                   "zip",
                   "DT",
                   "farver",
                   "fgsea",
                   "ggdendro",
                   # "ggplot2",
                   "gridExtra",
                   "knitr",
                   "MASS",
                   "packrat",
                   "pheatmap",
                   "plotly",
                   "RColorBrewer",
                   "readxl",
                   # "shiny",
                   # "shinydashboard",
                   # "shinyFiles",
                   # "shinythemes",
                   # "shinyWidgets",
                   "BiocManager")

if (length(setdiff(cran_packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(cran_packages, rownames(installed.packages())), update = T, ask = F) }
#
# # 2. Bioconductor packages----
bio_packages <- c("DESeq2",
                  "DEGseq",
                  "GOSemSim",
                  "ChIPseeker",
                  # "TxDb.Mmusculus.UCSC.mm10.knownGene",
                  "DSS",
                  "units",
                  "fgsea",
                  # "org.Mm.eg.db",
                  # "dada2"
                  "BiocInstaller")

if (length(setdiff(bio_packages, rownames(installed.packages()))) > 0) {
  BiocManager::install(setdiff(bio_packages, rownames(installed.packages())), update = T, ask = F)}

