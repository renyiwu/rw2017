# This app is used for analyzing RNA-seq and DNA methyl-seq data by Konglab at Rutgers University 
# Authors: R Wu and D Sargsyan




try(source("VERSION.R"), silent = TRUE)
LOCAL_VERSION <- UPDATE_VERSION
try(source(UPDATE_URL), silent = TRUE)
# try(source(paste0("https://raw.githubusercontent.com/", UPDATE_REPO, "/master/VERSION.R")), silent = TRUE)
REMOTE_VERSION <- UPDATE_VERSION



### changelog 

# Feb 6, 2021
# add descriptions for items in DNA methyl-seq analysis panel.

## 2-5-2021
# added support for replicates in DNA methyl-seq analysis
# R Wu


# May 6, 2019 R Wu.  added scripts for first time installation of needed packages.
# v1.1 4-17-2019
# 
 
 
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
                   "stats",
                   "units",
                   "VennDiagram",
                   "zip",
                   "tidyverse",
                   "shinyBS",
                   # "NewPackage", # uncomment this line and change the name in duoble quotes to add more package. Add more lines if desire
                   "BiocManager")

if (length(setdiff(cran_packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(cran_packages, rownames(installed.packages()))) 
  }
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
                  "fgsea",
                  "org.Mm.eg.db",
                  # "NewPackage", # uncomment this line and change the name in duoble quotes to add more package. Add more lines if desire
                  # "dada2",
                  "units")

if (length(setdiff(bio_packages, rownames(installed.packages()))) > 0) {
  BiocManager::install(setdiff(bio_packages, rownames(installed.packages())),
                       update = F)
  }

 

options(stringsAsFactors = FALSE)
options(shiny.maxRequestSize = 60*1024^2) 
# options(repos = BiocInstaller::biocinstallRepos())

library(shinydashboard)
library(DT)
library(shiny)
library(shinyBS)
library(shinythemes)
library(shinyFiles)
library(shinyWidgets)
library(data.table)
library(ggdendro)
library(ggplot2)
library(gridExtra)
library(knitr)
library(MASS)
library(pheatmap)
library(plotly)
library(RColorBrewer)
library(readxl)
library(stringr)
library(tibble)
library(VennDiagram)
library(zip)

# library(BiocInstaller) # no need in R 3.6
library(DESeq2)
library(DEGseq)
library(BiocParallel)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(DSS)

# load dplyr last to avoid "select" being masked from biomaRt or plotly
library(tidyverse)

# Determine the OS
sysOS <- Sys.info()[['sysname']]

# Specify root folder----
if (sysOS == "Linux") {
  volumes <- c("Home" = "/home/administrator/",
               "Data Storage" = "/datastorage/")
} else if (sysOS == "Windows") {
  volumes <- c("Local Disk (C:)" = "C:/")
}

# source functions
source("source/plots.R")
source("source/data_prep.R")

ui <- dashboardPage(dashboardHeader(title = "NGS Pipeline"),
                    dashboardSidebar(sidebarMenu(menuItem(text = "Introduction", 
                                                          tabName = "introduction", 
                                                          icon = icon("dashboard")),
                                                 menuItem(text = "RNA-seq Analysis", 
                                                          tabName = "rna-seq_analysis", 
                                                          icon = icon("th")),
                                                 menuItem(text = "DNA MethylSeq Analysis",
                                                          tabName = "dna-seq_analysis",
                                                          icon = icon("th")),
                                                 menuItem(text = "DNA vs RNA",
                                                          tabName = "dna_vs_rna",
                                                          icon = icon("th")))),
                    dashboardBody(tabItems(
                      
                      ## ---------- Introduction ------------
                      tabItem(tabName = "introduction",
                              h2("Introduction to the App",
                                 style = "padding-left: 1em"),
                              h4("Current local version is: ", LOCAL_VERSION,
                                 ". Latest remote version is: ", REMOTE_VERSION,
                                 ". In Rstudio, click Git then Pull to update.",
                                 style = "padding-left: 2em"), # by RW
                              # h3("Current version of this App was updated on ", update_date, ". In Rstudio, click Git then Pull for new version.",
                              #   style = "padding-left: 2em"), # by RW
                              p("This interactive web application (NGS pipeline) is developed in R with Shiny to 
                                (1) Perform explorative analysis on RNA-seq and DNA methylation data
                                (2) Conduct RNA differential expression (DE) analysis with ",
                                a("DEGseq",
                                  href = "https://bioconductor.org/packages/release/bioc/html/DEGseq.html"),
                                " and ",
                                a("DESeq2 ",
                                  href = "http://bioconductor.org/packages/release/bioc/html/DESeq2.html"),
                                " based on the provided input data (raw count table and experimental design table). 
                                (3) Condcut DNA differential methylation ratio (DMR) analysis with ",
                                a("DSS",
                                  href = "http://bioconductor.org/packages/release/bioc/html/DSS.html"),
                                "based on annotated DNA methylation ratio data.",
                                style = "padding-left: 2em"),
                              
                              h3("1. RNA-seq Analysis", 
                                 style = "padding-left: 1em"),
                              
                              h4("1.1 Read In Data",
                                 style = "padding-left: 2em"), 
                              p("This tab read in the data for the explorative analysis and DE analysis tabs. First read in RNA-seq count table, a tab or comma seperated file with which the first row will be deleted. 
                                Then read in the coresponding design table, a template design_table.csv file can be ",
                                a("downloaded here.", 
                                  href = "https://drive.google.com/uc?export=download&id=1iWeDU8JK5mZxpp6fd9ipf9JXuKKztv7Y"),
                                "To tailor the design table based on your count table, do not change the first two column names \'Sample Name\' and \'Sample Label\'.
                                In the \'Sample Name\' column, enter the sample names from your count table. In the \'Sample Label\' column,
                                 give unique sample label for each sample (no duplicates in this column). 
                                \'Sample Name\' is used to mapping design table with count table, \'Sample Label\' is used as choices for user to select samples and to show in the plots. 
                                User can think it as renaming each Sample Name with Sample Label.
                                After read in the design table, the app will match design table with count table using \'Sample Name\', 
                                a warning message box will show and tell user which row in design table is not found in the count table, if there is any."
                                ,style = "padding-left: 3em"),
                              
                              h4("1.2 Exploratory Analysis",
                                 style = "padding-left: 2em"), 
                              p("This tab gives some explorative plots for the samples user selected: a histogram of the log-transformed count to look at the distribution of the count 
                                and check number genes with zero count; a boxplot of the log-transformed count to compare the spread of each sample;
                                a heatmap of hierarchy clustering on samples and a PCA plot to show the sample-to-sample distance and help spot potential clustering.",
                                style = "padding-left: 3em"),
                              
                              h4("1.3 DE Analysis: DEGseq",
                                 style = "padding-left: 2em"), 
                              p("This tab is developed for different gene expression analysis using", 
                                a("DEGseq",
                                  href = "https://bioconductor.org/packages/release/bioc/html/DEGseq.html"),
                                  ". This package only takes design without replications. Using DEGexp function, 
                                it performs a pairwise comparison to identify differentially expressed genes and return a result 
                                table. Based on the result table, along with the adjusted p-value threshold and absolute log foldchage threshold provided by the user, the app also returns an MA plot.",
                                style = "padding-left: 3em"),
                              
                              h4("1.4 DE Analysis: DESeq2",
                                 style = "padding-left: 2em"), 
                              p("This tab is developed for different gene expression analysis using", 
                                a("DESeq2 ",
                                  href = "http://bioconductor.org/packages/release/bioc/html/DESeq2.html"),
                                ". It uses DESeqDataSetFromMatrix function to generate a DESeqDataSet object which is then used in function DESeq to estimate size factors,
                                dispersion, and then a Negative Binomila GLM, with which the Wald test is perfomred to calcuated the p value. (Note, experiments without 
                                replicates do not allow for estimation of the dispersion of counts around the expected value for each group, which is critical for differential expression analysis. 
                                Analysis without replicates was deprecated in v1.20 and is no longer supported since v1.22.) User then could extract a result table from a DESeq analysis which gives
                                base means across samples, log2 fold changes, standard errors, test statistics, p-values and adjusted p-values. Based on the extracted result table, along with the 
                                adjusted p-value threshold and absolute log foldchage threshold provided by the user, the app also returns an MA plot.",
                                style = "padding-left: 3em"),
                              
                              h4("1.5 Change In Gene Expression",
                                 style = "padding-left: 2em"), 
                              p("This tab takes in contrast 1 and contrast 2 result tables (.csv file) from Differential Gene Expression Analysis and returns two venn diagrams and a two-column heatmap 
                                of change in gene expression of contrast 1 and contrast 2. User has to select the gene name column, the expression difference (log foldchnage) column, the p-value column, and enter
                                the absolute log foldchage threshold and p-value threshold, seperately for each table. The App filters out the significant genes based on the provided threshold seperately for each 
                                table and then inner join the two tables based on gene name column to plot the venn diagrams and the heatmap. The app will give warning message if there is no gene left after filtering.",
                                style = "padding-left: 3em"),
                              h3("2. DNA Methylseq Analysis", 
                                 style = "padding-left: 1em"),
                              
                              h4("2.1 Read In Data",
                                 style = "padding-left: 2em"),
                              p("This tab takes in annotated DNA methylation table with N(reads) and X(read coverage), either tab or comma sepearted. 
                                Besides the reads and read coverage column, the input table should also include gene name, CpG, region column. The input table will be used in the following Exploratory Analysis and the DMR Analysis: DSS. ",
                                style = "padding-left: 3em"),
                              
                              h4("2.2 Exploratory Analysis",
                                 style = "padding-left: 2em"),
                              p("The tab gives 3 descripative plots of the DNA table. An annotation pie chart takes in the region column selected by the user, and plots the portion/percentage of each annotated region, 
                                Some regions are combined to give a clear plot. A histogram of CpG facet by region takes in the region and CpG column selected by the user, and plots the distribution of CpG numbers for each region. The last plot, methylation percentage
                                by region, takes in the region, the CpG column, the N and X columns seleted by the user, then plot the average methylation ratio by region, colored by treatment.",
                                style = "padding-left: 3em"),
                              
                              h4("2.3 DMR Analysis: DSS",
                                 style = "padding-left: 2em"),
                              p("This tab uses DSS package to do pairwise difference-in-methylation-ratio analysis (DMR). User please refer to the instruction in the Data Preparation section to construct BSseqData object, which is later used in 
                                function DMLtest to perform statistical test of differntially methylated loci (DML).",
                                style = "padding-left: 3em"),
                              
                              h4("2.4 Change In Methylation Ratio",
                                 style = "padding-left: 2em"),
                              p("This tab takes in contrast 1 and contrast 2 result tables (.csv file) from DMR Analysis: DSS. Like Change In Gene Expression tab, this tab returns two venn diagrams and a two-column heatmap of change in methylation ratio difference
                                of contrast 1 and contrast 2. User need to select the gene identifier column (chr_start, the first column in DSS result table), the DNA methylation ratio difference column, and the p-value column, and enter the absolute 
                                methylation ratio difference threshold and p-value threshold, seperately for each table. The App filters out the significant genes based on the provided threshold seperately for each 
                                table and then inner join the two tables based on gene identifier column to plot the venn diagrams and the heatmap. User could rename the region column and choose the region they want to plot. The inner joined 
                                table will change accordingly. Please select at least one region. The app will give a warning message if there is no gene left after filtering.",
                                style = "padding-left: 3em"),
                              
                              h3("3. RNA vs DNA Analysis", 
                                 style = "padding-left: 1em"),
                              h4("3.1 Read In Data",
                                 style = "padding-left: 2em"),
                              p("This tab takes in an RNA DE analysis result table and a DNA DMR analysis result table for a starburst plot. Note, there is no threshold filtering for the two table, all the genes in the two table will be plotted after inner joining. 
                                So user might want to read in the tables that are already filtered to only significant genes or other criteria of interest.",
                                style = "padding-left: 3em"),
                              
                              h4("3.2 Plot Data",
                                 style = "padding-left: 2em"),
                              p("This tab gives a starburst plot of inner joined RNA and DNA table by gene name, colored by region. User need to select the gene name column and expression difference column for RNA table, and 
                                 gene name column, methylation ratio column and region column for DNA table. To make the color legend clearer in the plot, after inner joinning two table, user could rename the region and filter the region of insterest to plot. Please select at least one region.
                                The app will give a warning message if there is no gene left to plot.",
                                style = "padding-left: 3em"),
                              ),
                      
 
                    ## ---------------- RNA seq analysis -----------------
                    tabItem(tabName = "rna-seq_analysis",
                            tabsetPanel(
                              ## ------- read in data -------------
                              tabPanel("Read In Data", 
                                       wellPanel(
                                         fluidRow(
                                           column(5,
                                                  textInput(inputId = "project_name",
                                                            label = "Enter project name below:",
                                                            value = "Enter project name")),
                                           column(4,
                                                  bsButton("q1", label = "", icon = icon("question"),style = "info", size = "extra-small"),
                                                  bsTooltip(id = "q1", title = "Will be used in the filenames of the downloaded files", placement = "right"))
                                           ),
                                         fluidRow(
                                           column(5,
                                                  fileInput(inputId = "rna_count",
                                                            label = "Select Count Table File")),
                                           column(4,
                                                  bsButton("q2", label = "", icon = icon("question"),style = "info", size = "extra-small"),
                                                  bsTooltip(id = "q2", title = "Refer to Introduction 1.1, take in tab or comma sepearted RNA count table, and delete the first row", placement = "right"))
                                           ),
                                         fluidRow(
                                           column(5,
                                                  fileInput(inputId = "rna_info",
                                                            label = "Select Design File",
                                                            accept = c(".csv"))),
                                           column(4,
                                                  bsButton("q3", label = "", icon = icon("question"),style = "info", size = "extra-small"),
                                                  bsTooltip(id = "q3", title = "Refer to Introduction 1.1 to downloaded and fill the design_table.csv. Warning if value of Sample Name is not found in count table", placement = "right"))
                                           )),
                                       
                                       
                                       wellPanel(
                                         tags$h4("View Count Table:"),
                                         DT::dataTableOutput(outputId = 'display_rna_count')),
                                       wellPanel(
                                         tags$h4("View Design Table:"),
                                         DT::dataTableOutput(outputId = 'display_rna_info'))),
                              
                              
                              ## ------- exploratory analysis ---------
                              tabPanel("Exploratory Analysis", 
                                       wellPanel(
                                         tags$h3("Data Preparation"),
                                         fluidRow(
                                           column(2,
                                                  selectInput(inputId = "gene_col_selected_expl",
                                                              label = "Select gene column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(5,
                                                  selectInput(inputId = "sample_label_selected_expl",
                                                              label = "Select samples (same order as you want them in the plots):",
                                                              choices = NULL,
                                                              multiple = TRUE,
                                                              selected = NULL)),
                                           column(3,
                                                  selectInput(inputId = "covariate_selected_expl",
                                                              label = "Group by:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL),
                                                  bsTooltip(id = "covariate_selected_expl", title = "Applied to the fill color in boxplot and PCA plot"))),
                                         fluidRow(
                                           column(2,
                                                  actionButton(inputId = "generate_result_expl",
                                                               label = "Generate Result")),
                                           column(2,
                                                  downloadButton(outputId = "download_expl",
                                                                 label = "Download Results"))
                                         )),
                                       wellPanel(
                                         tags$h3("Count Histogram"),
                                         plotlyOutput(outputId = "display_count_hist")),
                                       
                                       wellPanel(
                                         tags$h3("Between-sample Distribution: Boxplot"),
                                         fluidRow(
                                           column(6,
                                                  plotlyOutput(outputId = "display_boxplot")))),
                                       wellPanel(
                                         fluidRow(
                                           column(6,
                                                  wellPanel(
                                                    tags$h3("Heatmap Sampel-to-sample Distance"),
                                                    fluidRow(
                                                      plotOutput(outputId = "display_heatmap_expl")
                                                    ))),
                                           column(6,
                                                  wellPanel(
                                                    tags$h3("PCA Plot"),
                                                    fluidRow(plotlyOutput(outputId = "display_pca_expl"))))))),
                              
                              
                              ## ------------- DEGseq --------------
                              tabPanel("DE Analysis: DEGseq",
                                       # Data preparation
                                       wellPanel(tags$h3("Data Preparation"),
                                                 fluidRow(
                                                   column(2,
                                                          selectInput(inputId = "gene_col_selected_degseq",
                                                                      label = "Select gene column",
                                                                      choices = NULL,
                                                                      multiple = FALSE,
                                                                      selected = NULL)),
                                                   column(5,
                                                          selectInput(inputId = "sample_label_selected_degseq",
                                                                      label = "Select 2 Samples in order: trt1(denominator) trt2(numerator)",
                                                                      choices = NULL,
                                                                      multiple = TRUE,
                                                                      selected = NULL),
                                                          bsTooltip(id = "sample_label_selected_degseq", title = "Refer to Introduction 1.3")),
                                                   column(4,
                                                          textInput(inputId = "row_sum",
                                                                    label = "Remove if sum across 2 samples is <",
                                                                    value = "10"))),
                                                 fluidRow(
                                                   column(2,
                                                          actionButton(inputId = "trim",
                                                                       label = "Generate Count Table"))),
                                                 fluidRow(
                                                   column(6,
                                                          tags$h4("Trimmed Count Table Summary:"),
                                                          verbatimTextOutput("trimmed_left"),
                                                          verbatimTextOutput("display_trimmed_ct")))),
                                       
                                       ## DE analysis
                                       wellPanel(
                                         tags$h3("DE Analysis"),
                                         fluidRow(
                                           column(3,
                                                  numericInput(inputId = "degseq_q_value",
                                                               label = "Set q-value(Storey et al. 2003) <",
                                                               value = 0.01,
                                                               min = 0,
                                                               max = 1,
                                                               step = 0.01),
                                                  bsTooltip(id = "degseq_q_value", title = "Applied to MA plot and the Signature column of the result table, showing signaficant or not. All genes from trimmed table will be included in the result table")),
                                           column(3,
                                                  numericInput(inputId = "degseq_fold_change",
                                                              label = "Set Obs(log foldchange) >=",
                                                              value = 1,
                                                              min = 0),
                                                  bsTooltip(id = "degseq_fold_change", title = "Applied only to MA plot together with the q-value")),
                                           column(3,
                                                  actionButton(inputId = "run_DEGexp",
                                                               label = "RUN DEGexp")),
                                           column(3,
                                                  downloadButton(outputId = "download_degseq",
                                                                 label = "Download Results"),
                                                  bsTooltip(id = "download_degseq", title = "Including result table and MA plot"))),
                                         wellPanel(tags$h4("Result Table:"),
                                           DT::dataTableOutput(outputId = "result_table1_degseq")),
                                         
                                         wellPanel(tags$h4("MA Plot:"),
                                                   fluidRow(column(3,
                                                                   textInput(inputId = "ma1_title",
                                                                             label = "Enter MA plot title below:",
                                                                             value = ""),
                                                                   bsTooltip(id = "ma1_title", title = "Plot will update automatically"))),
                                                   fluidRow(column(6,
                                                                   plotOutput(outputId = "ma1_degseq")),
                                                            column(6,
                                                                   tableOutput(outputId = "degseq_sign_number")))))),
                              
                              
                              ## --------------- DEseq2 ----------------
                              tabPanel("DE Analysis: DEseq2",
                                       
                                       # Data Preparation
                                       wellPanel(
                                         tags$h3("Data Preparation"),
                                         fluidRow(
                                           column(2,
                                                  selectInput(inputId = "gene_col_selected_deseq2",
                                                              label = "Select gene column",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  selectInput(inputId = "sample_label_selected_deseq2",
                                                              label = "Select samples for DE analysis",
                                                              choices = NULL,
                                                              multiple = TRUE,
                                                              selected = NULL),
                                                  bsTooltip(id = "sample_label_selected_deseq2", title = "Refer to Introduction 1.4")),
                                           column(3,
                                                  textInput(inputId = "row_sum_deseq2",
                                                            label = "Remove if total across samples is <",
                                                            value = "10")),
                                           column(2,
                                                  selectInput(inputId = "covariate_selected_deseq2",
                                                              label = "Select one covariate",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))),
                                         fluidRow(
                                           column(2,
                                                  actionButton(inputId = "generate_dds_from_matrix",
                                                               label = "Generate DESeqDataSet From Matrix"))),
                                         fluidRow(
                                           column(6,
                                                  tags$h4("DESeqDataSet Summary:"),
                                                  verbatimTextOutput(outputId = "trim_left_number"),
                                                  verbatimTextOutput("display_formula")))),
                                       
                                       # DE analysis
                                       wellPanel(
                                         tags$h3("DE Analysis"),
                                         fluidRow(column(4,
                                                         numericInput(inputId = "deseq2_p_value",
                                                                     label = "FDR adjusted p-value <",
                                                                     value = 0.01,
                                                                     min = 0,
                                                                     max = 1,
                                                                     step = 0.01),
                                                         bsTooltip(id = "deseq2_p_value", title = "Applied only to MA plot. All genes from trimmed table will be included in the result table")),
                                                  column(4,
                                                         numericInput(inputId = "deseq2_fold_change",
                                                                     label = "Set Obs(log foldchange) >=",
                                                                     value = 1,
                                                                     min = 0),
                                                         bsTooltip(id = "deseq2_fold_change", title = "Applied only to MA plot")),
                                                  column(3,
                                                         actionButton(inputId = "run_deseq",
                                                                      label = "Run DESeq"))),
                                         fluidRow(column(4,
                                                         selectInput(inputId = "comp1",
                                                                     label = "Select numerator level for the fold change",
                                                                     choices = NULL,
                                                                     multiple = FALSE,
                                                                     selected = NULL)),
                                                  column(4,
                                                         selectInput(inputId = "comp2",
                                                                     label = "Select denominator leve for the fold changel",
                                                                     choices = NULL,
                                                                     multiple = FALSE,
                                                                     selected = NULL)),
                                                  column(2,
                                                         actionButton(inputId = "extract_results_deseq2",
                                                                      label = "Extract Results")),
                                                  column(2,
                                                         downloadButton(outputId = "download_deseq2",
                                                                        label = "Download Results"),
                                                         bsTooltip(id = "download_deseq2", title = "Including result table, MA plot and count plot if generated"))),
                                         fluidRow(verbatimTextOutput(outputId = "deseq2_results_dir")),
                                         
                                         wellPanel(tags$h4("Result Table:"),
                                                   DT::dataTableOutput(outputId = "display_dtf_res_contrast")),
                                         wellPanel(tags$h4("MA Plot:"),
                                                   fluidRow(column(3,
                                                                   textInput(inputId = "enter_deseq2_ma_title",
                                                                             label = "Enter MA plot title below:",
                                                                             value = ""),
                                                                   bsTooltip(id = "enter_deseq2_ma_title", title = "Plot will update automatically"))),
                                                   fluidRow(column(6,
                                                                   plotOutput(outputId = "display_deseq2_ma")),
                                                            column(6,
                                                                   tableOutput(outputId = "deseq2_sign_number")))),
                                         wellPanel(tags$h4("Plot Counts"),
                                                   fluidRow(column(3,
                                                                   textInput(inputId = "enter_gene_name",
                                                                             label = "Enter one gene name:")),
                                                            column(3,
                                                                   selectInput(inputId = "intgroup",
                                                                               label = "Select one covariate :",
                                                                               choices = NULL,
                                                                               multiple = TRUE,
                                                                               selected = NULL)),
                                                            column(2,
                                                                   actionButton(inputId = "generate_count_plot",
                                                                                label = "Generate Count Plot"))),
                                                   fluidRow(column(6,
                                                                   plotOutput(outputId = "deseq2_count_plot")))))),
                              ## ----------- change in gene expression analysis----------------
                              tabPanel("Change In Gene Expression",
                                       wellPanel(
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "rna_contrast_1",
                                                            label = "Select result table from DE analysis for contrast 1",
                                                            accept = c(".csv")))),
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "rna_contrast_2",
                                                            label = "Select result table from DE analysis for contrast 2",
                                                            accept = c(".csv"))))),
                                       wellPanel(
                                         tags$h4("View Contrast 1 Table:"),
                                         DT::dataTableOutput(outputId = 'display_rna_contrast_1')),
                                       wellPanel(
                                         tags$h4("View Contrast 2 Table:"),
                                         DT::dataTableOutput(outputId = 'display_rna_contrast_2')),
                                       
                                       wellPanel(
                                         tags$h3("Data Preparation"),
                                         
                                         hr(),
                                         tags$h4("From contrast 1:"),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "gene_col_selected_rna_contrast_1",
                                                              label = "Select gene column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))
                                           ),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "logfold_change_selected_contrast_1",
                                                              label = "Select Logfold change column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  numericInput(inputId = "logfold_change_thresh_change_in_gene_expr_contrast1",
                                                              label = "Set Obs(log foldchange) >=",
                                                              min = 0,
                                                              value = 1))
                                         ),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "p_value_selected_rna_contrast_1",
                                                              label = "Select P value column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  numericInput(inputId = "p_thresh_change_in_gene_expr_contrast1",
                                                              label = "Set P value <",
                                                              min = 0, 
                                                              max = 1,
                                                              value = 0.05,
                                                              step = 0.01))
                                         ),
                                         hr(),
                                         tags$h4("From contrast 2:"),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "gene_col_selected_rna_contrast_2",
                                                              label = "Select gene column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))
                                           
                                        ),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "logfold_change_selected_contrast_2",
                                                              label = "Select Logfold change column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  numericInput(inputId = "logfold_change_thresh_change_in_gene_expr_contrast2",
                                                              label = "Set Obs(log foldchange) >=",
                                                              min = 0,
                                                              value = 1))
                                         ),
                                        fluidRow(
                                          column(4,
                                                  selectInput(inputId = "p_value_selected_rna_contrast_2",
                                                              label = "Select P value column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  numericInput(inputId = "p_thresh_change_in_gene_expr_contrast2",
                                                              label = "Set P value <",
                                                              min = 0, 
                                                              max = 1,
                                                              value = 0.05,
                                                              step = 0.01))
                                        ),
                                        hr(),
                                         fluidRow(
                                           column(3,
                                                  actionButton(inputId = "inner_join_2contrasts_rna",
                                                               label = "Inner Join Contrast 1 and Contrast 2 table"))),
                                         fluidRow(
                                           column(7,
                                                  tags$h4("Inner Joined Table Summary:"),
                                                  verbatimTextOutput("display_summary_inner_joined_2contrasts_rna")))),
                                       wellPanel(
                                         tags$h3("Change-in-gene-expression Heatmap"),
                                         fluidRow(
                                           column(3,
                                                  textInput(inputId = "change_in_gene_expr_plot_title",
                                                            label = "Enter plot title",
                                                            value = "Change in Gene Expression"),
                                                  bsTooltip(id = "change_in_gene_expr_plot_title", title = "Plot will update automatically")),
                                           
                                           column(2,
                                                  textInput(inputId = "change_in_gene_expr_plot_x_text_1",
                                                            label = "Enter x axis label 1",
                                                            value = "trt_2 vs trt_1"),
                                                  bsTooltip(id = "change_in_gene_expr_plot_x_text_1", title = "Plot will update automatically")),
                                           column(2,
                                                  textInput(inputId = "change_in_gene_expr_plot_x_text_2",
                                                            label = "Enter x axis label 2",
                                                            value = "trt_3 vs trt_2"),
                                                  bsTooltip(id = "change_in_gene_expr_plot_x_text_2", title = "Plot will update automatically")),
                                           column(2,
                                                  actionButton(inputId = "generate_change_in_gene_expr_plot",
                                                               label = "Generate Plots"),
                                                  bsTooltip(id = "generate_change_in_gene_expr_plot", title = "From filtered and inner-joined table created above")),
                                           column(2,
                                                  downloadButton(outputId = "download_change_in_gene_expr_plot",
                                                                 label = "Download Results"),
                                                  bsTooltip(id = "download_change_in_gene_expr_plot", title = "Including inner-joined table, venn diagrams and heatmap", trigger = "hover"))),
                                         hr(),
                                         fluidRow(
                                           column(5,
                                                  plotOutput(outputId = "display_venn_diagram1_rna")),
                                           column(5,
                                                  plotOutput(outputId = "display_venn_diagram2_rna"))),
                                         hr(),
                                         fluidRow(
                                           column(7,
                                                  plotlyOutput(outputId = "change_in_gene_expr_plot")))))
                              # tabItem "rna-seq_analysis" ends
                              )),
                    
                    ## ---------------- DNA methylseq analysis ----------------
                    tabItem(tabName = "dna-seq_analysis",
                            tabsetPanel(
                              ## ------------ read in data -----------
                              tabPanel("Read In Data",
                                       wellPanel(
                                         fluidRow(
                                           column(5,
                                                  textInput(inputId = "project_name_dna",
                                                            label = "Enter project name below:",
                                                            value = "Enter project name")),
                                           column(4,
                                                  bsButton("prj_name_dna_button", label = "", icon = icon("question"),style = "info", size = "extra-small"),
                                                  bsTooltip(id = "prj_name_dna_button", title = "Will be used as base filenames of the downloaded files", placement = "right"))
                                           ),
                                         fluidRow(
                                           column(5,
                                                  fileInput(inputId = "dna_table",
                                                            label = "Select Annotated DNA file")),
                                           column(4,
                                                  bsButton("dna_table_button", label = "", icon = icon("question"),style = "info", size = "extra-small"),
                                                  bsTooltip(id = "dna_table_button", title = "Refer to Introduction 2.1, take in tab or comma sepearted DNA methylation table with N(reads) and X(read coverage)", placement = "right"))
                                           )),
                                       wellPanel(
                                         tags$h4("View DNA Table:"),
                                         DT::dataTableOutput(outputId = 'display_dna_table')
                                       )),
                              ## ----- exploratory analysis -----------
                              tabPanel("Exploratory Analysis",
                                       wellPanel(tags$h3("Data Preparation"),
                                                 fluidRow(column(3,
                                                                 selectInput(inputId = "dna_sample_cols1",
                                                                             label = "Select sample columns for dtN",
                                                                             choices = NULL,
                                                                             multiple = TRUE,
                                                                             selected = NULL)),
                                                          column(3,
                                                                 selectInput(inputId = "dna_sample_cols2",
                                                                             label = "Select sample columns for dtX",
                                                                             choices = NULL,
                                                                             multiple = TRUE,
                                                                             selected = NULL),
                                                                 bsTooltip(id = "dna_sample_cols2", title = "Same order as dtN")),
                                                          
                                                          column(5,
                                                                 textInput(inputId = "sample_name",
                                                                           label = "Enter sample name in order, seperated by comma without tab",
                                                                           value = "name1,name2,name.."))),
                                                 fluidRow(column(3,
                                                                 selectInput(inputId = "dna_feature_col",
                                                                             label = "Select region column",
                                                                             choices = NULL,
                                                                             multiple = FALSE,
                                                                             selected = NULL)),
                                                          column(3,
                                                                 selectInput(inputId = "dna_cpg_col",
                                                                             label = "Select CpG column",
                                                                             choices = NULL,
                                                                             multiple = FALSE,
                                                                             selected = NULL))),
                                                 fluidRow(column(2,
                                                                 actionButton(inputId = "generate_expl_dna",
                                                                              label = "Generate Results"))),
                                                 verbatimTextOutput("test_dna")),
                                       
                                       wellPanel(tags$h3("Annotation Percentage(%) Pie Chart"),
                                                 plotOutput(outputId = "anno_by_reg")),
                                       wellPanel(tags$h3("CpG Histogram by Region"),
                                                 plotOutput(outputId = "cpg_hist")),
                                       wellPanel(tags$h3("Methyl% by region"),
                                                 plotOutput(outputId = "meth_by_region"))),
                              
                              ## ------------ DSS ---------------
                              tabPanel("DMR Analysis: DSS",
                                       wellPanel(
                                         tags$h3("Data Preparation"),
                                         tags$h4("Instruction"),
                                         p("To make BSseq Data for pairwise comparison, select the sample columns for each comparison, following the order 
                                           (1)N: Read coverage of the position from BS-seq data; 
                                           (2)X: Number of reads showing methylation of the position.",
                                           style = "padding-left: 0em"),
                                         p("For example, comparison1 (Con_N, Con_X), comparison2 (Exp_N, Exp_X). 
                                           Corresponding name can be Con and  Exp. For more information, please refer to  ",
                                           a("DSS. ", href = "https://bioconductor.org/packages/release/bioc/manuals/DSS/man/DSS.pdf"),
                                           style = "padding-left: 0em"),
                                         
                                         
                                         tags$h4("Construct BSseqData"),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "dss_gene",
                                                              label = "Select the gene name column (geneid, etc)",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  selectInput(inputId = "dss_chr",
                                                              label = "Select the chromosome column (chr, seqname, etc)",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "dss_start",
                                                              label = "Select the start position column (start, etc)",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  selectInput(inputId = "dss_region",
                                                              label = "Select the region (feature) column (annotation, etc)",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "dss_comp1",
                                                              label = "Select N and X column(s) for group 1 (for each sample, select N FIRST then X)",
                                                              choices = NULL,
                                                              multiple = TRUE,
                                                              selected = NULL)),
                                           column(4,
                                                  textInput(inputId = "dss_comp1_name",
                                                            label = "Enter a name for group 1 (eg, treated)",
                                                            value = ""))
                                         ),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "dss_comp2",
                                                              label = "Select N and X column(s) for group 2 (for each sample, select N FIRST then X)",
                                                              choices = NULL,
                                                              multiple = TRUE,
                                                              selected = NULL)),
                                           column(4,
                                                  textInput(inputId = "dss_comp2_name",
                                                            label = "Enter a name for group 2 (eg, control)",
                                                            value = ""))
                                         ),
                                         tags$h4("Set Parameters for DML test. Do not change these settings unless you know what you are doing."),
                                         fluidRow(
                                           column(3,
                                                  checkboxInput(inputId = "dss_wo_rep",
                                                                label = "No replicates (equal.disp=True in DMLtest)",
                                                                value = FALSE, width = NULL)),
                                                  # selectInput(inputId = "dss_wo_rep",
                                                  #             label = "Sample without replicates? (equal.disp in DMLtest)",
                                                  #             choices = c("yes" = TRUE,"no" = FALSE),
                                                  #             multiple = FALSE,
                                                  #             selected = TRUE)),
                                           column(3,
                                                  checkboxInput(inputId = "dss_smoothing",
                                                                label = "Smoothing",
                                                                value = TRUE, width = NULL)),
                                                  # selectInput(inputId = "dss_smoothing",
                                                  #             label = "Smoothing.",
                                                  #             choices = c(TRUE,FALSE),
                                                  #             multiple = FALSE,
                                                  #             selected = TRUE)),
                                           column(3,
                                                  textInput(inputId = "dss_smoothing_span",
                                                            label = "smoothing.span",
                                                            value = 500))),
                                         tags$h4("Set threshold"),
                                         fluidRow(
                                           column(3,
                                                  selectInput(inputId = "dss_pvalue",
                                                              label = "P value cutoff (select 1 to keep all)",
                                                              choices = c(0.01, 0.05, 0.1, 0.25, 1),
                                                              selected = 0.01)),
                                           column(3,
                                                  selectInput(inputId = "dss_methyl_diff",
                                                              label = "Obs(diff methyl ratio) select 0 to keep all",
                                                              choices = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5),
                                                              selected = 0.1)),
                                           column(2,
                                                  actionButton(inputId = "generate_dss",
                                                               label = "Generate result")),
                                           column(2,
                                                  downloadButton(outputId = "download_dss",
                                                                 label = "Download result table")))),
                                         wellPanel(
                                         tags$h3("Result of DML test"),
                                         DT::dataTableOutput(outputId = "dss_dml_tb"))),
                              
                              ## ------------ change in methyl ratio analysis ---------------
                              tabPanel("Change In Methylation Ratio",
                                       wellPanel(
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "dna_contrast_1",
                                                            label = "Select result table from DSS analysis for contrast 1",
                                                            accept = c(".csv")))),
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "dna_contrast_2",
                                                            label = "Select result table from DSS analysis for contrast 2",
                                                            accept = c(".csv"))))),
                                       wellPanel(
                                         tags$h4("View Contrast 1 Table:"),
                                         DT::dataTableOutput(outputId = 'display_contrast_1')),
                                       wellPanel(
                                         tags$h4("View Contrast 2 Table:"),
                                         DT::dataTableOutput(outputId = 'display_contrast_2')),
                                       
                                       wellPanel(
                                         tags$h3("Data Preparation"),
                                         hr(),
                                         tags$h4("From contrast 1:"),
                                         
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "gene_col_selected_dna_contrast_1",
                                                              label = "Select gene identifier column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL),
                                                  bsTooltip(id = "gene_col_selected_dna_contrast_1",
                                                            title = "The unique column: chr_start")),
                                           column(4,
                                                  selectInput(inputId = "gene_name_col_selected_dna_contrast_1",
                                                              label = "Select gene name column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL),
                                                  bsTooltip(id = "gene_name_col_selected_dna_contrast_1",
                                                            title = "The gene name column")),
                                           column(4,
                                                  selectInput(inputId = "region_col_selected_dna_contrast_1",
                                                              label = "Select region column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))
                                           ),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "methyl_diff_dna_contrast_1",
                                                              label = "Select DNA methylation ratio difference column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  numericInput(inputId = "change_in_dmr_diff_thresh1",
                                                              label = "Set Obs(log foldchange) >=",
                                                              min = 0,
                                                              value = 0))
                                         ),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "pval_col_selected_dna_contrast_1",
                                                              label = "Select P value column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  numericInput(inputId = "change_in_dmr_p_thresh1",
                                                               label = "Set P value <",
                                                               value = 0.01,
                                                               min = 0,
                                                               max = 1))
                                           ),
                                         
                                         
                                         hr(),
                                         tags$h4("From contrast 2:"),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "gene_col_selected_dna_contrast_2",
                                                              label = "Select gene identifier column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL),
                                                  bsTooltip(id = "gene_col_selected_dna_contrast_2",
                                                            title = "The unique column: chr_start")),
                                           column(4,
                                                  selectInput(inputId = "region_col_selected_dna_contrast_2",
                                                              label = "Select region column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))
                                           ),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "methyl_diff_dna_contrast_2",
                                                              label = "Select DNA methylation ratio difference column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  numericInput(inputId = "change_in_dmr_diff_thresh2",
                                                              label = "Set Obs(log foldchange) >=",
                                                              min = 0,
                                                              value = 0))
                                         ),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "pval_col_selected_dna_contrast_2",
                                                              label = "Select P value column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  numericInput(inputId = "change_in_dmr_p_thresh2",
                                                               label = "Set P value <",
                                                               value = 0.01,
                                                               min = 0,
                                                               max = 1))
                                           ),
                                         
                                         hr(),
                                         fluidRow(
                                           column(3,
                                                  actionButton(inputId = "inner_join_2contrasts",
                                                               label = "Inner Join Contrast 1 and Contrast 2 table"))
                                         ),
                                         fluidRow(
                                           column(9,
                                                  tags$h4("Inner Joined Table Summary:"),
                                                  verbatimTextOutput("display_summary_inner_joined_2contrasts"))),
                                         hr(),
                                         fluidRow(
                                           column(4,
                                                  textInput(inputId = "region_name_contian_dna",
                                                            label = "(Optional) Rename region that contain:")),
                                           column(4,
                                                  textInput(inputId = "region_new_name_dna",
                                                            label = "As:")),
                                           column(2,
                                                  actionButton(inputId = "rename_region_dna",
                                                               label = "Apply"),
                                                  bsTooltip(id = "rename_region_dna", title = "Rename one at a time, case sensitive, hit Inner Join button to start over"))
                                         ),
                                         hr(),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "region_selcted_dna",
                                                              label = "Include region:",
                                                              choices = NULL,
                                                              multiple = TRUE,
                                                              selected = NULL)),
                                           column(2,
                                                  actionButton(inputId = "filter_region_dna",
                                                               label = "Apply"),
                                                  bsTooltip(id = "filter_region_dna", title = "Filter the inner joined table by region, hit Inner Join button to start over"))
                                         )),
                                       wellPanel(
                                         tags$h3("Change-in-Methyl-Ratio Heatmap"),
                                         fluidRow(
                                           column(3,
                                                  textInput(inputId = "change_in_methyl_plot_title",
                                                            label = "Enter plot title",
                                                            value = "Change in Methylation Ratio"),
                                                  bsTooltip(id = "change_in_methyl_plot_title", title = "Plot will update automatically")),
                                           
                                           column(2,
                                                  textInput(inputId = "change_in_methyl_plot_x_text_1",
                                                            label = "x axis text 1",
                                                            value = "trt_2 vs trt_1"),
                                                  bsTooltip(id = "change_in_methyl_plot_x_text_1", title = "Plot will update automatically")),
                                           column(2,
                                                  textInput(inputId = "change_in_methyl_plot_x_text_2",
                                                            label = "x axis text 2",
                                                            value = "trt_3 vs trt_2"),
                                                  bsTooltip(id = "change_in_methyl_plot_x_text_2", title = "Plot will update automatically")),
                                           column(2,
                                                  actionButton(inputId = "generate_change_in_methyl_plot",
                                                               label = "Generate Plots"),
                                                  bsTooltip(id = "generate_change_in_methyl_plot", title = "From filtered and inner-joined table created above")),
                                           column(2,
                                                  downloadButton(outputId = "download_change_in_methyl_plot",
                                                                 label = "Download Results"),
                                                  bsTooltip(id = "download_change_in_methyl_plot", title = "Including inner-joined table, venn diagrams and heatmap", trigger = "hover"))),
                                         hr(),
                                         fluidRow(
                                           column(5,
                                                  plotOutput(outputId = "display_venn_diagram1_dna")),
                                           column(5,
                                                  plotOutput(outputId = "display_venn_diagram2_dna"))),
                                         hr(),
                                         fluidRow(
                                           column(7,
                                                  plotlyOutput(outputId = "change_in_methyl_plot")))))
                              # tabItem DNA end
                              )),
                    
                    ## --------------- DNA vs RNA -----------------
                    tabItem(tabName = "dna_vs_rna",
                            tabsetPanel(
                              tabPanel("Read In data",
                                       wellPanel(
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "rna_sig",
                                                            label = "Select RNA File of Significant Genes",
                                                            accept = c(".csv")))),
                                         fluidRow(
                                           column(6,
                                                  fileInput(inputId = "dna_sig",
                                                            label = "Select DNA File of Significant Genes",
                                                            accept = c(".csv"))))),
                                       wellPanel(
                                         tags$h4("View RNA Table:"),
                                         DT::dataTableOutput(outputId = 'display_rna_sig')),
                                       wellPanel(
                                         tags$h4("View DNA Table:"),
                                         DT::dataTableOutput(outputId = 'display_dna_sig'))),
                              tabPanel("Plot Data",
                                       wellPanel(
                                         tags$h3("Data Preparation"),
                                         hr(),
                                         tags$h4("From RNA table:"),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "rna_gene_col_selected_rna_vs_dna",
                                                              label = "Select gene column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  selectInput(inputId = "logfoldchange_col_selected_rna_vs_dna",
                                                              label = "Select expression difference column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))),
                                         hr(),
                                         tags$h4("From DNA table:"),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "dna_gene_col_selected_rna_vs_dna",
                                                              label = "Select gene column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           column(4,
                                                  selectInput(inputId = "methyl_diff_col_selected_rna_vs_dna",
                                                              label = "Select methylation ratio difference column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL)),
                                           
                                           column(4,
                                                  selectInput(inputId = "region_col_selected_rna_vs_dna",
                                                              label = "Select region column:",
                                                              choices = NULL,
                                                              multiple = FALSE,
                                                              selected = NULL))),
                                         hr(),
                                         fluidRow(
                                           column(2,
                                                  actionButton(inputId = "inner_join_rna_vs_dna",
                                                               label = "Inner Join RNA and DNA table"))),
                                         fluidRow(
                                           column(8,
                                                  tags$h4("Inner Joined Table Summary:"),
                                                  verbatimTextOutput("display_summary_inner_joined_rna_dna"))),
                                         hr(),
                                         fluidRow(
                                           column(4,
                                                  textInput(inputId = "region_name_contian_rna_vs_dna",
                                                            label = "(Optional) Rename region that contain:")),
                                           column(4,
                                                  textInput(inputId = "region_new_name_rna_vs_dna",
                                                            label = "As:")),
                                           column(2,
                                                  actionButton(inputId = "rename_region_rna_vs_dna",
                                                               label = "Apply"),
                                                  bsTooltip(id = "rename_region_rna_vs_dna", title = "Rename one at a time, case sensitive, hit Inner Join button to start over"))),
                                         hr(),
                                         fluidRow(
                                           column(4,
                                                  selectInput(inputId = "region_selcted_dna_rna",
                                                              label = "Include region:",
                                                              choices = NULL,
                                                              multiple = TRUE,
                                                              selected = NULL)),
                                           column(4,
                                                  actionButton(inputId = "filter_region_dna_rna",
                                                               label = "Apply"),
                                                  bsTooltip(id = "filter_region_dna_rna", title = "Filter the inner joined table by region, hit Inner Join button to start over")
                                                  )
                                                  )),
                                       wellPanel(
                                         tags$h3("Starburst Plot"),
                                         fluidRow(
                                           column(3,
                                                  numericInput(inputId = "starburst_rna_hline",
                                                               label = "Enter RNA expression diff threshold",
                                                               value = 0.1),
                                                  bsTooltip(id = "starburst_rna_hline", title = "Hit generate button to update")),
                                           column(3,
                                                  numericInput(inputId = "starburst_dna_vline",
                                                               label = "Enter DNA methylation % diff threshold",
                                                               value = 10),
                                                  bsTooltip(id = "starburst_dna_vline", title = "Hit generate button to update")),
                                           column(3,
                                                  textInput(inputId = "starburst_title",
                                                            label = "Enter plot title",
                                                            value = "DNA vs RNA: starburst plot"),
                                                  bsTooltip(id = "starburst_title", title = "Plot will update automatrically")),
                                           column(2,
                                                  actionButton(inputId = "generate_starburst_plot_rna_vs_dna",
                                                               label = "Generate Plot")),
                                           column(2,
                                                  downloadButton(outputId = "download_starburst_plot_and_joined_table",
                                                                 label = "Download Results"),
                                                  bsTooltip(id = "download_starburst_plot_and_joined_table",
                                                            title = "Including inner-joined table and starburst plot"))),
                                         hr(),
                                         fluidRow(
                                           column(7,
                                                  plotlyOutput(outputId = "starburst_plot")))))
                            # tabItem DNA vs RNA end
                            ))
                    )))

###############################################################################

## ------------- sssssSERVERrrrrr -------------------
server <- function(input, output, session) {
  rv <- reactiveValues()
  shinyDirChoose(input = input, 
                 id = "directory", 
                 roots = volumes, 
                 session = session, 
                 restrictions = system.file(package = "base"))
  
  ## ------------  RNA seq: Read in data and update SelectInput -------------
  
  observeEvent(input$rna_count,{
    infile_rna_count <- input$rna_count
    tb <- fread(infile_rna_count$datapath, 
                skip = 1)
    rv$ct <- as.data.frame(tb)
    
    output$display_rna_count = DT::renderDataTable({
      DT::datatable(rv$ct,
                    options = list(scrollX = TRUE))})
  })
  
  observeEvent(input$rna_info,{
    infile_rna_info <- input$rna_info
    tb <- fread(infile_rna_info$datapath, 
                header = TRUE)
    tb <- as.data.frame(tb)
    rv$info <- tb
    
    output$display_rna_info = DT::renderDataTable({
      DT::datatable(rv$info,
                    options = list(scrollX = TRUE))})
    
    if (any(!(rv$info$`Sample Name` %in% colnames(rv$ct)))) {
      showNotification(paste("Design table row", paste(which(!(rv$info$`Sample Name` %in% colnames(rv$ct))), collapse = ","), "\'sample name is not found in count table"),
                      type = "error",
                      duration = NULL)
    } else {
      
    updateSelectInput(session, 
                      inputId = "gene_col_selected_expl", 
                      label = "Select gene column", 
                      choices = colnames(rv$ct), 
                      selected = NULL)
    updateSelectInput(session, 
                      inputId = "sample_label_selected_expl", 
                      label = "Select samples (same order as in histogram and boxplot)", 
                      choices = rv$info$`Sample Label`, 
                      selected = NULL)
    updateSelectInput(session, 
                      inputId = "covariate_selected_expl",
                      label = "Group by (select one):",
                      choices = colnames(rv$info)[c(-1,-2)], 
                      selected = NULL)
    
    updateSelectInput(session, 
                      inputId = "gene_col_selected_degseq", 
                      label = "Select gene column", 
                      choices = colnames(rv$ct), 
                      selected = NULL)
    updateSelectInput(session, 
                      inputId = "sample_label_selected_degseq", 
                      label = "Select 3 samples in order: trt1 trt2 trt3", 
                      choices = rv$info$`Sample Label`, 
                      selected = NULL)
    
    updateSelectInput(session, 
                      inputId = "gene_col_selected_deseq2", 
                      label = "Select gene column", 
                      choices = colnames(rv$ct), 
                      selected = NULL)
    updateSelectInput(session, 
                      inputId = "sample_label_selected_deseq2", 
                      label = "Select samples for DE analysis", 
                      choices = rv$info$`Sample Label`, 
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "covariate_selected_deseq2",
                      label = "Select one covariate",
                      choices = colnames(rv$info)[c(-1,-2)], 
                      selected = NULL)
    }
    
  })
  
  ## ---------------------- Exploratory Analysis -------------
  observeEvent(input$generate_result_expl,{
    
    withProgress(message = "Generating Result",
                 value = 0,
                 expr = {
                   # prepare selected count table
                   incProgress(0.2, detail = "Preparing Data")
                   rv$ct[is.na(rv$ct)] <- 0
                   dt <- try(
                     match.ct.info(rv$ct, 
                                   rv$info, 
                                   input$sample_label_selected_expl, 
                                   input$gene_col_selected_expl), 
                     silent = T)
                   if (class(dt) == "try-error") {
                     showNotification(dt[1],
                                      type = "error",
                                      duration = NULL)
                     
                   }else{
                     # selected design table
                     selected_info <- rv$info[match(input$sample_label_selected_expl, 
                                                    rv$info$`Sample Label`), 
                                              c("Sample Label", 
                                                input$covariate_selected_expl)]
                     colnames(selected_info) <- c("sample", "covariate")
                     # log transform
                     dt_log2 <- dt
                     dt_log2[, 2:ncol(dt_log2)] <- log2(dt_log2[, 2:ncol(dt_log2)] + 1)
                     dt_log2_long <- melt(dt_log2,
                                          id.vars = 1,
                                          measure.vars = 2:ncol(dt_log2),
                                          variable.name = "sample",
                                          value.name = "count")
                     # merged long 
                     dt1 <- merge(dt_log2_long, selected_info, by = "sample")
                     
                     # histogram
                     rv$count_hist <- ggplot(data = dt1,
                                             aes(x = count)) +
                       geom_histogram() +
                       xlab("") +
                       ylab("log2(count + 1)") + 
                       facet_wrap( ~ sample)
                     
                     p1 <- ggplotly(rv$count_hist)
                     output$display_count_hist <- renderPlotly({print(p1)})
                     
                     # boxplot
                     rv$boxplot <- ggplot(data = dt1,
                                          aes(x = sample,
                                              y = count,
                                              fill = covariate )) + 
                       geom_boxplot() +
                       xlab("") +
                       ylab("log2(count + 1)") + 
                       labs(fill = "legend")
                     
                     p2 <- ggplotly(rv$boxplot)
                     output$display_boxplot <- renderPlotly({print(p2)})
                     
                     # heatmap s2s
                     tmp <- dt_log2[, -1]
                     sample_dist <- as.matrix(dist(t(as.matrix(tmp))))

                     rv$expl_heatmap <- pheatmap(sample_dist,
                                                 col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255))

                     output$display_heatmap_expl <- renderPlot({print(rv$expl_heatmap)})

                     # pca
                     matrix_pca <- t(dt_log2[, -1])
                     matrix_pca <- matrix_pca[, apply(matrix_pca, 2, var) != 0]
                     res_pca <- prcomp(matrix_pca,
                                       center = TRUE,
                                       scale. = TRUE)
                     df_out <- as.data.frame(res_pca$x)

                     percentage <- round(res_pca$sdev^2 / sum(res_pca$sdev^2) * 100, 2)
                     percentage <- paste( colnames(df_out), "(", paste0( as.character(percentage), "%"), ")")

                     df_out$sample <- colnames(dt_log2)[-1]
                     df_out <- merge(df_out, selected_info, by = "sample")

                     rv$pca <- ggplot(data = df_out,
                                      aes(x = PC1,
                                          y = PC2)) +
                       geom_point(aes(fill = df_out[, ncol(df_out)]),
                                  shape = 21,
                                  size = 3,
                                  alpha = 0.5) +
                       xlab(percentage[1]) +
                       ylab(percentage[2]) +
                       geom_text(data = df_out,
                                 aes(x = PC1 + 20,
                                     y = PC2,
                                     label = df_out$sample),
                                 size = 2,
                                 hjust = 0.5) +
                       labs(fill = "legend")

                     p4 <- ggplotly(rv$pca)

                     output$display_pca_expl <- renderPlotly({print(p4)})
                     
                     # end of if else   
                   }
                   
                   
                   incProgress(0.8, detail = "Processing 100%")
                   Sys.sleep(1)
                 })
    
  })
  
  ## ---------------------- DEGseq  ----------------
  ## Data Preparation
  observeEvent(input$trim,{
    if (length(input$sample_label_selected_degseq) != 2) {
      showNotification("Please select two samples in order",
                       type = "error",
                       duration = NULL)
    }else{
      dt <- try(
        match.ct.info(rv$ct, 
                      rv$info, 
                      input$sample_label_selected_degseq, 
                      input$gene_col_selected_degseq), 
        silent = T)
      if (class(dt) == "try-error") {
        showNotification(dt[1],
                         type = "error",
                         duration = NULL)
        
      }else{
        # trim 
        rv$ct_degseq <- subset(dt,
                               rowSums(dt[, -1]) > as.numeric(input$row_sum))
        
        output$trimmed_left <- renderPrint({
          cat(
            paste(nrow(rv$ct_degseq),
                  "genes left, down from",
                  nrow(rv$ct),
                  "genes",
                  sep = " "))})
        output$display_trimmed_ct <- renderPrint({summary(rv$ct_degseq)})
      }                         
    }  
  })
  
  ## DE analysis
  observeEvent(input$run_DEGexp,{
    if (is.null(rv$ct_degseq)) {
      showNotification("Please create a count table in Data Preparation first",
                       type = "error",
                       duration = NULL)
    } else {
      withProgress(message = "Running DEGexp",
                 value = 0,
                 expr = {
                   # run DEGexp
                   incProgress(0.4, 
                               detail = "Processing trt2-trt1")
                   outputDir <-  file.path(tempdir())
                   rv$ct_degseq_matrix <- as.matrix(rv$ct_degseq)
                   DEGexp(geneExpMatrix1 = rv$ct_degseq_matrix,
                          geneCol1 = 1, 
                          expCol1 = 3, 
                          groupLabel1 = colnames(rv$ct_degseq_matrix)[3],
                          
                          geneExpMatrix2 = rv$ct_degseq_matrix,
                          geneCol2 = 1, 
                          expCol2 = 2,
                          groupLabel2 = colnames(rv$ct_degseq_matrix)[2],
                          
                          foldChange = as.numeric(input$degseq_fold_change),
                          
                          qValue = as.numeric(input$degseq_q_value),
                          thresholdKind = 5, 
                          rawCount = TRUE,
                          normalMethod = "none",
                          method = "MARS",
                          outputDir = outputDir)
                   
                   rv$table_trt2_trt1 <- fread(paste(outputDir,"\\output_score.txt", sep = ""))
                   
                   
                   incProgress(0.2, detail = "Displaying Result Tables")
                   
                   
                   rv$table_trt2_trt1_round <- rv$table_trt2_trt1
                   rv$table_trt2_trt1_round[, c(4:9)] <- round(rv$table_trt2_trt1_round[, c(4:9)], 2)
                   output$result_table1_degseq = DT::renderDataTable({
                     DT::datatable(rv$table_trt2_trt1_round,
                                   filter = "top",
                                   options = list(scrollX = TRUE,
                                                  columnDefs = list(list(searchable = FALSE, 
                                                                         targets = c(1:4,6:7)))))
                   })
                   
                   # MA plot
                   incProgress(0.2, detail = "Generating MA plot")
                   # add mu clr and pch
                   dt_ma <- as.tibble(rv$table_trt2_trt1) %>% 
                     mutate(mu = case_when(value1 != 0 & value2 != 0 ~ (log2(value1) + log2(value2))/2,
                                           value1 != 0 & value2 == 0 ~ (log2(value1) + log2(value2 + 1))/2,
                                           value1 == 0 & value2 != 0 ~ (log2(value1 + 1) + log2(value2))/2,
                                           TRUE ~ (log2(value1 + 1) + log2(value2 + 1))/2)) %>% 
                     mutate(clr = case_when( `q-value(Storey et al. 2003)` < input$degseq_q_value & 
                                               abs(`log2(Fold_change) normalized`) >= input$degseq_fold_change ~ "red",
                                             `q-value(Storey et al. 2003)` < input$degseq_q_value & 
                                               abs(`log2(Fold_change) normalized`) < input$degseq_fold_change ~ "purple",
                                             TRUE ~ "black")) %>% 
                     mutate(clr = factor(clr, levels = c("black", "purple", "red"))) %>% 
                     mutate(pch = case_when(`q-value(Storey et al. 2003)` < input$degseq_q_value & 
                                               abs(`log2(Fold_change) normalized`) >= input$degseq_fold_change ~ 4,
                                             `q-value(Storey et al. 2003)` < input$degseq_q_value & 
                                               abs(`log2(Fold_change) normalized`) < input$degseq_fold_change ~ 3,
                                             TRUE ~ 46)) %>% 
                     mutate(pch = factor(pch, levels = c(46, 3, 4))) %>% 
                     rename(q_value = `q-value(Storey et al. 2003)`, log_foldchange = `log2(Fold_change) normalized`) %>% 
                     select(mu, clr, pch, log_foldchange)
                     
                   
                   rv$ma1_degseq <- ma(dt.tb = dt_ma,
                                       q_value_thresh = input$degseq_q_value, 
                                       fold_change_thresh = input$degseq_fold_change,
                                       tt = input$ma1_title) 
                   
                   output$ma1_degseq <- renderPlot({print(rv$ma1_degseq)})
                   
                   result <- c("up-regulated DEG",
                                 "non-significant DEG",
                                 "down-regulated DEG")
                     number <- c(sum(rv$table_trt2_trt1$`q-value(Storey et al. 2003)` < as.numeric(input$degseq_q_value) & 
                                       rv$table_trt2_trt1$`log2(Fold_change) normalized` >= as.numeric(input$degseq_fold_change), 
                                     na.rm = TRUE),
                                 nrow(rv$table_trt2_trt1) - sum(rv$table_trt2_trt1$`q-value(Storey et al. 2003)` < as.numeric(input$degseq_q_value) & 
                                                                   abs(rv$table_trt2_trt1$`log2(Fold_change) normalized`) >= as.numeric(input$degseq_fold_change),
                                                                 na.rm = TRUE),
                                 sum(rv$table_trt2_trt1$`q-value(Storey et al. 2003)` < as.numeric(input$degseq_q_value) & 
                                       rv$table_trt2_trt1$`log2(Fold_change) normalized` <= -as.numeric(input$degseq_fold_change),
                                     na.rm = TRUE))
                     degseq_sign_number_tb <- data.frame(result, number)
                     output$degseq_sign_number <- renderTable(degseq_sign_number_tb)
                   
                   
                   
                   
                   incProgress(0.2, detail = "Processing 100%")
                   Sys.sleep(1)
                 })
    }
  })
  
  # renew MA plot title
  observeEvent(input$ma1_title,{
    if (!is.null(rv$ma1_degseq)) {
      rv$ma1_degseq <- rv$ma1_degseq +
        ggtitle(input$ma1_title)
      output$ma1_degseq <- renderPlot({print(rv$ma1_degseq)})
    }
  })
  
  
  ## ---------------------- DEseq2 -----------------------
  
      
  observeEvent(input$generate_dds_from_matrix,{  
    
      
      # tailor the count tb
      rv$sample_name_list_deseq2 <- rv$info$`Sample Name`[match(input$sample_label_selected_deseq2, 
                                                              rv$info$`Sample Label`)]
      rv$ct_deseq2 <- as.matrix(rv$ct[, rv$sample_name_list_deseq2])
      rv$ct_deseq2[is.na(rv$ct_deseq2)] <- 0
      rownames(rv$ct_deseq2) <- rv$ct[, input$gene_col_selected_deseq2]
      colnames(rv$ct_deseq2) <- input$sample_label_selected_deseq2
      
      # tailor the design tb
      if (length(input$covariate_selected_deseq2) != 1) {
          showNotification("Please select one covariate",
                         type = "error",
                         duration = NULL)
      }else{
       rv$info_deseq2 <- rv$info[match(input$sample_label_selected_deseq2,
                                    rv$info$`Sample Label`), 
                              c("Sample Label",
                                input$covariate_selected_deseq2)]
       colnames(rv$info_deseq2) <- c("sample", 
                                  input$covariate_selected_deseq2)
       rv$info_deseq2[-1] <- lapply(rv$info_deseq2[-1], factor) 
       
       # design formula
       rv$formula <- as.formula(paste("~",input$covariate_selected_deseq2))
       output$display_formula <- renderPrint({
       cat(paste0("The design formula: ", paste(as.character(rv$formula), collapse = " ")))})
      }
      
      
      
      
      rv$dds <- try(DESeqDataSetFromMatrix(countData = rv$ct_deseq2,
                                           colData = rv$info_deseq2,
                                           design = rv$formula),
                    silent = TRUE)
      if (class(rv$dds) == "try-error") {
        showNotification(rv$dds[1],
                         type = "error",
                         duration = NULL)
        output$trim_left_number <- renderPrint({cat(rv$dds)})
        
      }else{
        rv$dds_trimmed <- rv$dds[rowSums(counts(rv$dds)) >= as.numeric(input$row_sum_deseq2), ]
        output$trim_left_number <- renderPrint({
          cat(paste(nrow(counts(rv$dds_trimmed)), "genes left, down from", nrow(counts(rv$dds)), "genes"))})}
  })
  
  
  
  observeEvent(input$run_deseq,{
    
    if (!is.null(rv$dds_trimmed)) {
      snowparam <- SnowParam(workers = snowWorkers(), 
                             type = "SOCK")
      register(snowparam, 
               default = TRUE)
      
      withProgress(message = "Running DESeq", 
                   value = 0,
                   expr = {
                     incProgress(0.25, detail = "Processing 25%")
                     incProgress(0.25, detail = "Processing 50%")
                     
                     rv$dds_res <- try(
                       DESeq(rv$dds_trimmed, fitType = "local", parallel = FALSE),
                       silent = TRUE
                     )
                     
                     if (class(rv$dds_res) == "try-error") {
                       showNotification(rv$dds_res[1],
                                        type = "error",
                                        duration = NULL)
                     }else{
                       contrast <- levels(rv$info_deseq2[,2])
                       updateSelectInput(session,
                                         inputId = "comp1",
                                         label = "Select numerator level for the fold change",
                                         choices = contrast,
                                         selected = NULL)
                       updateSelectInput(session,
                                         inputId = "comp2",
                                         label = "Select denominator level for the fold change",
                                         choices = contrast,
                                         selected = NULL)
                     }
                     
                     incProgress(0.25, detail = "Processing 75%")
                     incProgress(0.25, detail = "Processing 100%")
                     Sys.sleep(1)
                   })
      
    }else{
      showNotification("Please generate dds object from matrix first",
                       type = "error",
                       duration = NULL)
    }  
  })
  
  observeEvent(input$extract_results_deseq2,{
    if (is.null(rv$dds_res)) {
      showNotification("Please Run DESeq first before you extract result",
                       type = "error",
                       duration = NULL)
    }else if (input$comp1 == input$comp2) {
      showNotification("Please select two different samples for comparison",
                       type = "error",
                       duration = NULL)
    }else{
      withProgress(message = "Extracting results with contrast", 
                   value = 0,
                   expr = {
                     # extract result
                     incProgress(0.2, detail = "Preparing result table")
                     
                     contrast <- c(input$covariate_selected_deseq2, input$comp1, input$comp2)
                     
                     rv$dds_res_contrast <- results(rv$dds_res, contrast = contrast)
                     rv$dtf_res_contrast <- data.frame(gene = rownames(rv$dds_res_contrast),
                                                       do.call("cbind", 
                                                               rv$dds_res_contrast@listData))
                     
                     rv$dtf_res_contrast_round <- rv$dtf_res_contrast
                     rv$dtf_res_contrast_round[, c(2:7)] <- round(rv$dtf_res_contrast_round[, c(2:7)], 2)

                     output$display_dtf_res_contrast <- DT::renderDataTable({
                       DT::datatable(rv$dtf_res_contrast_round,
                                     filter = "top",
                                     options = list(scrollX = TRUE,
                                                    columnDefs = list(list(searchable = FALSE,
                                                                           targets = c(1:2, 4:6)))))})
                     
                    
                     
                     # MA plot
                     incProgress(0.2, detail = "Preparing MA plot")
                     dt_ma <- rv$dtf_res_contrast %>%
                       mutate(mu = log(baseMean + 1)) %>%
                       mutate(log_foldchange = log2FoldChange) %>%
                       mutate(clr = case_when( padj < input$deseq2_p_value &
                                               abs(log_foldchange) >= input$deseq2_fold_change ~ "red",
                                               padj < input$deseq2_p_value &
                                               abs(log_foldchange) < input$deseq2_fold_change ~ "purple",
                                             TRUE ~ "black")) %>%
                       mutate(clr = factor(clr, levels = c("black", "purple", "red"))) %>%
                       mutate(pch = case_when(padj < input$deseq2_p_value &
                                               abs(log_foldchange) >= input$deseq2_fold_change ~ 4,
                                               padj < input$deseq2_p_value &
                                               abs(log_foldchange) < input$deseq2_fold_change ~ 3,
                                             TRUE ~ 46)) %>%
                       mutate(pch = factor(pch, levels = c(46, 3, 4))) %>%
                       select(mu, clr, pch, log_foldchange)

                     rv$deseq2_ma <- ma(dt_ma, input$deseq2_p_value, input$deseq2_fold_change, input$enter_deseq2_ma_title)

                     output$display_deseq2_ma <- renderPlot({
                       print(rv$deseq2_ma)
                     })

                     # sign gene number
                     result <- c("up-regulated DEG",
                                 "non-significant DEG",
                                 "down-regulated DEG")
                     number <- c(sum(rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) &
                                       rv$dtf_res_contrast$log2FoldChange >= as.numeric(input$deseq2_fold_change),
                                     na.rm = TRUE),
                                 nrow(rv$dtf_res_contrast) - sum(rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) &
                                                                   abs(rv$dtf_res_contrast$log2FoldChange) >= as.numeric(input$deseq2_fold_change),
                                                                 na.rm = TRUE),
                                 sum(rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) &
                                       rv$dtf_res_contrast$log2FoldChange <= -as.numeric(input$deseq2_fold_change),
                                     na.rm = TRUE))
                     deseq2_sign_number_tb <- data.frame(result, number)
                     output$deseq2_sign_number <- renderTable(deseq2_sign_number_tb)
                     
                     incProgress(0.6, detail = "Processing 100%")
                     Sys.sleep(2)
                   })
      # renew MA plot
      observeEvent(input$enter_deseq2_ma_title,{
        if (!is.null(rv$deseq2_ma)) {
          rv$deseq2_ma <- rv$deseq2_ma +
          ggtitle(input$enter_deseq2_ma_title)
        output$display_deseq2_ma <- renderPlot({
          print(rv$deseq2_ma)})
        }
      })
    }
  })
  
  ## generate count plot
  observeEvent(input$generate_count_plot,{
    if (is.null(input$enter_gene_name) | is.null(input$intgroup)) {
      showNotification("Please enter one gene name and Select one covariate",
                       type = "error",
                       duration = NULL)
    }else if (input$enter_gene_name %in% rownames(rv$dds_trimmed)  == FALSE ) {
      showNotification("Please enter the exact correct gene name",
                       type = "error",
                       duration = NULL)
    }else{
      dt <- plotCounts(rv$dds_res, gene = input$enter_gene_name,
                       intgroup = input$intgroup, 
                       returnData = TRUE)
      if (length(input$intgroup) > 1) {
        dt$group <- apply(dt[, -1], 
                          1, 
                          paste,
                          collapse = "_")
      }
      
      rv$count_plot <- ggplot(dt,
                              aes(x = dt[,ncol(dt)],
                                  y = count)) +
        geom_point(position = position_jitter(width = 0.1, height = 0)) +
        labs(y = "normalized count", 
             x = NULL)
      output$deseq2_count_plot  <- renderPlot({print(rv$count_plot)})
    }
  })
  
  ## ---------------------- RNA: download results ----------------
  output$download_expl <- downloadHandler(
    filename = function() {
      paste("exploratory_results", 
            ".zip", 
            sep = "")
    },
    content = function(fname) {
      
      file_path1 <- paste(input$project_name,
                          "_expl_count_hist",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path1,
           height = 6,
           width = 10,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$count_hist)
      dev.off()
      
      file_path2 <- paste(input$project_name,
                          "_expl_boxplot",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path2,
           height = 6,
           width = 6,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$boxplot)
      dev.off()
      
      file_path3 <- paste(input$project_name,
                          "_expl_pca",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path3,
           height = 6,
           width = 6,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$pca)
      dev.off()
      
      file_path4 <- paste(input$project_name,
                          "_expl_heatmap",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path4,
           height = 6,
           width = 6,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$expl_heatmap)
      dev.off()
      
      fs <- c(file_path1, file_path2, file_path3, file_path4)
      zip(zipfile = fname, files = fs)
    },
    contentType = "application/zip"
  )
  
  
  output$download_degseq <- downloadHandler(
    filename = function() {
      paste("degseq_results", ".zip", sep = "")
    },
    content = function(fname) {
      
      file_path1 <- paste(input$project_name,
                          "_degseq_ma1",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path1,
           height = 6,
           width = 6,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$ma1_degseq)
      dev.off()
      
      file_path2 <- paste(input$project_name,
                          "degseq_trt2_trt1.csv",
                          sep = "")
      write.csv(rv$table_trt2_trt1, file_path2, row.names = FALSE)
      
      fs <- c(file_path1,
              file_path2)
      zip(zipfile = fname, 
          files = fs)
    },
    contentType = "application/zip"
  )
  

  
  output$download_deseq2 <- downloadHandler(
    filename = function() {
      paste("deseq2_results", 
            ".zip", 
            sep = "")
    },
    content = function(fname) {
      
      file_path1 <- paste(input$project_name,
                          "_deseq2_result_table",
                          ".csv",
                          sep = "")
      write.csv(rv$dtf_res_contrast,
                file_path1,
                row.names = FALSE)
      
      
      file_path2 <- paste(input$project_name,
                          "_deseq2_ma",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path2,
           height = 6,
           width = 6,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$deseq2_ma)
      dev.off()
      
      if (!is.null(rv$count_plot)) {
        file_path4 <- paste(input$project_name,
                            "_deseq2_count_plot.tiff",
                            sep = "")
        tiff(filename = file_path4,
             height = 6,
             width = 6,
             units = 'in',
             res = 300,
             compression = "lzw+p")
        print(rv$count_plot)
        dev.off()
      }else{
        file_path4 <- NULL
      }
      
      fs <- c(file_path1, 
              file_path2,
              file_path4)
      zip(zipfile = fname,
          files = fs)
    },
    contentType = "application/zip"
  )
  
  ## ----------- change in gene expre analysis ------------
  observeEvent(input$rna_contrast_1,{
    infile_rna_contrast_1 <- input$rna_contrast_1
    tb <- read.table(infile_rna_contrast_1$datapath,
                     sep = ",",
                     header = T,
                     quote = "\"")
    rv$tb_rna_contrast_1 <- as_tibble(tb)
    output$display_rna_contrast_1 = DT::renderDataTable({
      DT::datatable(rv$tb_rna_contrast_1,
                    options = list(scrollX = TRUE))})
    updateSelectInput(session,
                      inputId = "gene_col_selected_rna_contrast_1",
                      label = "Select gene column:",
                      choices = colnames(rv$tb_rna_contrast_1),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "logfold_change_selected_contrast_1",
                      label = "Select Logfold change column:",
                      choices = colnames(rv$tb_rna_contrast_1),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "p_value_selected_rna_contrast_1",
                      label = "Select P value column:",
                      choices = colnames(rv$tb_rna_contrast_1),
                      selected = NULL)})
  
  observeEvent(input$rna_contrast_2,{
    infile_rna_contrast_2 <- input$rna_contrast_2
    tb <- read.table(infile_rna_contrast_2$datapath,
                     sep = ",",
                     header = T,
                     quote = "\"")
    rv$tb_rna_contrast_2 <- as_tibble(tb)
    output$display_rna_contrast_2 = DT::renderDataTable({
      DT::datatable(rv$tb_rna_contrast_2,
                    options = list(scrollX = TRUE))})
    updateSelectInput(session,
                      inputId = "gene_col_selected_rna_contrast_2",
                      label = "Select gene column:",
                      choices = colnames(rv$tb_rna_contrast_2),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "logfold_change_selected_contrast_2",
                      label = "Select Logfold change column:",
                      choices = colnames(rv$tb_rna_contrast_2),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "p_value_selected_rna_contrast_2",
                      label = "Select P value column:",
                      choices = colnames(rv$tb_rna_contrast_2),
                      selected = NULL)})
  
  observeEvent(input$inner_join_2contrasts_rna,{
    
    contrast1 <- try(
      rv$tb_rna_contrast_1 %>%
      filter( !!as.name(input$p_value_selected_rna_contrast_1) < input$p_thresh_change_in_gene_expr_contrast1 &
              abs(!!as.name(input$logfold_change_selected_contrast_1)) >= input$logfold_change_thresh_change_in_gene_expr_contrast1  ) %>% 
      select(input$gene_col_selected_rna_contrast_1,
             input$logfold_change_selected_contrast_1) %>%
      rename("gene" = input$gene_col_selected_rna_contrast_1,
             "log_foldchange_1" = input$logfold_change_selected_contrast_1),
      silent = TRUE
    )
    
    if (inherits(contrast1, "try-error")) {
      showNotification("Please choose the correct columns for contrast1 and try again",
                       type = "error",
                       duration = NULL)}  
    
    contrast2 <- try(
      rv$tb_rna_contrast_2 %>%
      filter( !!as.name(input$p_value_selected_rna_contrast_2) < input$p_thresh_change_in_gene_expr_contrast2 &
              abs(!!as.name(input$logfold_change_selected_contrast_2)) >= input$logfold_change_thresh_change_in_gene_expr_contrast2  ) %>% 
      select(input$gene_col_selected_rna_contrast_2,
             input$logfold_change_selected_contrast_2) %>%
      rename("gene" = input$gene_col_selected_rna_contrast_2,
             "log_foldchange_2" = input$logfold_change_selected_contrast_2),
      silent = TRUE
    )
    
    if (inherits(contrast2, "try-error")) {
      showNotification("Please choose the correct columns for contrast2 and try again",
                       type = "error",
                       duration = NULL)}  
    
    rv$joined_2contrasts_rna <- try(
      inner_join(contrast1, contrast2, by = "gene") %>%
        filter((log_foldchange_1*log_foldchange_2) < 0),
      silent = TRUE)
    
    if (inherits(rv$joined_2contrasts_rna, "try-error")) {
      showNotification(rv$joined_2contrasts_rna[1],
                       type = "error",
                       duration = NULL)

    }else{
      rv$rna_contrast1_up_num <-  nrow(filter(contrast1, log_foldchange_1 > 0))
      rv$rna_contrast1_down_num <-  nrow(filter(contrast1, log_foldchange_1 < 0))
      rv$rna_contrast2_up_num <-  nrow(filter(contrast2, log_foldchange_2 > 0))
      rv$rna_contrast2_down_num <-  nrow(filter(contrast2, log_foldchange_2 < 0))
      rv$rna_up_donw <- nrow(filter(rv$joined_2contrasts_rna, log_foldchange_1 > 0))
      rv$rna_down_up <- nrow(filter(rv$joined_2contrasts_rna, log_foldchange_1 < 0))
      
      output$display_summary_inner_joined_2contrasts_rna <- renderPrint({summary(rv$joined_2contrasts_rna)})
    }
    
    })
  
  
  observeEvent(input$generate_change_in_gene_expr_plot,{
    if (nrow(rv$joined_2contrasts_rna) == 0) {
      showNotification("There is no genes left to plot. Please reset your threshold.",
                       type = "error",
                       duration = NULL)
    } else {
      output$display_venn_diagram1_rna <- renderPlot({
      p1 <- draw.pairwise.venn(area1 = rv$rna_contrast1_up_num,
                               area2 = rv$rna_contrast2_down_num,
                               cross.area = rv$rna_up_donw,
                               scaled = TRUE,
                               col = c("green3", "firebrick"),
                               cex = rep(2, 3))
      rv$venn_diagram1_rna <- grid.arrange(gTree(children = p1), top = textGrob("Trt2-Trt1 Up Trt3-Trt2 Down",gp = gpar(fontsize = 18)))
    })
      
    output$display_venn_diagram2_rna <- renderPlot({
      p2 <- draw.pairwise.venn(area1 = rv$rna_contrast1_down_num,
                                       area2 = rv$rna_contrast2_up_num,
                                       cross.area = rv$rna_down_up,
                                       scaled = TRUE,
                                       col = c("firebrick", "green3"),
                                       cex = rep(2, 3))
      rv$venn_diagram2_rna <- grid.arrange(gTree(children = p2), top = textGrob("Trt2-Trt1 Down Trt3-Trt2 Up",gp = gpar(fontsize = 18)))
    })
    rv$change_in_gene_expr_heatmap <-  two_column_heatmap(rv$joined_2contrasts_rna,
                                                          gene,
                                                          log_foldchange_1,
                                                          log_foldchange_2,
                                                          input$change_in_gene_expr_plot_title,
                                                          input$change_in_gene_expr_plot_x_text_1,
                                                          input$change_in_gene_expr_plot_x_text_2)
    
    p <- ggplotly(rv$change_in_gene_expr_heatmap) %>% 
      layout(height = 800, width = 800)
    output$change_in_gene_expr_plot <- renderPlotly({print(p)})
    }
  })
  
  # update change_in_gene_expr_plot title and x tick
  observeEvent(input$change_in_gene_expr_plot_title,{
    if (!is.null(rv$change_in_gene_expr_heatmap)) {
      rv$change_in_gene_expr_heatmap <- rv$change_in_gene_expr_heatmap +
        ggtitle(input$change_in_gene_expr_plot_title)
      p <- ggplotly(rv$change_in_gene_expr_heatmap) %>% 
      layout(height = 800, width = 800)
    output$change_in_gene_expr_plot <- renderPlotly({print(p)})
    }
  })
  
  observeEvent(input$change_in_gene_expr_plot_x_text_1,{
    if (!is.null(rv$change_in_gene_expr_heatmap)) {
      rv$change_in_gene_expr_heatmap <- rv$change_in_gene_expr_heatmap +
        scale_x_discrete(labels = c(input$change_in_gene_expr_plot_x_text_1, input$change_in_gene_expr_plot_x_text_2))
      p <- ggplotly(rv$change_in_gene_expr_heatmap) %>% 
      layout(height = 800, width = 800)
    output$change_in_gene_expr_plot <- renderPlotly({print(p)})
    }
  })
  
  observeEvent(input$change_in_gene_expr_plot_x_text_2,{
    if (!is.null(rv$change_in_gene_expr_heatmap)) {
      rv$change_in_gene_expr_heatmap <- rv$change_in_gene_expr_heatmap +
        scale_x_discrete(labels = c(input$change_in_gene_expr_plot_x_text_1, input$change_in_gene_expr_plot_x_text_2))
      p <- ggplotly(rv$change_in_gene_expr_heatmap) %>% 
      layout(height = 800, width = 800)
    output$change_in_gene_expr_plot <- renderPlotly({print(p)})
    }
  })
  
  
  output$download_change_in_gene_expr_plot <- downloadHandler(
    filename = function() {
      paste("change_in_gene_expr_results", ".zip", sep = "")
    },
    content = function(fname) {
      
      file_path1 <- paste(input$project_name,
                          "_change_in_gene_expr_plot",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path1,
           height = 10,
           width = 10,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$change_in_gene_expr_heatmap)
      dev.off()
      
      file_path2 <- paste(input$project_name,
                          "_change_in_gene_expr_table",
                          ".csv",
                          sep = "")
      write.csv(rv$joined_2contrasts_rna, 
                file_path2, 
                row.names = FALSE)
      file_path3 <- paste(input$project_name,
                          "_rna_venn_up_donw.tiff",
                          sep = "")
      tiff(filename = file_path3,
           height = 4,
           width = 5,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      grid.draw(rv$venn_diagram1_rna)
      dev.off()
      
      file_path4 <- paste(input$project_name,
                          "_rna_venn_donw_up.tiff",
                          sep = "")
      tiff(filename = file_path4,
           height = 4,
           width = 5,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      grid.draw(rv$venn_diagram2_rna)
      dev.off()
      
      fs <- c(file_path1, 
              file_path2,
              file_path3,
              file_path4)
      
      zip(zipfile = fname, 
          files = fs)
    },
    contentType = "application/zip"
  )
  
  ## ------------  DNA mettyl: Read in data and update selectinput------------
  
  observeEvent(input$dna_table,{
    infile_dna_table <- input$dna_table
    tb <- fread(infile_dna_table$datapath)
    rv$dna_table <- as.tibble(tb)
    
    output$display_dna_table = DT::renderDataTable({
      DT::datatable(rv$dna_table,
                    options = list(scrollX = TRUE))})
    
    
    updateSelectInput(
      session,
      inputId = "dna_sample_cols1",
      label = "Select sample columns for dtN",
      choices = colnames(rv$dna_table),
      selected = NULL
    )
    updateSelectInput(
      session,
      inputId = "dna_sample_cols2",
      label = "Select sample columns for dtX (same order as dtN)",
      choices = colnames(rv$dna_table),
      selected = NULL
    )
    
    updateSelectInput(
      session,
      inputId = "dna_heat_N",
      label = "Select sample columns for dtN",
      choices = colnames(rv$dna_table),
      selected = NULL
    )
    updateSelectInput(
      session,
      inputId = "dna_heat_X",
      label = "Select sample columns for dtX (same order as dtN)",
      choices = colnames(rv$dna_table),
      selected = NULL
    )
    
    
    updateSelectInput(
      session,
      inputId = "dss_chr",
      label = "Select the chromosome column  (chr, seqname, etc)",
      choices = colnames(rv$dna_table),
      selected = NULL
    )
    updateSelectInput(
      session,
      inputId = "dss_start",
      label = "Select the start position column (start, etc)",
      choices = colnames(rv$dna_table),
      selected = NULL
    )
    updateSelectInput(
      session,
      inputId = "dss_region",
      label = "Select the region (feature) column (annotation, etc)",
      choices = colnames(rv$dna_table),
      selected = NULL
    )
    updateSelectInput(
      session,
      inputId = "dss_gene",
      label = "Select the gene name column (geneid, etc)",
      choices = colnames(rv$dna_table),
      selected = NULL
    )
    updateSelectInput(
      session,
      inputId = "dss_comp1",
      label = "Select N and X column(s) for group 1 (for each sample, select N FIRST then X)",
      choices = colnames(rv$dna_table),
      selected = NULL
    )
    updateSelectInput(
      session,
      inputId = "dss_comp2",
      label = "Select N and X column(s) for group 2 (for each sample, select N FIRST then X)",
      choices = colnames(rv$dna_table),
      selected = NULL
    )
    updateSelectInput(
      session,
      inputId = "dna_feature_col",
      label = "Select region column",
      choices = colnames(rv$dna_table),
      selected = NULL
    )
    updateSelectInput(
      session,
      inputId = "dna_cpg_col",
      label = "Select CpG column",
      choices = colnames(rv$dna_table),
      selected = NULL
    )
    
  })
  
  
  
                              
  
  ## ----------------------- Explotary analysis ------------------------
  observeEvent(input$generate_expl_dna,{
    
    if (length(input$dna_sample_cols1) != length(input$dna_sample_cols2) ) {
      showNotification("Number of samples selected for dtN is not equal to number of samples selected for dtX.",
                       type = "error",
                       duration = NULL)
    } else if (length(input$dna_sample_cols1) != length(unlist(strsplit(input$sample_name, ","))) ) {
      showNotification("Number of names entered is not equal to number of samples selected.",
                       type = "error",
                       duration = NULL)
    } else if (!is.numeric(pull(rv$dna_table[, input$dna_cpg_col]))) {
      showNotification("Selected CpG column is not numerical. Please select the right column.",
                       type = "error",
                       duration = NULL)
    } else if (!is.character(pull(rv$dna_table[, input$dna_feature_col]))) {
      showNotification("Selected region column is not character. Please select the right column.",
                       type = "error",
                       duration = NULL)
    } else {
      withProgress(message = "Generating Result",
                 value = 0,
                 expr = {
                   incProgress(0.3, detail = "Generate annotation percentage pie chart")
                   # anno pie chart (not related to samples)
                   t1 <- rv$dna_table %>%
                     rename(region = !!as.name(input$dna_feature_col)) %>%
                     mutate(region = replace(region,
                                             str_detect(region, "Intro"), "Intron")) %>%
                     mutate(region = replace(region,
                                             str_detect(region, "Exon"), "Exon")) %>%
                     count(region) %>%
                     rename(Frequency = n) %>%
                     mutate(Frequency = Frequency / sum(Frequency) * 100) %>%
                     mutate(region = fct_reorder(region, Frequency))
                   
                   rv$anno_by_reg_plot <-
                     ggplot(t1, aes(x = "", y = Frequency, fill = region)) +
                     geom_bar(width = 1,
                              stat = "identity",
                              color = "black") +
                     coord_polar("y",
                                 start = 0,
                                 direction = 1) +
                     scale_x_discrete("") +
                     ggtitle("Annotation by Region (%)")
                   
                   output$anno_by_reg <- renderPlot({
                     print(rv$anno_by_reg_plot)
                   })
                   
                   incProgress(0.3, detail = "Generate CpG histogram")
                   # cpg hist (not related to samples)
                   t2 <- rv$dna_table %>%
                     select(
                       cpg = !!as.name(input$dna_cpg_col),
                       region = !!as.name(input$dna_feature_col)
                     ) %>%
                     mutate(
                       region = case_when(
                         str_detect(region, "Intr") | str_detect(region, "Exon") ~ "Body",
                         str_detect(region, "3' UTR") ~ "3' UTR",
                         str_detect(region, "5' UTR") ~ "5' UTR",
                         str_detect(region, "Dist") |
                           str_detect(region, "Down") ~ "Downstream",
                         str_detect(region, "Promoter") ~ "Promoter"
                       )
                     )
                   
                   rv$cpg_hist <- ggplot(t2,
                                         aes(x = cpg)) +
                     facet_wrap(~ region, scales = "free_y") +
                     geom_histogram(color = "black",
                                    fill = "grey",
                                    binwidth = 5) +
                     scale_x_continuous(name = "Number of CpG") +
                     scale_y_continuous(name = "Frequency") +
                     ggtitle("Distribution of CpG by Region")
                   
                   output$cpg_hist <- renderPlot({
                     print(rv$cpg_hist)
                   })
                   
                   
                   incProgress(0.2, detail = "Generate average methy% by region plot")
                   # avergae methy% by region (relate to samples selected)
                   # dtN
                   t3 <- rv$dna_table %>%
                     select(input$dna_sample_cols1)
                   t3 <- apply(t3,
                               2,
                               function(col_t3) {
                                 col_t3[col_t3 == 0] <- NA
                                 return(col_t3)
                               })
                   
                   #dtX
                   t4 <- rv$dna_table %>%
                     select(input$dna_sample_cols2)
                   t4 <- apply(t4,
                               2,
                               function(col_t4) {
                                 col_t4[col_t4 == 0] <- NA
                                 return(col_t4)
                               })
                   # pct
                   pct <- as.tibble(as.matrix(t4) / as.matrix(t3))
                   names(pct) <- unlist(strsplit(input$sample_name, ","))
                   
                   t5 <- rv$dna_table %>%
                     select(region = !!as.name(input$dna_feature_col)) %>%
                     mutate(
                       region = case_when(
                         str_detect(region, "Intr") | str_detect(region, "Exon") ~ "Body",
                         str_detect(region, "3' UTR") ~ "3' UTR",
                         str_detect(region, "5' UTR") ~ "5' UTR",
                         str_detect(region, "Dist") |
                           str_detect(region, "Down") ~ "Downstream",
                         str_detect(region, "Promoter") ~ "Promoter"
                       )
                     )
                   
                   pct_tb <- cbind(t5, pct)
                   
                   pct_tb <- pct_tb %>%
                     drop_na() %>%
                     group_by(region) %>%
                     summarise_all(mean) %>%
                     gather(key = "trt", value = "Methylation (%)",-region)
                   
                   
                   rv$meth_by_region <- ggplot(pct_tb,
                                               aes(
                                                 x = region,
                                                 y = 100 * `Methylation (%)`,
                                                 group = trt,
                                                 fill = trt
                                               )) +
                     geom_bar(position = position_dodge(),
                              stat = "identity") +
                     scale_x_discrete("Region") +
                     scale_y_continuous("Average Methylation (%)",
                                        limits = c(0, 100)) +
                     scale_fill_discrete("Sample") +
                     theme(plot.title = element_text(hjust = 0.5),
                           axis.text.x = element_text(angle = 45,
                                                      hjust = 1))
                   
                   output$meth_by_region <- renderPlot({
                     print(rv$meth_by_region)
                   })
                   
                   incProgress(0.2, detail = "Generate average methy% by region plot")
                   Sys.sleep(1)
                 }) # end of withProgress
    }
    
    
  })
  
  
  
  ## ----------------------- DSS -------------  ## UI at line 860
  fdtdss <- reactive({
    
    comp1 <- input$dss_comp1 ## list of columns for group a from user input. N X N X ...
    comp2 <- input$dss_comp2 ## list of columns for group b from user input. N X N X ...
    name1 <- input$dss_comp1_name ## group name for group a from user input.
    name2 <- input$dss_comp2_name ## group name for group b from user input.
    
    
    
    rv$dss_w_gene <- rv$dna_table %>% 
      select(input$dss_comp1, input$dss_comp2, input$dss_chr, input$dss_start, input$dss_gene, input$dss_region) %>% 
      drop_na()
    
    ### below added by R Wu #####################
   names(rv$dss_w_gene)[(length(rv$dss_w_gene)-3):length(rv$dss_w_gene)] <- c("geneChr", "start", "gene", "region")
    
    # construct N data.frame instead of only two. 
    # first need to identify how many items in 
    
    # rv$dss_w_gene is the 
    comp1_samples <- substr(comp1, 1, nchar(comp1)-2)
    comp1_samples <- comp1_samples[!duplicated(comp1_samples)]
    
    comp2_samples <- substr(comp2, 1, nchar(comp2)-2)
    comp2_samples <- comp2_samples[!duplicated(comp2_samples)]
    
    combined_samples <- c(comp1_samples, comp2_samples)
    
    sample_list1 <- vector('list', length(comp1))
    i = 1
    while (i < length(comp1)) {
      df_temp <- data.frame(chr = rv$dss_w_gene$geneChr,
                            pos = rv$dss_w_gene$start,
                            N = rv$dss_w_gene[[comp1[i]]],
                            X = rv$dss_w_gene[[comp1[i+1]]]
      )
      sample_list1[[i]] <- df_temp
      i <- i + 2
    }
    sample_list1[sapply(sample_list1, is.null)] <- NULL
    
    sample_list2 <- vector('list', length(comp2))
    i = 1
    while (i < length(comp2)) {
      df_temp <- data.frame(chr = rv$dss_w_gene$geneChr,
                            pos = rv$dss_w_gene$start,
                            N = rv$dss_w_gene[[comp2[i]]],
                            X = rv$dss_w_gene[[comp2[i+1]]]
      )
      sample_list2[[i]] <- df_temp
      i <- i + 2
    }
    sample_list2[sapply(sample_list2, is.null)] <- NULL
    
    bsdata <- makeBSseqData(c(sample_list1, sample_list2), combined_samples)
    name1 <- comp1_samples
    name2 <- comp2_samples
    
#### end addition by R Wu ####
    
## below original code ### 
    
    # names(rv$dss_w_gene) <- c("sample1_N", "sample1_X", "sample2_N", "sample2_X", "geneChr", "start", "gene", "region")
    # df1 <- data.frame(chr = rv$dss_w_gene$geneChr,
    #                   pos = rv$dss_w_gene$start,
    #                   N = rv$dss_w_gene$sample1_N,
    #                   X = rv$dss_w_gene$sample1_X)
    # df2 <- data.frame(chr = rv$dss_w_gene$geneChr,
    #                   pos = rv$dss_w_gene$start,
    #                   N = rv$dss_w_gene$sample2_N,
    #                   X = rv$dss_w_gene$sample2_X)
    # frames <- list(df1, df2)
    # names(frames) <- c(name1, name2)
    # bsdata <- makeBSseqData(frames, names(frames))
    
    # perform DML test, differntially methylated loci (DML) for two group comparisons of bisulfite sequencing (BS-seq)
    # without replicates, must set equal.disp=T
    
##### end commented original code ## 
    
    dml <- DMLtest(bsdata, 
                   group1 = name1,
                   group2 = name2,
                   # equal.disp = TRUE,
                   equal.disp = input$dss_wo_rep, ## R Wu
                   smoothing = TRUE, 
                   smoothing.span = 500)
    
    ## below added by R Wu
    # change names to group names. 
    names(dml)[3:4] <- c(paste0(input$dss_comp1_name, ".mu"), paste0(input$dss_comp2_name, ".mu")) ## w
    ## end addition by R Wu
    dml
  })
  
  observeEvent(input$generate_dss, {
    withProgress(message = "Generating Result table",
                 value = 0,
                 expr = {
                   
                   incProgress(0.2, detail = "Smoothing")
                   incProgress(0.3, detail = "Estimating dispersion for each CpG site, this will take a while...")
                   dml <- fdtdss()
                   
                   dml <- as.tibble(dml) %>% 
                     unite("chr_start", chr, pos)
                   
                   # map gene and feature/region column from input table to dml result table 
                   gene_region <- rv$dss_w_gene %>% 
                     mutate(start = as.character(start)) %>% 
                     unite("chr_start", geneChr, start) %>% 
                     select(chr_start, region, gene)
                     
                   rv$dml_gene_mapped_back <- left_join(dml, gene_region, by = "chr_start")
                   
                   
                   output$dss_dml_tb <- DT::renderDataTable({
                     datatable(rv$dml_gene_mapped_back,
                               options = list(scrollX = TRUE))
                   })
                   
                   incProgress(0.5, detail = "Estimating dispersion for each CpG site, this will take a while...")
                   Sys.sleep(1)
                   })
  })
  
  
  output$download_dss <- downloadHandler(
    filename = function() {
      paste(input$project_name_dna,
            "_dss_result_table",
            ".csv",
            sep = "")
    },
    content = function(file) {
      write.csv(rv$dml_gene_mapped_back,
                file = file,
                row.names = FALSE)
    }
  )

  
  # # This should work but doesn't
  # # Source: https://yihui.shinyapps.io/DT-info/
  
  
  ## --------------- Change-in-Methyl Heatmap -------------
  observeEvent(input$dna_contrast_1,{
    infile_dna_contrast_1 <- input$dna_contrast_1
    tb <- read.table(infile_dna_contrast_1$datapath,
                     sep = ",",
                     header = T,
                     quote = "\"")
    rv$tb_dna_contrast_1 <- as_tibble(tb)
    output$display_contrast_1 = DT::renderDataTable({
      DT::datatable(rv$tb_dna_contrast_1,
                    options = list(scrollX = TRUE))})
    updateSelectInput(session,
                      inputId = "gene_col_selected_dna_contrast_1",
                      label = "Select gene identifier column:",
                      choices = colnames(rv$tb_dna_contrast_1),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "gene_name_col_selected_dna_contrast_1",
                      label = "Select gene name column:",
                      choices = colnames(rv$tb_dna_contrast_1),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "region_col_selected_dna_contrast_1",
                      label = "Select region column:",
                      choices = colnames(rv$tb_dna_contrast_1),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "methyl_diff_dna_contrast_1",
                      label = "Select DNA methylation ratio difference column:",
                      choices = colnames(rv$tb_dna_contrast_1),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "pval_col_selected_dna_contrast_1",
                      label = "Select P value column",
                      choices = colnames(rv$tb_dna_contrast_1))})
  
  observeEvent(input$dna_contrast_2,{
    infile_dna_contrast_2 <- input$dna_contrast_2
    tb <- read.table(infile_dna_contrast_2$datapath,
                     sep = ",",
                     header = T,
                     quote = "\"")
    rv$tb_dna_contrast_2 <- as_tibble(tb)
    output$display_contrast_2 = DT::renderDataTable({
      DT::datatable(rv$tb_dna_contrast_2,
                    options = list(scrollX = TRUE))})
    updateSelectInput(session,
                      inputId = "gene_col_selected_dna_contrast_2",
                      label = "Select gene identifier column:",
                      choices = colnames(rv$tb_dna_contrast_2),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "region_col_selected_dna_contrast_2",
                      label = "Select region column:",
                      choices = colnames(rv$tb_dna_contrast_2),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "methyl_diff_dna_contrast_2",
                      label = "Select DNA methylation ratio difference column:",
                      choices = colnames(rv$tb_dna_contrast_2),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "pval_col_selected_dna_contrast_2",
                      label = "Select P value column",
                      choices = colnames(rv$tb_dna_contrast_2))})
  
  observeEvent(input$inner_join_2contrasts,{
    rv$contrast1 <- try(rv$tb_dna_contrast_1 %>%
                       filter(!!as.name(input$pval_col_selected_dna_contrast_1) < input$change_in_dmr_p_thresh1 &
                              abs(!!as.name(input$methyl_diff_dna_contrast_1)) >=  input$change_in_dmr_diff_thresh1) %>% 
                       select(input$gene_col_selected_dna_contrast_1,
                              input$gene_name_col_selected_dna_contrast_1,
                              input$methyl_diff_dna_contrast_1,
                              input$region_col_selected_dna_contrast_1) %>%
                       rename("gene" = input$gene_col_selected_dna_contrast_1,
                              "gene_name" = input$gene_name_col_selected_dna_contrast_1,
                              "methyl_diff_1" = input$methyl_diff_dna_contrast_1,
                              "region" = input$region_col_selected_dna_contrast_1),
                     silent = TRUE)
    
    if (inherits(rv$contrast1, "try-error")) {
      showNotification("Please choose the correct columns for contrast1 and try again",
                       type = "error",
                       duration = NULL)}  
    
    rv$contrast2 <- try(rv$tb_dna_contrast_2  %>%
                       filter(!!as.name(input$pval_col_selected_dna_contrast_2) < input$change_in_dmr_p_thresh2 &
                              abs(!!as.name(input$methyl_diff_dna_contrast_2)) >=  input$change_in_dmr_diff_thresh2) %>% 
                       select(input$gene_col_selected_dna_contrast_2,
                              input$methyl_diff_dna_contrast_2,
                              input$region_col_selected_dna_contrast_2) %>%
                       rename("gene" = input$gene_col_selected_dna_contrast_2,
                              "methyl_diff_2" = input$methyl_diff_dna_contrast_2,
                              "region2" = input$region_col_selected_dna_contrast_2),
                     silent = TRUE)
    
    if (inherits(rv$contrast2, "try-error")) {
      showNotification("Please choose the correct columns for contrast2 and try again",
                       type = "error",
                       duration = NULL)}  
    
    
    rv$joined_2contrasts <- try(
      inner_join(rv$contrast1, rv$contrast2, by = "gene") %>%
        filter((methyl_diff_1*methyl_diff_2) < 0) %>% 
        select(gene_name, methyl_diff_1, methyl_diff_2, region) %>% 
        rename(gene = gene_name) %>% 
        mutate(region = as.factor(region)),
      silent = T)

    if (inherits(rv$joined_2contrasts, "try-error")) {
      showNotification(rv$joined_2contrasts[1],
                       type = "error",
                       duration = NULL)

    }else{
      
      output$display_summary_inner_joined_2contrasts <- renderPrint({summary(rv$joined_2contrasts)})
      
      updateSelectInput(session,
                        inputId = "region_selcted_dna",
                        label = "Select region to plot:",
                        choices =  levels(rv$joined_2contrasts$region),
                        selected = NULL)
    }
    
    })

  
  observeEvent(input$rename_region_dna,{
    
    rv$contrast1 <- rv$contrast1 %>% 
      mutate(region = as.character(region)) %>% 
      mutate(region = replace(region, str_detect(region, input$region_name_contian_dna), input$region_new_name_dna))
    
    rv$contrast2 <- rv$contrast2 %>% 
      mutate(region2 = as.character(region2)) %>% 
      mutate(region2 = replace(region2, str_detect(region2, input$region_name_contian_dna), input$region_new_name_dna))
    
    rv$joined_2contrasts <- rv$joined_2contrasts %>%
      mutate(region = as.character(region)) %>% 
      mutate(region = replace(region, str_detect(region, input$region_name_contian_dna), input$region_new_name_dna)) %>% 
      mutate(region = as.factor(region))
    
    
    updateSelectInput(session,
                      inputId = "region_selcted_dna",
                      label = "Select region to plot:",
                      choices =  levels(rv$joined_2contrasts$region),
                      selected = NULL)
    
    output$display_summary_inner_joined_2contrasts <- renderPrint({summary(rv$joined_2contrasts)})
  })
  
  
  observeEvent(input$filter_region_dna,{
    rv$contrast1 <- rv$contrast1 %>%
      filter(region %in% input$region_selcted_dna)

    rv$contrast2 <- rv$contrast2 %>%
      filter(region2 %in% input$region_selcted_dna)
    
    rv$joined_2contrasts <- rv$joined_2contrasts %>% 
      filter(region %in% input$region_selcted_dna) %>% 
      droplevels()
    
  output$display_summary_inner_joined_2contrasts <- renderPrint({summary(rv$joined_2contrasts)})
  })
  
  observeEvent(input$generate_change_in_methyl_plot,{
    if (nrow(rv$joined_2contrasts) == 0) {
      showNotification("There is no genes left to plot. Please reset your threshold.",
                       type = "error",
                       duration = NULL)
    } else {
      
      rv$dna_contrast1_up_num <-  nrow(filter(rv$contrast1, methyl_diff_1 > 0))
      rv$dna_contrast1_down_num <-  nrow(filter(rv$contrast1, methyl_diff_1 < 0))
      rv$dna_contrast2_up_num <-  nrow(filter(rv$contrast2, methyl_diff_2 > 0))
      rv$dna_contrast2_down_num <-  nrow(filter(rv$contrast2, methyl_diff_2 < 0))
      rv$dna_up_donw <- nrow(filter(rv$joined_2contrasts, methyl_diff_1 > 0))
      rv$dna_down_up <- nrow(filter(rv$joined_2contrasts, methyl_diff_1 < 0))
      
      output$display_venn_diagram1_dna <- renderPlot({
      p1 <- draw.pairwise.venn(area1 = rv$dna_contrast1_up_num,
                               area2 = rv$dna_contrast2_down_num,
                               cross.area = rv$dna_up_donw,
                               scaled = TRUE,
                               col = c("green3", "firebrick"),
                               cex = rep(2, 3))
      rv$venn_diagram1_dna <- grid.arrange(gTree(children = p1), top = textGrob("Trt2-Trt1 Up Trt3-Trt2 Down",gp = gpar(fontsize = 18)))
    })
      
    output$display_venn_diagram2_dna <- renderPlot({
      p2 <- draw.pairwise.venn(area1 = rv$dna_contrast1_down_num,
                                       area2 = rv$dna_contrast2_up_num,
                                       cross.area = rv$dna_down_up,
                                       scaled = TRUE,
                                       col = c("firebrick", "green3"),
                                       cex = rep(2, 3))
      rv$venn_diagram2_dna <- grid.arrange(gTree(children = p2), top = textGrob("Trt2-Trt1 Down Trt3-Trt2 Up",gp = gpar(fontsize = 18)))
    })
    
    rv$change_in_methyl_heatmap <-  two_column_heatmap(rv$joined_2contrasts,
                                                             gene,
                                                             methyl_diff_1,
                                                             methyl_diff_2,
                                                             input$change_in_methyl_plot_title,
                                                             input$change_in_methyl_plot_x_text_1,
                                                             input$change_in_methyl_plot_x_text_2)
    
    p <- ggplotly(rv$change_in_methyl_heatmap) %>% 
      layout(height = 800, width = 800)
    output$change_in_methyl_plot <- renderPlotly({print(p)})
    }
  })
  
  # update change_in_methyl_heatmap title and x tick
  observeEvent(input$change_in_methyl_plot_title,{
    if (!is.null(rv$change_in_methyl_heatmap)) {
      rv$change_in_methyl_heatmap <- rv$change_in_methyl_heatmap +
        ggtitle(input$change_in_methyl_plot_title)
      p <- ggplotly(rv$change_in_methyl_heatmap) %>% 
      layout(height = 800, width = 800)
    output$change_in_methyl_plot <- renderPlotly({print(p)})
    }
  })
  
  observeEvent(input$change_in_methyl_plot_x_text_1,{
    if (!is.null(rv$change_in_methyl_heatmap)) {
      rv$change_in_methyl_heatmap <- rv$change_in_methyl_heatmap +
        scale_x_discrete(labels = c(input$change_in_methyl_plot_x_text_1, input$change_in_methyl_plot_x_text_2))
      p <- ggplotly(rv$change_in_methyl_heatmap) %>% 
      layout(height = 800, width = 800)
    output$change_in_methyl_plot <- renderPlotly({print(p)})
    }
  })
  
  observeEvent(input$change_in_methyl_plot_x_text_2,{
    if (!is.null(rv$change_in_methyl_heatmap)) {
      rv$change_in_methyl_heatmap <- rv$change_in_methyl_heatmap +
        scale_x_discrete(labels = c(input$change_in_methyl_plot_x_text_1, input$change_in_methyl_plot_x_text_2))
      p <- ggplotly(rv$change_in_methyl_heatmap) %>% 
      layout(height = 800, width = 800)
    output$change_in_methyl_plot <- renderPlotly({print(p)})
    }
  })
  
  output$download_change_in_methyl_plot <- downloadHandler(
    filename = function() {
      paste("change_in_methyl_results", ".zip", sep = "")
    },
    content = function(fname) {
      
      file_path1 <- paste(input$project_name,
                          "_change_in_methyl_plot",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path1,
           height = 10,
           width = 10,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$change_in_methyl_heatmap)
      dev.off()
      
      file_path2 <- paste(input$project_name,
                          "_change_in_methyl_table",
                          ".csv",
                          sep = "")
      write.csv(rv$joined_2contrasts, 
                file_path2, 
                row.names = FALSE)
      
      file_path3 <- paste(input$project_name,
                          "_dna_venn_up_donw.tiff",
                          sep = "")
      tiff(filename = file_path3,
           height = 4,
           width = 5,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      grid.draw(rv$venn_diagram1_dna)
      dev.off()
      
      file_path4 <- paste(input$project_name,
                          "_dna_venn_donw_up.tiff",
                          sep = "")
      tiff(filename = file_path4,
           height = 4,
           width = 5,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      grid.draw(rv$venn_diagram2_dna)
      dev.off()
      
      fs <- c(file_path1,
              file_path2,
              file_path3,
              file_path4)
      zip(zipfile = fname, 
          files = fs)
    },
    contentType = "application/zip"
  )
  
  ## -------------  RNA vs DNA: Read in data and update SelectInput ------------
  
  observeEvent(input$rna_sig,{
    infile_rna_sig <- input$rna_sig
    tb <- read.table(infile_rna_sig$datapath,
                     sep = ",",
                     header = T,
                     quote = "\"")
    rv$tb_rna_sig <- as_tibble(tb)
    output$display_rna_sig = DT::renderDataTable({
      DT::datatable(rv$tb_rna_sig,
                    options = list(scrollX = TRUE))})
    updateSelectInput(session,
                      inputId = "rna_gene_col_selected_rna_vs_dna",
                      label = "From RNA table select gene column:",
                      choices = colnames(rv$tb_rna_sig),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "logfoldchange_col_selected_rna_vs_dna",
                      label = "Select RNA expression difference column:",
                      choices = colnames(rv$tb_rna_sig),
                      selected = NULL)})
  
  observeEvent(input$dna_sig,{
    infile_dna_sig <- input$dna_sig
    tb <- read.table(infile_dna_sig$datapath,
                     sep = ",",
                     header = T,
                     quote = "\"")
    rv$tb_dna_sig <- as_tibble(tb)
    output$display_dna_sig = DT::renderDataTable({
      DT::datatable(rv$tb_dna_sig,
                    options = list(scrollX = TRUE))})
    updateSelectInput(session,
                      inputId = "dna_gene_col_selected_rna_vs_dna",
                      label = "From DNA table select gene column:",
                      choices = colnames(rv$tb_dna_sig),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "methyl_diff_col_selected_rna_vs_dna",
                      label = "Select DNA methylation ratio difference column:",
                      choices = colnames(rv$tb_dna_sig),
                      selected = NULL)
    updateSelectInput(session,
                      inputId = "region_col_selected_rna_vs_dna",
                      label = "Select region column:",
                      choices = colnames(rv$tb_dna_sig),
                      selected = NULL)})
  
  observeEvent(input$inner_join_rna_vs_dna,{
    rna_sig <- rv$tb_rna_sig %>%
      select(input$rna_gene_col_selected_rna_vs_dna,
             input$logfoldchange_col_selected_rna_vs_dna) %>%
      rename("gene" = input$rna_gene_col_selected_rna_vs_dna,
             "RNA_exp_diff" = input$logfoldchange_col_selected_rna_vs_dna)
    
    dna_sig <- rv$tb_dna_sig  %>%
      select(input$dna_gene_col_selected_rna_vs_dna,
             input$methyl_diff_col_selected_rna_vs_dna,
             input$region_col_selected_rna_vs_dna)  %>%
      rename("gene" = input$dna_gene_col_selected_rna_vs_dna,
             "DNA_methyl_diff" = input$methyl_diff_col_selected_rna_vs_dna,
             "region" = input$region_col_selected_rna_vs_dna)
    
    rv$joined_rna_dna <- try(
      inner_join(rna_sig, dna_sig, by = "gene"),
      silent = T)
    
    if (inherits(rv$joined_rna_dna, "try-error")) {
      showNotification("Please choose the correct columns and try again",
                       type = "error",
                       duration = NULL)

    }else{
      rv$joined_rna_dna %>%
        mutate(region = as.factor(region)) -> rv$joined_rna_dna
      output$display_summary_inner_joined_rna_dna <- renderPrint({summary(rv$joined_rna_dna)})
      
      
      updateSelectInput(session,
                        inputId = "region_selcted_dna_rna",
                        label = "Select region to plot:",
                        choices =  levels(rv$joined_rna_dna$region),
                        selected = NULL)
    }
    })
  
  
  
  observeEvent(input$rename_region_rna_vs_dna,{
    rv$joined_rna_dna <- rv$joined_rna_dna %>%
      mutate(region = as.character(region)) %>%
      mutate(region = replace(
        region,
        str_detect(region,
                   input$region_name_contian_rna_vs_dna),
        input$region_new_name_rna_vs_dna
      )) %>%
      mutate(region =  as.factor(region))
    
    output$display_summary_inner_joined_rna_dna <- renderPrint({summary(rv$joined_rna_dna)})
    
    updateSelectInput(session,
                      inputId = "region_selcted_dna_rna",
                      label = "Select region to plot:",
                      choices = levels(rv$joined_rna_dna$region),
                      selected = NULL)
  })
  
  observeEvent(input$filter_region_dna_rna,{
    
    rv$joined_rna_dna <- rv$joined_rna_dna %>%
      filter(region %in% input$region_selcted_dna_rna) %>% 
      droplevels()
    
    output$display_summary_inner_joined_rna_dna <- renderPrint({summary(rv$joined_rna_dna)})
      
  })
  
  observeEvent(input$generate_starburst_plot_rna_vs_dna,{
    if (nrow(rv$joined_rna_dna) == 0) {
      showNotification("There is no genes left to plot",
                       type = "error",
                       duration = NULL)
    }else {
      rv$starburst_plot <- starburst(rv$joined_rna_dna,
                                   gene,
                                   RNA_exp_diff,
                                   DNA_methyl_diff,
                                   region,
                                   input$starburst_rna_hline,
                                   input$starburst_dna_vline,
                                   input$starburst_title)
    p1 <- ggplotly(rv$starburst_plot,
                   tooltip = c("text", "x", "y")) %>%
      layout(height = 800, width = 800)
    output$starburst_plot <- renderPlotly({print(p1)})
    }
  })
  
  
  observeEvent(input$input$starburst_title,{
    if (!is.null(rv$starburst_plot)) {
      rv$starburst_plot <- rv$starburst_plot +
        ggtitle(input$starburst_title)

      p1 <- ggplotly(rv$starburst_plot,
                   tooltip = c("text", "x", "y")) %>%
      layout(height = 800, width = 800)
      output$starburst_plot <- renderPlotly({print(p1)})
    }
  })
  
  
  
  output$download_starburst_plot_and_joined_table <- downloadHandler(
    filename = function() {
      paste("rna_vs_dna_results", ".zip", sep = "")
    },
    content = function(fname) {
      
      file_path1 <- paste(input$project_name,
                          "_rna_vs_dna_starburst_plot",
                          ".tiff",
                          sep = "")
      tiff(filename = file_path1,
           height = 10,
           width = 10,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$starburst_plot)
      dev.off()
      
      file_path2 <- paste(input$project_name,
                          "_rna_vs_dna_joined_table",
                          ".csv",
                          sep = "")
      write.csv(rv$joined_rna_dna, 
                file_path2, 
                row.names = FALSE)
      
      fs <- c(file_path1,
              file_path2)
      zip(zipfile = fname, 
          files = fs)
    },
    contentType = "application/zip"
  )
}

shinyApp(ui, server)
