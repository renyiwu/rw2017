## annotate DESeq2 output with gene description



# install.packages("BiocManager")
# BiocManager::install("biomaRt")
library("biomaRt")

dt <- read.table("data/UVB_SKIN/UVB_CON/results/25wksUvbTum-25wksConSki.csv",
                 sep = "\t",
                 row.names = 1,
                 header = T )

dt <- read.table("data/UVB_SKIN/UVB_CON/results/top-50_25wksUvbTum-25wksConSki.csv",
                 sep = ",",
                 row.names = 1,
                 header = T )


mart = useEnsembl('ENSEMBL_MART_ENSEMBL')
mart2 <- listDatasets(mart)
mart2$dataset
# mmusculus_gene_ensembl

ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
# IDs <- c("BRCA2","BRAF")

genedesc <- getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = rownames(dt), mart = ensembl)
genedesc
dt$geneid <- rownames(dt)
genedesc$external_gene_name
m <- merge(dt, genedesc, by.x = "geneid", by.y = "external_gene_name", all.x = T, all.y = F)

nrow(dt)
# 10402

nrow(genedesc)
# 10051

nrow(m)
# 10417

test <- c('.name.1.','name.2','.name.3.')
gsub('^\\.|\\.$', '', test)
# [1] "name.1" "name.2" "name.3"

m$description <- gsub(' \\[.*', '', m$description)


write.table(m, "data/UVB_SKIN/UVB_CON/results/top-50_25wksUvbTum-25wksConSki_description.csv",
            sep = "\t",
            row.names = F)

