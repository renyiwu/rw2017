require(data.table)
dt1 <- fread("data/canonical pathways comparison in DSS 8wks top hits.csv")
dt1
# Table1---change parameters as needed-
t1 <- dt1[1:34, c(1, 3)]
tmp <- strsplit(t1$Molecules, split = ",")
tmp

out <- list()
for (i in 1:length(tmp)) {
  out[[i]] <- data.table(pathway = rep(t1$`Ingenuity Canonical Pathways`[i], 
                                       length(tmp[[i]])),
                         gene = tmp[[i]])
}
out

tt1 <- do.call("rbind", out)
tt1
# dir.create(path="data/temp")
# write.csv(tt1,
#           file = "data/temp/8wk_aomdss_pathways.csv",
#           row.names = FALSE)
##################STOP HERE#####
# Table 2----
# t2 <- dt1[1:23, c(7, 11)]
# tmp <- strsplit(t2$Molecules, split = ",")
# tmp
# 
# out <- list()
# for (i in 1:length(tmp)) {
#   out[[i]] <- data.table(pathway = rep(t2$V7[i], 
#                                        length(tmp[[i]])),
#                          gene = tmp[[i]])
# }
# out
# 
# tt2 <- do.call("rbind", out)
# tt2

# Table3
# t3 <- data.table(gene = dt1$V13,
                 # found = TRUE)
#gene name, change lenth as neede)
t3 <- data.table(gene = dt1[1:26,gene],
                 found = TRUE)
t3

# Combine Table 1 and Table 3----
tt1.match <- subset(tt1,
                    tt1$gene %in% t3$gene)
tt1.match

length(unique(tt1$gene))
length(unique(t3$gene))

t3$gene[(t3$gene %in% tt1$gene)]

tt1.merged <- unique(merge(tt1, t3, by = "gene", all.y = TRUE))

setkey(tt1.merged, gene)
tt1.merged[, n := 1:.N,
           by = gene]

dt1.out <- dcast.data.table(tt1.merged,
                            gene ~ n,
                            value.var = "pathway")
dt1.out

# # Combine Table 2 and Table 3----
# length(unique(tt2$gene))
# length(unique(t3$gene))
# 
# t3$gene[(t3$gene %in% tt2$gene)]
# 
# tt2.merged <- unique(merge(tt2, t3, by = "gene", all.y = TRUE))
# 
# setkey(tt2.merged, gene)
# tt2.merged[, n := 1:.N,
#            by = gene]
# 
# dt2.out <- dcast.data.table(tt2.merged,
#                             gene ~ n,
#                             value.var = "pathway")
# dt2.out
# 
# # Save as CSV files----Change file name as needes---
write.csv(dt1.out,
          file = "data/temp/8 weeks dss genes and pathways.csv",
          row.names = FALSE)
# 
# write.csv(dt2.out,
#           file = "tmp/dt2.out.csv",
#           row.names = FALSE)
