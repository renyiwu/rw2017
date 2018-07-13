require(data.table)
dt1 <- data.frame(fread("data/yue/combraw.count"))
rownames(dt1) <- dt1$Geneid
dt2 <- dt1
dt2[,1:5] <- NULL
colnames(dt2) <- c("length", "C1", "C2", "C26", "C29", "C42", "C46", "C20", "C14", "C15", "C19", "C34", "C40", "C33", "C36", "C54", "C60", "C55", "C59")
dtt <- t(dt2)

require(cummeRbund)
PCAplot(diffGenes2, replicates=T)

#8-3-2017
dt1 <- read.table("data/John/samples_18_fpkm_all.csv", header = T, sep = "\t", row.names = "X")
write.table(dt1[,1, drop =F], file = "data/John/18wks_fpkm/C20.txt", sep = "\t", quote = F, row.names = T, col.names = F)
#and finish other columns.

dt2 <- read.table("data/combined_CpG_all_John.csv", header = T, sep = "\t")#, row.names = "Geneid")
dt22 <- dt2
dt2[,5:16] <- NULL
write.table(dt2[,c(1:4,5:6)], file = "data/John/18wks_dmr/C14_dmr.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(dt2[,c(1:4,7:8)], file = "data/John/18wks_dmr/C20_dmr.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(dt2[,c(1:4,9:10)], file = "data/John/18wks_dmr/C34_dmr.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(dt2[,c(1:4,11:12)], file = "data/John/18wks_dmr/C40_dmr.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(dt2[,c(1:4,13:14)], file = "data/John/18wks_dmr/C54_dmr.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(dt2[,c(1:4,15:16)], file = "data/John/18wks_dmr/C60_dmr.txt", sep = "\t", quote = F, row.names = F, col.names = T)
