#by Davit, 7-10-2017

dt1 <- read.table("data/methyl-results.head.csv",
           sep = "\t",
           header = T)
dt1

ndx <- seq(from = 5, 
    to = ncol(dt1),
    by = 2)
ndx
out <- list()
for (i in 1:length(ndx)) {
  out[[i]] <- dt1[, ndx[i] + 1]/dt1[, ndx[i]]
}
out<- do.call("cbind",
              out)
colnames(out) <- gsub(x = colnames(dt1)[ndx],
                      pattern = ".N",
                      replacement = ".cpt")
out

dt1 <- data.frame(dt1, out)
dt1
write.table(dt1, "methyl-results.head.mu.csv", sep='\t', quote=F, row.names=F)
# can we just append these colunms to the end (right) of the table in the original file? save as a new file...