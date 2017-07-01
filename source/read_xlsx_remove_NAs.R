require(xlsx)
m <- read.xlsx(file = "data/Methyl-seq results_John_18wks only_Renyi2.xlsx",1)
m
require(data.table)
m <- read.csv(file = "data/Methyl-seq results_John_18wks only_Renyi.csv")
m1 <- na.omit(m)
## Alternative ways:
m2 <- m[complete.cases(m),] #Remove Nas in all columns.
m3 <- m[complete.cases(m[,1:8]),] #Only remove rows with NAs in few columns.

class(m)
#Save as Excel sheet
# write.xlsx(m1,file = "data/temp.xlsx")
#Save as csv file.
write.csv(m1, file = "data/methyl-seq 18weeks Renyi.csv")
