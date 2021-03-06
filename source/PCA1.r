#PCA (principal component analysis)
#R Wu. 7-16-2017
# install.packages("ggfortify")
# install.packages("tidyr")
library(ggfortify)
df <- iris[c(1, 2, 3, 4)]
df1 <- iris[,1:4]
df
iris
class(iris)
df1
class(df1)
prcomp(df)
autoplot(prcomp(df))
autoplot(prcomp(df), data = iris, colour = 'Species')
autoplot(prcomp(df), data = iris, colour = 'Species', label = TRUE, label.size = 3)
autoplot(prcomp(df), data = iris, colour = 'Species', shape = FALSE, label.size = 3)
autoplot(prcomp(df), data = iris, colour = 'Species', loadings = TRUE)
autoplot(prcomp(df), data = iris, colour = 'Species',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)
autoplot(prcomp(df), scale = 0)
iris
#



#start my own project here
#This is just an example. Not real study design.
require(data.table)
dt1 <- fread("data/combined_CpG_all_John.csv")
dt1 <- fread("data/david/mes13_featurecounts_122017.csv")
dt1 <- na.omit(dt1)
dt2 <- dt1[, c(1, 6:13)]

# class(dt5)
ncol(dt5)
colnames(dt2) <- C(geneid, length, paste("DMR",1:123098, sep = "")
#or

colnames(dt5)[1:ncol(dt5)] <- paste("DMR",1:ncol(dt5), sep = "")
dt6 <- data.frame(dt5)
dt7 <- dt6
dt7$group <- c(rep("group1", 4), rep("group2", 4), rep("group3", 4))
dt7$sample <- rownames(dt7)
prcomp(dt6)
autoplot(prcomp(dt6), data = dt7, colour = 'sample')
#or
autoplot(prcomp(dt6), data = dt7, colour = 'group')

# dt2 <- dt1[dt1$CpG>=12 &
#              dt1$Control_18.mu<0.75 &
#              dt1$AOMDSS_18.mu<0.75 &
#              dt1$AOMDSSCur_18.mu<0.75,]
dt2 <- dt1[dt1$CpG>=12,]
dt3 <- dt2[,6:8]
dt3

