#Hierarchical clustering
#Renyi Wu 6/5/2017
library(ggplot2)
require(data.table)
##read file
dt <- fread("data/Methyl-seq-short.csv")
dt
## remove rows with NAs
dt2 <- na.omit(dt)
## Alternative ways:
# dt1 <- dt[complete.cases(dt),] #Take into account all columns
# dt1 <- dt[complete.cases(dt[,1:2]),] #Only remove rows with NAs in fisrt 2 columns.
dt <- na.omit(dt)
## transpose data table, rows <-> columns
dtt <- t(dt)
## caculating distance
ds <- dist(dtt)
##Hierarchical clustering
h <- hclust(ds)
dev.off()
plot(h, xlab = "sample", asp = 2)
# Parameters for plot
plot(h, hang = -1, xlab = "sample", ylab = "distance", cex = 1, lwd = 1, col = "black")
rect.hclust(h, 2)

#Other ways to show full lenth  of y axis.
hc <- h
plot(hc)
plot(hc, ylim = c(0,1))
plot(as.dendrogram(hc))
plot(as.dendrogram(hc), ylim = c(0,20),xlim = c(1,8), xlab = "sample", ylab = "distance")
# Define nodePar
nodePar <- list(lab.cex = 1.2, pch = c(NA, 19), 
                cex = 1.2, col = "black")
nodePar
# Define edgePar
edgePar <- list(col = 1, lwd = 2)
#Convert to dendrogram
hcd <- as.dendrogram(hc)
hcd
# Customized plot; remove labels
plot(hcd, ylim = c(0,15), xlab = NA, ylab = NA, nodePar = nodePar,edgePar = edgePar, horiz = TRUE, asp = 1.5)

# plot(table(rpois(100, 5)), type = "h", col = "red", lwd = 10,
#      main = "rpois(100, lambda = 5)")
# plot(table(rpois(100, 5)))
# Below are reference online
# n3 <- subset(dt3, 
#              substr(dt3$'18ad',1,1) == "0")
# 
# # dt1 <- fread("data/avg_methyl.csv")
# # unique(substr(dt1$feature, 1, 5))
# # 
# # pr <- subset(dt1, 
# #              substr(dt1$feature, 1, 8) == "Promoter")
# 
# x <- matrix(rnorm(100), nrow = 5)
# x
# dist(x)
# dist(x, diag = TRUE)
# dist(x, upper = TRUE)
# m <- as.matrix(dist(x))
# d <- as.dist(m)
# stopifnot(d == dist(x))
# 
# d1 <- dist(iris[1:5, -5])
# d1
# h = hclust(d1)
# plot(h)
# 
# x1 <- t(x)
# x1
