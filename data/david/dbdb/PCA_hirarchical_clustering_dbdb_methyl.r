###groups
# 16 db/m	16db/db	21 db/m	21 db/db
# DD1	DD3	DD5	DD7
# DD2	DD4	DD6	DD8
#               DD9	DD10

#1 Hirarchical clustering
library(ggplot2)


dt <- read.table("data/david/dbdb/combined_fr5_DD.csv",
                  header = T, sep = "\t") #, row.names = "Geneid")
dt1 <- na.omit(dt) # Remove rows that have at least one NA.
dt1 <- dt1[-c(1:4, 13:14)] # keep only DD01 ~ DD08
colnames(dt1) <- c("16wks_con_0", "16wks_con_1", "16wks_db_0", "16wks_db_1", "21wks_con_0", "21wks_con_1", "21wks_db_0", "21wks_db_1")
dtt1 <- t(dt1)
ds <- dist(dtt1)
h <- hclust(ds)
plot(h)
plot(h, xlab = "sample", asp = 2, ylab = "distance", main = "db/db model")
# Parameters for plotting
# plot(h, hang = -1, xlab = "sample", ylab = "distance", cex = 1.2, lwd = 2, col = "black")
# #rect.hclust(h, 5)
# dev.off()



#2 PCA analysis
#
library(ggfortify)
dt <- read.table("data/david/dbdb/combined_fr5_DD.csv",
                 header = T, sep = "\t") #, row.names = "Geneid")
dt1 <- na.omit(dt) # Remove rows that have at least one NA.
dt1 <- dt1[-c(1:4, 13:14)] # keep only DD01 ~ DD08
colnames(dt1) <- c("16wks_con_0", "16wks_con_1", "16wks_db_0", "16wks_db_1", "21wks_con_0", "21wks_con_1", "21wks_db_0", "21wks_db_1")
df1 <- data.frame(t(dt1))
df2 <- df1
df2$Sample <- rownames(df2)
df2$Group <- c("16 weeks Control","16 weeks Control", "16 weeks Diabetic", "16 weeks Diabetic", "21 weeks Control", "21 weeks Control", "21 weeks Diabetic", "21 weeks Diabetic")
summary(prcomp(df1))
# Importance of components:
#         PC1    PC2    PC3    PC4    PC5    PC6    PC7       PC8
# Standard deviation     8.5450 7.5686 7.4643 7.4094 7.2836 7.2309 7.1860 9.132e-14
# Proportion of Variance 0.1835 0.1440 0.1400 0.1380 0.1333 0.1314 0.1298 0.000e+00
# Cumulative Proportion  0.1835 0.3275 0.4675 0.6055 0.7388 0.8702 1.0000 1.000e+00
autoplot(prcomp(df1), data = df2, colour = 'Sample')
autoplot(prcomp(df1), data = df2, colour = 'Group')
autoplot(prcomp(df1), data = df2, colour = 'Group', label = TRUE, label.size = 3)
autoplot(prcomp(df1), data = df2, colour = 'Group', shape = FALSE, label.size = 3) #Labels from rownames(dtp2)
autoplot(prcomp(df1), data = df2, colour = 'Group', loadings = TRUE)
dev.off()

 