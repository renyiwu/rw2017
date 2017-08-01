#Draw boxplot
#R Wu. 2017-7-15
require(data.table)
dt1 <- fread("data/methyl-seq 18weeks Renyi.csv")
dt1 <- na.omit(dt1)
# dt2 <- dt1[dt1$CpG>=12 &
#              dt1$Control_18.mu<0.75 &
#              dt1$AOMDSS_18.mu<0.75 &
#              dt1$AOMDSSCur_18.mu<0.75,]
dt2 <- dt1[dt1$CpG>=12,]
dt3 <- dt2[,6:8]
#backup parameters
#par.b <- par()
#restore par.
#par(par.b)
boxplot(dt3,border = c("black"),
        at = c(1,3,5), # To make a gaps between groups.This is optional.
        names = c("Control", "AOM+DSS", "AOM+DSS+Cur"),
        col = c("red", "green", "blue"),
        ylim = c(0,1.3),par(lwd = 2, cex.lab = 1.5, cex.axis = 1.4),
        ylab = "Methylation ratio")
pas <- par()  #back current setting
tiff(filename = "data/boxplot/18wks.tiff",
    width = 400, height = 300)
png(filename = "data/boxplot/18wks.png",
    width = 400, height = 300, bg="transparent")
dev.off()

#
t1 <- t.test(dt3$Control_18.mu,dt3$AOMDSS_18.mu)
t2 <- t.test(dt3$AOMDSSCur_18.mu,dt3$AOMDSS_18.mu)
t1
# > t1
# 
# 	Welch Two Sample t-test
# 
# data:  dt3$Control_18.mu and dt3$AOMDSS_18.mu
# t = 0.22725, df = 91995, p-value = 0.8202
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.003226004  0.004072190
# sample estimates:
# mean of x mean of y 
# 0.1635548 0.1631317 
t2
# > t2
# 
# 	Welch Two Sample t-test
# 
# data:  dt3$AOMDSSCur_18.mu and dt3$AOMDSS_18.mu
# t = -0.20378, df = 92013, p-value = 0.8385
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.004016106  0.003259652
# sample estimates:
# mean of x mean of y 
# 0.1627534 0.1631317 
t1p <- t.test(dt3$Control_18.mu,dt3$AOMDSS_18.mu, paired = T)
# data:  dt3$Control_18.mu and dt3$AOMDSS_18.mu
# t = 3.1617, df = 46014, p-value = 0.00157
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.0001608044 0.0006853819
# sample estimates:
# mean of the differences 
#            0.0004230932
t2p <- t.test(dt3$AOMDSSCur_18.mu,dt3$AOMDSS_18.mu, paired = T)
# > t2p
# 
# 	Paired t-test
# 
# data:  dt3$AOMDSSCur_18.mu and dt3$AOMDSS_18.mu
# t = -3.0407, df = 46014, p-value = 0.002362
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.0006220337 -0.0001344204
# sample estimates:
# mean of the differences 
#            -0.000378227 
