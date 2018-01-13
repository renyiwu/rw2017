#Mann Whitney U test, aka rank sum test
#R Wu. 2017-9-5

require(data.table)
dt1 <- fread("data/methyl-seq 18weeks Renyi.csv")
dt1 <- na.omit(dt1)
# dt2 <- dt1[dt1$CpG>=12 &
#              dt1$Control_18.mu<0.75 &
#              dt1$AOMDSS_18.mu<0.75 &
#              dt1$AOMDSSCur_18.mu<0.75,]
dt2 <- dt1[dt1$CpG>=12,] #46015 DMRs remaining
dt21 <- (dt2[,6:8])

write.table(dt1[,1:8], file = "data/methy-seq 18 weeks-renyi-mu.csv", sep = ",", quote = F, row.names = F)


#################
##
#####################
#
####
dt3 <- fread("data/methy-seq 18 weeks-renyi-mu.csv")
dt13 <- dt3[dt3$CpG>=11,6:8]
wilcox.test(dt13$Control_18.mu - dt13$AOMDSS_18.mu, conf.int = TRUE) #or use this. CpG>=11, P=0.4652
wilcox.test(dt13$AOMDSSCur_18.mu - dt13$AOMDSS_18.mu, conf.int = TRUE) #or use this. CpG>=9, P=0.1924


# > wilcox.test(dt13$AOMDSSCur_18.mu - dt13$AOMDSS_18.mu, conf.int = TRUE) #or use this
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  dt13$AOMDSSCur_18.mu - dt13$AOMDSS_18.mu
# V = 568970000, p-value = 0.1924
# alternative hypothesis: true location is not equal to 0
# 95 percent confidence interval:
#   -1.262519e-05  6.628249e-06
# sample estimates:
#   (pseudo)median 
# 1.566505e-05 
# 
# > dt13 <- dt3[dt3$CpG>=11,6:8]
# > wilcox.test(dt13$Control_18.mu - dt13$AOMDSS_18.mu, conf.int = TRUE) #or use this.
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  dt13$Control_18.mu - dt13$AOMDSS_18.mu
# V = 321650000, p-value = 0.4652
# alternative hypothesis: true location is not equal to 0
# 95 percent confidence interval:
#   -2.186580e-05  4.389699e-05
# sample estimates:
#   (pseudo)median 
# 2.633525e-05 

dt211 <-data.table((format(round(dt21,2),nsmall=2, justify = "none")))
dt211$Control_18.mu <- as.numeric(dt211$Control_18.mu)
dt211$AOMDSS_18.mu <- as.numeric(dt211$AOMDSS_18.mu)
dt211$AOMDSSCur_18.mu <- as.numeric(dt211$AOMDSSCur_18.mu)
                                  
wilcox.test(dt211$Control_18.mu, dt211$AOMDSS_18.mu,paired = T, conf.int = TRUE)
class(dt211)
class(dt211$Control_18.mu)

dt3 <- dt211[!rowSums(dt211 >= 0.1)>=3,] #34504 remaining
dt33 <- dt211[rowSums(dt211 < 0.1)==0,] #11511 remaining, all values greater or equal to 0.1
#Run test
wilcox.test(dt3$Control_18.mu, dt3$AOMDSS_18.mu, paired = T, conf.int = TRUE) #use this
wilcox.test(dt13$Control_18.mu - dt13$AOMDSS_18.mu, conf.int = TRUE) #or use this

#
wilcox.test(dt13$AOMDSSCur_18.mu - dt13$AOMDSS_18.mu, conf.int = TRUE) #or use this
wilcox.test(dt3$Control_18.mu - dt3$AOMDSSCur_18.mu, conf.int = TRUE) #or use this

#Or melt data first
##Mann Whitney U test, aka Wilcoxon test.
#Melt data. Pick only control and treatment groups.
dt3m12 <- melt.data.table(dt3[,1:2],
                         variable.name = "Treatment",
                         value.name = "Methyl.ratio")
dt3m23 <- melt.data.table(dt3[,2:3],
                        variable.name = "Treatment",
                        value.name = "Methyl.ratio")
#Transform group values to factors. Seems optional??
dt3m12$Treatment <- factor(dt3m12$Treatment,
                         levels = c("AOMDSS_18.mu",
                                    "Control_18.mu"))

dt3m23$Treatment <- factor(dt3m23$Treatment,
                           levels = c("AOMDSSCur_18.mu",
                                      "AOMDSS_18.mu"))
dt3m23$Treatment <- factor(dt3m23$Treatment,
                           levels = c("AOMDSS_18.mu",
                                      "AOMDSSCur_18.mu"))
#Run test
wilcox.test(Methyl.ratio ~ Treatment, data = dt3m12, paired = TRUE)
wilcox.test(Methyl.ratio ~ Treatment, data = dt3m23, paired = TRUE)

wilcox.test(rnorm(500), rnorm(500), conf.int = TRUE)
r1 <-rnorm(5)
r2 <-rnorm(5)
wilcox.test(r1,r2, conf.int = TRUE)
##End mann Whitney test
plot(rnorm(100,100,1))
##
#t test. Just for example.
#Data here are not normal distributed (bell shap) so should not be used for t test before any kind of transformations. (e.g. log transform)
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
t1p <- t.test(dt4$Control_18.mu,dt3$AOMDSS_18.mu, paired = T)
# data:  dt3$Control_18.mu and dt3$AOMDSS_18.mu
# t = 3.1617, df = 46014, p-value = 0.00157
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.0001608044 0.0006853819
# sample estimates:
# mean of the differences 
#            0.0004230932
t2p <- t.test(dt4$AOMDSSCur_18.mu,dt3$AOMDSS_18.mu, paired = T)
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
