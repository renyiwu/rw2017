#Draw boxplot
#R Wu. 2017-7-15
#Davit 2017-9-8
require(data.table)
dt1 <- fread("data/methyl-seq 18weeks Renyi.csv")
dt1 <- na.omit(dt1)
# dt2 <- dt1[dt1$CpG>=12 &
#              dt1$Control_18.mu<0.75 &
#              dt1$AOMDSS_18.mu<0.75 &
#              dt1$AOMDSSCur_18.mu<0.75,]
dt2 <- dt1[dt1$CpG>=12,] #46015 DMRs remaining
dt3 <- dt2[,6:8]

dt4 <- log2(dt3)

#backup parameters
#par.b <- par()
#restore par.
#par(par.b)
boxplot(dt3,
        log = "y")
boxplot(dt4,
        border = c("black"),
        at = c(1,3,5), # To make a gaps between groups.This is optional.
        names = c("Control", "AOM+DSS", "AOM+DSS+Cur"),
        col = c("red", "green", "blue"),
        # ylim = c(0,1.3),
        log = "y",
        par(lwd = 2, 
            cex.lab = 1.5,
            cex.axis = 1.4),
        ylab = "Methylation ratio")
pas <- par()  #back current setting
tiff(filename = "data/boxplot/18wks.tiff",
    width = 400, height = 300)
png(filename = "data/boxplot/18wks.png",
    width = 400, height = 300, bg="transparent")
dev.off()

dt4.1 <- melt.data.table(dt4,
                        variable.name = "Treatment",
                        value.name = "log2.meth.rat")
dt4.1$Treatment <- factor(dt4.1$Treatment,
                          levels = c("AOMDSS_18.mu",
                                     "Control_18.mu",
                                     "AOMDSSCur_18.mu"))

m1 <- lm(log2.meth.rat ~ Treatment,
         data = dt4.1)
m2 <- aov(log2.meth.rat ~ Treatment,
         data = dt4.1)
summary(m1)
summary(m2)
# > summary(m1)
# 
# Call:
#   lm(formula = log2.meth.rat ~ Treatment, data = dt4.1)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.4523 -0.5244 -0.2470  0.4526  1.3978 
# 
# Coefficients:
#   Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)              -1.381092   0.003292 -419.518  < 2e-16 ***
#   TreatmentControl_18.mu   -0.018142   0.004656   -3.897 9.75e-05 ***
#   TreatmentAOMDSSCur_18.mu -0.017892   0.004656   -3.843 0.000122 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7062 on 138042 degrees of freedom
# Multiple R-squared:  0.0001447,	Adjusted R-squared:  0.0001302 
# F-statistic: 9.986 on 2 and 138042 DF,  p-value: 4.61e-05
anova(m1)

require(multcomp)
m1.1 <- glht(m1,
             linfct = mcp(Treatment = "Dunnett"))
summary(m1.1)
# > summary(m1.1)
# 
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Dunnett Contrasts
# 
# 
# Fit: lm(formula = log2.meth.rat ~ Treatment, data = dt4.1)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)    
# Control_18.mu - AOMDSS_18.mu == 0   -0.018142   0.004656  -3.897 0.000193 ***
#   AOMDSSCur_18.mu - AOMDSS_18.mu == 0 -0.017892   0.004656  -3.843 0.000241 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)
#

# Boxplot
require(ggplot2)
p1 <- ggplot(data = dt4.1) +
  geom_boxplot(aes(x = Treatment,
                   y = log2.meth.rat,
                   fill = Treatment)) +
  scale_x_discrete("Treatment") + 
  scale_y_continuous("Log2(Methylation Ratio)") + 
  ggtitle("Change in Methylation") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")
print(p1)











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
