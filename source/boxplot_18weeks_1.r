#Draw boxplot
#R Wu. 2017-7-15
require(data.table)
dt1 <- fread("data/methyl-seq 18weeks Renyi.csv")
dt1 <- na.omit(dt1)
# dt2 <- dt1[dt1$CpG>=12 &
#              dt1$Control_18.mu<0.75 &
#              dt1$AOMDSS_18.mu<0.75 &
#              dt1$AOMDSSCur_18.mu<0.75,]
dt2 <- dt1[dt1$CpG>=12,] #46015 DMRs remaining
dt3 <- dt2[,6:8]

dt4 <- log10(dt3)

#backup parameters
#par.b <- par()
#restore par.
#par(par.b)
boxplot(dt4,
        border = c("black"),
        at = c(1,3,5), # To make a gaps between groups.This is optional.
        names = c("Control", "AOM+DSS", "AOM+DSS+Cur."),
        col = c("red", "green", "blue"),
        ylim = c(-3.5,0.5),par(lwd = 2, cex.lab = 1.5, cex.axis = 1.4),
        ylab = expression(log[10]~Methylation~ratio))
pas <- par()  #back current setting
tiff(filename = "data/boxplot/18wks_log10.tiff",
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
summary(m1)
anova(m1)

require(multcomp)
m1.1 <- glht(m1,
             linfct = mcp(Treatment = "Dunnett"))
summary(m1.1)

# Boxplot
require(ggplot2)
p1 <- ggplot(data = dt4.1) +
  geom_boxplot(aes(x = Treatment,
                   y = log2.meth.rat)) +
  geom_point(aes(x = Treatment,
                 y = log2.meth.rat),
             size = 1,
             alpha = 0.6,
             position = position_dodge(0.3)) + 
  scale_x_discrete("Treatment") + 
  scale_y_continuous("Readout") + 
  ggtitle("Boxplot") +
  facet_wrap(~ read,
             nrow = 1) +


  geom_line(aes(x = trt,
                y = readout,
                group = id,
                colour = id),
            size = 2,
            alpha = 0.6,
            position = position_dodge(0.3)) + 
  guides(colour = guide_legend(title = "ID",
                               title.position = "top",
                               ncol = 1)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "left",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
print(p1)
aov()











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

t1p <- t.test(dt4$Control_18.mu,dt4$AOMDSS_18.mu, paired = T)
# > t1p
# 
# Paired t-test
# 
# data:  dt4$Control_18.mu and dt4$AOMDSS_18.mu
# t = -15.993, df = 46014, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.02036552 -0.01591876
# sample estimates:
#   mean of the differences 
# -0.01814214 

t2p <- t.test(dt4$AOMDSSCur_18.mu,dt4$AOMDSS_18.mu, paired = T)

# > t2p
# 
# Paired t-test
# 
# data:  dt4$AOMDSSCur_18.mu and dt4$AOMDSS_18.mu
# t = -15.797, df = 46014, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.02011219 -0.01567226
# sample estimates:
#   mean of the differences 
# -0.01789223 