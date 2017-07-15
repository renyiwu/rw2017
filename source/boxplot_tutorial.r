## boxplot on a formula:
boxplot(count ~ spray, data = InsectSprays, col = "lightgray")
# *add* notches (somewhat funny here):
boxplot(count ~ spray, data = InsectSprays,
        notch = TRUE, add = TRUE, col = "blue")

boxplot(decrease ~ treatment, data = OrchardSprays,
        log = "y", col = "bisque")

rb <- boxplot(decrease ~ treatment, data = OrchardSprays, col = "bisque")
title("Comparing boxplot()s and non-robust mean +/- SD")

mn.t <- tapply(OrchardSprays$decrease, OrchardSprays$treatment, mean)
sd.t <- tapply(OrchardSprays$decrease, OrchardSprays$treatment, sd)
xi <- 0.3 + seq(rb$n)
points(xi, mn.t, col = "orange", pch = 18)
arrows(xi, mn.t - sd.t, xi, mn.t + sd.t,
       code = 3, col = "pink", angle = 75, length = .1)

## boxplot on a matrix:
mat <- cbind(Uni05 = (1:100)/21, Norm = rnorm(100),
             `5T` = rt(100, df = 5), Gam2 = rgamma(100, shape = 2))
boxplot(mat) # directly, calling boxplot.matrix()

## boxplot on a data frame:
df. <- as.data.frame(mat)
par(las = 1) # all axis labels horizontal
boxplot(df., main = "boxplot(*, horizontal = TRUE)", horizontal = TRUE)

## Using 'at = ' and adding boxplots -- example idea by Roger Bivand :
boxplot(len ~ dose, data = ToothGrowth,
        boxwex = 0.25, at = 1:3 - 0.2,
        subset = supp == "VC", col = "yellow",
        main = "Guinea Pigs' Tooth Growth",
        xlab = "Vitamin C dose mg",
        ylab = "tooth length",
        xlim = c(0.5, 3.5), ylim = c(0, 35), yaxs = "i")
boxplot(len ~ dose, data = ToothGrowth, add = TRUE,
        boxwex = 0.25, at = 1:3 + 0.2,
        subset = supp == "OJ", col = "orange")
legend(2, 9, c("Ascorbic acid", "Orange juice"),
       fill = c("yellow", "orange"))

## With less effort (slightly different) using factor *interaction*:
boxplot(len ~ dose:supp, data = ToothGrowth,
        boxwex = 0.5, col = c("orange", "yellow"),
        main = "Guinea Pigs' Tooth Growth",
        xlab = "Vitamin C dose mg", ylab = "tooth length",
        sep = ":", lex.order = TRUE, ylim = c(0, 35), yaxs = "i")

## more examples in  help(bxp)


###demo 2 Box-plot with R – Tutorial
#
#Yesterday I wanted to create a box-plot for a small dataset 
# to see the evolution of 3 stations through a 3 days period. 
# I like box-plots very much because I think they are one of 
# the clearest ways of showing trend in your data. 
# R is extremely good for this type of plot and, 
# for this reason, I decided to add a post on my blog 
# to show how to create a box-plot, but also 
# because I want to use my own blog to 
# help me remember pieces of code 
# that I might want to use in the future but that I tend to forget.
# For this example I first created a dummy dataset using the 
# function rnorm() which generates random normal-distributed sequences.
# This function requires 3 arguments, 
# the number of samples to create, t
# he mean and the standard deviation of the distribution, for example: 
rnorm(n=100,mean=3,sd=1)
# This generates 100 numbers (floats to be exact), 
# which have mean equal to 3 and standard deviation equal to 1.
# To generate my dataset I used the following line of code:
data<-data.frame(Stat11=rnorm(100,mean=3,sd=2),
Stat21=rnorm(100,mean=4,sd=1),
Stat31=rnorm(100,mean=6,sd=0.5),
Stat41=rnorm(100,mean=10,sd=0.5),
Stat12=rnorm(100,mean=4,sd=2),
Stat22=rnorm(100,mean=4.5,sd=2),
Stat32=rnorm(100,mean=7,sd=0.5),
Stat42=rnorm(100,mean=8,sd=3),
Stat13=rnorm(100,mean=6,sd=0.5),
Stat23=rnorm(100,mean=5,sd=3),
Stat33=rnorm(100,mean=8,sd=0.2),
Stat43=rnorm(100,mean=4,sd=4))
 
# As I mentioned before, this should represent 4 stations
# for which the measure were replicated in 3 successive days.
# # Now, for the creation of the box-plot the simplest function 
# is boxplot() and can be simply called by adding the 
# name of the dataset as only argument:
boxplot(data)
boxplot(data, las = 2)
boxplot(data, las = 2,
        names = c("Station 1","Station 2","Station 3","Station 4",
                   "Station 1","Station 2","Station 3","Station 4",
                   "Station 1","Station 2","Station 3","Station 4"))
#
#If the names are too long and they do not fit into the plot’s window 
#you can increase it by using the option par:
boxplot(data, las = 2,
        par(mar = c(12, 5, 4, 2)+ 0.1),
        names = c("Station 1","Station 2","Station 3","Station 4",
                   "Station 1","Station 2","Station 3","Station 4",
                   "Station 1","Station 2","Station 3","Station 4"))
#
#Now I want to group the 4 stations so that the division in 3 successive days is clearer.
boxplot(data, las = 2,at =c(1,2,3,4, 6,7,8,9, 11,12,13,14),
           par(mar = c(12, 5, 4, 2)+ 0.1),
        names = c("Station 1","Station 2","Station 3","Station 4",
                   "Station 1","Station 2","Station 3","Station 4",
                   "Station 1","Station 2","Station 3","Station 4"))
#
#Add colors
boxplot(data, las = 2,at =c(1,2,3,4, 6,7,8,9, 11,12,13,14),
           par(mar = c(12, 5, 4, 2)+ 0.1),
     col = c("red","sienna","palevioletred1","royalblue2","red","sienna","palevioletred1", 
"royalblue2","red","sienna","palevioletred1","royalblue2"),
        names = c("Station 1","Station 2","Station 3","Station 4",
                   "Station 1","Station 2","Station 3","Station 4",
                   "Station 1","Station 2","Station 3","Station 4"))
#
#Add labels
boxplot(data, las = 2,ylab = "y abc", xlab = "hahah",
        at =c(1,2,3,4, 6,7,8,9, 11,12,13,14),
           par(mar = c(12, 5, 4, 2)+ 0.1),
     col = c("red","sienna","palevioletred1","royalblue2","red","sienna","palevioletred1", 
"royalblue2","red","sienna","palevioletred1","royalblue2"),
        names = c("Station 1","Station 2","Station 3","Station 4",
                   "Station 1","Station 2","Station 3","Station 4",
                   "Station 1","Station 2","Station 3","Station 4"))
#
#Add labels 2
boxplot(data, las = 2,ylab = "y abc", xlab = "hahah",
        at =c(1,2,3,4, 6,7,8,9, 11,12,13,14),
           par(mar = c(12, 5, 4, 2)+ 0.1),
     col = c("red","sienna","palevioletred1","royalblue2","red","sienna","palevioletred1", 
"royalblue2","red","sienna","palevioletred1","royalblue2"),
        names = c("Station 1","Station 2","Station 3","Station 4",
                   "Station 1","Station 2","Station 3","Station 4",
                   "Station 1","Station 2","Station 3","Station 4"))
#Add labels 2
mtext("x label here", side = 1, line = 5)
#Add labels 2
mtext("y label here", side = 2, line = 2)
