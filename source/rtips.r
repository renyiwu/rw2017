# ONLINE RESOURCES
# Programming with R.  https://swcarpentry.github.io/r-novice-inflammation/
# Data wrangling in R. http://clayford.github.io/dwir.html


# 1, sort a data.frame based on columns. (https://stackoverflow.com/questions/1296646/how-to-sort-a-dataframe-by-columns)
dd <- data.frame(b = factor(c("Hi", "Med", "Hi", "Low"), 
      levels = c("Low", "Med", "Hi"), ordered = TRUE),
      x = c("A", "D", "A", "C"), y = c(8, 3, 9, 9),
      z = c(1, 1, 1, 2))
dd
#     b x y z
# 1  Hi A 8 1
# 2 Med D 3 1
# 3  Hi A 9 1
# 4 Low C 9 2


# Do this:
dd[with(dd, order(-z, b)), ]
# or this, (better)
dd[ order(-dd[,4], dd[,1]), ]

# 2. data.table tips.
library("data.table")
load("~/Downloads/ghData/datasets_L07.Rda")
allStocks
nrow(allStocks)
class(allStocks)
allStocksDT <-data.table(allStocks)
class(allStocksDT)
str(allStocksDT) # structure of x
allStocksDT[1:5][Stock == "bbby"]
allStocksDT[1:5, "Stock"]
allStocksDT[1:5, .(Stock)]
# Finally we can use the by argument to do calculations by group. Here we
# calculate mean and SD of Open by levels of Stock:
allStocksDT[,.(meanOpen = mean(Open), sdOpen = sd(Open)), by = .(Stock)]

#The function is ":=" (read "colon equals"). It 
# updates or adds column(s) by reference. That is, it makes no copies of any
# part of memory at all. This can be very efficient for large data sets.
allStocksDT[, Day := weekdays(Date)]
# Remove the column we created:
allStocksDT[, Day := NULL]


