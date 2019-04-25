# data.table practice
# R Wu
# April 2019

# https://cran.r-project.org/web/packages/data.table/vignettes/
#
library("data.table")



#### 1. INtroduction
# https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html

input <- if (file.exists("flights14.csv")) {
  "flights14.csv"
} else {
  "https://raw.githubusercontent.com/Rdatatable/data.table/master/vignettes/flights14.csv"
} 
flights <- fread(input)
flights

ans <- flights[origin == "JFK" & month == 7L]
head(ans)

ans <- flights[1:2]
ans


ans <- flights[order(origin, -dest)]
head(ans)


ans <- flights[, .(arr_delay)]
head(ans)


ans <- flights[, .(delay_arr = arr_delay, delay_dep = dep_delay)]
head(ans)


ans <- flights[, .(total_delay = sum( (arr_delay + dep_delay) < 0 ), arr_delay = sum( arr_delay < 0), dep_delay = sum( dep_delay < 0)), c("month", "day")]
ans


ans <- flights[origin == "JFK" & month == 6L, length(dest)]
ans


ans <- flights[origin == "JFK" & month == 6L, .N]
ans



select_cols = c("arr_delay", "dep_delay")
flights[ , ..select_cols]


ans <- flights[, .(total_flights = .N), by = .(origin)]
ans



ans <- flights[carrier == "AA", .N, .(origin,dest)]
head(ans)


ans <- flights[carrier == "AA",
               .(mean(arr_delay), mean(dep_delay)),
               by = .(origin, dest, month)]
ans



ans <- flights[carrier == "AA",
               .(mean(arr_delay), mean(dep_delay)),
               keyby = .(origin, dest, month)]
ans



ans <- flights[carrier == "AA", .N, by = .(origin, dest)]
ans <- ans[order(origin, -dest)]
head(ans)


ans <- flights[carrier == "AA", .N, by = .(origin, dest)][order(origin, -dest)]
head(ans, 10)


ans <- flights[, .N, .(dep_delay>0, arr_delay>0)]
ans


flights[carrier == "AA", ## Only on trips with carrier "AA"
        lapply(.SD, mean), ## compute the mean
        by = .(origin, dest, month), ## for every 'origin,dest,month'
        .SDcols = c("arr_delay", "dep_delay")] ## for just those specified in .SDcols



DT = data.table(
  ID = c("b","b","b","a","a","c"),
  a = 1:6,
  b = 7:12,
  c = 13:18
)
DT

DT[, .(val = c(a,b)), by = ID]
# 
# ID val
# 1:  b   1
# 2:  b   2
# 3:  b   3
# 4:  b   7
# 5:  b   8
# 6:  b   9
# 7:  a   4
# 8:  a   5
# 9:  a  10
# 10:  a  11
# 11:  c   6
# 12:  c  12



DT[, .(val = list(c(a,b))), by = ID]
# ID val
# 1: b 1,2,3,7,8,9
# 2: a 4, 5,10,11
# 3: c 6,12



DT[, print(c(a,b)), by = ID]
# [1] 1 2 3 7 8 9
# [1] 4 5 10 11
# [1] 6 12
# Empty data.table (0 rows and 1 cols): ID
## (2) and
DT[, print(list(c(a,b))), by = ID]
# [[1]]
# [1] 1 2 3 7 8 9
# 
#[[1]]
# [1] 4 5 10 11
# [[1]]
# [1] 6 12
# Empty data.table (0 rows and 1 cols): ID


 ##### 2. Reference semantics
# https://cran.r-project.org/web/packages/data.table/vignettes/datatable-reference-semantics.html

download.file("https://raw.githubusercontent.com/Rdatatable/data.table/master/vignettes/flights14.csv", "prac/data.table/flights14.csv")

flights <- fread("prac/data.table/flights14.csv")
flights


flights[, `:=`(speed = distance / (air_time/60), # speed in mph (mi/h)
               delay = arr_delay + dep_delay)]   # delay in minutes
head(flights)

## alternatively, using the 'LHS := RHS' form
# flights[, c("speed", "delay") := list(distance/(air_time/60), arr_delay + dep_delay)]



# get all 'hours' in flights
flights[, sort(unique(hour))]
#  [1]  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24


# subassign by reference
flights[hour == 24L, hour := 0L][]   # L indicates integer (instead of double float)

flights[, c("delay") := NULL]
head(flights)

flights[, speed := NULL]



flights[, max_speed := max(speed), by = .(origin, dest)]
head(flights)



in_cols  = c("dep_delay", "arr_delay")
out_cols = c("max_dep_delay", "max_arr_delay")
flights[, c(out_cols) := lapply(.SD, max), by = month, .SDcols = in_cols]
head(flights)

flights[, c("speed", "max_speed", "max_dep_delay", "max_arr_delay") := NULL]
head(flights)


####     ## 3. Keys and fast binary search based subset#################################################### 3.
# https://cran.r-project.org/web/packages/data.table/vignettes/datatable-keys-fast-subset.html

flights <- fread("prac/data.table/flights14.csv")


setkey(flights, origin)
head(flights)


flights[.("JFK")]


## alternatively
# flights[J("JFK")] (or) 
# flights[list("JFK")]


flights[.("JFK", "LGA")]  # wrong. "LGA" was ignored.

flights[c("JFK", "LGA")]    ## same as flights[.(c("JFK", "LGA"))]

key(flights)

# remove keys
setkey(flights, NULL)


setkey(flights, origin, dest)
head(flights)


flights[.("JFK", "MIA")]


key(flights)
# [1] "origin" "dest"


flights[.("JFK")] ## or in this case simply flights["JFK"], for convenience

flights[.(unique(origin), "MIA")] 
# Here 3 elements were provided to the first key matching, but only 1 ("MIA") was provided for
# the second key matching. Here "MIA" is automatically recycled to fit the length of unique(origin), which is 3.
# Below works the same  
# flights[.(unique(origin), rep("MIA",3))]



key(flights)
# [1] "origin" "dest"

flights[.("LGA", "TPA"), .(arr_delay)]
flights[.("LGA", "TPA"), "arr_delay", with = FALSE]

flights[.("LGA", "TPA"), .(arr_delay)][order(-arr_delay)]

flights[.("LGA", "TPA"), max(arr_delay)]


flights[, sort(unique(hour))]


setkey(flights, hour)
key(flights)
# [1] "hour"
flights[.(24), hour := 0L]
key(flights)
# "hour"


setkey(flights, origin, dest)
key(flights)
# [1] "origin" "dest"

ans <- flights["JFK", max(dep_delay), keyby = month]

head(ans)

key(ans)


flights[.("JFK", "MIA"), mult = "first"]
#    year month day dep_delay arr_delay carrier origin dest air_time distance hour
# 1: 2014     1   1         6         3      AA    JFK  MIA      157     1089    5
# Good for selecting only one DMR for each gene. # methyl-seq 

flights[.(c("LGA", "JFK", "EWR"), "XNA"), mult = "last"]

flights[.(c("LGA", "JFK", "EWR"), "XNA"), mult = "last", nomatch = NULL]


# key by origin,dest columns
flights[.("JFK", "MIA")]
 
  
  flights[origin == "JFK" & dest == "MIA"]
  
  set.seed(2L)
  N = 9e7L
  DT = data.table(x = sample(letters, N, TRUE),
                  y = sample(1000L, N, TRUE),
                  val = runif(N))
  print(object.size(DT), units = "Mb")
  
  key(DT)

  t1 <- system.time(ans1 <- DT[x == "g" & y == 877L])
  t1
  
  
  head(ans1)
  
  dim(ans1)
  
  
  
  setkeyv(DT, c("x", "y"))
  key(DT)

  t2 <- system.time(ans2 <- DT[.("g", 877L)])
  t2 
  
  dim(ans2)
  
  identical(ans1$val, ans2$val)
  
  
  ############## 4. Secondary indices and auto indexing #################################################  4
  # https://cran.r-project.org/web/packages/data.table/vignettes/datatable-secondary-indices-and-auto-indexing.html
  
  flights <- fread("prac/data.table/flights14.csv")
  
  head(flights)
  dim(flights)
  
  names(attributes(flights))
  
  setindex(flights, origin)
  
  indices(flights)
  
  setindex(flights, origin, dest)
  indices(flights)
  
  
  
  # ############ 5. Efficient reshaping using data.tables   ############################################## 5
  # https://cran.r-project.org/web/packages/data.table/vignettes/datatable-reshape.html
  
  s1 <- "family_id age_mother dob_child1 dob_child2 dob_child3
1         30 1998-11-26 2000-01-29         NA
2         27 1996-06-22         NA         NA
3         26 2002-07-11 2004-04-05 2007-09-02
4         32 2004-10-10 2009-08-27 2012-07-21
5         29 2000-12-05 2005-02-28         NA"
  DT <- fread(s1)
DT  
str(DT)


DT.m1 = melt(DT, id.vars = c("family_id", "age_mother"),
             measure.vars = c("dob_child1", "dob_child2", "dob_child3"))
DT.m1
str(DT.m1)

DT.m1 = melt(DT, measure.vars = c("dob_child1", "dob_child2", "dob_child3"),
             variable.name = "child", value.name = "dob")
DT.m1



dcast(DT.m1, family_id + age_mother ~ child, value.var = "dob")

DT
str(DT)
# Classes ‘data.table’ and 'data.frame':	5 obs. of  5 variables:

class(DT)
# [1] "data.table" "data.frame"

# Below expressions return a vector
DT$family_id
DT[, family_id]
DT[["family_id"]]
DT[[family_id]] # wrong. variable family_id does not exist
# [1] 1 2 3 4 5

# Below expressions return a data table
DT[, "family_id"]
DT[, .(family_id)]
#    family_id
# 1:         1
# 2:         2
# 3:         3
# 4:         4
# 5:         5



DF = data.frame(x = 1:3, y = 4:6, z = 7:9)
DF
#   x y z
# 1 1 4 7
# 2 2 5 8
# 3 3 6 9
DF[ , c("y", "z")]
#   y z
# 1 4 7
# 2 5 8
# 3 6 9

DT = data.table(DF)
DT[ , c(y, z)]  # j is elevated first then returns a vector of 6 elements, 3 from each
# [1] 4 5 6 7 8 9
DT[, .(y, z)]  # j is a list of two elements, each one is a vector and has 3 elements.

l1 <- list(c = 1:10, d = sample(letters, 10))
l1
str(l1)

l2 <- list(1:10, sample(letters, 10))
l2
str(l2)

l3 <- list(1:3, l1, l2)
l4 <- list(l1, l3)
l4
