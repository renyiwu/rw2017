dt <- fread("data/RNA_JB6_Oct2017/RW_all_primary.dedup_hisat2_new.csv")
library("data.table")


head(dt)
class(dt)
colnames(dt)
dt2 <- dt[,.(Geneid, Length, RW01.dedup.bam, RW07.dedup.bam)]
length(dt2)/4
head(dt2)
iris
cars
dt <- as.data.table(iris)
head(dt)
dt[, head(.SD, 5), by = Species, .SDcols = c("Sepal.Length", "Petal.Length" )]




class(dt)




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
