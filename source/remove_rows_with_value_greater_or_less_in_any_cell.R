

require(raster)
r <- raster(ncols=2, nrows=5)
r <- rbind(c(1,2,3), c(4,5,6), c(7,8,9), c(10,11,12))
values(r) <- 1:10
as.matrix(r)
rowSums(r)
colSums(r)
mr <- as.matrix(r)
mr[!rowSums(mr < 4),]
#
#############################################################
#
#Delete rows in R if a cell contains a value larger than x  #
#
#############################################################
# I want to delete all rows containing a value larger than 7 in a cell in an arbitrary column, either across all columns or across specific columns.
# 
# a <- c(3,6,99,7,8,9)
# b <- c(99,6,3,4,5,6)
# c <- c(2,5,6,7,8,3)
# df <- data.frame (a,b,c)
# 
# a  b c
# 1  3 99 2
# 2  6  6 5
# 3 99  3 6
# 4  7  4 7
# 5  8  5 8
# 6  9  6 3
# V1: I want to delete all rows containing values larger than 7, regardless of the column.
# 
# # result V1
# a  b c
# 2  6  6 5
# 4  7  4 7
# V2: I want to delete all rows containing values larger than 7 in column b and c
# 
# # result V2
# a  b c
# 2  6  6 5
# 3 99  3 6
# 4  7  4 7
# 6  9  6 3
# There are plenty of similar problems on SOF, but I couldn't find a solution to this problem. So far I can only find rows that include 7using res <- df[rowSums(df != 7) < ncol(df), ].

# rowSums of the logical matrix df > 7 gives the number of 'TRUE' per each row. We get '0' if there are no 'TRUE' for that particular row. By negating the results, '0' will change to 'TRUE", and all other values not equal to 0 will be FALSE. This can be used for subsetting.
# 
# df[!rowSums(df >7),]
# #  a b c
# #2 6 6 5
# #4 7 4 7
# For the 'V2', we use the same principle except that we are getting the logical matrix on a subset of 'df'. ie. selecting only the second and third columns.
# 
# df[!rowSums(df[-1] >7),]
# #   a b c
# #2  6 6 5
# #3 99 3 6
# #4  7 4 7
# #6  9 6 3
