dt12 <- read.table("data/UVB_SKIN/May/fpkm_to_check_per_group_12-genes.csv",
                   sep = "\t",
                   header = F
                   )

dt12 <- t(dt12)
sink("~/UVB_SKIN/renyi_dedup_methylseq_counts_02032018/grep.sh")
for (i in 1:nrow(dt12)) {
m <- strsplit(dt12[i, 3], "")
lm <- length(m[[1]])
n <- strsplit(dt12[i, 4], "")

l <- lt <- 1
while ( m[[1]][lt] == n[[1]][lt] ) {
  l <- lt
  lt <- lt+1
}
l <- lm - l

str1 <- substring(dt12[i, 3], 1, length(m[[1]]) - l)

a <- data.frame(matrix(NA, nrow = 2*l-1, ncol = l))
for (j in 1:(2*l-1)) {
  for (k in 1:l) {
    if (j == 1) {
      if (k < l) {
        a[j,k] <- m[[1]][lm-l+k]
      } else {
        a[j,k] <- paste0("[", as.numeric(m[[1]][lm-l+k]), "-9]")
      }
    } else if (j < l){
      if (k < l-j+1){
        a[j,k] <- m[[1]][lm-l+k]
      } else if (k == l-j+1) {
        a[j,k] <- paste0("[", as.numeric(m[[1]][lm-l+k])+1, "-9]")
      } else if (k > l-j+1) {
        a[j,k] <- "[0-9]"
      }
    } else if (j == l) {
      if (k == 1){
        if (as.numeric(m[[1]][lm-l+k]) < as.numeric(n[[1]][lm-l+k])-1) {
          a[j,k] <- paste0("[", as.numeric(m[[1]][lm-l+k])+1, "-", as.numeric(n[[1]][lm-l+k])-1, "]")}
        else {
          a[j,k] <- "X"
        }
      } else {
        a[j,k] <- "[0-9]"
      }
    } else if (j == 2*l-1) {
      if (k < l){
        a[j,k] <- n[[1]][lm-l+k]
      } else if (k == l) {
        a[j,k] <- paste0("[0-", as.numeric(n[[1]][lm-l+k]), "]")
      }
    } else if (j > l) {
      if (k < j-l+1){
        a[j,k] <- n[[1]][lm-l+k]
      } else if (k == j-l+1) {
        if (as.numeric(n[[1]][lm-l+k]) > 0) {
          a[j,k] <- paste0("[0-", as.numeric(n[[1]][lm-l+k])-1, "]")
        } else {
          a[j,k] <- paste0("X") # ("[0-", as.numeric(n[[1]][lm-l+k]), "]")
        }
      } else if (k > j-l+1) {
        a[j,k] <- "[0-9]"
      }
    }
  }
}

str2 <- "("
for (r in 1:(nrow(a)-1)) {
  s1 <- paste0(a[r, ], collapse = "")
  str2 <- paste0(str2, s1, "|")
}
sr <- paste0(a[nrow(a), ], collapse = "")

str2 <- paste0(str2, sr, ")")

str0 <- dt12[i, 2]

str <- paste0("grep -P \"", str0, "\\t", str1, str2, '\" ', "r3_816.tsv ", "> 816/", dt12[i, 1], "_r3.tsv")
# strall <- paste0(strall, str, "\n")
cat(str, "\n")
}
sink()

#output
# grep -P "chr14\t46391(57[9-9]|5[8-9][0-9]|[6-8][0-9][0-9]|9[0-3][0-9]|94[0-9])" r3_816.tsv > 816/Bmp4_r3.tsv
# grep -P "chr11\t40754(08[7-9]|0[9-9][0-9]|X[0-9][0-9]|1[0-8][0-9]|19[0-4])" r3_816.tsv > 816/Ccng1_r3.tsv
# grep -P "chr3\t89416(22[2-9]|2[3-9][0-9]|X[0-9][0-9]|3[0-5][0-9]|36[0-2])" r3_816.tsv > 816/Cks1b_r3.tsv
# grep -P "chr2\t154568(10[0-9]|1[1-9][0-9]|[2-2][0-9][0-9]|3[0-2][0-9]|33[0-4])" r3_816.tsv > 816/E2f1_r3.tsv
# grep -P "chr4\t136178(47[6-9]|4[8-9][0-9]|[5-5][0-9][0-9]|6[0-8][0-9]|69[0-2])" r3_816.tsv > 816/E2f2_r3.tsv
# grep -P "chr13\t97242(37[9-9]|3[8-9][0-9]|[4-4][0-9][0-9]|5[0-5][0-9]|56[0-4])" r3_816.tsv > 816/Enc1_r3.tsv
# grep -P "chr3\t122246(53[6-9]|5[4-9][0-9]|X[0-9][0-9]|6[0-4][0-9]|65[0-7])" r3_816.tsv > 816/Gclm_r3.tsv
# grep -P "chr7\t125551(22[3-9]|2[3-9][0-9]|X[0-9][0-9]|3[0-7][0-9]|38[0-3])" r3_816.tsv > 816/Il4ra_r3.tsv
# grep -P "chr3\t5166(264[7-9]|26[5-9][0-9]|2[7-9][0-9][0-9]|X[0-9][0-9][0-9]|3X[0-9][0-9]|30[0-2][0-9]|303[0-9])" r3_816.tsv > 816/Mgst2_r3.tsv
# grep -P "chr9\t99141(72[6-9]|7[3-9][0-9]|[8-8][0-9][0-9]|9[0-4][0-9]|95[0-3])" r3_816.tsv > 816/Pik3cb_r3.tsv
# grep -P "chr4\t14970(262[0-9]|26[3-9][0-9]|2[7-9][0-9][0-9]|X[0-9][0-9][0-9]|3X[0-9][0-9]|30[0-3][0-9]|304[0-7])" r3_816.tsv > 816/Pik3cd_r3.tsv
# grep -P "chr6\t8619(687[4-9]|68[8-9][0-9]|6[9-9][0-9][0-9]|X[0-9][0-9][0-9]|7X[0-9][0-9]|70[0-7][0-9]|708[0-0])" r3_816.tsv > 816/Tgfa_r3.tsv
# grep -P "chr9\t116172(62[7-9]|6[3-9][0-9]|X[0-9][0-9]|7[0-8][0-9]|79[0-8])" r3_816.tsv > 816/Tgfbr2_r3.tsv
# grep -P "chr17\t46033(45[4-9]|4[6-9][0-9]|[5-5][0-9][0-9]|6[0-6][0-9]|67[0-7])" r3_816.tsv > 816/Vegfa_r3.tsv

get_group_ave <- function(infile) {
#infile <- "data/UVB_SKIN/Methyl/816/Enc1_r3h.tsv"
dt.g <- read.table(infile,
                 sep = "\t",
                 header = T)
dt.g.ave <- data.frame(row.names = dt.g$start,
                       week_02_control = round(apply(dt.g[5:6], 1, mean) ,2),
                       week_15_control = round(apply(dt.g[7:8], 1, mean) ,2),
                       week_25_control = round(apply(dt.g[9:10], 1, mean) ,2),
                       week_25_whole_s = round(apply(dt.g[11:12], 1, mean) ,2),
                       week_02_uvb = round(apply(dt.g[13:14], 1, mean) ,2),
                       week_15_uvb = round(apply(dt.g[15:16], 1, mean) ,2),
                       week_25_uvb = round(apply(dt.g[17:18], 1, mean) ,2),
                       week_25_tumor = round(apply(dt.g[19:20], 1, mean) ,2)
                       )
outdir <- dirname(infile)
basename <- basename(infile)
outfile <- paste0(outdir, "/group/", basename)
write.table(t(dt.g.ave), outfile, sep = "\t", col.names = NA, row.names = T)
}

for ( files in Sys.glob(file.path("data/UVB_SKIN/Methyl/816", "*.tsv")) ) {
  get_group_ave(files)
#  get_group_ave("data/UVB_SKIN/Methyl/816/Gclm_r3h.tsv")
}

