setwd("~/Documents/klaus/mouse_human_liver_tumor_comparison_060214")

## read in raw microarray data and aggregate
raw <- read.delim("OneChannel_Allgenes_res.txt", as.is = T, fill = T)
raw <- raw[, c(3, 24:35)]
exprs <- aggregate(data.matrix(raw[, -1]) ~ raw$EntrezGene, FUN = "median")
rownames(exprs) <- exprs[, 1]
exprs <- data.matrix(exprs[, -1])

## write out processed raw data
write.csv(exprs, file = "mouse_exprs.csv")
