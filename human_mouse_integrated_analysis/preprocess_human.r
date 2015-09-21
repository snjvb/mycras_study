setwd("~/Documents/klaus/mouse_human_liver_tumor_comparison_060214")

## Dataset sourced from: 
## http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22058

## Main publication:
## http://www.ncbi.nlm.nih.gov/pubmed/20739924

## read in raw microarray intensities and format column names to 
## reflect tissue origin
raw <- read.delim("GSE22058-GPL6793_series_matrix.txt", as.is = T, fill = T, 
    row.names = 1)
sample.desc <- apply(raw, 2, function(col) {
    temp <- strsplit(col[1], split = " ")[[1]]
    gsub("-", ".", temp[length(temp)])
})
colnames(raw) <- paste(colnames(raw), sample.desc, sep = ".")
raw <- raw[-1, ]

## read in microarray meta data
meta.data <- read.delim("GPL6793.txt", header = F, as.is = T, fill = T, 
    row.names = 1)[, 1:2]
names(meta.data) <- c("entrez.id", "gene.symbol")
meta.data <- meta.data[rownames(raw), ]

## aggregate microarray intensities based on entrez.id
raw2 <- aggregate(data.matrix(raw) ~ meta.data$entrez.id, FUN = "median")
raw2 <- subset(raw2, !is.na(raw2[, "meta.data$entrez.id"]))
rownames(raw2) <- raw2[, "meta.data$entrez.id"]
raw2 <- raw2[, -1]

## log2 transform raw microarray intensities
exprs <- log2(data.matrix(raw2))

## write out processed raw data
write.csv(exprs, file = "human_exprs.csv")
