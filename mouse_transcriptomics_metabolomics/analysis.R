setwd("~/Documents/klaus/metabolomics_analysis_041414/analysis")

library(cluster)
library(gplots)
library(KEGGREST)
library(limma)
library(plyr)
library(RColorBrewer)

getGenesForPathway <- function(pathway) {
    record <- tryCatch({
        keggGet(pathway)[[1]]
    }, error = function(err) {
        NULL
    })

    record$GENE[seq(1, length(record$GENE), 2)]
}

getCompoundsForPathway <- function(pathway) {
    record <- tryCatch({
        keggGet(pathway)[[1]]
    }, error = function(err) {
        NULL
    })

    names(record$COMPOUND)
}

getGeneInfo <- function(gene, species) {
    print(paste("Getting info for gene ID: ", gene, sep=""))
    response <- c(EntrezID = gene, Name = NULL, Defintion = NULL)

    info <- tryCatch({
        keggGet(paste(species, gene, sep=":"))[[1]]
    }, error=function(err) {
        NULL
    })

    response["Name"] <- strsplit(info$NAME, ",")[[1]][1]
    response["Definition"] <- info$DEFINITION
    response
}

fillMissingValues <- function(data) {
    exprs <- data

    # Identify columns that contain missing values
    na.cols <- which(apply(exprs, 2, function(col) {
      any(is.na(col))
    }))

    # For cells with missing values, fill with minimum value for the metabolite
    for (col in na.cols) {
      col.data <- exprs[, col]
      col.data[is.na(col.data)] <- min(col.data, na.rm=T)

      exprs[, col] <- col.data
    }

    exprs
}

flipBranch <- function(dend) {
    d <- dend
    d[[1]] <- dend[[2]]
    d[[2]] <- dend[[1]]

    d
}

drawHeatMap <- function(data, groups, ...) {
    d <- data.matrix(data)

    # hc <- as.dendrogram(agnes(as.dist(1 - cor(d, method = "pearson")), 
    #         method = "average"))
    hr <- as.dendrogram(agnes(as.dist(1 - cor(t(d), method = "pearson")), 
            method = "average"))

    unique.groups <- unique(groups)
    col.cols <- as.character(mapvalues(groups, from = unique.groups, 
        to = brewer.pal(8, "Accent")[1:length(unique.groups)]))

    heatmap.2(data.matrix(d), 
        Rowv = hr, scale = "row", 
        trace = "none", col = bluered(25), keysize = 0.7, 
        cexRow = 1.7, cexCol = 2,
        ColSideColors = as.character(col.cols), ...)
}

pathways <- c("mmu00030", "mmu00471", "mmu00620", "mmu00020", "mmu00010")
genes <- sapply(pathways, getGenesForPathway)
compounds <- unique(c("C00236", "C00024", "C00026", "C00020", "C00002", "C00158", "C00354", "C00085", "C00122", "C00031", "C00103", "C00092", "C00186", "C00497", "C00006", "C00036", "C00074", "C00022", "C01151", "C00117", "C00042", "C00091", "C00993", "C00025", "C00064", "C00041", "C00003", "C00004", "C00068", 
    unlist(sapply(pathways, getCompoundsForPathway))))

# read in raw microarray data and aggregate
tx <- read.delim("OneChannel_Allgenes_res.txt", header=T, sep="\t", as.is=T, fill=T)
tx <- tx[, c(3, 24:35)]
tx.exprs <- aggregate(data.matrix(tx[, -1]) ~ tx$EntrezGene, FUN="median")
rownames(tx.exprs) <- tx.exprs[, 1]
tx.exprs <- data.matrix(tx.exprs[, -1])

# differential expression analysis
tx.conds <- factor(rep(c("Control", "Myc", "Ras"), each=4), 
    levels = c("Control", "Myc", "Ras"))
design <- model.matrix(~0 + tx.conds)
colnames(design) <- gsub("tx.conds", "", colnames(design))
fit <- lmFit(tx.exprs, design = design)
cont.matrix <- makeContrasts(MycEffect = Myc - Control, 
    RasEffect = Ras - Control, Ras_vs_Myc = Ras - Myc, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

tx.res <- data.frame(
    MycEffect.logFC = fit2$coefficients[, "MycEffect"],
    MycEffect.rawP = fit2$p.value[, "MycEffect"],
    MycEffect.fdr = p.adjust(fit2$p.value[, "MycEffect"], "fdr"),
    RasEffect.logFC = fit2$coefficients[, "RasEffect"],
    RasEffect.rawP = fit2$p.value[, "RasEffect"],
    RasEffect.fdr = p.adjust(fit2$p.value[, "RasEffect"], "fdr"),
    Ras_vs_Myc.logFC = fit2$coefficients[, "Ras_vs_Myc"],
    Ras_vs_Myc.rawP = fit2$p.value[, "Ras_vs_Myc"],
    Ras_vs_Myc.fdr = p.adjust(fit2$p.value[, "Ras_vs_Myc"], "fdr")
)

tx.res$EntrezID <- rownames(tx.res)
rownames(tx.res) <- NULL

tx.res <- subset(tx.res, tx.res$EntrezID %in% unlist(genes) & 
    tx.res$RasEffect.fdr < 0.05 & tx.res$Ras_vs_Myc.fdr < 0.05)
info <- as.data.frame(t(sapply(tx.res$EntrezID, getGeneInfo, "mmu")))
tx.res <- merge(tx.res, info)
tx.res <- tx.res[, c(1, 11, 12, 2:10)]

# transcriptomics heatmap
pdf("tx_heatmap.pdf", 11, 11)
tx.exprs2 <- tx.exprs[rownames(tx.exprs) %in% tx.res$EntrezID, ]
rownames(tx.exprs2) <- sapply(rownames(tx.exprs2), function(gene.id) {
    tx.res[tx.res$EntrezID == gene.id, "Name"]
})
tx.groups <- gsub("(.*)\\..*", "\\1", colnames(tx.exprs2))
## hack to convert exprs colnames
colnames(tx.exprs2) <- c(
    paste(rep("co", 4), 1:4, sep = "_"), 
    paste(rep("MT", 4), 1:4, sep = "_"), 
    paste(rep("RT", 4), 1:4, sep = "_"))
drawHeatMap(tx.exprs2, tx.groups, margins = c(15, 15))
dev.off()

# write out data
write.table(tx.res, file = "genes_of_interest.txt", sep = "\t", quote = F, 
    row.names = F)

# read in metabolomics data
mx <- read.table("20140410_metabolomics_complete.csv", header = T, sep = ",", 
    as.is = T, fill = T, row.names = 1)[
        -c(3:4, 8, 14:15, 16:21, 30:34, 35, 36:40, 46), ]
#[-c(3:4, 8, 14:15, 16:21, 30:34, 35, 36:40, 46), ] ## Leave out MT and RT_2w
#[-c(3:4, 8, 14:15, 22:28, 30:34, 35, 36:40, 46), ] ## Leave out RR and RT_2w

mx <- fillMissingValues(mx)

# read in compound name map
mx.names <- read.csv("name_map.csv", header = T , as.is = T, fill = T)
colnames(mx) <- mx.names$KEGG
mx.exprs <- t(mx)

# metabolomics heatmap
# heatmap
pdf("mx_heatmap.pdf", 11, 11)
mx.exprs2 <- mx.exprs[rownames(mx.exprs) %in% unlist(compounds), ]
rownames(mx.exprs2) <- sapply(rownames(mx.exprs2), function(kegg.id) {
    mx.names[mx.names$KEGG == kegg.id, "Match"]
})
mx.groups <- gsub("(.*)_.*", "\\1", colnames(mx.exprs2))
drawHeatMap(mx.exprs2, mx.groups, margins = c(15, 25))
dev.off()

source("~/Documents/Development/Repositories/r_tools/PCA.r")
pdf("mx_pca_analysis.pdf", 11, 11)
PCAAnalysis(mx.exprs2, mx.groups)
dev.off()
