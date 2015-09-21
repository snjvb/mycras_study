setwd("~/Documents/klaus/mouse_human_liver_tumor_comparison_060214")

library(FactoMineR)
library(gplots)
library(KEGGREST)
library(limma)
library(RColorBrewer)
library(sva)

## prepare orthologs dataframe (sourced from MGI Informatics)
orthologs <- read.delim("HOM_MouseHumanSequence.rpt", as.is = T, 
    fill = T)[, c(1, 2, 4, 5)]
orthologs <- split(orthologs[, -2], orthologs$Common.Organism.Name)
names(orthologs) <- c("human", "mouse")
orthologs <- merge(orthologs$mouse, orthologs$human, by = "HomoloGene.ID")
names(orthologs) <- c("homologene.id", "mouse.symbol", "mouse.entrez.id", 
    "human.symbol", "human.entrez.id")

## extract mouse genes of interest - belonging to specified metabolic
## pathways
pathways <- c("hsa00030", "hsa00471", "hsa00620", "hsa00020", "hsa00010")
genes <- sapply(pathways, function(pathway) {
    record <- tryCatch({
        keggGet(pathway)[[1]]
    }, error = function(err) {
        NULL
    })

    record$GENE[seq(1, length(record$GENE), 2)]
})

orthologs2 <- orthologs[orthologs$human.entrez.id %in% unlist(genes), ]
## manually correct ACAT2, ENO1 and GAPDH mappings
orthologs2 <- orthologs2[-which(orthologs2$mouse.entrez.id 
    %in% c(224530, 433182, 100042025)), ]

## read in human microarray data
hsa <- read.csv("human_exprs.csv", as.is = T, fill = T, row.names = 1)

## differential expression analysis for human microarray
hsa2 <- hsa
colnames(hsa2) <- gsub("GSM[0-9]{6}\\.", "", colnames(hsa2))

design <- model.matrix(~0 + factor(colnames(hsa2), 
    levels = c("non.tumor", "tumor")))
colnames(design) <- c("non.tumor", "tumor")

fit <- lmFit(hsa2, design = design)
cont.matrix <- makeContrasts(effect = tumor - non.tumor, levels = design)
fit2 <- contrasts.fit(fit, contrasts = cont.matrix)
fit2 <- eBayes(fit2)

res <- data.frame(
    logFC = fit2$coefficients[, "effect"],
    rawP = fit2$p.value[, "effect"],
    fdr = p.adjust(fit2$p.value[, "effect"], "BH")
)

degs <- rownames(res)[res$fdr < 0.01]

## read in mouse microarray data
mmu <- read.csv("mouse_exprs.csv", as.is = T, fill = T, row.names = 1)

## identify genes of interest that may not be represented in the mouse and
## human microarrays, and remove them
absent <- which(!orthologs2$human.entrez.id %in% rownames(hsa) | 
    !orthologs2$mouse.entrez.id %in% rownames(mmu))
orthologs3 <- orthologs2[-absent, ]

## only extract significantly DEG genes for subsequent analysis
orthologs3 <- orthologs3[orthologs3$human.entrez.id %in% degs, ]

## merge human and mouse microarray data based on gene orthology
hsa <- hsa[rownames(hsa) %in% orthologs3$human.entrez.id, ]
hsa$homologene <- sapply(rownames(hsa), function(id) {
    orthologs3[orthologs3$human.entrez.id == id, "homologene.id"]
})

mmu <- mmu[rownames(mmu) %in% orthologs3$mouse.entrez.id, ]
mmu$homologene <- sapply(rownames(mmu), function(id) {
    orthologs3[orthologs3$mouse.entrez.id == id, "homologene.id"]
})

exprs <- merge(hsa, mmu, by = "homologene")
rownames(exprs) <- exprs$homologene
exprs <- exprs[, -1]

## perform pre-normalization PCA analysis
extractGroup <- function(sample) {
    temp <- gsub("\\.[0-9]+", "", sample)
    gsub("GSM[0-9]{6}\\.", "", temp)
}
samples <- data.frame(id = colnames(exprs), 
    group = sapply(colnames(exprs), extractGroup), row.names = NULL)
groups <- samples$group
pca.exprs <- cbind(groups, as.data.frame(t(exprs)))
pca.res <- PCA(pca.exprs, quali.sup = 1, ncp = 5, graph = F)

concat <- cbind.data.frame(pca.exprs[, 1], pca.res$ind$coord)
ellipse.coord <- coord.ellipse(concat, level.conf = 0.99, bary = T)

png("combined_microarray_pre_norm_pca.png", 1024, 1024, "px")
plot.PCA(pca.res, axes=c(1, 2), choix = "ind", habillage = 1, 
    ellipse=ellipse.coord, label = "quali")
dev.off()

## correct for batch effects in merged dataset
batch <- sapply(groups, function(group) {
    if (grepl("tumor", group)) {
        "human"
    } else {
        "mouse"
    }
})

tissue <- sapply(groups, function(group) {
    if (grepl("Control", group) || grepl("non.tumor", group)) {
        "control"
    } else {
        "tumor"
    }
})

exprs2 <- ComBat(exprs, batch = batch, mod = model.matrix(~as.factor(tissue)), 
    numCovs = NULL, par.prior = TRUE, prior.plots = FALSE)

## perform post-normalization PCA analysis
pca.exprs2 <- cbind(groups, as.data.frame(t(exprs2)))
pca.res2 <- PCA(pca.exprs2, quali.sup = 1, ncp = 5, graph = F)

concat2 <- cbind.data.frame(pca.exprs2[, 1], pca.res2$ind$coord)
ellipse.coord2 <- coord.ellipse(concat2, level.conf = 0.99, bary = T)

png("combined_microarray_post_norm_pca.png", 1024, 1024, "px")
plot.PCA(pca.res2, axes=c(1, 2), choix = "ind", habillage = 1, 
    ellipse=ellipse.coord2, label = "quali")
dev.off()

## Update gene and column labels
rownames(exprs2) <- sapply(rownames(exprs2), function(id) {
    orthologs3[orthologs3$homologene.id == id, "human.symbol"]
}, USE.NAMES = F)

groups2 <- factor(groups)
levels(groups2) <- c("Mouse Non-Tumor", "Mouse MYC-Overexpressing Tumor", 
    "Human Non-Tumor", "Mouse RAS-Overexpressing Tumor", "Human Tumor")
colnames(exprs2) <- groups2

## plot a heatmap of the adjusted microarray expressions
hc <- agnes(as.dist(1 - cor(exprs2, method = "pearson")), 
    method = "average")
hr <- agnes(as.dist(1 - cor(t(exprs2), method = "pearson")), 
    method = "average")

## Generate color map for samples
groups.cols <- data.frame(group = names(table(groups)), 
    color = brewer.pal(length(table(groups)), "Accent"))
col.cols <- as.character(sapply(groups, function(group) {
    groups.cols[groups.cols$group == group, "color"]
}))

png("combined_microarray_sample_clustering.png", 2048, 640, "px")
plot(hc, axes = F, ann = F, ylab = "", sub = "", main = "", hang = -1)
dev.off()

png("combined_microarray_adj_exprs_heatmap.png", 1024, 1024, "px")
heatmap.2(data.matrix(exprs2), margins = c(2, 10), Colv = as.dendrogram(hc), 
    Rowv = as.dendrogram(hr), scale = "row", trace = "none", 
    col = bluered(25), keysize = 1, ColSideColors = as.character(col.cols), 
    labCol = rep("", ncol(exprs2)))
dev.off()
