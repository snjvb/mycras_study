setwd("~/Documents/klaus/mouse_human_liver_tumor_comparison_060214")

library(limma)

data <- read.csv("human_exprs.csv", as.is = T, fill = T, check.names = F, 
    row.names = 1)
groups <- gsub(".+?\\.(.*)", "\\1", colnames(data))

diffEq <- function(exprs, groups) {
    design <- model.matrix(~0 + groups)
    colnames(design) <- gsub("groups", "", colnames(design))

    fit <- lmFit(exprs, design = design)
    contrasts <- makeContrasts(effect = tumor - non.tumor, levels = design)
    fit2 <- contrasts.fit(fit, contrasts = contrasts)
    fit2 <- eBayes(fit2)

    res <- data.frame(
        logFC = fit2$coefficients[, "effect"],
        rawP = fit2$p.value[, "effect"],
        fdr = p.adjust(fit2$p.value[, "effect"], "fdr"),
        holm = p.adjust(fit2$p.value[, "effect"], "holm")
    )

    res[order(res$rawP), ]
}

res <- diffEq(data, groups)

source("~/Documents/r_tools/PathwayAnalysis.R")
ora.analysis <- ORAnalysis(rownames(res), res$holm)
write.csv(ora.analysis, file = "human_exprs_ora_analysis.csv", row.names = F)

pathways <- c("hsa00030", "hsa00471", "hsa00620", "hsa00020", "hsa00010")
ora.pathways <- gsub("(.*)\\|.*", "\\1", ora.analysis$Pathway)
ora.analysis.subset <- ora.analysis[ora.pathways%in% pathways, ]
write.csv(ora.analysis.subset, file = "human_exprs_ora_analysis_subset.csv", 
    row.names = F)
