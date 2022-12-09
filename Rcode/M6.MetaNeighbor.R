
DataPath <- "/home/yh/2-nano/backup/scrnaseq_resub/data"
setwd(DataPath)
library(dplyr)
library(Seurat)
library(patchwork)

pbmc1 <- readRDS("pbmc_rep1.rds")
pbmc2 <- readRDS("pbmc_rep2.rds")

ifnb.list <- list(n1 = pbmc1, n2 = pbmc2)

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nrow(x))
})
features <- VariableFeatures(ifnb.list$n1)[1:2000]
pbmc.combined <- readRDS("pbmc_integ.rds")
DimPlot(pbmc.combined, reduction = "umap", label = TRUE, repel = TRUE)


## MetaNeighborUS
library(SummarizedExperiment)

ricem <- GetAssayData(object = pbmc.combined, slot = "counts", assay = "RNA")
study_id <- pbmc.combined$orig.ident
cell_type <- dplyr::recode(pbmc.combined$seurat_clusters,
                           "0"="germ",
                           "1"="epi",
                           "2"="epi",
                           "3"="ov",
                           "4"="med",
                           "5"="sti",
                           "6"="med",
                           "7"="med",
                           "8"="cha",
                           "9"="epi",
                           "10"="ov")
cell_type <- paste(study_id, cell_type)
data.se <- SummarizedExperiment(list(gene_matrix=as.matrix(ricem, arabm)))
data.se$study_id <- study_id
data.se$cell_type <- cell_type


library(MetaNeighbor)
mn_data = data.se
var_genes = features
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = mn_data,
                             study_id = mn_data$study_id,
                             cell_type = mn_data$cell_type)
top_hits = topHitsByStudy(celltype_NV, threshold = 0.8)


dat <- celltype_NV
rownames(dat) <- colnames(dat)
annotation <- data.frame(Var1 = sub("[|].*", "", rownames(dat)))
annotation$Var1 <- sub("ricerep", "S", annotation$Var1)
rownames(annotation) <- colnames(dat)
annotation2 <- data.frame(Var2 =  sub(".* ", "C", rownames(dat)))
rownames(annotation2) <- colnames(dat)
annotation2$Var2 <- factor(annotation2$Var2, levels = c("Csti", "Cepi", "Cmed", "Ccha", "Cov", "Cgerm", "C11", "C12", "C13"))
pheatmap::pheatmap(dat, show_colnames = FALSE, annotation_col = annotation, annotation_row = annotation2, clustering_method="ward.D2")


