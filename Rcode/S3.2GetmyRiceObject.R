
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

pbmc.combined <- readRDS("pbmc_integ.rds")
features <- VariableFeatures(ifnb.list$n1)[1:2000]
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)

pbmc.combined <- RunPCA(pbmc.combined, npcs = 60, verbose = FALSE, features = features)
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:60)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.4)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:60, min.dist = 0.16)
DimPlot(pbmc.combined, reduction = "umap", label = TRUE, repel = TRUE)


pbmc <- CreateSeuratObject(counts = pbmc.combined@assays$RNA@counts, project = "myrice", min.cells = 0, min.features = 0)
pbmc$clusters <- pbmc.combined$seurat_clusters
pbmc$clusters <- dplyr::recode(pbmc$clusters,
                               "0"="nu",
                               "1"="ep",
                               "2"="ep",
                               "3"="ov",
                               "4"="ow",
                               "5"="sti",
                               "6"="ow",
                               "7"="ow",
                               "8"="ch",
                               "9"="ep",
                               "10"="ov",
                               "11"="sta",
                               "12"="c12",
                               "13"="c13")
saveRDS(pbmc, "mypbmc.rds")

