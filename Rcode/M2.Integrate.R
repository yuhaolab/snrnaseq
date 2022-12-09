
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

features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 25000)
pbmc.anchors <- FindIntegrationAnchors(object.list = ifnb.list, reference = 1, anchor.features = features)
pbmc.combined <- IntegrateData(anchorset = pbmc.anchors)
DefaultAssay(pbmc.combined) <- "integrated"

features <- VariableFeatures(ifnb.list$n1)[1:2000]
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)

pbmc.combined <- RunPCA(pbmc.combined, npcs = 60, verbose = FALSE, features = features)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:60, min.dist = 0.16)
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:60)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.4)
DimPlot(pbmc.combined, reduction = "umap", label = TRUE, repel = TRUE)
saveRDS(pbmc.combined, file = "pbmc_integ.rds")

## get makers
pbmc.markers <- FindAllMarkers(pbmc.combined, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 1, assay = "RNA")
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

anno <- read.delim("IRGSP-1.0_representative_annotation_2021-11-11.tsv", header = T, stringsAsFactors = F)
anno <- anno[match(pbmc.markers$gene, anno$Locus_ID),]
markers <- cbind(pbmc.markers, anno)
write.table(markers, file = "allmarkers.xls", sep = "\t", row.names = F, col.names = T, quote = F)




