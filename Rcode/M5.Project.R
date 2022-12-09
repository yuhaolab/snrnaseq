
DataPath <- "/home/yh/2-nano/backup/scrnaseq_resub/data"
setwd(DataPath)
library(dplyr)
library(Seurat)
library(patchwork)

pbmc.combined <- readRDS("pbmc_integ.rds")
DimPlot(pbmc.combined, reduction = "umap", label = TRUE, repel = TRUE)

pancreas.integrated <- pbmc.combined
pancreas.integrated$celltype <- pancreas.integrated$seurat_clusters

## get query 1n --------------------------------------------------------------------------------------------------------------
pancreas.query <- readRDS("pbmc_1n.rds")

pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query, dims = 1:60, reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$celltype, dims = 1:60)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:60, reduction = "pca", return.model = TRUE, min.dist = 0.16)
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.integrated, query = pancreas.query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")
pancreas.query$predicted.celltype <- factor(pancreas.query$predicted.celltype, levels = as.character(0:13))

## plot
library(ggplot2)
library(RColorBrewer)
library(scales)

cols <- brewer.pal(11, "Spectral")
cols2 <- brewer.pal(11, "RdBu")
cols <- c(cols2[1], cols, cols2[11])
pal <- colorRampPalette(cols)
cols <- pal(14)
cols[8] <- "#E5DAAE"
cols <- cols[c(2:14, 1)]
names(cols) <- as.character(c(1, 2, 4, 5, 6, 7, 9, 0, 3, 8, 10, 11, 12, 13))
cols[c("7", "9")] <- cols[c("9", "7")]
cols <- cols[as.character(0:13)]

p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 6, cols = cols, pt.size = 1,
              label.color = "black", label.box = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 6, cols = cols, pt.size = 1,
              label.color = "black", label.box = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

table(pancreas.query$predicted.celltype)
table(pancreas.query$predicted.celltype)/ncol(pancreas.query)


## plot hist -----------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)

dat0 <- data.frame()
for (i in unique(predictions$predicted.id)) {
  dat <- data.frame(score = predictions$prediction.score.max[predictions$predicted.id == as.numeric(i)],
                    cluster = rep(i, times = sum(predictions$predicted.id == as.numeric(i))))
  dat0 <- rbind(dat0, dat)
}
dat0$predict <- dat0$cluster
dat0$cluster[dat0$cluster != "0"] <- "other"

group_color <- c("#e41a1c", "#377eb8")
ggplot(dat0, aes(x = score, color = cluster, fill = cluster)) + 
  geom_histogram(binwidth = 0.025, alpha = 0.5, color = "white", boundary = -0.5, position="dodge") + 
  scale_color_manual(values = group_color) +
  scale_fill_manual(values = group_color) +
  theme_classic() + expand_limits(x=c(0,1)) +
  scale_x_continuous(breaks=seq(0, 1, 0.1), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) # save 4*8


## plot combined -----------------------------------------------------------------------------------------------------------------------------
## sub cluster
pbmc1 <- readRDS("pbmc_rep1.rds")
pbmc2 <- readRDS("pbmc_rep2.rds")

pbmc_sub <- subset(pancreas.query, subset = predicted.celltype == "0")
pbmc_sub

pbmc1$s1 <- pbmc.combined$seurat_clusters[match(colnames(pbmc1), colnames(pbmc.combined))]
pbmc_sub1 <- subset(pbmc1, subset = s1 == "0")
pbmc2$s1 <- pbmc.combined$seurat_clusters[match(colnames(pbmc2), colnames(pbmc.combined))]
pbmc_sub2 <- subset(pbmc2, subset = s1 == "0")

ifnb.list <- list(n2_1 = pbmc_sub1, n2_2 = pbmc_sub2, n1 = pbmc_sub)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reference = 1)

all_features <- lapply(ifnb.list, row.names) %>% Reduce(c, .)
all_features <- unique(all_features)

immune.combined <- IntegrateData(anchorset = immune.anchors, features.to.integrate = all_features)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 40, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# DimPlot(immune.combined)
DimPlot(immune.combined, group.by = 'orig.ident', shuffle = TRUE)


cell.n.type <- immune.combined$orig.ident
immune.combined$cell.n.type <- cell.n.type
Idents(immune.combined) <- "cell.n.type"
DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = 0.5) +
  scale_color_brewer(palette = "Spectral")

ncols <- cols[c("0", "3", "4")]
names(ncols) <- c("ricerep1", "ricerep2", "rice1n")
DimPlot(immune.combined, label = FALSE, pt.size = 1.6, cols = ncols,
        label.color = "black", label.box = TRUE)


## get DE gene ----------------------------------------------------------------------------------------------------------------------------
DimPlot(immune.combined, reduction = "umap", group.by = "cell.n.type")

cell.n.type <- immune.combined$cell.n.type
cell.n.type[cell.n.type != "rice1n"] <- "rice2n"
immune.combined$cell.n.type <- cell.n.type
Idents(immune.combined) <- "cell.n.type"
DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = 0.5)

x <- c("rice1n", "rice2n")
de.markers <- FindMarkers(immune.combined, ident.1 = x[1], ident.2 = x[2], logfc.threshold = 0, min.pct = 0)
de.markers$gene <- rownames(de.markers)

anno <- read.delim("IRGSP-1.0_representative_annotation_2021-11-11.tsv", header = T, stringsAsFactors = F)
anno <- anno[match(de.markers$gene, anno$Locus_ID),]
markers <- cbind(de.markers, anno)
names(markers)[3:4] <- x
write.table(markers, file = paste0(x[1], ".vs.", x[2], ".xls"), sep = "\t", row.names = F, col.names = T, quote = F)




