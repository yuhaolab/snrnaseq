
DataPath <- "/home/yh/2-nano/backup/scrnaseq_resub/data"
setwd(DataPath)
library(BSgenome)
library(rtracklayer)
library(GenomicFeatures)
load("PtMtGene.RData")

geopath <- "/home/yh/2-nano/backup/scrnaseq/geodata"
setwd(geopath)
library(data.table)
library(dplyr)
library(Seurat)
library(patchwork)

min.cells = 1

pbmc.data <- fread("GSE185068/GSM5605641_S1-1.counts.tsv")
cc <- fread("GSE185068/GSM5605641_S1-1.cellname.list.txt")
pbmc.data <- data.frame(pbmc.data)
rownames(pbmc.data) <- pbmc.data$gene
pbmc.data <- pbmc.data[,-1]
colnames(pbmc.data) <- cc$CellName[match(colnames(pbmc.data), cc$CellIndex)]
pbmc.s1.1 <- CreateSeuratObject(counts = pbmc.data, project = "rices1.1", min.cells = min.cells, min.features = 0)
median(colSums(pbmc.data != 0))
pbmc.s1.1

pbmc.data <- fread("GSE185068/GSM5605642_S1-2.counts.tsv")
cc <- fread("GSE185068/GSM5605642_S1-2.cellname.list.txt")
pbmc.data <- data.frame(pbmc.data)
rownames(pbmc.data) <- pbmc.data$gene
pbmc.data <- pbmc.data[,-1]
colnames(pbmc.data) <- cc$CellName[match(colnames(pbmc.data), cc$CellIndex)]
pbmc.s1.2 <- CreateSeuratObject(counts = pbmc.data, project = "rices1.2", min.cells = min.cells, min.features = 0)
median(colSums(pbmc.data != 0))
pbmc.s1.2

pbmc.data <- fread("GSE185068/GSM5605643_S2-1.counts.tsv")
cc <- fread("GSE185068/GSM5605643_S2-1.cellname.list.txt")
pbmc.data <- data.frame(pbmc.data)
rownames(pbmc.data) <- pbmc.data$gene
pbmc.data <- pbmc.data[,-1]
colnames(pbmc.data) <- cc$CellName[match(colnames(pbmc.data), cc$CellIndex)]
pbmc.s2.1 <- CreateSeuratObject(counts = pbmc.data, project = "rices2.1", min.cells = min.cells, min.features = 0)
median(colSums(pbmc.data != 0))
pbmc.s2.1

pbmc.data <- fread("GSE185068/GSM5605644_S2-2.counts.tsv")
cc <- fread("GSE185068/GSM5605644_S2-2.cellname.list.txt")
pbmc.data <- data.frame(pbmc.data)
rownames(pbmc.data) <- pbmc.data$gene
pbmc.data <- pbmc.data[,-1]
colnames(pbmc.data) <- cc$CellName[match(colnames(pbmc.data), cc$CellIndex)]
pbmc.s2.2 <- CreateSeuratObject(counts = pbmc.data, project = "rices2.2", min.cells = min.cells, min.features = 0)
median(colSums(pbmc.data != 0))
pbmc.s2.2

pbmc.data <- fread("GSE185068/GSM5605645_L.counts.tsv")
cc <- fread("GSE185068/GSM5605645_L.cellname.list.txt")
pbmc.data <- data.frame(pbmc.data)
rownames(pbmc.data) <- pbmc.data$gene
pbmc.data <- pbmc.data[,-1]
colnames(pbmc.data) <- cc$CellName[match(colnames(pbmc.data), cc$CellIndex)]
pbmc.le <- CreateSeuratObject(counts = pbmc.data, project = "rices.le", min.cells = min.cells, min.features = 0)
median(colSums(pbmc.data != 0))
pbmc.le

pbmc <- merge(pbmc.s1.1, y = c(pbmc.s1.2, pbmc.s2.1, pbmc.s2.2, pbmc.le))
median(pbmc$nFeature_RNA)
mean(pbmc$nCount_RNA)
pbmc

id <- rownames(x = pbmc)
pt <- pt[is.element(pt, id)]
mt <- mt[is.element(mt, id)]
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = mt)
pbmc[["percent.pt"]] <- PercentageFeatureSet(pbmc, features = pt)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.pt", "percent.mt"), ncol = 4)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pbmc0 <- pbmc
pbmc0 <- subset(pbmc0, subset = nFeature_RNA >= 200)
id <- !is.element(rownames(pbmc0), mt) & rowSums(pbmc0@assays$RNA@counts != 0) > 2
pbmc0 <- pbmc0[id,]
pbmc0

id <- is.element(rownames(pbmc), rownames(pbmc0))
pbmc <- pbmc[id,]
pbmc ## object with exactly expressed genes and cells as discribed in the original published paper
mean(pbmc$nCount_RNA)

pbmc <- subset(pbmc, subset = nCount_RNA >= 1000) ## get high quality cells
pbmc

setwd(DataPath)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA", "percent.mt"))
pbmc <- RunPCA(pbmc, npcs = 60, features = VariableFeatures(object = pbmc))

# DimHeatmap(pbmc, dims = 1:20, cells = 500, balanced = TRUE)
# DimHeatmap(pbmc, dims = 21:40, cells = 500, balanced = TRUE)

# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# JackStrawPlot(pbmc, dims = 1:20)


pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:15, min.dist = 0.1)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.8)

DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident")

# pbmc <- RunTSNE(pbmc, dims = 1:10)
# DimPlot(pbmc, reduction = "tsne", label = TRUE, repel = TRUE)
# DimPlot(pbmc, reduction = "tsne", label = TRUE, repel = TRUE, group.by = "orig.ident")

DimPlot(pbmc, reduction = "umap", label = TRUE,label.size = 6,label.color = "black",label.box = TRUE)

## plot raw cluster markers in the original paper
FeaturePlot(pbmc, features = "Os05g0469600") #0
FeaturePlot(pbmc, features = "Os05g0495300") #2
FeaturePlot(pbmc, features = "Os02g0203700") #10
FeaturePlot(pbmc, features = "Os05g0349800") #3
FeaturePlot(pbmc, features = "Os10g0329400") #4
FeaturePlot(pbmc, features = "Os01g0660200") #6

FeaturePlot(pbmc, features = "Os09g0551600") #1
FeaturePlot(pbmc, features = "Os04g0209200") #11
FeaturePlot(pbmc, features = "Os04g0627000") #spikelet

FeaturePlot(pbmc, features = "Os03g0764100") #5
FeaturePlot(pbmc, features = "Os10g0471100") #9
FeaturePlot(pbmc, features = "Os02g0700100") #13

FeaturePlot(pbmc, features = "Os12g0555000") #8
FeaturePlot(pbmc, features = "Os03g0149200") #7

saveRDS(pbmc, file = "pbmc1000.15.rds")
