
DataPath <- "/home/yh/2-nano/backup/scrnaseq_resub/data"
setwd(DataPath)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

pbmc.data <- Read10X(data.dir = "arab_stage5/", gene.column = 1)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "arab", min.cells = 0)
pbmc
median(pbmc$nFeature_RNA)
pbmc <- subset(pbmc, subset = nCount_RNA > 1000)

ptmt <- c(grep("ATM", rownames(pbmc)), grep("ATC", rownames(pbmc)))
pbmc <- pbmc[-ptmt,]
sum(pbmc$nFeature_RNA >= 1000)
sum(pbmc$nCount_RNA > 1000)

id <- rowSums(pbmc@assays$RNA@counts != 0) > 4
pbmc <- pbmc[id,]
pbmc ## object with exactly the same number of genes and cells as described in the original published paper

pbmc <- SCTransform(pbmc, verbose = FALSE, variable.features.n = 3000)

pbmc <- RunPCA(pbmc, npcs=9)
pbmc <- RunUMAP(pbmc, dims = 1:9, reduction = "pca", n.neighbors = 50, min.dist = 0.01, umap.method = "uwot")
pbmc <- FindNeighbors(pbmc, dims = 1:9)
pbmc <- FindClusters(pbmc, resolution = 2)
DimPlot(pbmc, reduction = "umap", label = TRUE)


## asigned raw celltype discribed in the published paper
celltype <- pbmc$seurat_clusters
celltype <- dplyr::recode(celltype,
                          "1"="epidermis", "5"="epidermis", "17"="epidermis", "24"="epidermis",
                          "2"="pro-cambium", "13"="pro-cambium", "14"="pro-cambium", "22"="pro-cambium",
                          "3"="cortex", "10"="cortex", "12"="cortex",
                          "11"="dividingCells", "15"="dividingCells", "20"="dividingCells", "25"="dividingCells",
                          "4"="floralMeristem", "9"="floralMeristem",
                          "0"="S-phaseCells",
                          "7"="dividingCells",
                          "6"="mesophyll",
                          "18"="xylemParenchyma", "21"="xylemParenchyma",
                          "8"="epidermis", "19"="epidermis",
                          "16"="epidermis/dividing",
                          "23"="epidermis/dividing",
                          "26"="phloem")

pbmc$celltype <- celltype
DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "celltype")


## asign combined cell type
celltype <- pbmc$seurat_clusters
celltype <- dplyr::recode(celltype,
                          "1"="epidermis", "5"="epidermis", "17"="epidermis", "24"="epidermis",
                          "2"="vis", "13"="vis", "14"="vis", "22"="vis",
                          "3"="cortex", "10"="cortex", "12"="cortex",
                          "11"="dividingCells", "15"="dividingCells", "20"="dividingCells", "25"="dividingCells",
                          "4"="floralMeristem", "9"="floralMeristem",
                          "0"="dividingCells",
                          "7"="dividingCells",
                          "6"="cortex",
                          "18"="vis", "21"="vis",
                          "8"="epidermis", "19"="epidermis",
                          "16"="dividingCells",
                          "23"="dividingCells",
                          "26"="vis")
pbmc$celltype <- celltype
DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "celltype")

arabm <- GetAssayData(object = pbmc, slot = "counts")
arabcluster <- paste0("A", pbmc$celltype)


## get rice
pbmc.combined <- readRDS("pbmc_integ.rds")
DimPlot(pbmc.combined, reduction = "umap", label = TRUE, repel = TRUE)

ricem <- GetAssayData(object = pbmc.combined, slot = "counts", assay = "RNA")
ricecluster <- paste0("R", pbmc.combined$seurat_clusters)


## combine -----------------------------------------------------------------------------------------------------------------
library(SummarizedExperiment)
geneid <- read.delim("arab2rice_gene52.txt", header = T, stringsAsFactors = F)
id <- grep("[.]1", geneid$Transcript.stable.ID)
geneid <- geneid[id,]
id <- grep("-01", geneid$Transcript.stable.ID.1)
geneid <- geneid[id,]

geneid <- data.frame(arab = geneid$Gene.stable.ID, rice = geneid$Gene.stable.ID.1)
geneid <- unique(geneid)
id <- is.element(geneid$arab, rownames(arabm)) & is.element(geneid$rice, rownames(ricem))
geneid <- geneid[id,]

riceid <- match(geneid$rice, rownames(ricem))
ricem <- ricem[riceid,]
arabid <- match(geneid$arab, rownames(arabm))
arabm <- arabm[arabid,]

idname <- paste(geneid$arab, geneid$rice, sep = "_")
rownames(ricem) <- idname
rownames(arabm) <- idname

colnames(ricem) <- paste0("rice", 1:ncol(ricem))
colnames(arabm) <- paste0("arab", 1:ncol(arabm))

study_id <- rep(c("Rice", "Arab"), times = c(ncol(ricem), ncol(arabm)))
cell_type <- c(ricecluster, arabcluster)

data.se <- SummarizedExperiment(list(gene_matrix=as.matrix(cbind(ricem, arabm))))
data.se$study_id <- study_id
data.se$cell_type <- cell_type

saveRDS(data.se, file = "combine.se.rds")














