
DataPath <- "/home/yh/2-nano/backup/scrnaseq_resub/data"
setwd(DataPath)
library(dplyr)
library(Seurat)
library(patchwork)
load("PtMtGene.RData")

pbmc.data <- Read10X(data.dir = "filtered_feature_bc_matrix_rep1/", gene.column = 1)
pbmc1 <- CreateSeuratObject(counts = pbmc.data, project = "ricerep1", min.cells = 5, min.features = 200)
pbmc1

pbmc.data <- Read10X(data.dir = "filtered_feature_bc_matrix_rep2/", gene.column = 1)
pbmc2 <- CreateSeuratObject(counts = pbmc.data, project = "ricerep2", min.cells = 5, min.features = 200)
pbmc2

pbmc.data <- Read10X(data.dir = "filtered_feature_bc_matrix_1n/", gene.column = 1)
pbmc1n <- CreateSeuratObject(counts = pbmc.data, project = "rice1n", min.cells = 5, min.features = 200)
pbmc1n

## filter pbmc1 -----------------------------------------------------------------------------------------------------------------------------
pbmc <- pbmc1
## remove doublet
doublet1 <- read.table("doublet1.1.xls", header = T, stringsAsFactors = F)
doublet2 <- read.table("doublet1.2.xls", header = T, stringsAsFactors = F)
doublet <- data.frame(cellname = colnames(pbmc), doublet = "Singlet")

id1 <- doublet1$doublet[match(doublet$cellname, doublet1$cellname)] == "Doublet"
id1[is.na(id1)] <- FALSE
id2 <- doublet2$doublet[match(doublet$cellname, doublet2$cellname)] == "Doublet"
id2[is.na(id2)] <- FALSE

doublet$doublet[id1 | id2] <- "Doublet"
sum(doublet$cellname == colnames(pbmc))
pbmc[["doublet_filter"]] <- doublet$doublet
pbmc <- subset(pbmc, subset = doublet_filter == "Singlet")
pbmc

## other filter
id <- rownames(x = pbmc)
pt <- pt[is.element(pt, id)]
mt <- mt[is.element(mt, id)]
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = mt)
pbmc[["percent.pt"]] <- PercentageFeatureSet(pbmc, features = pt)

pbmc <- subset(pbmc, subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.pt < 5 & percent.mt < 10)
pbmc

pbmc <- RenameCells(object = pbmc, add.cell.id = "rep1")
saveRDS(pbmc, file = "pbmc_rep1.rds")

## umi density filter --------------------------------------------------------------------------------------------------------------------------
rep1.nc <- pbmc1$nCount_RNA
plot(density(rep1.nc))

rep2.nc <- pbmc2$nCount_RNA
plot(density(rep2.nc))

d1 <- density(rep1.nc)
d2 <- density(rep2.nc)

dat <- data.frame(x = c(d1$x, d2$x), y = c(d1$y, d2$y), repn = rep(c("rep1", "rep2"), times = c(length(d1$x), length(d2$x))))
library(ggplot2)
s1d <- as.numeric(quantile(rep1.nc, 0.1))
p1 <- ggplot(dat, aes(x=x, y=y, color=repn)) + geom_line() + theme_classic() + xlab("Total mumber of UMIs") + ylab("Density") +
  geom_vline(xintercept = s1d, colour = "black", linetype = "dashed") +
  scale_x_continuous(expand = c(0,0), limits = c(min(dat$x), max(dat$x)*1.1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(dat$y)*1.1)) +
  annotate("text", x = s1d + 0.2*max(dat$x), y = max(dat$y)*1.05, label= paste0("Thresh = ", round(s1d, 2)), color = "black")
print(p1)

dat <- data.frame(genenum = c(pbmc1$nFeature_RNA, pbmc2$nFeature_RNA), group = c(rep("rep1", ncol(pbmc1)), rep("rep2", ncol(pbmc2))))
p1 <- ggplot(dat, aes(x = genenum, fill = group)) +
  geom_histogram(colour = "black", bins = 30, position = "dodge") +
  theme_bw()
p1

sum(pbmc2$nCount_RNA > s1d)

## rep2 --------------------------------------------------------------------------------------------------------------------------------------
pbmc <- pbmc2

## remove doublet
doublet <- read.table("doublet2.xls", header = T, stringsAsFactors = F)
rownames(doublet) <- doublet$cellname
doublet <- doublet[colnames(pbmc),]

sum(doublet$cellname == colnames(pbmc))
pbmc[["doublet_filter"]] <- doublet$doublet
pbmc <- subset(pbmc, subset = doublet_filter == "Singlet")
pbmc

id <- rownames(x = pbmc)
pt <- pt[is.element(pt, id)]
mt <- mt[is.element(mt, id)]
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = mt)
pbmc[["percent.pt"]] <- PercentageFeatureSet(pbmc, features = pt)

pbmc <- subset(pbmc, subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.pt < 5 & percent.mt < 10)
pbmc

pbmc <- subset(pbmc, subset = nCount_RNA > s1d)
pbmc

pbmc <- RenameCells(object = pbmc, add.cell.id = "rep2")
saveRDS(pbmc, file = "pbmc_rep2.rds")


## 1n ------------------------------------------------------------------------------------------------------------------------------------------
pbmc <- pbmc1n
## remove doublet
doublet <- read.table("doublet.1n.xls", header = T, stringsAsFactors = F)
pbmc[["doublet_filter"]] <- doublet$doublet
pbmc <- subset(pbmc, subset = doublet_filter == "Singlet")
pbmc

## other filter
id <- rownames(x = pbmc)
pt <- pt[is.element(pt, id)]
mt <- mt[is.element(mt, id)]
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = mt)
pbmc[["percent.pt"]] <- PercentageFeatureSet(pbmc, features = pt)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 500 & percent.pt < 5 & percent.mt < 10)
pbmc

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

saveRDS(pbmc, file = "pbmc_1n.rds")















