
DataPath <- "/home/yh/2-nano/backup/scrnaseq_resub/data"
setwd(DataPath)
library(dplyr)
library(Seurat)
library(patchwork)
pancreas.integrated <- readRDS("pbmc1000.15.rds")
DimPlot(pancreas.integrated, reduction = "umap", label = TRUE, label.size = 6, label.color = "black", label.box = TRUE)

celltype <- as.character(pancreas.integrated$seurat_clusters)
celltype <- dplyr::recode(celltype,
                          "0"="ra",
                          "1"="sp",
                          "2"="lf",
                          "3"="ra",
                          "4"="me",
                          "5"="lf",
                          "6"="ra",
                          "7"="ra",
                          "8"="ra",
                          "9"="ra",
                          "10"="sp",
                          "11"="me",
                          "12"="lf",
                          "13"="me",
                          "14"="lf",
                          "15"="me",
                          "16"="lf")
pancreas.integrated$celltype <- factor(as.character(celltype), levels = unique(as.character(celltype)))
DimPlot(pancreas.integrated, reduction = "umap", label = TRUE, repel = TRUE, group.by = "celltype")

# pbmc <- RunTSNE(pancreas.integrated, dims = 1:15)
# DimPlot(pbmc, reduction = "tsne", label = TRUE, repel = TRUE, group.by = "celltype")


pancreas.query <- readRDS("mypbmc.rds")
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query, dims = 1:15, reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$celltype, dims = 1:15)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:15, reduction = "pca", return.model = TRUE, min.dist = 0.1)
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.integrated, query = pancreas.query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

pancreas.query$predicted.celltype <- factor(pancreas.query$predicted.celltype, levels = unique(as.character(celltype)))
DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = F, label.size = 6, pt.size = 0.4, label.color = "black")
table(pancreas.query$predicted.celltype)


## plot -----------------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(scales)

group.colors <- brewer.pal(4, "Paired")
names(group.colors) <- c("sp", "me", "ra", "lf")

DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = F, label.size = 6, pt.size = 0.4,
        label.color = "black", cols = group.colors) ## save 6*6

p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = F, label.size = 6, pt.size = 0.4, order = "sp",
              label.color = "black", cols = group.colors) + ggtitle("Reference annotations")

p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = F, label.size = 6, pt.size = 0.4, order = "sp",
              label.color = "black", cols = group.colors) + ggtitle("Query transferred labels")
p1 + p2 ## save 6*12

pancreas.query$seurat_celltype <- pancreas.query$clusters
table(paste(pancreas.query$predicted.celltype, pancreas.query$seurat_celltype))

cols <- brewer.pal(9, "Spectral")
names(cols) <- c("nu", "ov", "ch", "ow", "ep", "sti", "sta", "c12", "c13")
DimPlot(pancreas.query, reduction = "ref.umap", group.by = "seurat_celltype", label = F, order = "ov",
        label.size = 6, pt.size = 0.4, label.color = "black", cols = cols) ## save 6*6.5


hist(pancreas.query$predicted.celltype.score)

dat <- data.frame(celltype = pancreas.query$clusters, predtype = pancreas.query$predicted.celltype, auc = pancreas.query$predicted.celltype.score)
library(ggplot2)
dat$predtype <- factor(dat$predtype, levels = c("sp", "ra", "me", "lf"))
group.colors <- brewer.pal(4, "Paired")
names(group.colors) <- c("sp", "me", "ra", "lf")
ggplot(dat, aes(x = auc, color = predtype, fill = predtype)) + geom_density(alpha = 0.2) + theme_classic() + scale_color_manual(values = group.colors) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values = group.colors)


## sk plot --------------------------------------------------------------------------------------------------------------------------------------
dat <- data.frame(celltype = pancreas.query$clusters, predtype = pancreas.query$predicted.celltype, auc = pancreas.query$predicted.celltype.score)
dat <- dat[dat$auc >= 0.95,]
indr <- as.character(unique(dat$celltype))
ind <- as.character(unique(dat$predtype))
links <- data.frame()
for (i in 1:length(indr)) {
  for (j in 1:length(ind)) {
    linkij <- data.frame(source = ind[j], target = indr[i], value = sum(dat$celltype == indr[i] & dat$predtype == ind[j]))
    links <- rbind(links, linkij)
  }
}

id <- c("c12", "ch", "c13", "nu", "ov", "ep", "ow", "sti", "sta", "ra", "sp")
nodes = data.frame("name" = id)
links$source <- match(links$source, nodes$name) - 1
links$target <- match(links$target, nodes$name) - 1
links$group <- "Pos"

## plot
library(networkD3)
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              fontSize= 12, nodeWidth = 30,
              iterations = 0) ## save 600*600


my_color <- 'd3.scaleOrdinal() 
.domain(["Pos","sti","ep","ov","nu","ch","13","c12","ow","sta","ra","sp"]) 
.range(["lightblue","#99003F","#CC344D","#EC6245","#FA9A57","#FDCC7A","#E5DAAE","#FEF0A7","#574C9C","#053061","#004C99","#00CCCC"])'
## plot
sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target", LinkGroup = "group",
              Value = "value", NodeID = "name", colourScale=my_color, fontSize= 16, nodeWidth = 30, nodePadding = 12,
              iterations = 0, sinksRight = F) ## save 600*600


