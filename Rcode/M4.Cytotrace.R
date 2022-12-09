
DataPath <- "/home/yh/2-nano/backup/scrnaseq_resub/data"
setwd(DataPath)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(monocle3)
library(ggplot2)
library(dplyr)

cds <- readRDS("cds_comb.rds")
cds <- preprocess_cds(cds, num_dim = 60)
cds <- align_cds(cds, alignment_group = "orig.ident")
cds <- reduce_dimension(cds)
plot_cells(cds, group_label_size = 6, cell_size = 0.5, color_cells_by = "seurat_clusters")

## get partition
cds <- cluster_cells(cds)
cds@colData$part <- partitions(cds, reduction_method = "UMAP")
plot_cells(cds, color_cells_by = "partition", label_cell_groups = FALSE)


## get cytotrace
exp_10x <- exprs(cds)
exp_10x <- as.matrix(exp_10x)
exp_10x <- data.frame(exp_10x)

library(CytoTRACE)
results <- CytoTRACE(exp_10x, enableFast = FALSE)
save(results, file = "alltrace.RData")
# load("alltrace.RData")

cytops <- results$CytoTRACE
names(cytops) <- sub("[.]", "-", names(cytops))
cytops <- 1-cytops
cds@colData$cytotrace = cytops


## sub1 ---------------------------------------------------------------------------------------------------------------------------
## plot sub traj 124567
cds_sub <- cds
valid_cells <- row.names(subset(pData(cds_sub), is.element(seurat_clusters, c(1,2,4,5,6,7))))
cds_sub <- cds_sub[,valid_cells]

cds_sub <- preprocess_cds(cds_sub, num_dim = 50)
cds_sub <- align_cds(cds_sub, alignment_group = "orig.ident")
cds_sub <- reduce_dimension(cds_sub)
plot_cells(cds_sub, group_label_size = 6, cell_size = 0.8, color_cells_by = "seurat_clusters")


## plot traj
cds_sub <- cluster_cells(cds_sub)
plot_cells(cds_sub, color_cells_by = "partition")
cds_sub <- learn_graph(cds_sub, use_partition = F, close_loop = F, learn_graph_control = list(minimal_branch_len = 15))
plot_cells(cds_sub, label_principal_points = T)

cds_sub <- order_cells(cds_sub, root_pr_nodes = "Y_44")
plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=4,
           cell_size = 1,
           label_principal_points = TRUE,
           trajectory_graph_segment_size = 1) ## save 6*7

## plot
## plot combine
p2 <- plot_cells(cds_sub,
                 color_cells_by = "cytotrace",
                 show_trajectory_graph = T,
                 # trajectory_graph_color = "grey40",
                 trajectory_graph_segment_size = 1,
                 label_cell_groups=F,
                 label_leaves=F,
                 label_branch_points=F,
                 label_principal_points = T,
                 graph_label_size=4,
                 cell_size = 0.8) +
  viridis::scale_color_viridis(discrete = F, option="C")

p3 <- plot_cells(cds_sub,
                 color_cells_by = "pseudotime",
                 show_trajectory_graph = T,
                 # trajectory_graph_color = "grey40",
                 trajectory_graph_segment_size = 1,
                 label_cell_groups=F,
                 label_leaves=F,
                 label_branch_points=F,
                 label_principal_points = T,
                 graph_label_size=4,
                 cell_size = 0.8)

library(patchwork)
p2+p3


## plot cor
groupid <- list(l721 = c("Y_44", "Y_43"),
                l746 = c("Y_44", "Y_34"),
                l745 = c("Y_44", "Y_74"))

monops <-  pseudotime(cds_sub, reduction_method = "UMAP")
corrdata <- data.frame(cyto = cytops[names(monops)], mono = monops, cluster = cds_sub@colData$seurat_clusters)
corrdata <- corrdata[corrdata$mono != Inf,]
cor(corrdata$cyto, corrdata$mono)

lapply(groupid, function(x, corrdata){
  id <- choose_graph_segments(cds_sub, starting_pr_node = x[1], ending_pr_nodes = x[2], return_list = T)
  id <- is.element(rownames(corrdata), id$cells)
  dat <- data.frame(cyto = corrdata$cyto[id], mono = corrdata$mono[id], cols = corrdata$cluster[id])
  ggplot(dat, aes(x = mono, y = cyto)) + geom_point(shape = 20, colour = "grey") + 
    geom_smooth(method = "loess", span = 2, n = 100, fill = "lightblue") + theme_classic() + 
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(limits = c(0,1), expand = c(0,0))
  ggsave(paste0("l", paste(x, collapse = ""), ".pdf"),  width = 5, height = 4)
}, corrdata)


## sub2 -------------------------------------------------------------------------------------------------------------------------------------
## plot sub traj 03810
cds_sub <- cds
valid_cells <- row.names(subset(pData(cds_sub), is.element(seurat_clusters, c(0,3,8,10)) & part == 2))
cds_sub <- cds_sub[,valid_cells]

cds_sub <- preprocess_cds(cds_sub , num_dim = 20)
cds_sub <- align_cds(cds_sub, alignment_group = "orig.ident")
cds_sub <- reduce_dimension(cds_sub)
plot_cells(cds_sub, group_label_size = 6, cell_size = 0.8, color_cells_by = "seurat_clusters")


## plot traj
cds_sub <- cluster_cells(cds_sub)
plot_cells(cds_sub, color_cells_by = "partition")
cds_sub <- learn_graph(cds_sub, use_partition = FALSE, close_loop = F, learn_graph_control = list(rann.k = 25, minimal_branch_len = 15, ncenter = 400))

cds_sub <- order_cells(cds_sub, root_pr_nodes = "Y_137")
plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=4,
           cell_size = 1,
           label_principal_points = TRUE,
           trajectory_graph_segment_size = 1)


## plot combine
p2 <- plot_cells(cds_sub,
                 color_cells_by = "cytotrace",
                 show_trajectory_graph = T,
                 # trajectory_graph_color = "grey40",
                 trajectory_graph_segment_size = 1,
                 label_cell_groups=F,
                 label_leaves=F,
                 label_branch_points=F,
                 graph_label_size=4,
                 label_principal_points = TRUE,
                 cell_size = 0.8) +
  viridis::scale_color_viridis(discrete = F, option="C")

p3 <- plot_cells(cds_sub,
                 color_cells_by = "pseudotime",
                 show_trajectory_graph = T,
                 # trajectory_graph_color = "grey40",
                 trajectory_graph_segment_size = 1,
                 label_cell_groups=F,
                 label_leaves=F,
                 label_branch_points=F,
                 graph_label_size=4,
                 label_principal_points = TRUE,
                 cell_size = 0.8)

library(patchwork)
p2+p3

## plot cor
groupid <- list(l108 = c("Y_137", "Y_5"),
                l100 = c("Y_137", "Y_111"),
                l103 = c("Y_144", "Y_21"))

monops <-  pseudotime(cds_sub, reduction_method = "UMAP")
corrdata <- data.frame(cyto = cytops[names(monops)], mono = monops, cluster = cds_sub@colData$seurat_clusters)
corrdata <- corrdata[corrdata$mono != Inf,]
cor(corrdata$cyto, corrdata$mono)


lapply(groupid, function(x, corrdata){
  id <- choose_graph_segments(cds_sub, starting_pr_node = x[1], ending_pr_nodes = x[2], return_list = T)
  id <- is.element(rownames(corrdata), id$cells)
  dat <- data.frame(cyto = corrdata$cyto[id], mono = corrdata$mono[id], cols = corrdata$cluster[id])
  ggplot(dat, aes(x = mono, y = cyto)) + geom_point(shape = 20, colour = "grey") + 
    geom_smooth(method = "loess", span = 1, n = 200, fill = "lightblue") + theme_classic() + 
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(limits = c(0,1), expand = c(0,0))
  ggsave(paste0("l", paste(x, collapse = ""), ".pdf"),  width = 5, height = 4)
}, corrdata)










