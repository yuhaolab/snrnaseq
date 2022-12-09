
DataPath <- "/home/yh/2-nano/backup/scrnaseq_resub/data"
setwd(DataPath)
library(Seurat)
library(monocle3)

pbmc <- readRDS("pbmc_integ.rds")
counts <- GetAssayData(object = pbmc, assay = "RNA", slot = "counts")
cell.metadata <- pbmc[[]]
feature.metadata <- pbmc[["RNA"]][[]]
feature.metadata$gene_short_name <- rownames(x = feature.metadata)

pd <- new(Class = "AnnotatedDataFrame", data = cell.metadata)
fd <- new(Class = "AnnotatedDataFrame", data = feature.metadata)

cds <- new_cell_data_set(counts,
                         cell_metadata = cell.metadata,
                         gene_metadata = feature.metadata)

saveRDS(cds, file = "cds_comb.rds")

## process
cds <- preprocess_cds(cds, num_dim = 60)
cds <- align_cds(cds, alignment_group = "orig.ident")
cds <- reduce_dimension(cds)

plot_cells(cds, group_label_size = 6, cell_size = 0.5, color_cells_by = "seurat_clusters")
plot_cells(cds, group_label_size = 6, cell_size = 0.5, color_cells_by = "orig.ident")

## get partition
cds <- cluster_cells(cds)
cds@colData$part <- partitions(cds, reduction_method = "UMAP")
plot_cells(cds, color_cells_by = "partition", label_cell_groups = FALSE)


## plot color
library(ggplot2)
library(RColorBrewer)
library(scales)
## cell type
colData(cds)$assigned_cell_type <- as.character(colData(cds)$seurat_clusters)
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
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
                                                 "10"="ov",
                                                 "11"="sta",
                                                 "12"="C12",
                                                 "13"="C13")
cols <- brewer.pal(6, "Set2")
names(cols) <- c("cha", "epi", "germ", "med", "ov", "sti")
cols2 <- brewer.pal(4, "Paired")[c(1,2,4)]
names(cols2) <- c("sta", "C13", "C12")
cols <- c(cols, cols2)
cols <- cols[c("sti", "epi", "med", "cha", "ov", "germ", "sta", "C12", "C13")]
plot_cells(cds, label_cell_groups = F, cell_size = 0.8, color_cells_by = "assigned_cell_type") +
  scale_color_manual(values = cols)


## get colore
cols <- brewer.pal(11, "Spectral")
cols2 <- brewer.pal(11, "RdBu")
cols <- c(cols2[1], cols, cols2[11])
pal <- colorRampPalette(cols)
cols <- pal(14)
cols[8] <- "#E5DAAE"
cols <- cols[c(2:14, 1)]
names(cols) <- as.character(c(1, 2, 4, 5, 6, 7, 9, 0, 3, 8, 10, 11, 12, 13))
cols[c("7", "9")] <- cols[c("9", "7")]

plot_cells(cds,
           color_cells_by = "seurat_clusters",
           show_trajectory_graph = F,
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=4,
           cell_size = 1,
           trajectory_graph_segment_size = 1) +
  scale_color_manual(values = cols)


## sub1 ---------------------------------------------------------------------------------------------------------------------------
## plot sub traj 124567
cds_sub <- cds
valid_cells <- row.names(subset(pData(cds_sub), is.element(seurat_clusters, c(1,2,4,5,6,7))))
cds_sub <- cds_sub[,valid_cells]

cds_sub <- preprocess_cds(cds_sub, num_dim = 50)
cds_sub <- align_cds(cds_sub, alignment_group = "orig.ident")
cds_sub <- reduce_dimension(cds_sub)
plot_cells(cds_sub, group_label_size = 6, cell_size = 0.8, color_cells_by = "seurat_clusters")

## plot color
cols1 <- cols[as.character(c(1,2,4,5,6,7))]
plot_cells(cds_sub,
           show_trajectory_graph = F,
           color_cells_by = "seurat_clusters",
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=4,
           cell_size = 1,
           trajectory_graph_segment_size = 1) +
  scale_color_manual(values = cols1) ## save 6*7.2

## plot traj
cds_sub <- cluster_cells(cds_sub)
plot_cells(cds_sub, color_cells_by = "partition")
cds_sub <- learn_graph(cds_sub, use_partition = F, close_loop = F, learn_graph_control = list(minimal_branch_len = 15))
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

cols1 <- cols[as.character(c(0,3,8,10))]
plot_cells(cds_sub,
           show_trajectory_graph = F,
           color_cells_by = "seurat_clusters",
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=4,
           cell_size = 1,
           trajectory_graph_segment_size = 1) +
  scale_color_manual(values = cols1) ## save 6*7.2

cds_sub <- order_cells(cds_sub, root_pr_nodes = "Y_137")
plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=4,
           cell_size = 1,
           label_principal_points = TRUE,
           trajectory_graph_segment_size = 1) ## save 6*7


## graph-based DE genes --------------------------------------------------------------------------------------------------------------------------------
## DE test
# ## sub
cds_sub <- cds
valid_cells <- row.names(subset(pData(cds_sub), is.element(seurat_clusters, c(0,1,2,3,4,5,6,7,8,9,10))))
cds <- cds_sub[,valid_cells]
plot_cells(cds, group_label_size = 6, cell_size = 0.8, color_cells_by = "seurat_clusters")

## DE test
cds <- cluster_cells(cds)
cds <- learn_graph(cds, close_loop = F, learn_graph_control = list(rann.k = 25, minimal_branch_len = 15, ncenter = 200))
plot_cells(cds, label_principal_points = TRUE)

subset_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=32)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))



