
DataPath <- "/home/yh/2-nano/backup/scrnaseq_resub/data"
setwd(DataPath)
library(SummarizedExperiment)
library(MetaNeighbor)

mn_data <- readRDS("combine.se.rds")
mn_data$cell_type <- dplyr::recode(mn_data$cell_type,
                                   "R0"="Rgerm",
                                   "R1"="Repi",
                                   "R2"="Repi",
                                   "R3"="Rov",
                                   "R4"="Rmed",
                                   "R5"="Rsti",
                                   "R6"="Rmed",
                                   "R7"="Rmed",
                                   "R8"="Rcha",
                                   "R9"="Repi",
                                   "R10"="Rov")
var_genes = variableGenes(dat = mn_data, exp_labels = mn_data$study_id)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = mn_data,
                             study_id = mn_data$study_id,
                             cell_type = mn_data$cell_type)
top_hits = topHitsByStudy(celltype_NV, threshold = 0.6)


## plot
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
library(RColorBrewer)
idname <- sub("Rice[|]", "", rownames(celltype_NV))
idname <- sub("Arab[|]", "", idname)
rownames(celltype_NV) <- idname
colnames(celltype_NV) <- idname
gplots::heatmap.2(celltype_NV,
                  margins=c(9,9),
                  keysize=1,
                  key.xlab="AUROC",
                  key.title="NULL",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 1.5,
                  cexCol = 1.5) ## save 8*8

