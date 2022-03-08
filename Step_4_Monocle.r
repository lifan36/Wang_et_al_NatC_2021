# Load data
lg343 <- readRDS('LG343_integrated_ready_4_DEGs.rds')
lg345 <- readRDS('LG345_integrated_ready_4_DEGs.rds')
micro_343 <- subset(lg343, subset = celltype == 'microglia')
micro_345 <- subset(lg345, subset = celltype == 'microglia')
rm(lg343)
rm(lg345)
project_list <- c(micro_343, micro_345)
all.gene <- rownames(micro_343)
all.gene <- subset(all.gene, subset = all.gene %in% rownames(micro_345))
anchors <- FindIntegrationAnchors(object.list = project_list, dims = 1:20)
micro <- IntegrateData(anchorset = anchors, dims = 1:20, features.to.integrate = all.gene)
micro <- subset(micro, subset = Condition == 'IKKbCtrl;P301S+' | Condition == 'IKKbKO;P301S+' | Condition == 'IKKbWT;P301S+'| Condition == 'IKKbCA;P301S+' | Condition == 'IKKbWT' | Condition == 'IKKbCtrl')
# Remove outlier
micro <- subset(micro, subset = orig.ident != 'IKKbWT;P301S+_2')
# rm(terc_146, terc_147, terc_f1, terc_f2, project_list, anchors)

micro@meta.data$genotype <- plyr::mapvalues(x = micro@meta.data$Condition,
                                            from = c("IKKbKO;P301S+", "IKKbCA;P301S+", "IKKbWT", "IKKbCtrl", "IKKbCtrl;P301S+", "IKKbWT;P301S+"),
                                            to = c("IKKbKO;P301S+", "IKKbCA;P301S+", "Ctrl", "Ctrl", "Ctrl;P301S+", "Ctrl;P301S+"))

## Subclustering of Microglia
micro <- subset(ikk, subset = celltype == 'microglia')
micro <- merge(micro_343, micro_345)
# rm(micro_343)
# rm(micro_345)
## Without any TAU
# micro <- subset(micro, subset = Condition == 'IKKbCA' | Condition == 'IKKbWT')
micro <- subset(micro, subset = Condition == 'IKKbCtrl;P301S+' | Condition == 'IKKbKO;P301S+')


DefaultAssay(micro) <- 'RNA'
micro <- NormalizeData(micro, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(micro)
micro <- ScaleData(micro, features = all.genes)
micro <- FindVariableFeatures(object = micro)
micro <- RunPCA(micro, features = VariableFeatures(object = micro))
DimPlot(micro, reduction = 'pca', group.by = 'Condition')
micro <- FindNeighbors(micro, dims = 1:20)
micro <- FindClusters(micro, resolution = 1.1)
micro <- RunUMAP(micro, dims = 1: 20, perplexity = 1)
Idents(micro) <- 'seurat_clusters'
DimPlot(micro, pt.size = 1, label = T, label.size = 8)
DimPlot(micro,split.by = 'Condition', pt.size = 1, label = T, label.size = 8)


ikk_seurat_micro_0 <- FindMarkers(micro, ident.1 = 0, test.use = 'wilcox')
ikk_seurat_micro_1 <- FindMarkers(micro, ident.1 = 1, test.use = 'wilcox')
ikk_seurat_micro_2 <- FindMarkers(micro, ident.1 = 2, test.use = 'wilcox')
write.csv(ikk_seurat_micro_0, file = 'ikk_seurat_micro_0.csv')
write.csv(ikk_seurat_micro_1, file = 'ikk_seurat_micro_1.csv')
write.csv(ikk_seurat_micro_2, file = 'ikk_seurat_micro_2.csv')

## Calculate the numbers of dimensions to be used in the following analysis.
micro <- JackStraw(micro, num.replicate = 100)
DimPlot(micro, pt.size = 1, ncol = 5, split.by = 'orig.ident', label = T, label.size = 5)
micro <- ScoreJackStraw(micro, dims = 1:20)
JackStrawPlot(micro, dims = 1:15)

# Convert seurat object to monocle3 object
DefaultAssay(micro) <- 'RNA'
cds <- as.cell_data_set(micro)
## Return to this point to remove non-microglial cells (cells too far from the main clusters)
cds <- choose_cells(cds)
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 9)
cds <- reduce_dimension(cds, reduction_method = 'UMAP')
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(cds)
cds <- cluster_cells(cds, reduction_method = 'UMAP', k = 9)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE, cell_size = 1, reduction_method = 'UMAP', group_label_size = 8) + theme(legend.position = 'right')
p2 <- plot_cells(cds, show_trajectory_graph = FALSE, cell_size = 1, reduction_method = 'UMAP', group_cells_by = 'partition', labels_per_group = 1, group_label_size = 4) + theme(legend.position = 'right')
p3 <- plot_cells(cds, show_trajectory_graph = FALSE, cell_size = 1, reduction_method = 'UMAP', color_cells_by = 'genotype', labels_per_group = 1, group_label_size = 4) + theme(legend.position = 'right')
p4 <- plot_cells(cds, show_trajectory_graph = FALSE, cell_size = 1, reduction_method = 'UMAP', color_cells_by = 'orig.ident', labels_per_group = 1, group_label_size = 4) + theme(legend.position = 'right')
wrap_plots(p1, p2, p3, p4)
cds_temp <- as.Seurat(cds, assay = 'RNA')
cds_temp <- subset(cds_temp, subset = monocle3_clusters != 6& monocle3_clusters != 7)
cds <- as.cell_data_set(cds_temp)
rm(cds_temp)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(cds)

# Plot monocle3 clusters
p1 <- plot_cells(cds, show_trajectory_graph = FALSE, cell_size = 1, reduction_method = 'UMAP', group_label_size = 8) + theme(legend.position = 'right')
p2 <- plot_cells(cds, show_trajectory_graph = FALSE, cell_size = 1, reduction_method = 'UMAP', group_cells_by = 'partition', labels_per_group = 1, group_label_size = 0) + theme(legend.position = 'right')
p3 <- plot_cells(cds, show_trajectory_graph = FALSE, cell_size = 1, reduction_method = 'UMAP', color_cells_by = 'Condition', labels_per_group = 0, group_label_size = 0) + theme(legend.position = 'right')
wrap_plots(p1, p2, p3)
plot_cells(cds, show_trajectory_graph = FALSE, cell_size = 0.5, reduction_method = 'UMAP', group_cells_by = 'partition',color_cells_by = 'Condition', labels_per_group = 1, group_label_size = 0) + 
  theme(legend.position = 'right') + facet_wrap(~Condition)  + theme(legend.position = 'right')

##
## Return to 'choose_cell' function to remove non-microglia cells
##

# Plot the composition of clusters
wt_mata <- cds@colData$Condition == 'IKKbWT'
wt_mata <- cds@clusters@listData$UMAP$clusters[wt_mata]
g3tercko_mata <- cds@colData$Condition == 'IKKbWT;P301S+'
g3tercko_mata <- cds@clusters@listData$UMAP$clusters[g3tercko_mata]
ikk_mata <- cds@colData$Condition == 'IKKbCA'
ikk_mata <- cds@clusters@listData$UMAP$clusters[ikk_mata]
ikkca_mata <- cds@colData$Condition == 'IKKbCA;P301S+'
ikkca_mata <- cds@clusters@listData$UMAP$clusters[ikkca_mata]
meta_frame <- as.data.frame(summary(wt_mata))
meta_frame$ ko <- summary(g3tercko_mata)
meta_frame$ikk <- summary(ikk_mata) 
meta_frame$ikk_ca <- summary(ikkca_mata)
colnames(meta_frame) <- c('dap12_wt', 'dap12_ko', 'ikk_wt', 'ikk_ca')
meta_frame$celltype <- rownames(meta_frame)
ratio_frame <- meta_frame
ikk_total <- sum(meta_frame$ikk_wt)
dap12_wt <- sum(meta_frame$dap12_wt)
dap12_ko <- sum(meta_frame$dap12_ko)
ikk_ca <- sum(meta_frame$ikk_ca)
ratio_frame$ikk_wt <- ratio_frame$ikk_wt / ikk_total
ratio_frame$dap12_wt <- ratio_frame$dap12_wt / dap12_wt
ratio_frame$dap12_ko <- ratio_frame$dap12_ko /dap12_ko
ratio_frame$ikk_ca <- ratio_frame$ikk_ca /ikk_ca
colnames(meta_frame) <- c('IKKbWT', 'IKKbWT;P301S+', 'IKKbCA', 'IKKbCA;P301S+', 'celltype')
colnames(ratio_frame) <- c('IKKbWT', 'IKKbWT;P301S+', 'IKKbCA', 'IKKbCA;P301S+', 'celltype')
bar_frame <- melt(meta_frame[,c('celltype', 'IKKbWT', 'IKKbWT;P301S+', 'IKKbCA', 'IKKbCA;P301S+')], id.vars = 1)
ratio_bar_frame <- melt(ratio_frame[,c('celltype', 'IKKbWT', 'IKKbWT;P301S+', 'IKKbCA', 'IKKbCA;P301S+')], id.vars = 1)
colnames(ratio_bar_frame) <- c('subcluster', 'variable', 'ratio')

cell_bar <- ggplot(ratio_bar_frame[order(ratio_bar_frame$subcluster),], aes(x = subcluster, y = ratio)) + 
  geom_bar(aes(fill = variable), stat="identity", position=position_dodge()) + theme(axis.text = element_text(size = 8))
cell_bar + scale_fill_brewer(palette="Paired") + theme_minimal() + theme(axis.text.x = element_text(angle = -45, hjust =0), axis.text = element_text(size =12)) + scale_x_discrete(name = 'Interneuron Subclusters', limits = c('1','2','3','4'))

###Next monocle steps: generate pseudotime 
cds <- learn_graph(cds, use_partition = F, learn_graph_control = list(minimal_branch_len = 10, orthogonal_proj_tip = F, geodesic_distance_ratio = 2))
p7 <- plot_cells(cds, label_groups_by_cluster = F, label_leaves = F, color_cells_by = 'genotype',label_branch_points = FALSE, cell_size = 1, group_label_size = 5) + theme(legend.position = 'right')

cds <- order_cells(cds)
p5 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5,
                 cell_size = 1)
p6 <- plot_cells(cds, show_trajectory_graph = T, cell_size = 1, reduction_method = 'UMAP', color_cells_by = 'genotype', trajectory_graph_segment_size = 0.75, graph_label_size = 3)
# Differnt cobinaitons of plot
wrap_plots(p1, p3, p5, p7)
wrap_plots(p3, p4)
wrap_plots(p3, p4, p5, p6)