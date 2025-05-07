## monocle3 trajectory ---------------------------------------------------------------------------
library(monocle3)
combined.2b.mn3 <- combined.2b %>% 
  RunUMAP(reduction = "harmony", dims = dimensions, reduction.name = "UMAP")
DimPlot(combined.2b.mn3, reduction = "UMAP")
# convert to monocle3 object
combined.2b.mono3cds <- as.cell_data_set(combined.2b.mn3)
# cluster partitions
combined.2b.mono3cds <- monocle3::cluster_cells(combined.2b.mono3cds, reduction_method = "UMAP", k=10)
combined.2b.mono3cds <- monocle3::learn_graph(combined.2b.mono3cds, use_partition = TRUE)
combined.2b.mono3cds <- monocle3::order_cells(combined.2b.mono3cds, reduction_method = "UMAP")
p1 <- plot_cells(combined.2b.mono3cds, show_trajectory_graph = TRUE)
p2 <- plot_cells(combined.2b.mono3cds, color_cells_by = "partition", show_trajectory_graph = FALSE, label_cell_groups = TRUE)
wrap_plots(p1, p2)

monocle3::plot_cells(combined.2b.mono3cds, 
                     color_cells_by = "pseudotime", 
                     label_branch_points = F, 
                     label_cell_groups = F,
                     label_leaves = F,
                     label_groups_by_cluster = F,
                     label_roots = F)
ggsave("plots/pseudotime.trajectory.umap.svg", width = 4, height = 3)
saveRDS(combined.2b.mono3cds, "rds/combined.2b.mono3cds.rds")

monocle3::plot_cells(combined.2b.mono3cds, color_cells_by = "cell.type")
monocle3::plot_cells(combined.2b.mono3cds, color_cells_by = "partition")

# add pseudotime value to the orginal seurat.

pseudotime(combined.2b.mono3cds)
pseudo.time.df <- data.frame(row.names = names(pseudotime(combined.2b.mono3cds)), pseudotime(combined.2b.mono3cds))
head(pseudo.time.df)
names(pseudo.time.df) <- "pseudotime"
combined.2b %>% AddMetaData(pseudo.time.df) -> combined.2b
saveRDS(combined.2b, "rds/combined.2b.rds")
FeaturePlot(combined.2b, "pseudotime", split.by = "genotype",  reduction = "umap.harmony")
