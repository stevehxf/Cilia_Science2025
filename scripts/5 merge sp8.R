

# aggregate endogenous and transgentic Sp8
combined.2b.sp8merge <- combined.2b

combined.2b.sp8merge[["RNA"]]$counts["Sp8", ] <- 
  combined.2b.sp8merge[["RNA"]]$counts["Sp8", ] +
  combined.2b.sp8merge[["RNA"]]$counts["SP8-IRES2-EGFP-SV40polyA", ]



# remove original SP8-IRES2-EGFP-SV40polyA
keep <- rownames(combined.2b.sp8merge) %in% "SP8-IRES2-EGFP-SV40polyA"
keep <- rownames(combined.2b.sp8merge)[!keep]
combined.2b.sp8merge <- subset(combined.2b.sp8merge, features = keep)
rownames(combined.2b.sp8merge) %in% "SP8-IRES2-EGFP-SV40polyA" %>% any()

# renormalization
combined.2b.sp8merge[["RNA"]] <- split(combined.2b.sp8merge[["RNA"]], f = combined.2b.sp8merge$orig.ident)

combined.2b.sp8merge <- combined.2b.sp8merge %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = "percent.mt") %>% 
  RunPCA(features = VariableFeatures(combined.2b.sp8merge))


dimensions <- 1:15
resolutions <- 0.1

# clustering and dimension reduction
combined.2b.sp8merge <- combined.2b.sp8merge %>% 
  FindNeighbors(dims = dimensions)  %>% 
  FindClusters(resolution = resolutions) %>% 
  RunUMAP(dims = dimensions)
DimPlot(combined.2b.sp8merge, group.by = "cell.type")

# Batch correction
combined.2b.sp8merge <- combined.2b.sp8merge %>% 
  IntegrateLayers(
    method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    verbose = FALSE
  )

combined.2b.sp8merge <- combined.2b.sp8merge %>% 
  FindNeighbors(reduction = "harmony", dims = dimensions) %>% 
  FindClusters(resolution = resolutions, cluster.name = "harmony_clusters") %>% 
  RunUMAP(reduction = "harmony", dims = dimensions, reduction.name = "umap.harmony")


p1 <- DimPlot(
  combined.2b.sp8merge,
  reduction = "umap.harmony",
  group.by = c("cell.type"),
  split.by = "orig.ident",
  cols = my36colors,
  combine = FALSE, label.size = 2
)
p1

combined.2b.sp8merge <- JoinLayers(combined.2b.sp8merge)

FeaturePlot(combined.2b.sp8merge, "Sp8", split.by = "genotype",  reduction = "umap.harmony", max.cutoff = 2)
saveRDS(combined.2b.sp8merge, "rds/combined.2b.sp8merge.rds")
