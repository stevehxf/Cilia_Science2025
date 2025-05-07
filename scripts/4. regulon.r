setwd("~/2025_Yinwen_SP8OE")
## load the dataset
combined.2b <- readRDS("rds/combined.2b.rds")
dir.create("SCENIC")

# Analyze regulon in seurat -----
# load AUC scores per cell
AUCmatrix <- read.csv(file = "data/Regulon/merge_sce_SCENIC.csv", row.names = "Cell")
rownames(AUCmatrix) 
RegulonName_AUC <- colnames(AUCmatrix)
# Remove dots ONLY at the end of strings
RegulonName_AUC <- str_remove(RegulonName_AUC, "\\.+$")
# add auc_ prefix
colnames(AUCmatrix) <- paste0("auc_", RegulonName_AUC)

## anlysis AUC with seurat  ----
#setwd("/home/steveexternal/ginsc")
combined.2b.AUC <- AddMetaData(combined.2b, AUCmatrix)
combined.2b.AUC[['AUC']] <- CreateAssayObject(counts = t(AUCmatrix))
combined.2b.AUC[["RNA"]] 

#### rds ----
saveRDS(combined.2b.AUC, "rds/combined.2b.AUC.rds")

#### regulon auc visualization-------------
TF <- c("Sp5", "Sp8", "Sp2", "Foxj1", "Trp73", "Myb", "Rfx2", "Rfx3", "Rfx7", "Zfp423", "E2f7")
idx <- TF %in% RegulonName_AUC
auc.name <- paste0("auc_", TF)

VlnPlot(combined.2b.AUC, features = auc.name[idx], group.by = "cell.type", split.by = "genotype", stack = T, flip = T)
VlnPlot(combined.2b.AUC, "auc_Sp2")

VlnPlot(combined.2b.AUC, "SP8-IRES2-EGFP-SV40polyA", group.by = "cell.type", split.by = "genotype")

FeaturePlot(
  combined.2b.AUC,
  reduction = "umap.harmony",
  split.by = "genotype",
  features = "auc_Sp2",
  max.cutoff = 0.05)

FeaturePlot(
  combined.2b.AUC,
  reduction = "umap.harmony",
  split.by = "genotype",
  features = "auc_Myb")

FeaturePlot(
  combined.2b.AUC,
  reduction = "umap.harmony",
  split.by = "genotype",
  features = "auc_Rfx2",
  max.cutoff = 0.05)

FeaturePlot(
  combined.2b.AUC,
  reduction = "umap.harmony",
  split.by = "genotype",
  features = "auc_E2f7",
  max.cutoff = 0.03)

