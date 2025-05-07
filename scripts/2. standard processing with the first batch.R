# Load 10x Data and Create Seurat Objects
wt1.data <- Read10X_h5("data/2nd_batch/WT_Ys1/outs/filtered_feature_bc_matrix.h5")
wt2.data <- Read10X_h5("data/2nd_batch/WT_Ys2/outs/filtered_feature_bc_matrix.h5")
mut1.data <- Read10X_h5("data/2nd_batch/Mut_Ys1/outs/filtered_feature_bc_matrix.h5")
mut2.data <- Read10X_h5("data/2nd_batch/Mut_Ys2/outs/filtered_feature_bc_matrix.h5")
wt3.data <- Read10X_h5("data/WT_E8.75_Outs/Ctrl1/count/sample_feature_bc_matrix.h5")
wt4.data <- Read10X_h5("data/WT_E8.75_Outs/Ctrl2/count/sample_feature_bc_matrix.h5")
wt5.data <- Read10X_h5("data/WT_E8.75_Outs/Ctrl3/count/sample_feature_bc_matrix.h5")
wt3.data <- wt3.data[["Gene Expression"]]
wt4.data <- wt4.data[["Gene Expression"]]
wt5.data <- wt5.data[["Gene Expression"]]
class(wt1.data)
class(wt3.data)

wt1 <- CreateSeuratObject(counts = wt1.data, min.cells = 3, min.features = 200, project = "wt-1")
wt2 <- CreateSeuratObject(counts = wt2.data, min.cells = 3, min.features = 200, project = "wt-2")
mut1 <- CreateSeuratObject(counts = mut1.data, min.cells = 3, min.features = 200, project = "mut-1")
mut2 <- CreateSeuratObject(counts = mut2.data, min.cells = 3, min.features = 200, project = "mut-2")
wt3 <- CreateSeuratObject(counts = wt3.data, min.cells = 3, min.features = 200, project = "wt-3")
wt4 <- CreateSeuratObject(counts = wt4.data, min.cells = 3, min.features = 200, project = "wt-4")
wt5 <- CreateSeuratObject(counts = wt5.data, min.cells = 3, min.features = 200, project = "wt-5")

wt1.data <- NULL
wt2.data <- NULL
mut1.data <- NULL
mut2.data <- NULL
wt3.data <- NULL
wt4.data <- NULL
wt5.data <- NULL

# Combine datasets and add mitochondria percentage to meta
combined.2b <- merge(wt1, y = c(wt2, mut1, mut2, wt3, wt4, wt5))
combined.2b[["percent.mt"]] <- PercentageFeatureSet(combined.2b, pattern = "^mt-")

# Visualize QC metrics
VlnPlot(combined.2b, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p1 <- FeatureScatter(combined.2b, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(combined.2b, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1+p2

## quantile
quantile(combined.2b$nFeature_RNA, 0.02)
quantile(combined.2b$nFeature_RNA, 0.98)
quantile(combined.2b$nCount_RNA, 0.02)
quantile(combined.2b$nCount_RNA, 0.98)
## remove low quality 
combined.2b <- subset(combined.2b, subset = nCount_RNA >= 500 & nFeature_RNA >= 500 & nFeature_RNA < 8000 & percent.mt < 20)

# standard processing
combined.2b <- combined.2b %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = "percent.mt")

# PCA
combined.2b <- combined.2b %>% 
  RunPCA(features = VariableFeatures(combined.2b))

DimHeatmap(combined.2b, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(combined.2b, dims = 16:30, cells = 500, balanced = TRUE)
ElbowPlot(combined.2b, ndims = 30)

dimensions <- 1:15
resolutions <- 0.1

# clustering and dimension reduction
combined.2b <- combined.2b %>% 
  FindNeighbors(dims = dimensions) 
combined.2b <- combined.2b %>% 
  FindClusters(resolution = resolutions) 
combined.2b <- combined.2b %>% 
  RunUMAP(dims = dimensions)
DimPlot(combined.2b, reduction = "umap", cols = my36colors)
DimPlot(combined.2b, reduction = "umap", cols = my36colors, split.by = "orig.ident")
idx <- rownames(combined.2b) %>%  str_detect("IRES")
rownames(combined.2b)[idx]
# batch correction and clustering
combined.2b <- combined.2b %>%  
  IntegrateLayers(
    method = CCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.cca", k.anchor = 5,
    verbose = FALSE
  )

combined.2b <- combined.2b %>% 
  IntegrateLayers(
    method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    verbose = FALSE
  )

combined.2b <- combined.2b %>% 
  FindNeighbors(reduction = "integrated.cca", dims = dimensions)
combined.2b <- combined.2b %>% 
  FindClusters(resolution = resolutions, cluster.name = "cca_clusters")
combined.2b <- combined.2b %>% 
  FindNeighbors(reduction = "harmony", dims = dimensions)
combined.2b <- combined.2b %>% 
  FindClusters(resolution = resolutions, cluster.name = "harmony_clusters")
combined.2b <- combined.2b %>% 
  RunUMAP(reduction = "harmony", dims = dimensions, reduction.name = "umap.harmony")

p1 <- DimPlot(
  combined.2b,
  reduction = "umap.harmony",
  group.by = c("harmony_clusters"),
  split.by = "orig.ident",
  cols = my36colors,
  combine = FALSE, label.size = 2
)
p1

FeaturePlot(combined.2b, "Sp8", split.by = "orig.ident",  reduction = "umap.harmony")
DimPlot(combined.2b, reduction = "umap", cols = my36colors)
DimPlot(combined.2b, reduction = "umap", cols = my36colors, split.by = "orig.ident", group.by = "cca_clusters")
DimPlot(combined.2b, reduction = "umap", cols = my36colors, split.by = "orig.ident", group.by = "harmony_clusters")

# collapses the individual datasets together and recreates the original counts and data layers
combined.2b <- combined.2b %>% JoinLayers()

# group by genotype
idx <- combined.2b@meta.data$orig.ident %>% str_detect("wt")
combined.2b@meta.data$genotype <- "mut"
combined.2b@meta.data$genotype[idx] <- "wt"

# marker vis.
ExE.ectoderm <- c("Tfap2c","Rhox5","Epcam","Perp","P4ha2","Trap1a","Fabp3","Tex19.1","Elf5","Ass1")
PGC <- c("Pou5f1","Tfap2c","Rhox5","Epcam")
YsVE <- c("Epcam","Cdh1","Rhox5","Sox17","Afp","Hoxd1","Cldn2","Slc2a2","Ttr","Apoa1","Apoa2","Dab2","S100g","Apob","Apoa4")
Endothelium.progenitors <- c("Lyve1","Ramp2","Plvap","Cdh5", "Etv2","Swap70","Hhex","Clec14a")
Mesothelium <- c("Col1a1","Col3a1","Lum","Tdo2")
Erythroid <- c("Hspe1","Hba-a2","Alas2","Gypa","Hbb-y")

DimPlot(combined.2b, reduction = "umap.harmony", cols = my36colors, split.by = "genotype")

DefaultAssay(combined.2b)

DoHeatmap(combined.2b, features = Erythroid) #0,2
DoHeatmap(combined.2b, features = Mesothelium) #1
DoHeatmap(combined.2b, features = Endothelium.progenitors) #6
DoHeatmap(combined.2b, features = ExE.ectoderm) #3
DoHeatmap(combined.2b, features = YsVE) #5

## annotation

Idents(combined.2b) <- "harmony_clusters"
combined.2b <- RenameIdents(combined.2b, 
                            "0" = "Erythroid", 
                            "1" = "Mesothelium", 
                            "2" = "Erythroid_mut", 
                            "3" = "ExE.ectoderm",
                            "4" = "Erythroid", 
                            "5" = "YsVE", 
                            "6" = "Endothelium.progenitors")
combined.2b$cell.type <- Idents(combined.2b)
DimPlot(combined.2b, reduction = "umap.harmony", cols = my36colors, split.by = "genotype")


## cell cycle scoring ####
# function to convert human to mouse gene names

# extract cell cylce genes
g2m_genes <- read_csv("gene_sets/g2m_genes.csv")
g2m.genes<- g2m_genes$x
s_genes <- read_csv("gene_sets/s_genes.csv")
s.genes<- s_genes$x
# scoring
combined.2b <- CellCycleScoring(combined.2b, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# visual.

# Generate the violin plot
combined.2b@meta.data %>% head()
VlnPlot(combined.2b, 
        cols = my36colors, 
        group.by = "cell.type.genotype", 
        features = c("S.Score", "G2M.Score"),
        split.plot = F) &
  stat_summary(fun = mean, geom="point", size = 10, colour = "white", shape = 95)


## cell type numbers
combined.2b@meta.data$cell.type.genotype <- paste(combined.2b@meta.data$cell.type, combined.2b@meta.data$genotype, sep = ".")
table(combined.2b@meta.data$cell.type.genotype)

DoHeatmap(combined.2b, c(Erythroid, Mesothelium, Endothelium.progenitors, ExE.ectoderm, YsVE, PGC), group.by = "cell.type") 

## DEG (findmarkers) OE vs WT pair-wise comparison ------------------
DefaultAssay(combined.2b) 
## split idents
Idents(combined.2b) <- "genotype"
cell.type.list <-  SplitObject(combined.2b, split.by = "cell.type")
oeVSwt.DEGlist <- cell.type.list %>% lapply(FindMarkers, ident.1 = "mut", ident.2 = "wt", min.pct = 0.3)
saveRDS(oeVSwt.DEGlist, "rds/oeVSwt.DEGlist.rds")

## saveRDS two batches combined.
saveRDS(combined.2b, "rds/combined.2b.rds")
combined.2b <- readRDS("rds/combined.2b.rds")

#### Pseudobulk ----
# collapses the individual datasets together and recreates the original counts and data layers
combined.2b <- combined.2b %>% JoinLayers()
# Preparing the single-cell dataset for pseudobulk analysis
table(combined.2b@meta.data$cell.type, combined.2b@meta.data$orig.ident)
## make pseudo replicate for mutant to make n=2.
idx.mut.ExE <- combined.2b@meta.data$cell.type == "ExE.ectoderm" & combined.2b@meta.data$genotype == "mut"
combined.2b@meta.data$orig.ident.pseudo <- combined.2b@meta.data$orig.ident
total.num.ExE.mut <- combined.2b@meta.data$orig.ident.pseudo[idx.mut.ExE] %>% length()
rep.num <- 2
random_integer <- sample(1:total.num.ExE.mut, total.num.ExE.mut/rep.num)
combined.2b@meta.data$orig.ident.pseudo[idx.mut.ExE][random_integer] <- "mut-1"
DimPlot(combined.2b, reduction = "umap.harmony", cols = my36colors, 
        group.by = "cell.type", split.by = "orig.ident.pseudo")

## pseudobulk DEG
DEG.pseudobulk <- run_de(combined.2b,   
                         meta = meta,   
                         cell_type_col = "cell.type",  # Cluster identity  
                         label_col = "genotype",                # Experimental condition (OE vs ctrl)  
                         replicate_col = "orig.ident.pseudo",        # Sample ID  
                         de_family = 'pseudobulk',            # Pseudobulk mode  
                         de_method = "edgeR",                 # Use edgeR  
                         de_type = "QLF"                      # Quasi-Likelihood F-test
)

DEG.pseudobulk %>% write.csv("markers/DEG.pseudobulk.csv")
