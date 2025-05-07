## server setup ####
library(tidyverse)
library(Seurat)
#library(SoupX)
library(DoubletFinder)
#library(monocle)
#library(monocle3)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
#library(annotables)
library(enrichplot)
#library(DOSE)
library(RColorBrewer)
#library(velocyto.R)
#library(ggrepel)
library(SeuratWrappers)
library(cowplot)
library(patchwork)
library(viridis)
library(pheatmap)
library(SingleCellExperiment)
library(Libra) # pseudobulk
library(RcisTarget)
library(AUCell)
library(GENIE3)
library(SCENIC)
#devtools::load_all("C:/Users/zhoulab/AppData/Local/R/win-library/4.2/monocle/")

## in-house functions
#source("scripts/helper_function.r")

## colors of collection ------------------------------------------------------
my36colors <- c("#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87",
                "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658",
                "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398",
                "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B",
                "#712820", "#DCC1DD", "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963",
                "#968175")
my9colors <- c("#fdbf6f", # yellow
               "#a6cee3", # soft blue
               "#0096c7", # cyan
               "#b2df8a", # lightgreen
               "#33a02c", # darkgreen
               "#cab2d6", # grayish violet
               "#fb9a99", # pink
               "#e63946", # brightred
               "#6a3d9a"  # moderate violet
)
my5colors <- c("#ff595e", # red
               "#1982c4", # blue
               "#ffca3a", # yellow
               "#8ac926", # green
               "#6a4c93" # purple
)

#
grad <- colorRampPalette(RColorBrewer::brewer.pal(5,"RdPu"))(50) 
hmcols <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100)
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(100)


