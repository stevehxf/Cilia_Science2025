# Functional analysis option 2. use FindMarker DEGs ----
#
oeVSwt.DEGlist <- readRDS("rds/oeVSwt.DEGlist.rds")
##convert findmarker deg table into a tibble. And rank the genes by foldchange
oeVSwt.DEG.tb.list <- oeVSwt.DEGlist %>% 
  lapply(rownames_to_column,"gene") %>% 
  lapply(as_tibble) %>% 
  lapply(arrange, desc(avg_log2FC))
# symbols will be converted to ENTREZID.
# We will lose some genes here because not all IDs will be converted and genes will be sorted based on fold change
#keytypes(org.Mm.eg.db) 
## create an empty list to store the conversion results
deg.ID.list <- oeVSwt.DEG.tb.list %>% lapply(function(x) NULL)
deg.rank.list <- deg.ID.list

for (i in names(oeVSwt.DEG.tb.list)) {
  
  ids <- bitr(oeVSwt.DEG.tb.list[[i]]$gene, 
              fromType = c("ALIAS"),
              toType = c("ENTREZID"), 
              OrgDb = "org.Mm.eg.db")
  ## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
  ## so we need to remove duplicate IDs prior to assessing enriched GO terms
  non_duplicates <- which(duplicated(ids$ALIAS) == FALSE)
  
  ids <- ids[non_duplicates, ] %>% as.tibble()
  colnames(ids) <- c("gene","ENTREZID")
  
  ## Merge the Ensembl IDs with the results     
  print("Merging the Ensembl IDs with the results")
  deg.ID.list[[i]] <- right_join(x=oeVSwt.DEG.tb.list[[i]], y=ids, by = "gene") %>% as.data.frame()
  
  ## Remove any NA values
  print("Removing any NA values")
  deg.ID.list[[i]] <- base::subset(deg.ID.list[[i]], deg.ID.list[[i]]$ENTREZID != "NA")
  deg.ID.list[[i]] %>% head() %>% print()
  ## Remove any duplicates
  print("Removing any duplicates")
  deg.ID.list[[i]] <- deg.ID.list[[i]][which(duplicated(deg.ID.list[[i]]$ENTREZID) == F), ]
  deg.rank.list[[i]] <- deg.ID.list[[i]]$avg_log2FC
  names(deg.rank.list[[i]]) <- deg.ID.list[[i]]$ENTREZID
  deg.rank.list[[i]] <- sort(deg.rank.list[[i]], decreasing = TRUE)
}

## KEGG pathway gene set enrichment analysis
kegg.list <- deg.rank.list
for (i in names(kegg.list)){
  kegg.list[[i]] <- gseKEGG(geneList     = deg.rank.list[[i]],
                            organism     = 'mmu', # house mouse
                            minGSSize    = 120,
                            pvalueCutoff = 0.1,
                            verbose      = FALSE)
  write.csv(kegg.list[[i]]@result, file = paste0("gsea/kegg_",i,".csv"))
}

## GO GSEA
go.list <- deg.rank.list
for (i in names(go.list)){
  go.list[[i]] <- gseGO(geneList     = deg.rank.list[[i]],
                        OrgDb = org.Mm.eg.db, 
                        ont = 'BP', 
                        nPerm = 1000, 
                        minGSSize = 120, 
                        pvalueCutoff = 0.1,
                        verbose = FALSE) 
  write.csv(go.list[[i]]@result, file = paste0("gsea/GO_BP_",i,".csv"))
}

for (i in names(go.list)){
  go.list[[i]] <- gseGO(geneList     = deg.rank.list[[i]],
                        OrgDb = org.Mm.eg.db, 
                        ont = 'CC', 
                        nPerm = 1000, 
                        minGSSize = 120, 
                        pvalueCutoff = 0.1,
                        verbose = FALSE) 
  write.csv(go.list[[i]]@result, file = paste0("gsea/GO_CC_",i,".csv"))
}

for (i in names(go.list)){
  go.list[[i]] <- gseGO(geneList     = deg.rank.list[[i]],
                        OrgDb = org.Mm.eg.db, 
                        ont = 'MF', 
                        nPerm = 1000, 
                        minGSSize = 120, 
                        pvalueCutoff = 0.1,
                        verbose = FALSE) 
  write.csv(go.list[[i]]@result, file = paste0("gsea/GO_MF_",i,".csv"))
}
## custom gene set

## run GSEA on custome gene sets.####
# create custome gene set list
custom.genes <- read.csv("gene_sets/GSEA genelist_SP58boundCilia.csv")
#View(custom.genes)
custom.gene.list <- vector("list", ncol(custom.genes))
names(custom.gene.list) <- colnames(custom.genes)
for (name in names(custom.gene.list)) {
  custom.gene.list[[name]] <- data.frame(
    pathwayName = name,
    geneSymbol = custom.genes[[name]]
  )
  custom.gene.list[[name]] <- 
    custom.gene.list[[name]][custom.gene.list[[name]]$geneSymbol != "", ]
}
custom.geneset.df <- bind_rows(custom.gene.list) # convert the list into one dataframe.
# get ranked gene list.
deg.gene.symbol.list <- oeVSwt.DEG.tb.list %>% lapply(function(x) NULL)
for(i in names(oeVSwt.DEG.tb.list)){
  deg.gene.symbol.list[[i]] <- oeVSwt.DEG.tb.list[[i]]$avg_log2FC
  names(deg.gene.symbol.list[[i]]) <- oeVSwt.DEG.tb.list[[i]]$gene
  deg.gene.symbol.list[[i]] <- sort(deg.gene.symbol.list[[i]], decreasing = TRUE)
}

gsea.ChIP.list <- deg.gene.symbol.list %>% lapply(function(x) NULL)
# run GSEA
for (i in names(deg.gene.symbol.list)){
  
  gsea.ChIP.list[[i]] <- GSEA(geneList = deg.gene.symbol.list[[i]],
                              TERM2GENE = custom.geneset.df,
                              pvalueCutoff = 0.1) 
  
}

View(gsea.ChIP.list)

for (name in names(gsea.ChIP.list)) {
  # visualization and export
  if (nrow(gsea.ChIP.list[[name]])>0) {
    print("plotting...")
    plot <- gseaplot2(gsea.ChIP.list[[name]], geneSetID = 1:nrow(gsea.ChIP.list[[name]]))
    print("exporting...")
    path <- paste0("plots/pathways/gsea.",name,".pdf")
    path <- paste0("plots/pathways/gsea.",name,".svg")
    ggsave(path, plot = plot, width = 5, height = 5)
  } 
}



