## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  githubURL <- "https://github.com/feiyoung/FAST/blob/main/vignettes_data/seulist2_ID9_10.RDS?raw=true"
#  download.file(githubURL,"seulist2_ID9_10.RDS",mode='wb')
#  

## ----eval =FALSE--------------------------------------------------------------
#  dlpfc2 <- readRDS("./seulist2_ID9_10.RDS")
#  dlpfc <- dlpfc2[[1]]

## ----eval =FALSE--------------------------------------------------------------
#  library(ProFAST) # load the package of FAST method
#  library(PRECAST)
#  library(Seurat)

## ----eval =FALSE--------------------------------------------------------------
#  dlpfc ## a list including three Seurat object with default assay: RNA

## ----eval =FALSE--------------------------------------------------------------
#  print(row.names(dlpfc)[1:10])
#  count <- dlpfc[['RNA']]@counts
#  row.names(count) <- unname(transferGeneNames(row.names(count), now_name = "ensembl",
#                                                   to_name="symbol",
#                                                   species="Human", Method='eg.db'))
#  print(row.names(count)[1:10])
#  seu <- CreateSeuratObject(counts = count, meta.data = dlpfc@meta.data)
#  seu

## ----eval =FALSE--------------------------------------------------------------
#  ## Check the spatial coordinates: they are named as "row" and "col"!
#  print(head(seu@meta.data))

## ----eval =FALSE--------------------------------------------------------------
#  seu <- NormalizeData(seu)
#  seu <- FindVariableFeatures(seu)
#  print(seu)
#  print(seu[['RNA']]@var.features[1:10])

## ----eval = FALSE-------------------------------------------------------------
#  set.seed(2023)
#  seu <- DR.SC::FindSVGs(seu)
#  ### Check the results
#  print(seu[['RNA']]@var.features[1:10])

## ----eval =FALSE--------------------------------------------------------------
#  Adj_sp <- AddAdj(as.matrix(seu@meta.data[,c("row", "col")]), platform = "Visium")
#  ### set q= 15 here
#  set.seed(2023)
#  seu <- FAST_single(seu, Adj_sp=Adj_sp, q= 15, fit.model='poisson')
#  seu

## ----eval = FALSE-------------------------------------------------------------
#  Adj_sp <- AddAdj(as.matrix(seu@meta.data[,c("row", "col")]), platform = "Visium")
#  set.seed(2023)
#  seu <- FAST_single(dlpfc, Adj_sp=Adj_sp, q= 15, fit.model='gaussian')
#  ### Check the results
#  seu

## ----eval =FALSE--------------------------------------------------------------
#  ## Obtain the true labels
#  y <- seu$layer_guess_reordered
#  ### Evaluate the MacR2
#  Mac <-  get_r2_mcfadden(Embeddings(seu, reduction='fast'), y)
#  ### output them
#  print(paste0("MacFadden's R-square of FAST is ", round(Mac, 3)))

## ----eval =FALSE--------------------------------------------------------------
#  seu <- FindNeighbors(seu, reduction = 'fast')
#  seu <- FindClusters(seu, resolution = 0.4)
#  seu$fast.cluster <- seu$seurat_clusters
#  ARI.fast <- mclust::adjustedRandIndex(y, seu$fast.cluster)
#  print(paste0("ARI of PCA is ", round(ARI.fast, 3)))

## ----eval =FALSE--------------------------------------------------------------
#  seu <- ScaleData(seu)
#  seu <- RunPCA(seu, npcs=15, verbose=FALSE)
#  Mac.pca <- get_r2_mcfadden(Embeddings(seu, reduction='pca'), y)
#  print(paste0("MacFadden's R-square of PCA is ", round(Mac.pca, 3)))
#  set.seed(1)
#  seu <- FindNeighbors(seu, reduction = 'pca', graph.name ="pca.graph")
#  seu <- FindClusters(seu, resolution = 0.8,graph.name = 'pca.graph')
#  seu$pca.cluster <- seu$seurat_clusters
#  ARI.pca <- mclust::adjustedRandIndex(y, seu$pca.cluster)
#  print(paste0("ARI of PCA is ", round(ARI.pca, 3)))

## ----eval =FALSE--------------------------------------------------------------
#  cols_cluster <- chooseColors(palettes_name = "Nature 10", n_colors = 8, plot_colors = TRUE)

## ----eval =FALSE, fig.width=8, fig.height=3-----------------------------------
#  seu <- PRECAST::Add_embed(embed = as.matrix(seu@meta.data[,c("row", "col")]), seu, embed_name = 'Spatial')
#  seu
#  p1 <- DimPlot(seu, reduction = 'Spatial', group.by = 'pca.cluster',cols = cols_cluster, pt.size = 1.5)
#  p2 <- DimPlot(seu, reduction = 'Spatial', group.by = 'fast.cluster',cols = cols_cluster, pt.size = 1.5)
#  drawFigs(list(p1, p2),layout.dim = c(1,2) )

## ----eval =FALSE, fig.width=5, fig.height=4-----------------------------------
#  seu <- RunUMAP(seu, reduction = "fast", dims=1:15)
#  seu
#  DimPlot(seu, reduction='umap', group.by = "fast.cluster")

## ----eval =FALSE--------------------------------------------------------------
#  Idents(seu) <- seu$fast.cluster
#  dat_deg <- FindAllMarkers(seu)
#  library(dplyr)
#  n <- 5
#  dat_deg %>%
#      group_by(cluster) %>%
#      top_n(n = n, wt = avg_log2FC) -> top5
#  top5

## -----------------------------------------------------------------------------
sessionInfo()

