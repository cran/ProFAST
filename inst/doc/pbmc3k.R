## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  set.seed(2024) # set a random seed for reproducibility.
#  library(Seurat)
#  pbmc3k <- SeuratData::LoadData("pbmc3k")
#  ## filter the seurat_annotation is NA
#  idx <- which(!is.na(pbmc3k$seurat_annotations))
#  pbmc3k <- pbmc3k[,idx]
#  pbmc3k

## ----eval = FALSE-------------------------------------------------------------
#  library(ProFAST) # load the package of FAST method
#  library(Seurat)

## ----eval = FALSE-------------------------------------------------------------
#  pbmc3k <- NormalizeData(pbmc3k)

## ----eval = FALSE-------------------------------------------------------------
#  pbmc3k <- FindVariableFeatures(pbmc3k)

## ----eval = FALSE, fig.width= 6, fig.height= 4.5------------------------------
#  dat_cor <- diagnostic.cor.eigs(pbmc3k)
#  q_est <- attr(dat_cor, "q_est")
#  cat("q_est = ", q_est, '\n')

## ----eval = FALSE-------------------------------------------------------------
#  pbmc3k <- NCFM(pbmc3k, q = q_est)
#  pbmc3k

## ----eval = FALSE-------------------------------------------------------------
#  pbmc3k <- pdistance(pbmc3k, reduction = "ncfm")

## ----eval = FALSE-------------------------------------------------------------
#  print(table(pbmc3k$seurat_annotations))
#  Idents(pbmc3k) <- pbmc3k$seurat_annotations
#  df_sig_list <- find.signature.genes(pbmc3k)
#  str(df_sig_list)

## ----eval = FALSE-------------------------------------------------------------
#  dat <- get.top.signature.dat(df_sig_list, ntop = 5, expr.prop.cutoff = 0.1)
#  head(dat)

## ----eval = FALSE, fig.width=10,fig.height=7----------------------------------
#  pbmc3k <- coembedding_umap(
#    pbmc3k, reduction = "ncfm", reduction.name = "UMAP",
#    gene.set = unique(dat$gene))

## ----eval = FALSE, fig.width=8,fig.height=5-----------------------------------
#  ## choose beutifual colors
#  cols_cluster <- c("black", PRECAST::chooseColors(palettes_name = "Light 13", n_colors = 9, plot_colors = TRUE))
#  p1 <- coembed_plot(
#     pbmc3k, reduction = "UMAP",
#     gene_txtdata = subset(dat, label=='B'),
#     cols=cols_cluster,pt_text_size = 3)
#  p1

## ----eval = FALSE, fig.width=8,fig.height=5-----------------------------------
#  p2 <- coembed_plot(
#     pbmc3k, reduction = "UMAP",
#     gene_txtdata = dat, cols=cols_cluster,
#     pt_text_size = 3)
#  p2

## ----eval = FALSE, fig.width=7,fig.height=4-----------------------------------
#  cols_type <- cols_cluster[-1]
#  names(cols_type)<-  sort(levels(Idents(pbmc3k)))
#  DimPlot(pbmc3k, reduction = 'UMAP', cols=cols_type)

## ----eval = FALSE, fig.width=8,fig.height=3.6---------------------------------
#  FeaturePlot(pbmc3k, reduction = 'UMAP', features = c("CD79A", "VPREB3"))

## -----------------------------------------------------------------------------
sessionInfo()

