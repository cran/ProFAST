## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  set.seed(2024) # set a random seed for reproducibility.
#  library(ProFAST) # load the package of FAST method
#  data(CosMx_subset)
#  CosMx_subset

## ----eval = FALSE-------------------------------------------------------------
#  library(Seurat)

## ----eval = FALSE-------------------------------------------------------------
#  CosMx_subset <- NormalizeData(CosMx_subset)

## ----eval = FALSE-------------------------------------------------------------
#  CosMx_subset <- FindVariableFeatures(CosMx_subset)

## ----eval = FALSE, fig.width= 6, fig.height= 4.5------------------------------
#  dat_cor <- diagnostic.cor.eigs(CosMx_subset)
#  q_est <- attr(dat_cor, "q_est")
#  cat("q_est = ", q_est, '\n')

## ----eval = FALSE-------------------------------------------------------------
#  pos <- as.matrix(CosMx_subset@meta.data[,c("x", "y")]) # Extract the spatial coordinates
#  Adj_sp <- AddAdj(pos) ## calculate the adjacency matrix
#  CosMx_subset <- NCFM_fast(CosMx_subset, Adj_sp = Adj_sp, q = q_est)
#  CosMx_subset

## ----eval = FALSE-------------------------------------------------------------
#  CosMx_subset <- pdistance(CosMx_subset, reduction = "fast")

## ----eval = FALSE-------------------------------------------------------------
#  print(table(CosMx_subset$cell_type))
#  Idents(CosMx_subset) <- CosMx_subset$cell_type
#  df_sig_list <- find.signature.genes(CosMx_subset)
#  str(df_sig_list)

## ----eval = FALSE-------------------------------------------------------------
#  dat <- get.top.signature.dat(df_sig_list, ntop = 2, expr.prop.cutoff = 0.1)
#  head(dat)

## ----eval = FALSE, fig.width=10,fig.height=7----------------------------------
#  CosMx_subset <- coembedding_umap(
#    CosMx_subset, reduction = "fast", reduction.name = "UMAP",
#    gene.set = unique(dat$gene))

## ----eval = FALSE, fig.width=8,fig.height=5-----------------------------------
#  ## choose beutifual colors
#  cols_cluster <- c("black", PRECAST::chooseColors(palettes_name = "Blink 23", n_colors = 21, plot_colors = TRUE))
#  p1 <- coembed_plot(
#     CosMx_subset, reduction = "UMAP",
#     gene_txtdata = subset(dat, label=='tumor 5'),
#     cols=cols_cluster, pt_text_size = 3)
#  p1

## ----eval = FALSE, fig.width=9,fig.height=6-----------------------------------
#  p2 <- coembed_plot(
#     CosMx_subset, reduction = "UMAP",
#     gene_txtdata = dat, cols=cols_cluster,
#     pt_text_size = 3, alpha=0.2)
#  p2

## ----eval = FALSE, fig.width=9,fig.height=6-----------------------------------
#  cols_type <- cols_cluster[-1]
#  names(cols_type)<-  sort(levels(Idents(CosMx_subset)))
#  DimPlot(CosMx_subset, reduction = 'UMAP', cols=cols_type)

## ----eval = FALSE, fig.width=8,fig.height=3.6---------------------------------
#  FeaturePlot(CosMx_subset, reduction = 'UMAP', features = c("PSCA", "CEACAM6"))

## -----------------------------------------------------------------------------
sessionInfo()

