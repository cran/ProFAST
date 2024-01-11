## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  githubURL <- "https://github.com/feiyoung/ProFAST/blob/main/vignettes_data/simu3.rds?raw=true"
#  download.file(githubURL,"simu3.rds",mode='wb')
#  

## ----eval = FALSE-------------------------------------------------------------
#  load("simu3.rds")

## ----eval = FALSE-------------------------------------------------------------
#  library(ProFAST) # load the package of FAST method
#  library(PRECAST)
#  library(Seurat)

## ----eval = FALSE-------------------------------------------------------------
#  simu3 ## a list including three Seurat object with default assay: RNA

## ----eval= FALSE--------------------------------------------------------------
#  head(simu3[[1]])

## ----eval= FALSE--------------------------------------------------------------
#  row.names(simu3[[1]])[1:10]

## ----eval = FALSE-------------------------------------------------------------
#  ## Get the gene-by-spot read count matrices
#  countList <- lapply(simu3, function(x) x[["RNA"]]@counts)
#  
#  ## Check the spatial coordinates: Yes, they are named as "row" and "col"!
#  head(simu3[[1]]@meta.data)
#  
#  ## Get the meta data of each spot for each data batch
#  metadataList <- lapply(simu3, function(x) x@meta.data)
#  
#  
#  ## ensure the row.names of metadata in metaList are the same as that of colnames count matrix in countList
#  M <- length(countList)
#  for(r in 1:M){
#    row.names(metadataList[[r]]) <- colnames(countList[[r]])
#  }
#  
#  
#  ## Create the Seurat list  object
#  
#  seuList <- list()
#  for(r in 1:M){
#    seuList[[r]] <- CreateSeuratObject(counts = countList[[r]], meta.data=metadataList[[r]], project = "FASTsimu")
#  }
#  

## ----eval = FALSE-------------------------------------------------------------
#  
#  ## Create PRECASTObject
#  custom_genelist <- row.names(seuList[[1]])
#  set.seed(2023)
#  PRECASTObj <-  CreatePRECASTObject(seuList, customGenelist=custom_genelist)
#  
#  ## User can retain the raw seuList by the following commond.
#  ##  PRECASTObj <-  CreatePRECASTObject(seuList, customGenelist=row.names(seuList[[1]]), rawData.preserve = TRUE)
#  ## check the number of genes/features after filtering step
#  PRECASTObj@seulist

## ----eval = FALSE-------------------------------------------------------------
#  
#  ## seuList is null since the default value `rawData.preserve` is FALSE.
#  PRECASTObj@seuList
#  
#  ## Add adjacency matrix list for a PRECASTObj object to prepare for FAST model fitting.
#  PRECASTObj <- AddAdjList(PRECASTObj, platform = "ST")
#  
#  ## Add a model setting in advance for a PRECASTObj object: verbose =TRUE helps outputing the information in the algorithm;
#  PRECASTObj <- AddParSettingFAST(PRECASTObj, verbose=TRUE)
#  ## Check the parameters
#  PRECASTObj@parameterList

## ----eval = FALSE-------------------------------------------------------------
#  ### set q= 20 here
#  set.seed(2023)
#  PRECASTObj <- FAST(PRECASTObj, q=20)
#  ### Check the results
#  str(PRECASTObj@resList)

## ----eval =FALSE--------------------------------------------------------------
#  set.seed(2023)
#  PRECASTObj <- FAST(PRECASTObj, q=20, fit.model='gaussian')
#  ### Check the results
#  str(PRECASTObj@resList)

## ----eval = FALSE-------------------------------------------------------------
#  ## Obtain the true labels
#  yList <- lapply(PRECASTObj@seulist, function(x) x$Group)
#  ### Evaluate the MacR2
#  MacVec <- sapply(1:length(PRECASTObj@seulist), function(r) get_r2_mcfadden(PRECASTObj@resList$FAST$hV[[r]], yList[[r]]))
#  ### output them
#  print(MacVec)

## ----eval = FALSE-------------------------------------------------------------
#  PRECASTObj <- RunHarmonyLouvain(PRECASTObj, resolution = 0.4)
#  ARI_vec <- sapply(1:M, function(r) mclust::adjustedRandIndex(PRECASTObj@resList$Louvain$cluster[[r]], yList[[r]]))
#  print(ARI_vec)

## ----eval = FALSE-------------------------------------------------------------
#  str(PRECASTObj@resList$Harmony)
#  str(PRECASTObj@resList$Louvain)

## ----eval = FALSE-------------------------------------------------------------
#  seulist_HK <- NULL
#  seuInt <- IntegrateSRTData(PRECASTObj, seulist_HK=seulist_HK, Method = "HarmonyLouvain", seuList_raw=NULL, covariates_use=NULL, verbose=TRUE)
#  seuInt

## ----eval = FALSE-------------------------------------------------------------
#  cols_cluster <- chooseColors(palettes_name = "Nature 10", n_colors = 7, plot_colors = TRUE)

## ----eval = FALSE, fig.width=8, fig.height=6.3--------------------------------
#  p12 <- SpaPlot(seuInt, item = "cluster", batch = NULL, point_size = 1, cols = cols_cluster, combine = TRUE)
#  p12

## ----eval = FALSE, fig.width=9, fig.height=3.3--------------------------------
#  pList <- SpaPlot(seuInt, item = "cluster", title_name= 'Section',batch = NULL, point_size = 1, cols = cols_cluster, combine = FALSE)
#  drawFigs(pList, layout.dim = c(1, 3), common.legend = TRUE, legend.position = "right", align = "hv")

## ----eval = FALSE-------------------------------------------------------------
#  seuInt <- AddUMAP(seuInt, n_comp=3, reduction = 'harmony', assay = 'RNA')
#  seuInt

## ----eval = FALSE, fig.width=9, fig.height=3.3--------------------------------
#  p13 <- SpaPlot(seuInt, batch = NULL, item = "RGB_UMAP", point_size = 1, combine = FALSE, text_size = 15)
#  drawFigs(p13, layout.dim = c(1, 3), common.legend = TRUE, legend.position = "right", align = "hv")
#  

## ----eval = FALSE, fig.width=9, fig.height=3.3--------------------------------
#  seuInt <- AddTSNE(seuInt, n_comp = 2, reduction = 'harmony', assay = 'RNA')
#  p1 <- dimPlot(seuInt, item = "cluster", point_size = 0.5, font_family = "serif", cols = cols_cluster,
#      border_col = "gray10",  legend_pos = "right")  # Times New Roman
#  p2 <- dimPlot(seuInt, item = "batch", point_size = 0.5, font_family = "serif", legend_pos = "right")
#  
#  drawFigs(list(p1, p2), layout.dim = c(1, 2), legend.position = "right")

## ----eval = FALSE-------------------------------------------------------------
#  dat_deg <- FindAllMarkers(seuInt)
#  library(dplyr)
#  n <- 5
#  dat_deg %>%
#      group_by(cluster) %>%
#      top_n(n = n, wt = avg_log2FC) -> top10
#  seuInt <- ScaleData(seuInt)

## ----eval = FALSE, fig.width=9, fig.height=5.4--------------------------------
#  col_here <- c("#F2E6AB", "#9C0141")
#  library(ggplot2)
#  p1 <- DotPlot(seuInt, features=unname(top10$gene), cols=col_here, #  idents = ident_here,
#                col.min = -1, col.max = 1) + coord_flip()+ theme(axis.text.y = element_text(face = "italic"))+
#    ylab("Domain") + xlab(NULL) + theme(axis.text.x = element_text(size=12, angle = 25, hjust = 1, family='serif'),
#                                        axis.text.y = element_text(size=12, face= "italic", family='serif'))
#  p1

## -----------------------------------------------------------------------------
sessionInfo()

