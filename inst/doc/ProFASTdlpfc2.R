## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  githubURL <- "https://github.com/feiyoung/ProFAST/blob/main/vignettes_data/seulist2_ID9_10.RDS?raw=true"
#  download.file(githubURL,"seulist2_ID9_10.RDS",mode='wb')
#  

## ----eval = FALSE-------------------------------------------------------------
#  dlpfc2 <- readRDS("./seulist2_ID9_10.RDS")

## ----eval = FALSE-------------------------------------------------------------
#  library(ProFAST)
#  library(PRECAST)
#  library(Seurat)

## ----eval = FALSE-------------------------------------------------------------
#  dlpfc2 ## a list including three Seurat object with default assay: RNA

## ----eval= FALSE--------------------------------------------------------------
#  head(dlpfc2[[1]])

## ----eval= FALSE--------------------------------------------------------------
#  row.names(dlpfc2[[1]])[1:10]

## ----eval = FALSE-------------------------------------------------------------
#  ## Get the gene-by-spot read count matrices
#  countList <- lapply(dlpfc2, function(x) x[["RNA"]]@counts)
#  M <- length(countList)
#  ### transfer the Ensembl ID to symbol
#  for(r in 1:M){
#    row.names(countList[[r]]) <- transferGeneNames(row.names(countList[[r]]), now_name = "ensembl",
#                                                   to_name="symbol",
#                                                   species="Human", Method='eg.db')
#  }
#  
#  ## Check the spatial coordinates: Yes, they are named as "row" and "col"!
#  print(head(dlpfc2[[1]]@meta.data))
#  
#  ## Get the meta data of each spot for each data batch
#  metadataList <- lapply(dlpfc2, function(x) x@meta.data)
#  
#  
#  ## ensure the row.names of metadata in metaList are the same as that of colnames count matrix in countList
#  for(r in 1:M){
#    row.names(metadataList[[r]]) <- colnames(countList[[r]])
#  }
#  
#  
#  ## Create the Seurat list  object
#  
#  seuList <- list()
#  for(r in 1:M){
#    seuList[[r]] <- CreateSeuratObject(counts = countList[[r]], meta.data=metadataList[[r]], project = "ProFASTdlpfc2")
#  }
#  print(seuList)
#  row.names(seuList[[1]])[1:10]

## ----eval = FALSE-------------------------------------------------------------
#  ## Create PRECASTObject
#  set.seed(2023)
#  PRECASTObj <- CreatePRECASTObject(seuList, project = "ProFASTdlpfc2", gene.number = 2000, selectGenesMethod = "HVGs",
#      premin.spots = 20, premin.features = 20, postmin.spots = 1, postmin.features = 10)
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
#  ## Add adjacency matrix list for a PRECASTObj object to prepare for ProFAST model fitting.
#  PRECASTObj <- AddAdjList(PRECASTObj, platform = "Visium")
#  
#  ## Add a model setting in advance for a PRECASTObj object: verbose =TRUE helps outputing the information in the algorithm;
#  PRECASTObj <- AddParSettingProFAST(PRECASTObj, verbose=TRUE)
#  ## Check the parameters
#  PRECASTObj@parameterList

## ----eval = FALSE-------------------------------------------------------------
#  ### set q= 15 here
#  set.seed(2023)
#  PRECASTObj <- ProFAST(PRECASTObj, q=15, fit.model='gaussian')
#  
#  ### Check the results
#  str(PRECASTObj@resList)

## ----eval =FALSE--------------------------------------------------------------
#  set.seed(2023)
#  PRECASTObj <- ProFAST(PRECASTObj, q=15)
#  ### Check the results
#  str(PRECASTObj@resList)

## ----eval = FALSE-------------------------------------------------------------
#  ## Obtain the true labels
#  yList <- lapply(PRECASTObj@seulist, function(x) x$layer_guess_reordered)
#  ### Evaluate the MacR2
#  MacVec <- sapply(1:length(PRECASTObj@seulist), function(r) get_r2_mcfadden(PRECASTObj@resList$ProFAST$hV[[r]], yList[[r]]))
#  ### output them
#  print(MacVec)

## ----eval = FALSE-------------------------------------------------------------
#  PRECASTObj <- RunHarmonyLouvain(PRECASTObj, resolution = 0.3)
#  ARI_vec_louvain <- sapply(1:M, function(r) mclust::adjustedRandIndex(PRECASTObj@resList$Louvain$cluster[[r]], yList[[r]]))
#  print(ARI_vec_louvain)

## ----eval = FALSE-------------------------------------------------------------
#  PRECASTObj <- RuniSCMEB(PRECASTObj, seed=1)
#  str(PRECASTObj@resList$iSCMEB)
#  sapply(1:M, function(r) mclust::adjustedRandIndex(PRECASTObj@resList$iSCMEB$cluster[[r]], yList[[r]]))

## ----eval = FALSE-------------------------------------------------------------
#  str(PRECASTObj@resList$iSCMEB)

## ----eval = FALSE-------------------------------------------------------------
#  ## select the HK genes
#  HKgenes <- SelectHKgenes(seuList, species= "Human", HK.number=200)
#  seulist_HK <- lapply(seuList, function(x) x[HKgenes, ])
#  seulist_HK

## ----eval = FALSE-------------------------------------------------------------
#  seuInt <- IntegrateSRTData(PRECASTObj, seulist_HK=seulist_HK, Method = "iSC-MEB",
#                             seuList_raw=NULL, covariates_use=NULL, verbose=TRUE)
#  seuInt

## ----eval = FALSE-------------------------------------------------------------
#  cols_cluster <- chooseColors(palettes_name = "Nature 10", n_colors = 8, plot_colors = TRUE)

## ----eval = FALSE, fig.width=6.5, fig.height=3--------------------------------
#  p12 <- SpaPlot(seuInt, item = "cluster", batch = NULL, point_size = 1, cols = cols_cluster, combine = FALSE)
#  library(ggplot2)
#  p12 <- lapply(p12, function(x) x+ coord_flip() + scale_x_reverse())
#  drawFigs(p12, layout.dim = c(1,2), common.legend = T)

## ----eval = FALSE-------------------------------------------------------------
#  seuInt <- AddUMAP(seuInt, n_comp=3, reduction = 'iscmeb', assay = 'ProFAST')
#  seuInt

## ----eval = FALSE, fig.width=6, fig.height=3.4--------------------------------
#  p13 <- SpaPlot(seuInt, batch = NULL, item = "RGB_UMAP", point_size = 1, combine = FALSE, text_size = 15)
#  p13 <- lapply(p13, function(x) x+ coord_flip() + scale_x_reverse())
#  drawFigs(p13, layout.dim = c(1, 2), common.legend = TRUE, legend.position = "right", align = "hv")

## ----eval = FALSE, fig.width=8, fig.height=3----------------------------------
#  seuInt <- AddTSNE(seuInt, n_comp = 2, reduction = 'iscmeb', assay = 'ProFAST')
#  p1 <- dimPlot(seuInt, item = "cluster", point_size = 0.5, font_family = "serif", cols = cols_cluster,
#      border_col = "gray10",  legend_pos = "right")  # Times New Roman
#  p2 <- dimPlot(seuInt, item = "batch", point_size = 0.5, font_family = "serif", legend_pos = "right")
#  
#  drawFigs(list(p1, p2), layout.dim = c(1, 2), legend.position = "right", align = "hv")

## ----eval = FALSE-------------------------------------------------------------
#  dat_deg <- FindAllMarkers(seuInt)
#  library(dplyr)
#  n <- 5
#  dat_deg %>%
#      group_by(cluster) %>%
#      top_n(n = n, wt = avg_log2FC) -> top10
#  

## ----eval = FALSE, fig.width=8, fig.height=6----------------------------------
#  col_here <- c("#F2E6AB", "#9C0141")
#  library(ggplot2)
#  marker_Seurat <- readRDS("151507_DEGmarkerTop_5.rds")
#  marker_Seurat <- lapply(marker_Seurat, function(x) transferGeneNames(x, Method='eg.db'))
#  maker_genes <- unlist(lapply(marker_Seurat, toupper))
#  features <- maker_genes[!duplicated(maker_genes)]
#  p_dot <- DotPlot(seuInt, features=unname(features), cols=col_here, #  idents = ident_here,
#                col.min = -1, col.max = 1) + coord_flip()+ theme(axis.text.y = element_text(face = "italic"))+
#    ylab("Domain") + xlab(NULL) + theme(axis.text.x = element_text(size=12, angle = 25, hjust = 1, family='serif'),
#                                        axis.text.y = element_text(size=12, face= "italic", family='serif'))
#  p_dot
#  # table(unlist(yList), Idents(seuInt))

## ----eval = FALSE, fig.width=8, fig.height=6.5--------------------------------
#  seuInt2 <- RenameIdents(seuInt, '1' = 'Layer6', '2' = 'Layer2/3', '3'="Layer5", '4'="WM",
#                             '5'='Layer3/4', '6'= 'Layer1')
#  levels(seuInt2)  <- c("Layer1", "Layer2/3", "Layer3/4", "Layer5", "Layer6", "WM")
#  p_dot2 <- DotPlot(seuInt2, features=unname(features), cols=col_here, #  idents = ident_here,
#                col.min = -1, col.max = 1) + coord_flip()+ theme(axis.text.y = element_text(face = "italic"))+
#    ylab("Domain") + xlab(NULL) + theme(axis.text.x = element_text(size=12, angle = 25, hjust = 1, family='serif'),
#                                        axis.text.y = element_text(size=12, face= "italic", family='serif'))
#  p_dot2

## ----eval =FALSE, fig.width=8, fig.height=6.5---------------------------------
#  suppressPackageStartupMessages(library(CellChat))
#  data.input =  seuInt2[["ProFAST"]]@data-min(seuInt2[["ProFAST"]]@data)# normalized data matrix with a movement to make the value not smaller than zero.
#  meta = data.frame(labels=Idents(seuInt2)) # a dataframe with rownames containing cell mata data
#  row.names(meta) <- colnames(data.input)
#  head(meta)
#  # Prepare input data for CelChat analysis
#  dim(data.input)
#  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
#  cellchat <- addMeta(cellchat, meta = meta)
#  cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
#  levels(cellchat@idents) # show factor levels of the cell labels
#  groupSize <- as.numeric(table(cellchat@idents)) #

## ----eval =FALSE--------------------------------------------------------------
#  ## Set the ligand-receptor interaction database
#  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#  showDatabaseCategory(CellChatDB)
#  dplyr::glimpse(CellChatDB$interaction)

## ----eval =FALSE--------------------------------------------------------------
#  cellchat@DB <- CellChatDB
#  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#  future::plan("multiprocess", workers = 2) # do parallel
#  cellchat <- identifyOverExpressedGenes(cellchat)
#  cellchat <- identifyOverExpressedInteractions(cellchat)
#  # project gene expression data onto PPI network (optional)
#  cellchat <- projectData(cellchat, PPI.human)
#  ## Compute the communication probability and infer cellular communication network
#  cellchat <- computeCommunProb(cellchat)
#  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
#  cellchat <- filterCommunication(cellchat, min.cells = 10)
#  cellchat <- computeCommunProbPathway(cellchat)
#  cellchat <- aggregateNet(cellchat)
#  groupSize <- as.numeric(table(cellchat@idents))
#  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#  

## ----eval =FALSE--------------------------------------------------------------
#  het1 <- netVisual_heatmap(cellchat, signaling = NULL,measure = "count",
#                            font.size = 12, font.size.title = 14) # color.heatmap = "Reds"
#  
#  
#  het2 <- netVisual_heatmap(cellchat, signaling = NULL,measure = "weight", font.size = 12,
#                            font.size.title = 14)
#  het2

## -----------------------------------------------------------------------------
sessionInfo()

