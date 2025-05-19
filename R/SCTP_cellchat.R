SCTP_cellchat <- function(my_seurat, spatial_path,scale_factors='scalefactors_json.json', json_file){
  my_seurat <- st_dataset

  pseudotime0 <- my_seurat@meta.data$pseudotime
  spot_type <- my_seurat@meta.data$spot_type
  names(pseudotime0) <- colnames(my_seurat)
  times_sub <- pseudotime0[which(!is.infinite(pseudotime0))]
  times_cat <- rep("Root", length(times_sub)) #<1
  times_cat[findInterval(times_sub, c(quantile(times_sub, 0.15), quantile(times_sub, 0.60))) == 1] = "T2"
  times_cat[findInterval(times_sub, c(quantile(times_sub, 0.60), max(times_sub[is.finite(times_sub)]))) == 1] = "Invasion"

  pseudotime <- rep("Normal", dim(my_seurat)[2])
  pseudotime[which(colnames(my_seurat)%in%names(times_sub))]=times_cat
  pseudotime[which(spot_type=="Front")]="Front"
  my_seurat@meta.data[,"pseudotime_cat"] = pseudotime

  ###### Cell-cell communication ------

  data.input <- Seurat::GetAssayData(my_seurat, layer = "data", assay = "SCT")
  meta <- my_seurat@meta.data # manually create a dataframe consisting of the cell labels


  # spatial image information
  spatial.locs <- Seurat::GetTissueCoordinates(my_seurat, scale = NULL,
                                              cols = c("imagerow", "imagecol"))
  # Scale factors and spot diameters information
  scale.factors <- jsonlite::fromJSON(txt =
                                       file.path(json_file, scale_factors))
  scale.factors <- list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                       fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
  )

  cellchat <- createCellChat(object = data.input,
                             meta = meta,
                             group.by = "pseudotime_cat",
                             datatype = "spatial",
                             coordinates = spatial.locs,
                             scale.factors = scale.factors)

  CellChatDB <- CellChatDB.human
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB

  # set the used database in the object
  cellchat@DB <- CellChatDB.use

  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

  #future::plan("multisession", workers = 1)
  ## identify Over Expressed Genes
  cellchat <- identifyOverExpressedGenes(cellchat)
  # identify Over Expressed Interactions
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat,
                                type = "truncatedMean", trim = 0.1,
                                distance.use = TRUE,
                                scale.distance = 0.01)

  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  return(cellchat)


  }
