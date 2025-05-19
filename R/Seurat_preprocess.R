Seurat_preprocess <- function(counts,type="SC",
                              min.cells = 400, min.features = 0,
                              normalization.method = "LogNormalize", scale.factor = 10000,
                              selection.method = "vst", resolution = 0.6,
                              dims_Neighbors = 1:10, dims_TSNE = 1:10, dims_UMAP = 1:10,
                              verbose = FALSE, data_dir, st_file_name="filtered_feature_bc_matrix.h5"){

  if(type=="SC"){
  data <- CreateSeuratObject(counts = counts, min.cells = min.cells, min.features = min.features)
  data <- NormalizeData(object = data, normalization.method = normalization.method, scale.factor = scale.factor, verbose = verbose)
  data <- FindVariableFeatures(object = data, selection.method = selection.method, verbose = verbose)
  data <- ScaleData(object = data, verbose = verbose)
  data <- RunPCA(object = data, features = VariableFeatures(data), verbose = verbose)
  data <- FindNeighbors(object = data, dims = dims_Neighbors, verbose = verbose)
  data <- FindClusters( object = data, resolution = resolution, verbose = verbose)
  data <- RunTSNE(object = data, dims = dims_TSNE,check_duplicates = FALSE)
  data <- RunUMAP(object = data, dims = dims_UMAP, verbose = verbose)
  } else if(type=="ST"){
    data <- Load10X_Spatial(data.dir = data_dir,
                                filename = st_file_name)

    #data <- CreateSeuratObject(counts = counts, min.cells = min.cells, min.features = min.features)
    data <- SCTransform(object = data, verbose = verbose)
    #data@images$image <- Read10X_Image(image.dir=image_dir)
    #data@images$image@assay <- "RNA"
    #name_st="ST"
    #data@images$image@key <- tolower(paste0(name_st,"_"))

  } else{
    print("type should be SC for single-cell data or ST for spatial transcriptomic data")
  }

  return(data)

}
