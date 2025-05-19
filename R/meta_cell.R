#' Select meta cells from huge single cell dataset
#'
#' @param sc_data scRNA-seq data of the original large dataset
#' @param cell_type Clusters, each each cluster, the most representative cells will be selected
#' @returns selected meta cells
#' @examples
#' meta_cell(sc_dataset, Clusters)

meta_cell <- function(sc_data,  cell_types){

  if (class(sc_data) == "Seurat"){
    network_sc  <- as.matrix(sc_data@graphs$RNA_snn)
  }else{
    sc_data <- CreateSeuratObject(sc_data)
    sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", verbose = F)
    sc_data <- ScaleData(sc_data, verbose = F)
    sc_data <- RunPCA(sc_data, features = VariableFeatures(Seurat_tmp), verbose = F)
    sc_data <- FindNeighbors(sc_data, dims = 1:10, verbose = F)
    network_sc <- as.matrix(sc_data@graphs$RNA_snn)
  }
  cell_degree <- data.frame(snn_degree= colSums(network_sc), Cell_type=cell_types, Index=c(1:length(cell_types)))

  meta.cells <- cell_degree %>%
    group_by(Cell_type) %>%
    arrange(desc(snn_degree)) %>%
    dplyr::slice(1:150) # 5081 cells

  new_scdata <- sc_data[, meta.cells$Index]

  return(new_scdata)
}




