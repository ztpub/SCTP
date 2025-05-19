#' Run deep learning training process with python torch module
#'
#' @param sc_seurat seurat object of scRNA-seq data
#' @param st_seurat seurat object of ST data
#' @param st_dataset spatial transcriptomic data in matrix
#' @param bulk_dataset bulk gene expression matrix
#' @param pheno phenotype of bulk samples
#' @param out_dir the path in which the features and networks will be stored
#' @returns create networks and features and phenotype for model training
#' @examples
#' SCTP_prepare(sc_expr=sc_dataset, st_dataset = st_dataset, bulk_dataset = bulk_dataset, pheno=pheno, out_dir="OUTPUT")


SCTP_prepare <- function(sc_seurat, st_seurat, bulk_dataset, pheno, out_dir){
  # SC data
  Seurat_tmp <- FindVariableFeatures(sc_seurat, selection.method = "vst", verbose = F)
  Seurat_tmp <- ScaleData(Seurat_tmp, verbose = F)
  Seurat_tmp <- RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp), verbose = F)
  Seurat_tmp <- FindNeighbors(Seurat_tmp, dims = 1:10, verbose = F)
  network_sc  <- as.matrix(Seurat_tmp@graphs$RNA_snn)
  diag(network_sc) <- 0
  network_sc[which(network_sc != 0)] <- 1 # SC adjacency matrix

  sc_exprs <- as.matrix(Seurat_tmp@assays$RNA@counts)
  common <- intersect(rownames(bulk_dataset), rownames(sc_exprs))
  dataset0 <- as.matrix(cbind(bulk_dataset[common,], sc_exprs[common,]))
  dataset1 <- normalize.quantiles(dataset0)
  rownames(dataset1) <- rownames(dataset0)
  colnames(dataset1) <- colnames(dataset0)
  Expression_bulk <- dataset1[,1:ncol(bulk_dataset)]
  Expression_cell <- dataset1[,(ncol(bulk_dataset) + 1):ncol(dataset1)]
  X_sc <- cor(Expression_bulk, Expression_cell)
  dim(X_sc) # SC feature matrix

  # ST data
  if(is.null(st_dataset@assays$SCT) == FALSE){Seurat_tmp = st_seurat}else{
    Seurat_tmp <- SCTransform(st_seurat, verbose = FALSE)
  }
  my_coordinates = GetTissueCoordinates(Seurat_tmp)
  my_image = GetImage(Seurat_tmp , mode="raster")
  x=unlist(my_coordinates[,1])
  y=unlist(my_coordinates[,2])
  my_adjacency = calculate_adj_matrix(x=x, y=y, x_pixel=x, y_pixel=y, image = my_image)
  p=0.02 # node percentage
  #Find the l value given p
  l=search_l(p, my_adjacency, start=0.1, end=10, tol=0.01, max_run=100)
  adj_exp <- exp(-1*(my_adjacency^2)/(2*(l^2)))
  network_st <- adj_exp
  diag(network_st) <- 0
  network_st[which(network_st != 0)] <- 1 # ST adjacency matrix

  st_exprs <- as.matrix(Seurat_tmp@assays$RNA@counts)
  common <- intersect(rownames(bulk_dataset), rownames(st_exprs))
  dataset0 <- cbind(bulk_dataset[common,], st_exprs[common,])
  dataset1 <- normalize.quantiles(dataset0)
  rownames(dataset1) <- rownames(dataset0)
  colnames(dataset1) <- colnames(dataset0)
  Expression_bulk <- dataset1[,1:ncol(bulk_dataset)]
  Expression_spot <- dataset1[,(ncol(bulk_dataset) + 1):ncol(dataset1)]
  X_st <- cor(Expression_bulk, Expression_spot) # ST feature matrix

  write.csv(X_sc, file = paste0(out_dir, "/features_sc.csv"), row.names = FALSE)
  write.csv(network_sc, file = paste0(out_dir, "/network_sc.csv"), row.names = FALSE)
  write.csv(X_st, file = paste0(out_dir, "/features_st.csv"), row.names = FALSE)
  write.csv(network_st, file = paste0(out_dir, "/network_st.csv"), row.names = FALSE)
  write.csv(pheno, file=paste0(out_dir, "/pheno.csv"))

  print("Networks and features are well saved.")
  #return(list(network_sc= network_sc, network_st=network_st, X_sc=X_sc, X_st=X_st, Y=pheno))
}


