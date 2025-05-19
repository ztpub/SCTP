SCTP_mono <- function(st_dataset, bulk_dataset, pheno, out_dir, my_seed=12345, n.epoch=100){

  use_condaenv("py39", required = FALSE)
  script_path <- "inst/run_GATMLP.py"

  if(is.null(st_dataset@assays$SCT) == FALSE){Seurat_tmp = st_dataset}else{
    Seurat_tmp <- SCTransform(st_dataset, verbose = FALSE)
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

  st_exprs <- as.matrix(Seurat_tmp@assays$SCT@counts)
  common <- intersect(rownames(bulk_dataset), rownames(st_exprs))
  dataset0 <- cbind(bulk_dataset[common,], st_exprs[common,])
  dataset1 <- normalize.quantiles(dataset0)
  rownames(dataset1) <- rownames(dataset0)
  colnames(dataset1) <- colnames(dataset0)
  Expression_bulk <- dataset1[,1:ncol(bulk_dataset)]
  Expression_spot <- dataset1[,(ncol(bulk_dataset) + 1):ncol(dataset1)]
  X_st <- cor(Expression_bulk, Expression_spot) # ST feature matrix

  write.csv(X_st, file = paste0(out_dir, "/features_st.csv"), row.names = FALSE)
  write.csv(network_st, file = paste0(out_dir, "/network_st.csv"), row.names = FALSE)
  write.csv(pheno, file=paste0(out_dir, "/pheno.csv"))
  print("Network and features are well saved.")


  print("Training your model.")
  source_python(script_path)
  beta <- train_data_mono(work_dir = out_dir, seed=my_seed, num_epoch=n.epoch)

  return(list(beta=beta,
              st_expr=st_dataset@assays$RNA@counts) )

}
