#' Run deep learning training process with python torch module
#'
#' @param sc_expr scRNA-seq data
#' @param cell_types Clusters, each each cluster, the most representative cells will be selected
#' @param st_dataset spatial transcriptomic data in matrix
#' @param bulk_dataset bulk gene expression matrix
#' @param pheno phenotype of bulk samples
#' @param out_dir the path in which the features and networks will be stored
#' @returns estimated coefficient and the scRNA-seq and ST data used in the model
#' @examples
#' SCTP(sc_expr=sc_filtered, cell_types=filtered_annot$Cell_subtype, st_dataset = st_dataset, bulk_dataset = bulk_dataset, pheno=pheno, out_dir="OUTPUT")


SCTP <- function(sc_expr, cell_types, st_dataset, bulk_dataset, pheno, out_dir, my_seed=12345, n.epoch=100){

  use_condaenv("py39", required = FALSE)
  script_path <- "inst/run_GATMLP.py"

  k = ncol(sc_expr)
  l = ncol(st_dataset)
  if(k > (2*l)){
    print("Too large single-cell data, finding meta-cells")
    sc_dataset <- meta_cell(sc_expr, cell_types = cell_types)
  } else {
    sc_dataset <- sc_expr
  }

  print("Multi-scale networks construction")
  SCTP_prepare(sc_seurat = sc_dataset, st_seurat = st_dataset, bulk_dataset = bulk_dataset, pheno=pheno, out_dir=out_dir)

  print("Training your model. Please be patient, it may take time (1-19)")
  source_python(script_path)
  res <- train_data(work_dir = out_dir, seed=my_seed, num_epoch=n.epoch)

  kl_divergence <- res[[3]]
  gamma=0.99
  ratio_kl <- c()
  for(k in 1:(length(kl)-1)){
    ratio_kl[k] <- round(kl[k+1]/kl[k], 6)
  }

  opt_kl <- min(which(ratio_kl >= gamma))
  #opt_kl <- which.min(kl_divergence)
  beta_sc <- matrix(unlist(res[[1]]), ncol = 19)
  beta_st <- matrix(unlist(res[[2]]), ncol = 19)

  beta_sc_opt <- beta_sc[,opt_kl]
  beta_st_opt <- beta_st[,opt_kl]
  alpha=0.05*opt_kl
  #source_python(script_path)
  #res <- load_data(work_dir = out_dir)

  return(list(para=list(beta_sc=beta_sc_opt, beta_st=beta_st_opt, alpha = alpha),
              para_all=list(beta_sc_mat = beta_sc, beta_st_mat=beta_st, kl_divergence=kl_divergence),
              sc_expr=sc_dataset@assays$RNA@counts,
              st_expr=st_dataset@assays$RNA@counts) )
}









