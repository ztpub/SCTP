#' Predict malignancy from trained SCTP model
#'
#' @param SCTP_model trained SCTP model, should be result from SCTP() function
#' @param my_seurat input data that you want to predict,should be seurat object

#' @returns a seurat object with predicted malignancy
#' @examples
#' SCTP_pred(SCTP_model = res, my_seurat=st_dataset)

SCTP_pred <- function(SCTP_model, my_seurat){

  sc_data <- as.matrix(SCTP_model$sc_expr) #  SC data
  st_data <- as.matrix(SCTP_model$st_expr) #  ST data

  beta_sc <- SCTP_model$para$beta_sc
  beta_st <- SCTP_model$para$beta_st
  alpha<- SCTP_model$para$alpha

  pred_data <- my_seurat@assays$RNA@data
  print("Calculating correlations with SC DATA")
  common_sc <- intersect(rownames(pred_data), rownames(sc_data))
  dataset0 <- as.matrix(cbind(pred_data[common_sc,], sc_data[common_sc,]))
  dataset1 <- normalize.quantiles(dataset0)
  rownames(dataset1) <- rownames(dataset0)
  colnames(dataset1) <- colnames(dataset0)
  Expression_pred <- dataset1[,1:ncol(pred_data)]
  Expression_sc <- dataset1[,(ncol(pred_data) + 1):ncol(dataset1)]
  Corr_sc <- scale(cor(Expression_pred, Expression_sc))

  print("Calculating correlations with ST DATA")
  common_st <- intersect(rownames(pred_data), rownames(st_data))
  dataset0 <- as.matrix(cbind(pred_data[common_st,], st_data[common_st,]))
  dataset1 <- normalize.quantiles(dataset0)
  rownames(dataset1) <- rownames(dataset0)
  colnames(dataset1) <- colnames(dataset0)
  Expression_pred <- dataset1[,1:ncol(pred_data)]
  Expression_st <- dataset1[,(ncol(pred_data) + 1):ncol(dataset1)]
  Corr_st <- scale(cor(Expression_pred, Expression_st))

  print("Predicting malignancy for all cells or spots")
  weight = (1-alpha)/alpha
  pred <- unlist(1-sigmoid(Corr_sc%*%beta_sc+weight*Corr_st%*%beta_st))
  pred <- pred[,1]

  my_seurat@meta.data[,"malignancy"] = pred
  my_seurat@meta.data[,"condition"] = ifelse(pred>0.5, "Tumor", "Normal")

  return (my_seurat)
}



