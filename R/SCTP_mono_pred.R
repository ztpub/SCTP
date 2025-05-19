#' Predict malignancy from trained SCTP model with single modality
#'
#' @param SCTP_model trained SCTP model, should be result from SCTP_mono() function
#' @param my_seurat input data that you want to predict,should be seurat object

#' @returns a seurat object with predicted malignancy
#' @examples
#' SCTP_mono_pred(SCTP_model = res, my_seurat=st_dataset)

SCTP_mono_pred <- function(SCTP_model, my_seurat){

  st_data <- as.matrix(SCTP_model$st_expr) #  ST data

  beta <- SCTP_model$beta

  pred_data <- my_seurat@assays$RNA@data
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
  pred <- unlist(1-sigmoid(Corr_st%*%beta[,1]))
  pred <- pred[,1]

  my_seurat@meta.data[,"malignancy"] = pred
  my_seurat@meta.data[,"condition"] = ifelse(pred>0.5, "Tumor", "Normal")

  return (my_seurat)
}



