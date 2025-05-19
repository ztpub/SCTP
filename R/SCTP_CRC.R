#' Run deep learning training process with python torch module
#'
#' @param my_seurat input data that you want to predict,should be seurat object

#' @returns a seurat object with predicted malignancy
#' @examples
#' SCTP_CRC(my_seurat = st_dataset)

SCTP_CRC <- function(my_seurat){

    pred_data <- as.matrix(my_seurat[["RNA"]]$data)

    weight = (1-alpha)/alpha
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

    print("Predicting malignancy")
    pred <- unlist(sigmoid(Corr_sc%*%beta_sc+weight*Corr_st%*%beta_st))
    pred <- pred[,1]

    my_seurat@meta.data[,"malignancy"] = pred
    my_seurat@meta.data[,"condition"] = ifelse(pred>0.5, "Tumor", "Normal")

  return(my_seurat)

}
