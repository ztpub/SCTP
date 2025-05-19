#' Run deep learning training process with python torch module
#'
#' @param my_seurat input data that you want to predict,should be a seurat object

#' @returns a seurat object with estimated speudotime
#' @examples
#' SCTP_invasion(my_seurat = st_dataset)

SCTP_invasion <- function(my_seurat){

    if("malignancy" %in% names(my_seurat@meta.data)){
      print("Estimating tumor pseudotime")
      cds <- subset(my_seurat, subset= malignancy > 0.5) #only tumors
      cds <- as.cell_data_set(cds)
      cds <- cluster_cells(cds, resolution=1e-3)
      cds <- monocle3::learn_graph(cds, verbose = FALSE)
      cds <- order_cells(cds, root_cells = colnames(cds[,monocle3::clusters(cds) == 1]))

      times_sub <- pseudotime(cds, reduction_method = "UMAP")
      Front <- names(times_sub)[is.infinite(times_sub)]
      Tumor <- names(times_sub)[!is.infinite(times_sub)]
      spot_type <- rep("Normal", dim(my_seurat)[2])
      spot_type[which(colnames(my_seurat)%in%Front)] <- "Front"
      spot_type[which(colnames(my_seurat)%in%Tumor)] <- "Tumor"
      my_seurat@meta.data[,"spot_type"] = spot_type

      times_sub[which(is.infinite(times_sub))] <- - max(times_sub[is.finite(times_sub)])
      times_sub <- times_sub+max(times_sub[is.finite(times_sub)])

      pseudotime <- rep(Inf, dim(my_seurat)[2])
      pseudotime[which(colnames(my_seurat)%in%names(times_sub))]=-times_sub
      my_seurat@meta.data[,"pseudotime"] = pseudotime

      return(my_seurat)
    } else {
      stop("malignancy not found, run SCTP_CRC() first or provide one.")
    }
}
