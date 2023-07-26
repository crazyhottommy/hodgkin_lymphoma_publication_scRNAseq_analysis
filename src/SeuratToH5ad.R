SeuratToH5ad <- function(seurat_object, filename, assay = NULL, res = 1) {
  library(reticulate)
  
  if (!py_module_available("anndata") | !py_module_available("scanpy") | !py_module_available("igraph") | !py_module_available("louvain")) {
    stop("Please install the anndata python module")
  }
  ad <- import("anndata")
  sc <- import("scanpy")
  
  message(paste("Starting to fix the mess..."))
  
  raw <- seurat_object@assays$RNA@data
  if (assay == "RNA") {
    X <- seurat_object@assays$RNA@scale.data
  } else if (assay == "SCT") {
    X <- seurat_object@assays$SCT@scale.data
  } else {
    stop("Please select an existent assay")
  }
  
  cell_names <- colnames(x = X)
  gene_names <- rownames(x = X)
  raw <- as(object = raw, Class = "dgCMatrix")
  
  scipy <- import(module = 'scipy.sparse', convert = FALSE)
  sp_sparse_csc <- scipy$csc_matrix
  raw.rownames <- rownames(x = raw)
  raw <- sp_sparse_csc(
    tuple(np_array(raw@x), np_array(raw@i), np_array(raw@p)),
    shape = tuple(raw@Dim[1], raw@Dim[2])
  )
  
  raw <- raw$T
  raw <- dict(X = raw, var = dict(var_names = raw.rownames))
  
  X <- np_array(t(x = X))
  
  obsm <- list()
  for (dr in names(seurat_object@reductions)) {
    obsm[[paste0("X_",dr)]] <- np_array(Embeddings(
      object = seurat_object,
      reduction = dr
    ))
  }
  obsm <- dict(obsm)
  meta_data <- seurat_object@meta.data
  if ("nCount_RNA" %in% colnames(x = meta_data)) {
    colnames(x = meta_data) <- gsub(
      pattern = "nCount_RNA",
      replacement = "n_counts",
      x = colnames(x = meta_data)
    )
  }
  if ("nFeature_RNA" %in% colnames(x = meta_data)) {
    colnames(x = meta_data) <- gsub(
      pattern = "nFeature_RNA",
      replacement = "n_genes",
      x = colnames(x = meta_data)
    )
  }
  colnames(x = meta_data) <- gsub(
    pattern = "\\.",
    replacement = "_",
    x = colnames(x = meta_data)
  )
  
  anndata.object <- ad$AnnData(
    raw = raw,
    X = X,
    obs = meta_data,
    obsm = obsm
  )
  anndata.object$var_names <- gene_names
  anndata.object$obs_names <- cell_names
  
  #message(paste("Clustering for resolution:", res))
  #sc$pp$neighbors(anndata.object)
  #sc$tl$louvain(anndata.object, resolution=res, key_added = paste0("res", res))
  
  message(paste("Writing to h5ad file..."))
  anndata.object$write(filename)
  message(paste("Finished!!"))
}