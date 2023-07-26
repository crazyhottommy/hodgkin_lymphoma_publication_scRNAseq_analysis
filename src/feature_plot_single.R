library(Seurat)
FeaturePlotSingle<- function(obj, feature, metadata_column, ...){
  all_cells<- colnames(obj)
  groups<- levels(obj@meta.data[, metadata_column])
  
  minimal<- min(obj[['RNA']]@data[feature, ])
  maximal<- max(obj[['RNA']]@data[feature, ])
  ps<- list()
  for (group in groups) {
    subset_indx<- obj@meta.data[, metadata_column] == group
    subset_cells<- all_cells[subset_indx]
    p<- FeaturePlot(obj, features = feature, cells= subset_cells, ...) +
      scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
      ggtitle(group) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    ps[[group]]<- p
  }
  
  
  return(ps)
}

