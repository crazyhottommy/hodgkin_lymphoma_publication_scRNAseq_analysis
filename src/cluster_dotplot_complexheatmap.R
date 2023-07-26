library(Seurat)
library(ComplexHeatmap)
library(tidyverse)

GetMatrixFromSeurat<- function(obj, features, ...) {
  p<- DotPlot(object = obj, features = features, ...)
  df<- p$data
  
  exp_mat<-df %>% 
    select(-pct.exp, -avg.exp.scaled) %>%  
    arrange(id) %>%
    pivot_wider(names_from = id, values_from = avg.exp) %>% 
    as.data.frame() 
  
  row.names(exp_mat) <- exp_mat$features.plot  
  exp_mat <- exp_mat[,-1] %>% as.matrix()
  
  percent_mat<-df %>% 
    select(-avg.exp, -avg.exp.scaled) %>%  
    arrange(id) %>%
    pivot_wider(names_from = id, values_from = pct.exp) %>% 
    as.data.frame()
  
  row.names(percent_mat) <- percent_mat$features.plot  
  percent_mat <- percent_mat[,-1] %>% as.matrix()
  
  if (!identical(dim(exp_mat), dim(percent_mat))) {
    stop("the dimension of the two matrice should be the same!")
  }
  
  if(! all.equal(colnames(exp_mat), colnames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  
  if(! all.equal(rownames(exp_mat), rownames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  
  return(list(exp_mat = exp_mat, percent_mat = percent_mat))
  
}



MakeClusterDotPlot<- function(exp_mat, percent_mat, col_fun, legend_title = "expression",
                              column_title = "NK cells", row_km = 4, 
                              row_names_font_size = 5, column_ha,
                              legend_dot_size = unit(2,"mm"), ...){
  layer_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h, 
              gp = gpar(col = NA, fill = NA))
    grid.circle(x=x,y=y,r= sqrt(pindex(percent_mat, i, j)/100) * legend_dot_size,
                gp = gpar(fill = col_fun(pindex(exp_mat, i, j)), col = NA))}
  
  hp<- Heatmap(exp_mat,
               heatmap_legend_param=list(title= legend_title),
               column_title = column_title, 
               col=col_fun,
               rect_gp = gpar(type = "none"),
               layer_fun = layer_fun,
               row_names_gp = gpar(fontsize = row_names_font_size),
               row_km = row_km ,
               border = "black",
               top_annotation = column_ha,
               ...)
  # see https://github.com/jokergoo/ComplexHeatmap/issues/737 for retain k-means
  
  return(hp)
}


