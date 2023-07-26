library(Seurat)
library(aplot)
library(Polychrome)

cluster_dot_plot<- function(obj, features, annotation, color_scale_limits = c(-1.5, 2.5), 
                            annotation_colors = NULL){
  p<- DotPlot(object = obj, features = features)
  df<- p$data
  df<- left_join(df, annotation, by =c("id" = "seurat_clusters")) %>%
    mutate(id = paste0(annotation, "_", id))
  
  dot_plot <- df %>% 
    filter(!is.nan(avg.exp.scaled)) %>%
    ggplot(aes(x=id, y = features.plot, color = avg.exp.scaled, size = pct.exp)) + 
    geom_point() + 
    cowplot::theme_cowplot() + 
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab(NULL) +
    xlab(NULL) +
    theme(axis.ticks = element_blank()) +
    scale_color_gradientn(colours = viridis::viridis(20), limits = color_scale_limits , oob = scales::squish, name = 'scaled expression') +
    scale_y_discrete(position = "right")
  
  mat <- df %>% 
    select(-pct.exp, -avg.exp, -annotation) %>%  # drop unused columns to facilitate widening
    filter(!is.nan(avg.exp.scaled)) %>%
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
    data.frame() # make df as tibbles -> matrix annoying
  row.names(mat) <- mat$features.plot  # put gene in `row`
  mat <- mat[,-1] #drop gene column as now in rows
  mat<- mat[complete.cases(mat), ] # remove some rows with NaNs
  clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
  
  ggtree_plot <- ggtree::ggtree(clust)
  
  v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
  library(ggtree)
  ggtree_plot_col <- ggtree(v_clust) + layout_dendrogram()
  
  if (is.null(annotation_colors)) {
    labels<- ggplot(df, aes(id, y=1, fill= annotation)) + geom_tile() +
      scale_fill_brewer(palette = 'Paired',name="Cell Type") + 
      theme_void()
  } else {
    labels<- ggplot(df, aes(id, y=1, fill= annotation)) + geom_tile() +
      scale_fill_manual(values = annotation_colors ,name="Cell Type") + 
      theme_void()
  }
  
  dot_plot %>% 
    insert_left(ggtree_plot, width=.2) %>%
    insert_top(labels, height=.02) %>%
    insert_top(ggtree_plot_col, height=.1)
  
}