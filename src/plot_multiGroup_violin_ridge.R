library(ggplot2)
library(ggridges)

PlotMultiGroupViolin<- function(obj, features, assay = "RNA", slot = "data", cluster, groups = NULL, x = "tigl_id" , ncol =1){
  meta<- obj@meta.data %>%
    tibble::rownames_to_column(var = "cell")
  
  if (is.null(groups)){
    groups<- unique(meta$group)
  }
  
  cells<- meta %>%
    filter(clusters_anno == cluster, group %in% groups) %>%
    pull(cell)
  
  expression_df<- GetAssayData(obj, assay = assay, slot = slot)[features, cells, drop=FALSE] %>%
    as.matrix() %>% as.data.frame() %>%
    tibble::rownames_to_column(var= "gene") %>%
    gather(-1, key = "cell", value = "expression")
  
  expression_df<- left_join(expression_df, meta) 
  
  expression_df<- arrange(expression_df, group)
  
  expression_df$tigl_id<- factor(expression_df$tigl_id, levels = unique(expression_df$tigl_id))
  
  ## plotting
  g<- expression_df %>% 
    ggplot(aes(x = .data[[x]], y = expression)) +
    geom_violin(aes(color = group)) +
    geom_jitter(aes(color = group), width = 0.2, size = 0.5) +
    facet_wrap(~ gene, ncol = ncol) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab("") +
    ggtitle(cluster)
  return(g)
}




PlotMultiGroupRidges<- function(obj, features, assay = "RNA", slot = "data", cluster, groups = NULL, y = "tigl_id" , ncol =1){
  meta<- obj@meta.data %>%
    tibble::rownames_to_column(var = "cell")
  
  if (is.null(groups)){
    groups<- unique(meta$group)
  }
  
  cells<- meta %>%
    filter(clusters_anno == cluster, group %in% groups) %>%
    pull(cell)
  
  expression_df<- GetAssayData(obj, assay = assay, slot = slot)[features, cells, drop =FALSE] %>%
    as.matrix() %>% as.data.frame() %>%
    tibble::rownames_to_column(var= "gene") %>%
    gather(-1, key = "cell", value = "expression")
  
  expression_df<- left_join(expression_df, meta) 
  
  expression_df<- arrange(expression_df, group)
  
  expression_df$tigl_id<- factor(expression_df$tigl_id, levels = unique(expression_df$tigl_id))
  
  ## plotting
  g<- expression_df %>% 
    ggplot(aes(x = expression, y = .data[[y]])) +
    geom_density_ridges(aes(fill = group)) + 
    theme_bw(base_size = 14) +
    facet_wrap(~ gene, ncol = ncol)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab("") +
    ggtitle(cluster)
  return(g)
}