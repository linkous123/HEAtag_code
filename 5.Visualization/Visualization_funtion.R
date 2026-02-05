plot_Alluvial <- function(obj, left, right, left_select, right_select, my_color, labels) {
      require(dplyr)
      require(ggplot2)
      require(ggalluvial)
      require(ggforce)
      
      # Convert and relevel factors
      obj@meta.data <- obj@meta.data %>%
        mutate(
          !!left := factor(!!sym(left), levels = left_select),
          !!right := factor(!!sym(right), levels = right_select)
        )
      
      # Calculate proportions
      prop <- prop.table(table(obj@meta.data[[right]])) * 100
      prop_num <- paste0(round(prop, 2), "%")
      prop_num1 <- setNames(prop_num, right_select)
      
      # Update metadata
      obj@meta.data <- obj@meta.data %>%
        mutate(
          prop = plyr::mapvalues(!!sym(right), from = names(prop_num1), to = prop_num1),
          !!right := factor(paste0(!!sym(right), " (", prop, ")"),
                            levels = paste0(right_select, " (", prop_num, ")"))
        )
      
      # Prepare data for ggplot
      data <- obj@meta.data %>%
      count(!!sym(left), !!sym(right)) %>%
      gather_set_data(1:2) %>%
      mutate(x = factor(x, levels = unique(x)))  
      
      colnames(data) <- c("sample", "seurat_clusters", "n", "id", "x", "y")
      
      # Color assignments
      names <- c(left_select, levels(obj@meta.data[[right]]))
      color_assignments <- setNames(my_color, names)
      
      # Create plot
      ggplot(data, aes(x, id = id, split = y, value = n)) +
        geom_parallel_sets(aes(fill = seurat_clusters), alpha = 0.5, axis.width = 0.15) +
        geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
        geom_text(
          aes(y = n, split = y), 
          stat = 'parallel_sets_axes', 
          fontface = 'plain', 
          size = 5,
          hjust = c(rep(1, length(left_select)), rep(0, length(right_select))),
          nudge_x = c(rep(-0.1, length(left_select)), rep(0.1, length(right_select)))
        ) +
        scale_x_discrete(labels = labels) +
        scale_fill_manual(values = color_assignments) +
        theme_bw() +
        theme(
          legend.position = 'none',
          axis.title = element_blank(),
          axis.text.x = element_text(colour = 'black', size = 14),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 12)
        )
}

plot_genome_3d_scater_htodemux <- function(obj){
    data <- tibble(human = obj$human_umi, mouse= obj$mouse_umi, hamster = obj$hamster_umi, color=obj$HTODemux)

    data$color = fct_relevel(data$color, c('Doublet','Negative','ST-MEF','ST-K562','ST-CHO','ST-293T'))
    data$color = fct_relevel(data$color, c('Negative','Doublet','ST-293T','ST-CHO','ST-K562','ST-MEF'))
    
    names = c('Negative','Doublet','ST-293T','ST-CHO','ST-K562','ST-MEF')
    colors = c("#D9D9D9", "#737373","#8551A4", "#517939", "#4255A1","#D1292F")
    # c("#D25D19", "#E30B69", "#D1292F","#4255A1","#517939", "#8551A4")
    my_color = set_names(colors, names)

    range_all <- range(c(data$human, data$mouse, data$hamster))
    
    p = scatterplot3d(data$human, data$mouse, data$hamster, 
                color = adjustcolor(my_color[data$color], alpha.f = 0.7), pch = 16, type = "p", 
                xlab = "Human RNA UMIs", 
                ylab = "Mouse RNA UMIs", zlab = "Hamster RNA UMIs",
                grid=TRUE, box=FALSE,tick.marks=TRUE,
                angle= 30,
                     xlim = range_all, ylim = range_all, zlim = range_all)
    
        legend("topright", legend = c('Negative','Doublet','ST-293T','ST-CHO','ST-K562','ST-MEF'), col = my_color, 
           pch = 16,box.lty = 0, pt.cex = 1.5)

    return(p)
}

plot_umap = function(obj){
    # singlet
    umap = obj@reductions$umap.rpca@cell.embeddings %>%  
    as.data.frame() %>% 
    cbind(cell_type = obj@meta.data$HTODemux) 

    colnames(umap) = c("UMAP_1", "UMAP_2", "cell_type")

    umap$cell_type <- factor(umap$cell_type, levels = c('ST-MEF', 'ST-K562', 'ST-CHO', 'ST-293T'))

    allcolour1=c("#D1292F","#4255A1","#517939", "#8551A4", "#D25D19")

    p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  geom_point(size = 1 , alpha =0.5 ) +
    scale_color_manual(values = allcolour1)+
    labs(color = "Singlet") +
    theme(

      text = element_text(family = "FS Me Regular", size = 12),
      axis.text = element_text(family = "FS Me Regular", size = 12),
      legend.position = "bottom",
      legend.text = element_text(family = "FS Me Regular", size = 14),
      legend.title = element_text(family = "FS Me Regular", size = 14),
      legend.key = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.7)
    ) + guides(color = guide_legend(override.aes = list(size = 5)))+ 
    theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = 'white'))

    return(p)
}

plot_htodemux <- function(obj){
    hto_plot <- HTOHeatmap(obj, assay = "ST") +
            scale_fill_gradient(low = "#950095", high = "#FFFF00") +
            theme(
                text = element_text(family = "FS Me Regular", size = 12), 
                axis.title = element_text(size = 14), 
                axis.text = element_text(size = 14) 
            )
    return(hto_plot)
}

plot_cell_number_scale <- function(obj, cluster_col, sample_col, color) {
    
  table_clusters_by_samples <- obj@meta.data %>%
    dplyr::rename(cluster = !!sym(cluster_col), 
                  sampleK = !!sym(sample_col)) %>%
    group_by(cluster, sampleK) %>%
    summarize(count = n(), .groups = 'drop') %>%
    spread(sampleK, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    select(c('cluster', 'total_cell_count', everything())) %>%
    arrange(factor(cluster, levels =  levels(obj@meta.data[[cluster_col]])))
  
    print(knitr::kable(table_clusters_by_samples))

    temp_labels <- obj@meta.data %>%
      group_by(.data[[cluster_col]]) %>%
      tally() %>%
      dplyr::rename('cluster' = cluster_col)

    p1 <- table_clusters_by_samples %>%
      select(-c('total_cell_count')) %>%
      reshape2::melt(id.vars = 'cluster') %>%
      #mutate(cluster = factor(cluster, levels = levels(obj@meta.data[[cluster_col]]))) %>%
      ggplot(aes(cluster, value)) +
      geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
      geom_text(
        data = temp_labels,
        aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
        color = 'black', size = 3.0
      ) +
      scale_fill_manual(name = 'Sample', values = color) +
      scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
      coord_cartesian(clip = 'off') +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(t = 20, r = 0, b = 0, l = 15, unit = 'pt')
      )

    p2 <- table_clusters_by_samples %>%
      select(-c('total_cell_count')) %>%
      reshape2::melt(id.vars = 'cluster') %>%
      #mutate(cluster = factor(cluster, levels = levels(obj@meta.data[[cluster_col]]))) %>%
      ggplot(aes(cluster, value)) +
      geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
      geom_text(
        data = temp_labels, aes(x = cluster, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
        color = 'black', size = 3.0
      ) +
      scale_fill_manual(name = 'Sample', values = color) +
      scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
      coord_cartesian(clip = 'off') +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(t = 20, r = 0, b = 0, l = 15, unit = 'pt')
      )

    plot_list = list(p1, p2)
    return(plot_list)

}
