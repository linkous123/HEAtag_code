## This is only for displaying the drawing code. The specific parameters are subject to actual conditions. Repetitive types of images will no longer be shown.
## The visualizations that appear in the other sections are not repeated.
source("./Visualization_funtion.R")

library(Seurat)
library(ggplot2)
library(SCP)
library(gghalves)
library(dplyr)
library(tidyverse)

# fig2 b
plot_genome_3d_scater_htodemux(obj)

## fig2 c
plot_Alluvial(obj, 
              left = "species", 
              right = "MULTIseq", 
              left_select = c("mouse", "rat"), 
              right_select = c("Doublet", "Negative", "ST-mm", "ST-rat"), 
              my_color = c("#D1292F","#4255A1","#517939", "#8551A4", "#D1292F","#4255A1"),
              labels = c('Cell type', 'Demux'))

## fig2 d
plot_htodemux(obj1_htodemux)

## fig2 e
obj_clean_long <- obj_369_clean@meta.data[, c('Tag_3W', 'Tag_6W', "Tag_9W", 'Celltype')] %>%
  pivot_longer(cols = c('Tag_3W', 'Tag_6W', "Tag_9W",),
               names_to = 'Feature',
               values_to = 'Expression')
head(obj_clean_long)

p <- ggplot(data = obj_clean_long, aes(x = Celltype, y = log10(Expression), fill = Celltype)) +
  geom_violin(
    show.legend = FALSE, trim = TRUE, color = "transparent", alpha = 0.8
  ) +
  stat_summary(
    fun = median, geom = "crossbar", width = 0.2, color = "black", linewidth = 0.2, show.legend = FALSE
  )  +
  scale_color_manual(values = celltype_fresh_colors) +
  scale_fill_manual(values = celltype_fresh_colors) +
  # theme_bw() +
  labs(x = '', y = 'log10(Sampletag Count)') +  # 更新 y 轴标签
  theme(
      panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
      strip.background = element_blank(),  # 去除分组标题的背景
    strip.text = element_text(face = "bold", size = 17),
  panel.border = element_rect(color = "black", linewidth  = 0.2, fill = NA),
  plot.margin = margin(t = 0, r = 0, b = 0, l = 30, unit = 'pt')
  ) +
  facet_wrap(~ Feature, ncol = 4) +
  theme(aspect.ratio = 0.2/0.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13)  # 设置 y 轴标签的角度为 90 度
  ) +
  coord_cartesian(ylim = c(0, 5))  # 调整 y 轴范围，设定 max_limit

## fig2 f
plot_umap(obj)

## fig2 g
ht <- GroupHeatmap(
  srt = obj_369_clean,
  features = select_gene,
  group.by = "Celltype",
  heatmap_palette = "YlOrRd",
    show_row_names = TRUE, 
  add_dot = TRUE, add_reticle = F,row_names_side = "left", nlabel = 0,
)


## fig2 h
temp_labels <- obj@meta.data %>%
  group_by(sample) %>%
  tally()

p2 <- table_samples_by_clusters %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'sample') %>%
  mutate(sample = factor(sample, levels = levels(merge_obj@meta.data$sample))) %>%
  ggplot(aes(sample, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity', alpha = 0.8) +
  geom_text(
    data = temp_labels,
    aes(x = sample, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Demux', values = my_color) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 20, b = 0, l = 0, unit = 'pt'), 
      aspect.ratio = 1.218/1
  )

## fig2 i
p1 = ggplot(plot_line_df, aes(x = celltype, y = value, fill = celltype)) +
  geom_col() +
  geom_col_pattern(
    colour         = "black",      
    pattern_fill   = "white",       
    pattern_angle  = 45,            
    pattern_density= 0.01,           
    pattern_spacing= 0.05,          
    pattern_key_scale_factor = 0.6   
  ) +
  facet_wrap(~ experiment, ncol = 4) +   
coord_cartesian(ylim = c(0.95, 1)) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    breaks = seq(0.95, 1.00, by = 0.01)
  ) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Experiment",
       y = "Demux accuracy (Ratio)",
       fill = "Cell Type") +
  theme_classic() +
  scale_fill_manual(values = my_color2) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(face = "bold"),
      aspect.ratio = 0.8/1
)


## fig3 e
CellDimPlot(
  srt =obj, group.by = "Celltype", theme_use = "theme_blank",
  reduction = "umap.unintegrated", label = TRUE, label_repel = TRUE, label_segment_color = "transparent", label.size = 6, label.fg = "black", label.bg = "transparent",
    label_point_color = "transparent",xlab = "UMAP_1", ylab = "UMAP_2"
) + 
  theme(
    legend.position = "bottom",   
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 15),        
    legend.key.size = unit(1.5, "lines"),            
    legend.background = element_rect(fill = "transparent") 
  ) +
  guides(color = guide_legend(nrow = 5, override.aes = list(size = 5)))

## fig3 f
plot_cell_number_scale(obj, cluster_col = "Celltype", sample_col = "MULTIseq", color = c("#D9D9D9", "#737373", "#C6307C", "#4991C1"))

## fig3 h
df <- obj@meta.data %>%
  transmute(
    celltype = factor(Celltype),
    group    = factor(HTO_call, levels = c("label", "unlabel")),
    n_genes  = nFeature_RNA
  )

ggplot() +
  geom_half_violin(
    data = df %>% filter(group == "label"),
    aes(x = celltype, y = n_genes),
    side = "l", trim = FALSE, color = "black", fill = "#C6307C",
      width = 0.95 # 小提琴宽度
  ) +
  geom_half_violin(
    data = df %>% filter(group == "unlabel"),
    aes(x = celltype, y = n_genes),
    side = "r", trim = FALSE, color = "black", fill = "#D9D9D9",
      width = 0.95 # 小提琴宽度
  ) +
  labs(y = "No. of expressed genes", x = NULL, fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
options(repr.plot.width=10, repr.plot.height=5)


## fig3 i
avg_exp <- AverageExpression(object = obj, assays = "RNA", slot = "data", group.by = "HTO_call")

df = as.data.frame(avg_exp$RNA)

cor_test <- cor.test(df$label, df$unlabel, method = "pearson")
r_val  <- signif(cor_test$estimate,  3)
R2 = r_val^2

ggplot(df, aes(x = log1p(label), y = log1p(unlabel))) +
  
  geom_point(size = 3, alpha = 0.7, color = "#2C3E50") +
  
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "#E74C3C") +
  
  annotate("text",
           x    = Inf, y = Inf,
           hjust = 1.1, vjust = 1.1,
           label = paste0("R^2 = ", sprintf("%.2f", R2), "\n", "p < 2.2e-16"),
           size = 5) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(color = "black"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  
  labs(
    title = "",
    x     = "label",
    y     = "unlabel"
  )