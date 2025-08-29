########## library ###########
library(tidyverse)
library(Seurat)
library(dplyr)
library(pastecs)
library(rtracklayer)
library(qs)

source("./functions.R")

obj = qread("obj_369.qs")

obj_leydig = subset(obj_369, RNA_snn_res.0.4 %in% c("6", "8", "20"))

obj_leydig$stage = obj_leydig$freemuxlet

obj_leydig$stage = gsub("SNP-3W", "early", obj_leydig$stage)
obj_leydig$stage = gsub("SNP-6W", "early", obj_leydig$stage)
obj_leydig$stage = gsub("SNP-9W", "lata", obj_leydig$stage)

dge_vsm = FindMarkers(obj_leydig, group.by = "stage",
                        ident.1 = 'early', ident.2 = "lata")


dge_vsm_sig <- dge_vsm %>% subset(p_val_adj < 0.05)
dge_vsm_sig %>% head()

counts <- GetAssayData(obj_leydig, assay = "RNA", slot = "counts")
universe <- rownames(counts)[Matrix::rowSums(counts > 0) >= 3]

dge_vsm_sig <- dge_vsm %>% rownames_to_column("SYMBOL") %>% filter(p_val_adj < 0.05)
genes_up   <- dge_vsm_sig %>% filter(avg_log2FC > 0) %>% pull(SYMBOL) %>% unique()
genes_down <- dge_vsm_sig %>% filter(avg_log2FC < 0) %>% pull(SYMBOL) %>% unique()

library(Seurat)
library(tidyverse)
library(Matrix)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)

ego_up <- enrichGO(
  gene          = genes_up,
  universe      = universe,
  OrgDb         = org.Mm.eg.db,   
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
) %>% simplify(cutoff = 0.5, by = "p.adjust", select_fun = min)  

ego_down <- enrichGO(
  gene = genes_down, universe = universe, OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
  qvalueCutoff = 0.05, readable = TRUE
) %>% simplify(cutoff = 0.5, by = "p.adjust", select_fun = min)

p1 <- plot_go_lollipop_divergent(ego_up, ego_down, top_n = 10,
                                 term_text_size = 17,
                                 wrap_width = 40,
                                 colors = c("#006D2C","white","#CB181D"),
                                 title = "GO BP")
p1


ego_up <- enrichGO(
  gene          = genes_up,
  universe      = universe,
  OrgDb         = org.Mm.eg.db,  
  keyType       = "SYMBOL",
  ont           = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
) %>% simplify(cutoff = 0.5, by = "p.adjust", select_fun = min)  

ego_down <- enrichGO(
  gene = genes_down, universe = universe, OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH",
  qvalueCutoff = 0.05, readable = TRUE
) %>% simplify(cutoff = 0.5, by = "p.adjust", select_fun = min)


ego_up <- enrichGO(
  gene          = genes_up,
  universe      = universe,
  OrgDb         = org.Mm.eg.db, 
  keyType       = "SYMBOL",
  ont           = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
) %>% simplify(cutoff = 0.5, by = "p.adjust", select_fun = min) 

ego_down <- enrichGO(
  gene = genes_down, universe = universe, OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH",
  qvalueCutoff = 0.05, readable = TRUE
) %>% simplify(cutoff = 0.5, by = "p.adjust", select_fun = min)