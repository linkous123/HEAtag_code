########## library ###########
library(tidyverse)
library(Seurat)
library(dplyr)
library(pastecs)
library(rtracklayer)
library(qs)

source("./functions.R")

obj = qread("obj_1.qs")

obj = add_percent_spe(obj, c("mm", "rat"))

plot_spe_density(obj, c("percent.mm_genes", "percent.rat_genes"))


obj = subset(obj, species percent.mm_genes > 80 & percent.rat_genes > 80)

## Optionally
obj$mm.umi = obj$nCount_RNA * obj$percent.mm_genes / 100
obj$rat.umi = obj$nCount_RNA * obj$percent.rat_genes / 100
st_mat = obj@meta.data[, c("percent.mm_genes", "percent.rat_genes")] |> as.matrix() |> t()
HTODemux <- CreateSeuratObject(counts = st_mat,project = "agglutinin", assay = "ST")
HTODemux <- NormalizeData(HTODemux, assay = "ST", normalization.method = "CLR")
HTODemux <- HTODemux(HTODemux, assay = "ST", positive.quantile = 0.99)

# Re-analysis
obj = down_stream(obj)

DimPlot(obj, reduction = "umap.unintegrated", label = T)

## Select parameters

ElbowPlot(obj, ndims = 50)

obj <- FindNeighbors(obj, dims = 1:24, reduction = "pca")

res <- auto_cluster_resolution(
  obj,
  reduction = "pca",
  dims = 1:24,
  resolutions = seq(0.1, 1, 0.1),
  min_cluster_size = 10,
  max_cells_for_sil = 10000
)

obj <- FindClusters(obj, resolution =0.4, cluster.name = "unintegrated_clusters")

obj <- RunUMAP(obj, dims = 1:pcs, reduction = "pca", reduction.name = "umap.unintegrated")

options(repr.plot.width=10, repr.plot.height=10)
FeaturePlot(obj,
                       features = c("percent.mm_genes"),
                       raster=FALSE       
                  ,cols=c("grey","brown"),
                      reduction = "umap.unintegrated",
                  ) +
        #scale_color_viridis_c()+
        theme_bw()+
        theme(panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank())+
        xlab('UMAP_1')+
        ylab('UMAP_2')



qsave(obj, "obj_1.qs")
