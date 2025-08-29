rm(list = ls())

work_dir <- "dir"
set.seed(123)
setwd(work_dir)
getwd()

.libPaths()
Sys.info()["nodename"]
gc()


########## library ###########
library(tidyverse)
library(Seurat)
library(dplyr)
library(pastecs)
library(rtracklayer)
library(qs)

source("./functions.R")

obj_dir = "/Share/user/limaor/project/1sampeltag/20250721/mice_lable_unlable/lable_unlable/step3/filtered_feature_bc_matrix"

obj_matrix <- Read10X(obj_dir)
dim(obj_matrix)

CreateSeuratObject_min_cells =3
CreateSeuratObject_min_features = 200

obj <- CreateSeuratObject(counts = obj_matrix, project = "names"
                          , min.cells = CreateSeuratObject_min_cells
                              , min.features = CreateSeuratObject_min_features
)

obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

p0 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.rb"),
                ncol = 4, pt.size = 0)+geom_boxplot()
p0

## Optionally
mt_genes <- grep("^mt-", rownames(obj), ignore.case = TRUE, value = TRUE)
obj <- subset(
  obj,
  features = setdiff(rownames(obj), mt_genes)
)

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

## Select markers
seu = obj
all.markers = FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, 
                                          logfc.threshold = 0.25, group.by ="umap.unintegrated")

all.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20

p1 = DoHeatmap(seu, features = top20$gene, group.by = "umap.unintegrated") + ggtitle("") + NoLegend() +
    theme(
    plot.margin = unit(c(1, 3, 1, 1), "cm")
    )

options(repr.plot.width=5, repr.plot.height=5)
Nebulosa::plot_density(obj, features = c('Nr5a1'), joint =TRUE, reduction = "umap.unintegrated")

qsave(obj, "obj_1.qs")