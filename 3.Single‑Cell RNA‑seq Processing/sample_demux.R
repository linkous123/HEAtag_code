########## library ###########
library(tidyverse)
library(Seurat)
library(dplyr)
library(pastecs)
library(rtracklayer)
library(qs)

source("./functions.R")

obj = qread("obj_1.qs")

data_df <- read.csv("your_heattg.csv")
rownames(data_df) <- data_df$cellbarcode
data_df$cellbarcode <- NULL
data_df$nomatch <- NULL
data_df <- data_df[!(rownames(data_df) == "nomatch"), ]

colnames(data_df) <- c("rat", "mm")
colnames(data_df) <- paste0("ST_", colnames(data_df))

data_df_sub <- data_df[intersect(colnames(obj),rownames(data_df) ), ]
st_mat <- as.matrix(t(data_df_sub))

## demux
HTODemux <- CreateSeuratObject(counts = st_mat,project = "agglutinin", assay = "ST")
HTODemux <- NormalizeData(HTODemux, assay = "ST", normalization.method = "CLR")
HTODemux <- HTODemux(HTODemux, assay = "ST", positive.quantile = 0.99)

MULTIseq <- CreateSeuratObject(counts = st_mat,project = "agglutinin", assay = "ST")
MULTIseq <- NormalizeData(MULTIseq, assay = "ST", normalization.method = "CLR")
MULTIseq <- MULTIseqDemux(MULTIseq,
      assay = "ST",
      #quantile = 0.5,
      autoThresh = T,
      maxiter = 10,
      qrange = seq(from = 0.01, to = 0.99, by = 0.01),
      verbose = TRUE
)

library(hashDemux)
hashDemux <- CreateSeuratObject(counts = st_mat, assay="HTO")
hashDemux <- NormalizeData(hashDemux, assay = "HTO", normalization.method = "CLR",margin = 2)
hashDemux = clustering_based_demux(obj, assay = "HTO", nCores = 1)


