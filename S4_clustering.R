library(Seurat)
library(dplyr)
library(ggsci)
library(ggplot2)
library(patchwork)
library(monocle)
library(RColorBrewer)
library(plotrix)
library(ggforce)
library(ggpubr)
library(png)
library(clustree)
setwd("~/Desktop/Nymphaea_colorata/")
source("~/Desktop/Nymphaea_colorata/analysis/code/palette_84.R")

set.seed(123)

##### S4-1 carpels cluster #####
S4_1=Read10X(data.dir = "../S4-1/L7_carpels/")
position4 = read.table("../S4-1/L7_carpels/barcodes_pos.tsv",
                       header = F, row.names = 1, sep = "\t")
S4_1=CreateSeuratObject(counts = S4_1,project = "TenX_data",min.cells = 10) 
S4_1[['x']] = position4[1]
S4_1[['y']] = position4[2]
colnames(position4)=c("x","y")

S4_1 = NormalizeData(S4_1)
S4_1 = FindVariableFeatures(S4_1, nfeatures = 1000, selection.method = "vst")
S4_1 = ScaleData(S4_1)
S4_1 = RunPCA(S4_1,nps=50)
S4_1 = FindNeighbors(S4_1,dims = 1:10,k.param = 8)
S4_1 = FindClusters(S4_1, resolution = .9)

### 2nd round ###
marker_genes = FindAllMarkers(S4_1)
marker_genes = marker_genes[marker_genes$p_val_adj<0.05,]
marker_genes$cluster=as.integer(marker_genes$cluster)
VariableFeatures(S4_1)=marker_genes$gene

S4_1 = ScaleData(S4_1)
S4_1 = RunPCA(S4_1,nps=50)
S4_1 = FindNeighbors(S4_1, dims = 1:6,k.param = 15)
S4_1 = FindClusters(S4_1, resolution = .65)
S4_1_markers = FindAllMarkers(S4_1, min.pct = 0.2)



#### S4-1 all clusters ####
S4_1=Read10X(data.dir = "../S4-1/")
position7 = read.table("../S4-1/barcodes_pos.tsv",
                       header = F, row.names = 1, sep = "\t")
S4_1=CreateSeuratObject(counts = S4_1)
S4_1[['x']] = position7[1]
S4_1[['y']] = position7[2]
S4_1 = NormalizeData(S4_1)
S4_1 = FindVariableFeatures(S4_1, nfeatures = 2500, selection.method = "vst")
S4_1 = ScaleData(S4_1)
S4_1 = RunPCA(S4_1,nps=50)
S4_1 = FindNeighbors(S4_1,dims = 1:10,k.param = 8)
S4_1 = FindClusters(S4_1, resolution = 1.5)

# 2nd round
S4_1$avg_exp = colMeans(GetAssayData(object = S4_1))
median_avg = median(S4_1$avg_exp)
S4_1$exp_group = ifelse(S4_1$avg_exp > median_avg, "high", "low")
S4_1_high = subset(S4_1, exp_group == "high")
S4_1_low = subset(S4_1, exp_group == "low")

marker_genes_high = FindAllMarkers(S4_1_high)
marker_genes_high = marker_genes_high[marker_genes_high$p_val_adj<0.05,]
marker_genes_high$cluster=as.integer(marker_genes_high$cluster)

marker_genes_low = FindAllMarkers(S4_1_low)
marker_genes_high = marker_genes_high[marker_genes_high$p_val_adj<0.05,]
marker_genes_low$cluster=as.integer(marker_genes_low$cluster)

marker_genes = rbind(marker_genes_high,marker_genes_low)
VariableFeatures(S4_1)=marker_genes$gene

S4_1 = ScaleData(S4_1)
S4_1 = RunPCA(S4_1,nps=50)
S4_1 = FindNeighbors(S4_1, dims = 1:6,k.param = 12)
S4_1 = FindClusters(S4_1, resolution = .3)
S4_1_markers = FindAllMarkers(S4_1)
S4_1_top100 = S4_1_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)



##### S4-2 stamens cluster #####
S4_2=Read10X(data.dir = "../S4-2/L7_stamen/")
position5 = read.table("../S4-2/L7_stamen/barcodes_pos.tsv",
                       header = F, row.names = 1, sep = "\t")
S4_2=CreateSeuratObject(counts = S4_2)
S4_2[['x']] = position5[1]
S4_2[['y']] = position5[2]
colnames(position5)=c("x","y")

S4_2 = NormalizeData(S4_2)
S4_2 = FindVariableFeatures(S4_2, nfeatures = 1000, selection.method = "vst")
S4_2 = ScaleData(S4_2)
S4_2 = RunPCA(S4_2,nps=50)
S4_2 = FindNeighbors(S4_2)
S4_2 = FindClusters(S4_2, resolution = .2)

# 2nd round
marker_genes = FindAllMarkers(S4_2)
marker_genes = marker_genes[marker_genes$p_val_adj<0.05,]
marker_genes$cluster=as.integer(marker_genes$cluster)
VariableFeatures(S4_2)=marker_genes$gene

S4_2 = ScaleData(S4_2)
S4_2 = RunPCA(S4_2,nps=50)
S4_2 = FindNeighbors(S4_2, dims = 1:7,k.param = 15)
S4_2 = FindClusters(S4_2, resolution = .3)



##### S4-2 all clusters #####
S4_2=Read10X(data.dir = "../S4-2/")
position6 = read.table("../S4-2/barcodes_pos.tsv", header = F, row.names = 1, sep = "\t")
S4_2=CreateSeuratObject(counts = S4_2)
S4_2[['x']] = position6[1]
S4_2[['y']] = position6[2]
colnames(position6)=c("x","y")
S4_2_clu = read.table("../S4_2/clusters.txt",header=T,row.names = 1)
S4_2_clu$seurat_clusters[S4_2_clu$seurat_clusters == 8] =6
S4_2$seurat_clusters = as.factor(S4_2_clu$seurat_clusters)
Idents(S4_2) = factor(S4_2$seurat_clusters, levels = c("0","1","2","3","4","5","6","7"))
S4_2 = NormalizeData(S4_2)
S4_2 = FindVariableFeatures(S4_2, nfeatures = 4000, selection.method = "vst")
S4_2 = ScaleData(S4_2)
S4_2_markers = FindAllMarkers(S4_2, min.pct = 0.25)
