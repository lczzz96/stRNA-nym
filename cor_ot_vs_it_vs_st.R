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
library("gridExtra")
library(pheatmap)
library(reshape2)
library(edgeR)
library(ggrepel)
library(FNN)
library(tidyr)
library(matrixStats)

S3=Read10X(data.dir = "../S3/L7/")
position3 = read.table("../S3/L7/barcodes_pos.tsv",header = F, row.names = 1, sep = "\t")
S3=CreateSeuratObject(counts = S3) 
S3[['x']] = position3[1]
S3[['y']] = position3[2]
colnames(position3)=c("x","y")

S3$sample = c("all")
S3$cell = rownames(S3@meta.data)

sepal_cell = read.table("../S3/L7/L7_sepal_barcodes_pos.tsv")[,1]
S3$sample = ifelse(S3$cell %in% sepal_cell, "sepal", S3$sample)
petal_cell = read.table("../S3/L7/L7_petal_barcodes_pos.tsv")[,1]
S3$sample = ifelse(S3$cell %in% petal_cell, "petal", S3$sample)
stamen1_cell = read.table("../S3/L7/L7_stamen1_barcodes_pos.tsv")[,1]
S3$sample = ifelse(S3$cell %in% stamen1_cell, "stamen1", S3$sample)
stamen2_cell = read.table("../S3/L7/L7_stamen2_barcodes_pos.tsv")[,1]
S3$sample = ifelse(S3$cell %in% stamen2_cell, "stamen2", S3$sample)
stamen3_cell = read.table("../S3/L7/L7_stamen3_barcodes_pos.tsv")[,1]
S3$sample = ifelse(S3$cell %in% stamen3_cell, "stamen3", S3$sample)
stamen4_cell = read.table("../S3/L7/L7_stamen4_barcodes_pos.tsv")[,1]
S3$sample = ifelse(S3$cell %in% stamen4_cell, "stamen4", S3$sample)

S3 = subset(S3, subset = nFeature_RNA > 5 & nFeature_RNA < 3000)
S3 = NormalizeData(S3)
S3 = FindVariableFeatures(S3, nfeatures = 2000, selection.method = "vst")
all.genes = rownames(S3)
S3 = ScaleData(S3, features = all.genes)
S3 = RunPCA(S3, features = VariableFeatures(object = S3))
OT2S = subset(S3, sample != "all")


#### Stage3 OT vs It vs Stamens
sepal = subset(OT2S, sample == "sepal")
petal = subset(OT2S, sample == "petal")
stamen1 = subset(OT2S, sample == "stamen1")
stamen2 = subset(OT2S, sample == "stamen2")
stamen3 = subset(OT2S, sample == "stamen3")
stamen4 = subset(OT2S, sample == "stamen4")

sepal_pca = t(sepal@reductions$pca@cell.embeddings)
petal_pca = t(petal@reductions$pca@cell.embeddings)
stamen1_pca = t(stamen1@reductions$pca@cell.embeddings)
stamen2_pca = t(stamen2@reductions$pca@cell.embeddings)
stamen3_pca = t(stamen3@reductions$pca@cell.embeddings)
stamen4_pca = t(stamen4@reductions$pca@cell.embeddings)

sepal_pca1 = data.frame(rowMeans(sepal_pca))
petal_pca1 = data.frame(rowMeans(petal_pca))
stamen1_pca1 = data.frame(rowMeans(stamen1_pca))
stamen2_pca1 = data.frame(rowMeans(stamen2_pca))
stamen3_pca1 = data.frame(rowMeans(stamen3_pca))
stamen4_pca1 = data.frame(rowMeans(stamen4_pca))

all_pca = cbind(sepal_pca1,petal_pca1,stamen4_pca1,stamen3_pca1,stamen2_pca1,stamen1_pca1)
colnames(all_pca) = c("sepal", "petal", "stamen4", "stamen3", "stamen2","stamen1")

cor_matrix = cor(all_pca, method = "pearson")
pheatmap(cor_matrix, cluster_rows=F, cluster_cols=F, annotation_row=NULL, annotation_col=NULL)

cor_matrix_1 = as.data.frame(cor_matrix)
cor_matrix_1$CellType = rownames(cor_matrix_1) 
cor_matrix_m = melt(cor_matrix_1, id.vars = c("CellType"))
cor_matrix_m$CellType = factor(cor_matrix_m$CellType, levels = rev(c("stamen1", "stamen2","stamen3","stamen4",
                                                                     "petal","sepal")))
cor_matrix_m$variable = factor(cor_matrix_m$variable, levels = c("stamen1", "stamen2","stamen3","stamen4",
                                                                 "petal","sepal"))

ggplot(cor_matrix_m, aes(x=variable, y=CellType)) +
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low = "oldlace", high = "tomato", na.value = "whitesmoke") + 
  theme_bw() +
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 12, color = "black", face = "bold"),
        axis.text.x = element_text(margin = margin(t = 1, unit = "pt"), angle = 30, hjust = 1),
        axis.text.y = element_text(margin = margin(r = 0.1, unit = "pt")),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        legend.key.width = unit(0.3, "cm"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.title = element_text(size = 10, color = "black", face = "bold"),
        legend.text = element_text(size = 7),
        legend.justification = "top",            
        legend.margin = margin(0, 0, 0, 0),       
        legend.box.margin = margin(0, 0, 0, 0))    


#### ST exp count vs bulk exp count ####
# seapl to stamen raw ST RNA counts
sepal = subset(S3, sample == "sepal")
petal = subset(S3, sample == "petal")
stamen = subset(S3, sample %in% c("stamen1","stamen2","stamen3","stamen4"))

sepal_expr = as.matrix(t(GetAssayData(sepal, slot = "data")))
petal_expr = as.matrix(t(GetAssayData(petal, slot = "data")))
stamen_expr = as.matrix(t(GetAssayData(stamen, slot = "data")))

sepal_sum = Matrix::colSums(sepal_expr)
petal_sum = Matrix::colSums(petal_expr)
stamen_sum = Matrix::colSums(stamen_expr)

sepal_ST = data.frame(count = sepal_sum, row.names = colnames(sepal_expr))
sepal_ST = sepal_ST[order(rownames(sepal_ST)), , drop = FALSE]

petal_ST = data.frame(count = petal_sum, row.names = colnames(petal_expr))
petal_ST = petal_ST[order(rownames(petal_ST)), , drop = FALSE]

stamen_ST = data.frame(count = stamen_sum, row.names = colnames(stamen_expr))
stamen_ST = stamen_ST[order(rownames(stamen_ST)), , drop = FALSE]

#corrected gene count
Nym_gff = read.table("../Nymphaea_colorata.gff", sep = "\t")
colnames(Nym_gff) = c("chromsome","gene","start","end")
Nym_gff = Nym_gff[order(Nym_gff$gene), , drop = FALSE]
Nym_gff$length = Nym_gff$end - Nym_gff$start

sepal_ST_2 = log1p(sepal_ST*10^9/sum(sepal_ST)/Nym_gff$length)
petal_ST_2 = log1p(petal_ST*10^9/sum(petal_ST)/Nym_gff$length)
stamen_ST_2 = log1p(stamen_ST*10^9/sum(stamen_ST)/Nym_gff$length)

# seapl to stamen bulk RNA counts
sepal_bulk1 = read.table("../readcounts/L1MLC2503075-D1_Se_gene_abundances.tsv", sep = "\t", header = T)
sepal_bulk2 = read.table("../readcounts/L1MLC2503078-D2_Se_gene_abundances.tsv", sep = "\t", header = T)
sepal_bulk3 = read.table("../readcounts/L1MLC2503081-D3_Se_gene_abundances.tsv", sep = "\t", header = T)

petal_bulk1 = read.table("../readcounts/L1MLC2503076-D1_Pe_gene_abundances.tsv", sep = "\t", header = T)
petal_bulk2 = read.table("../readcounts/L1MLC2503079-D2_Pe_gene_abundances.tsv", sep = "\t", header = T)
petal_bulk3 = read.table("../readcounts/L1MLC2503082-D3_Pe_gene_abundances.tsv", sep = "\t", header = T)

stamen_bulk1 = read.table("../readcounts/L1MLC2503077-D1_St_gene_abundances.tsv", sep = "\t", header = T)
stamen_bulk2 = read.table("../readcounts/L1MLC2503080-D2_St_gene_abundances.tsv", sep = "\t", header = T)
stamen_bulk3 = read.table("../readcounts/L1MLC2503083-D3_St_gene_abundances.tsv", sep = "\t", header = T)

sepal_3rep = data.frame(sepal_rep1 = sepal_bulk1$FPKM,
                        sepal_rep2 = sepal_bulk2$FPKM,
                        sepal_rep3 = sepal_bulk3$FPKM,
                        row.names = sepal_bulk1$Gene.ID)
petal_3rep = data.frame(petal_rep1 = petal_bulk1$FPKM,
                        petal_rep2 = petal_bulk2$FPKM,
                        petal_rep3 = petal_bulk3$FPKM,
                        row.names = petal_bulk1$Gene.ID)
stamen_3rep = data.frame(stamen_rep1 = stamen_bulk1$FPKM,
                         stamen_rep2 = stamen_bulk2$FPKM,
                         stamen_rep3 = stamen_bulk3$FPKM,
                         row.names = stamen_bulk1$Gene.ID)

sepal_avg = log1p(rowMeans(sepal_3rep))
petal_avg = log1p(rowMeans(petal_3rep))
stamen_avg = log1p(rowMeans(stamen_3rep))

sepal_bulk = data.frame(count = sepal_avg, row.names = rownames(sepal_3rep))
sepal_bulk = sepal_bulk[order(rownames(sepal_bulk)), , drop = FALSE]

petal_bulk = data.frame(count = petal_avg, row.names = rownames(petal_3rep))
petal_bulk = petal_bulk[order(rownames(petal_bulk)), , drop = FALSE]

stamen_bulk = data.frame(count = stamen_avg, row.names = rownames(stamen_3rep))
stamen_bulk = stamen_bulk[order(rownames(stamen_bulk)), , drop = FALSE]

# ST vc bulk correlation
sepal_combined = cbind(sepal_ST_2, sepal_bulk)
colnames(sepal_combined) = c("ST_count", "bulk_count")
petal_combined = cbind(petal_ST_2, petal_bulk)
colnames(petal_combined) = c("ST_count", "bulk_count")
stamen_combined = cbind(stamen_ST_2, stamen_bulk)
colnames(stamen_combined) = c("ST_count", "bulk_count")

sepal_combined = subset(sepal_combined, ST_count > 0 & bulk_count > 0 )
petal_combined = subset(petal_combined, ST_count > 0 & bulk_count > 0 )
stamen_combined = subset(stamen_combined, ST_count > 0 & bulk_count > 0 )

ggplot(sepal_combined, aes(x = ST_count, y = bulk_count)) +
  geom_point(data = subset(sepal_combined, !target), 
             color = "gray", size = 0.2) +
  geom_point(data = subset(sepal_combined, target), 
             color = "red", size = 1.5) +
  geom_smooth(method = "rlm", se = FALSE, color = "blue", linewidth = 0.5)  +
  geom_text_repel(data = subset(sepal_combined, target),
                  aes(label = gene),
                  color = "black",
                  size = 3.5,
                  segment.color = "black",
                  segment.size = 0.5,
                  min.segment.length = 0.3,
                  box.padding = 0.5) +
  stat_cor(method = "spearman", label.x = 0, label.y = 12.5,size = 4) +
  labs(x = "ST_count", y = "bulk_count") +
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 12.5))+
  theme_bw()