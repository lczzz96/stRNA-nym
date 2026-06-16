library(Seurat)
library(dplyr)
library(ggsci)
library(ggplot2)
library(cowplot)
library(patchwork)
library(plotrix)
library(ggforce)
library(ggpubr)
library(monocle)
library(RColorBrewer)
library(plotrix)
library(png)
library(clustree)
library(pheatmap)
library(grid)
library(reshape2)

set.seed(777)

S1flower=Read10X(data.dir = "../S1_flower/L7_flower/")
position1 = read.table("../S1_flower/L7_flower/barcodes_pos.tsv",
                       header = F, row.names = 1, sep = "\t")
S1flower=CreateSeuratObject(counts = S1flower) 
S1flower[['x']] = position1[1]
S1flower[['y']] = position1[2]
S1flower$sample = c("stage1")
png1 = readPNG("../S1_flower/L7_flower/he_roi_small_gray.png")

S2=Read10X(data.dir = "../S2/L4/")
position2 = read.table("../S2/L4/barcodes_pos.tsv",
                       header = F, row.names = 1, sep = "\t")
S2=CreateSeuratObject(counts = S2) 
S2[['x']] = position2[1]
S2[['y']] = position2[2]
S2$sample = c("stage2")
png2 = readPNG("../S2/L4/he_roi_small_gray.png")

S3=Read10X(data.dir = "../S3/L4/")
position3 = read.table("../S3/L4/barcodes_pos.tsv",
                       header = F, row.names = 1, sep = "\t")
S3=CreateSeuratObject(counts = S3) 
S3[['x']] = position3[1]
S3[['y']] = position3[2]
S3$sample = c("stage3")
png3 = readPNG("../S3/L4/he_roi_small_gray.png")

colnames(position1) = c("x","y")
colnames(position2) = c("x","y")
colnames(position3) = c("x","y")

steel_clu1 = read.table("../S1_flower/L7_flower/S1L7_flower.map.7", header = T)
S1flower$steel_clusters = as.factor(steel_clu1$Cluster)
Idents(S1flower)=factor(steel_clu1[,4], levels=c(1:max(steel_clu1[,4])))

steel_clu2 = read.table("../S2/L4/S2L4.map.18", header = T)
S2$steel_clusters = as.factor(steel_clu2$Cluster)
Idents(S2)=factor(steel_clu2[,4], levels=c(1:max(steel_clu2[,4])))

steel_clu3 = read.table("../S3/L4/S3L4.map.70", header = F)
S3$steel_clusters = as.factor(steel_clu3$V4)
Idents(S3)=factor(steel_clu3[,4], levels=c(1:max(steel_clu3[,4])))


sub1 = subset(S1flower, idents = c(3,7))
sub1@meta.data = sub1@meta.data %>% mutate(sample_clusters = paste0("stage1_", as.character(steel_clusters)))
Idents(sub1) = factor(sub1@meta.data[,8], levels = c("stage1_3","stage1_7"))

sub2 = subset(S2, idents = c(7,14))
sub2@meta.data = sub2@meta.data %>% mutate(sample_clusters = paste0("stage2_", as.character(steel_clusters)))
Idents(sub2) = factor(sub2@meta.data[,8], levels = c("stage2_7","stage2_14"))

sub5 = subset(S3, idents = c(36,37,38,39))
sub5@meta.data = sub5@meta.data %>% mutate(sample_clusters = paste0("stage3_", as.character(steel_clusters)))
Idents(sub5) = factor(sub5@meta.data[,8], levels = c("stage3_39","stage3_38","stage3_37","stage3_36"))

sub10 = subset(S3, idents = c(36,37,38,54,4,65,66,67,69,70))
sub10@meta.data = sub10@meta.data %>% mutate(sample_clusters = paste0("stage3_", as.character(steel_clusters)))

sub11 = subset(S3, idents = c(34,35,36,37,38,39,40,41,42))
sub11@meta.data = sub11@meta.data %>% mutate(sample_clusters = paste0("stage3_", as.character(steel_clusters)))


#### S1/2/3 merge meristem monocle ####
set.seed(777)
mergelist = list(sub1,sub2,sub5)
MergeData = merge(mergelist[[1]],
                  y = c(mergelist[[2]], mergelist[[3]]))

mnc1 = as(as.matrix(MergeData@assays$RNA@counts), 'sparseMatrix')
pdata = new('AnnotatedDataFrame', data = MergeData@meta.data)
fdata = data.frame(gene_short_name = row.names(mnc1),row.names = row.names(mnc1))
fdata = new('AnnotatedDataFrame', data = fdata)

cds = newCellDataSet(mnc1,
                     featureData = fdata,
                     phenoData = pdata,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
cds = detectGenes(cds, min_expr = 0.1)
disp_table = dispersionTable(cds)
unsup_clustering_genes = subset(disp_table, mean_expression >= 0.05 & dispersion_empirical > 0)
dim(unsup_clustering_genes)
cds = setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)
cds = reduceDimension(cds, max_components = 2, num_dim = 20,
                      reduction_method = 'DDRTree',
                      residualModelFormulaStr = "~sample + num_genes_expressed",
                      verbose = T)
cds = orderCells(cds)

plot_cell_trajectory(cds, color_by = "sample_clusters", cell_size = 1) + scale_color_manual(values = palette_84)
plot_cell_trajectory(cds, color_by = "sample_clusters", cell_size = 1) + facet_wrap(~sample_clusters, nrow = 3) + scale_color_manual(values = palette_84)
plot_cell_trajectory(cds, color_by = "sample", cell_size = 1) 
plot_cell_trajectory(cds, color_by = "sample", cell_size = 1) + facet_wrap(~sample, nrow = 1)
plot_cell_trajectory(cds, color_by = "State", cell_size = 1)+ scale_color_manual(values = palette_84)
plot_cell_trajectory(cds, color_by = "State", cell_size = 1) + facet_wrap(~State, nrow = 2) + scale_color_manual(values = palette_84)
plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 1) + scale_color_gradient(low = "darkblue",high = "lightblue") 
plot_complex_cell_trajectory(cds, color_by = "State", cell_size = 1) + scale_color_manual(values = palette_84)


cds@phenoData@data$State[cds@phenoData@data$State == 5] = 1
plot_cell_trajectory(cds, color_by = "State", cell_size = 1.5)+ 
  scale_color_manual(values = palette_84) +
  theme_void() + theme(legend.position = "none")

stage1_state = cds@phenoData@data$State[cds@phenoData@data$sample == "stage1"]
stage1_cell = colnames(sub1)
S1_state = data.frame(cellid = stage1_cell, State = stage1_state)

stage2_state = cds@phenoData@data$State[cds@phenoData@data$sample == "stage2"]
stage2_cell = colnames(sub2)
S2_state = data.frame(cellid = stage2_cell, State = stage2_state)

stage3_state = cds@phenoData@data$State[cds@phenoData@data$sample == "stage3"]
stage3_cell = colnames(sub5)
S3_state = data.frame(cellid = stage3_cell, State = stage3_state)

# 查看State空间分布
stage1_state = cds@phenoData@data$State[cds@phenoData@data$sample == "stage1"]
stage1_cell = colnames(sub1)
S1_state = data.frame(cellid = stage1_cell, State = stage1_state)
position1$cellid = rownames(position1)
sub1_position = S1_state %>%
  left_join(position1 %>% dplyr::select(cellid, x, y), by = "cellid")
ggplot(sub1_position, aes(x = y, y = x, color=State)) +
  background_image(png1)+
  geom_point(shape = 19, size = 1) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_blank()) +
  scale_x_continuous(expand=c(0,0),limits = c(0,1000))+
  scale_y_continuous(expand=c(0,0),limits = c(0,1000))+
  labs(title = 'seurat_clusters') +
  scale_color_manual(values = palette_84) +
  guides(color = guide_legend(override.aes = list(size = 6)))

stage2_state = cds@phenoData@data$State[cds@phenoData@data$sample == "stage2"]
stage2_cell = colnames(sub2)
S2_state = data.frame(cellid = stage2_cell, State = stage2_state)
position2$cellid = rownames(position2)
sub2_position = S2_state %>%
  left_join(position2 %>% dplyr::select(cellid, x, y), by = "cellid")
ggplot(sub2_position, aes(x = 0.978*y, y = x - 1.5, color=State)) +
  background_image(png2)+
  geom_point(shape = 19, size = 0.8) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_blank()) +
  scale_x_continuous(expand=c(0,0),limits = c(0,1000))+
  scale_y_continuous(expand=c(0,0),limits = c(0,1000))+
  labs(title = 'seurat_clusters') +
  scale_color_manual(values = palette_84) +
  guides(color = guide_legend(override.aes = list(size = 6)))

stage3_state = cds@phenoData@data$State[cds@phenoData@data$sample == "stage3"]
stage3_cell = colnames(sub5)
S3_state = data.frame(cellid = stage3_cell, State = stage3_state)
position3$cellid = rownames(position3)
sub5_position = S3_state %>%
  left_join(position3 %>% dplyr::select(cellid, x, y), by = "cellid")
ggplot(sub5_position, aes(x = x, y = 1000-y, color=State)) +
  background_image(png3)+
  geom_point(shape = 20, size = 0.8) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_blank()) +
  scale_x_continuous(expand=c(0,0),limits = c(0,1000))+
  scale_y_continuous(expand=c(0,0),limits = c(0,1000))+
  labs(title = 'seurat_clusters') +
  scale_color_manual(values = palette_84) +
  guides(color = guide_legend(override.aes = list(size = 6)))


# 查看Pseudotime空间分布
time_sub1 = cds@phenoData@data[,c("x","y","Pseudotime")][cds@phenoData@data$sample == "stage1",]
time_S1 = S1flower@meta.data[,c("x","y")] %>%
  left_join(time_sub1 %>% select(x,y,Pseudotime), by = c("x","y"))
ggplot(time_S1,aes(x=y,y=x,color=Pseudotime))+
  geom_point(shape = 20, size = 5.5) + 
  theme_void() +
  scale_color_gradient(
    name = "Pseudotime",
    low = "darkblue", 
    high = "lightblue",  
    na.value = "lightgray" 
  ) +
  theme(legend.position = "none")

time_sub2 = cds@phenoData@data[,c("x","y","Pseudotime")][cds@phenoData@data$sample == "stage2",]
time_S2 = S2@meta.data[,c("x","y")] %>%
  left_join(time_sub2 %>% select(x,y,Pseudotime), by = c("x","y"))
ggplot(time_S2,aes(x=y,y=x,color=Pseudotime))+
  geom_point(shape = 20, size = 2) + 
  theme_void() +
  scale_color_gradient(
    name = "Pseudotime",
    low = "darkblue",  
    high = "lightblue",  
    na.value = "lightgray" 
  ) +
  theme(legend.position = "none")

time_sub3 = cds@phenoData@data[,c("x","y","Pseudotime")][cds@phenoData@data$sample == "stage3",]
time_S3 = S3@meta.data[,c("x","y")] %>%
  left_join(time_sub3 %>% select(x,y,Pseudotime), by = c("x","y"))
ggplot(time_S3,aes(x=x,y=1000-y,color=Pseudotime))+
  geom_point(shape = 20, size = 1) + 
  theme_void() +
  scale_color_gradient(
    name = "Pseudotime",
    low = "darkblue",  
    high = "lightblue",  
    na.value = "lightgray"  
  ) +
  theme(legend.position = "none")


# 定义pseudotime起点
pseudotime = cds@phenoData@data[,c("sample","State","Pseudotime")]
pseudotime$cellid = rownames(pseudotime)
pseudotime$State = as.numeric(pseudotime$State)

S1_state2_cell = pseudotime[pseudotime$sample == "stage1" & pseudotime$State == "2",]
median(S1_state2_cell$Pseudotime)
cds$Pseudotime_new = abs(cds$Pseudotime - 6.625372)

plot_cell_trajectory(cds, color_by = "Pseudotime_new", cell_size = 1) +
  scale_color_gradient(low = "darkblue",high = "lightblue")

newtime_sub1 = cds@phenoData@data[,c("x","y","Pseudotime_new")][cds@phenoData@data$sample == "stage1",]
newtime_S1 = S1flower@meta.data[,c("x","y")] %>%
  left_join(newtime_sub1 %>% select(x,y,Pseudotime_new), by = c("x","y"))
ggplot(newtime_S1,aes(x=y,y=x,color=Pseudotime_new))+
  geom_point(shape = 20, size = 5.5) + 
  theme_void() +
  scale_color_gradient(
    name = "Pseudotime",
    low = "darkblue",  
    high = "lightblue",  
    na.value = "lightgray" 
  ) +
  theme(legend.position = "none")

newtime_sub2 = cds@phenoData@data[,c("x","y","Pseudotime_new")][cds@phenoData@data$sample == "stage2",]
newtime_S2 = S2@meta.data[,c("x","y")] %>%
  left_join(newtime_sub2 %>% select(x,y,Pseudotime_new), by = c("x","y"))
ggplot(newtime_S2,aes(x=y,y=x,color=Pseudotime_new))+
  geom_point(shape = 20, size = 2) + 
  theme_void() +
  scale_color_gradient(
    name = "Pseudotime",
    low = "darkblue",  
    high = "lightblue", 
    na.value = "lightgray"  
  ) +
  theme(legend.position = "none")

newtime_sub3 = cds@phenoData@data[,c("x","y","Pseudotime_new")][cds@phenoData@data$sample == "stage3",]
newtime_S3 = S3@meta.data[,c("x","y")] %>%
  left_join(newtime_sub3 %>% select(x,y,Pseudotime_new), by = c("x","y"))
ggplot(newtime_S3,aes(x=x,y=1000-y,color=Pseudotime_new))+
  geom_point(shape = 20, size = 1) + 
  theme_void() +
  scale_color_gradient(
    name = "Pseudotime",
    low = "darkblue", 
    high = "lightblue", 
    na.value = "lightgray" 
  ) +
  theme(legend.position = "none")


#### S3 meristem monocle 1####
mnc3 = as(as.matrix(sub11@assays$RNA@counts), 'sparseMatrix')
pdata = new('AnnotatedDataFrame', data = sub11@meta.data)
fdata = data.frame(gene_short_name = row.names(mnc3),row.names = row.names(mnc3))
fdata = new('AnnotatedDataFrame', data = fdata)

cds3 = newCellDataSet(mnc3,
                     featureData = fdata,
                     phenoData = pdata,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())
cds3 = estimateSizeFactors(cds3)
cds3 = estimateDispersions(cds3)
cds3 = detectGenes(cds3, min_expr = 0.1)
disp_table = dispersionTable(cds3)
unsup_clustering_genes = subset(disp_table, mean_expression >= 0.4 & dispersion_empirical >= dispersion_fit )
dim(unsup_clustering_genes)
cds3 = setOrderingFilter(cds3, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds3)
cds3 = reduceDimension(cds3, max_components = 3, num_dim = 8, 
                       reduction_method = 'DDRTree',
                       verbose = T)
cds3 = orderCells(cds3)

plot_cell_trajectory(cds3, color_by = "State", cell_size = 0.5) + scale_color_manual(values = palette_84)
plot_cell_trajectory(cds3, color_by = "State") + facet_wrap(~State, nrow = 1) + scale_color_manual(values = palette_84)
plot_cell_trajectory(cds3, color_by = "steel_clusters", cell_size = 0.5)+ scale_color_manual(values = palette_84)
plot_cell_trajectory(cds3, color_by = "Pseudotime", cell_size = 0.5) + NoLegend() +


cds_tmp = cds3
cds_tmp@phenoData@data$State[cds3@phenoData@data$State == 7] = 4
plot_cell_trajectory(cds_tmp, color_by = "State", cell_size = 1) + 
  scale_color_manual(values = c("deepskyblue2","plum3","violetred2","wheat2","olivedrab2","chartreuse4"))+
  theme_void() + NoLegend()

ggplot(cds_tmp@phenoData@data, aes(x = x, y = 1000-y, color=State)) +
  background_image(png3)+
  geom_point(shape = 19, size = 1.4) +
  scale_x_continuous(expand=c(0,0),limits = c(0,1000))+
  scale_y_continuous(expand=c(0,0),limits = c(0,1000))+
  scale_color_manual(values = c("deepskyblue2","plum3","violetred2","wheat2","olivedrab2","chartreuse4")) +
  theme_void() + theme(legend.position = "none")


#### S3 meristem monocle 2####
mnc4 = as(as.matrix(sub10@assays$RNA@counts), 'sparseMatrix')
pdata = new('AnnotatedDataFrame', data = sub10@meta.data)
fdata = data.frame(gene_short_name = row.names(mnc4),row.names = row.names(mnc4))
fdata = new('AnnotatedDataFrame', data = fdata)

cds4 = newCellDataSet(mnc4,
                     featureData = fdata,
                     phenoData = pdata,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())
cds4 = estimateSizeFactors(cds4)
cds4 = estimateDispersions(cds4)
cds4 = detectGenes(cds4, min_expr = 0.1)
disp_table = dispersionTable(cds4)
unsup_clustering_genes = subset(disp_table, mean_expression >= 0.05 & dispersion_empirical > 0)
dim(unsup_clustering_genes)
cds4 = setOrderingFilter(cds4, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds4)
cds4 = reduceDimension(cds4, max_components = 8, num_dim = 20, reduction_method = 'DDRTree', verbose = T)
cds4 = orderCells(cds4)
cds4 = orderCells(cds4, root_state = 11)

plot_cell_trajectory(cds4, color_by = "sample_clusters", cell_size = 1) + scale_color_manual(values = palette_84)
plot_cell_trajectory(cds4, color_by = "sample", cell_size = 1) 
plot_cell_trajectory(cds4, color_by = "sample", cell_size = 1) + facet_wrap(~sample, nrow = 1)
plot_cell_trajectory(cds4, color_by = "State", cell_size = 1)+ scale_color_manual(values = palette_84)
plot_cell_trajectory(cds4, color_by = "State", cell_size = 1) + facet_wrap(~State, nrow = 2) + scale_color_manual(values = palette_84)
plot_cell_trajectory(cds4, color_by = "Pseudotime", cell_size = 1) 


## 查看State空间分布
ggplot(cds4@phenoData@data, aes(x = x, y = 1000-y, color=State)) +
  background_image(png3)+
  geom_point(shape = 19, size = .5) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_blank()) +
  scale_x_continuous(expand=c(0,0),limits = c(0,1000))+
  scale_y_continuous(expand=c(0,0),limits = c(0,1000))+
  scale_color_manual(values = palette_84) +
  guides(color = guide_legend(override.aes = list(size = 6)))
