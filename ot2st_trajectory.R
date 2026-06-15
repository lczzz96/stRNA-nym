library(Seurat)
library(dplyr)
library(ggsci)
library(ggplot2)
library(patchwork)
#devtools::load_all("/Library/Frameworks/R.framework/Versions/4.2/Resources/library/monocle")
library(monocle)
library(RColorBrewer)
library(plotrix)
library(ggforce)
library(ggpubr)
library(png)
library(clustree)
library(ggridges)
setwd("~/Desktop/Nymphaea_colorata/")
source("analysis/code/palette_84.R")
set.seed(1234)

#S3 sepaltostamen
S3=Read10X(data.dir = "../S3/L7/")
position3 = read.table("../S3/L7/barcodes_pos.tsv",header = F, row.names = 1, sep = "\t")
S3 = CreateSeuratObject(counts = data, project = "S3", min.cells = 3, min.features = 20)
S3[['x']] = position3[1]
S3[['y']] = position3[2]

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

S3 = subset(S3, sample!="all")
S3 = NormalizeData(S3)
S3 = FindVariableFeatures(S3, nfeatures = 2000, selection.method = "vst")
all.genes = rownames(S3)
S3 = ScaleData(S3, features = all.genes)
S3 = RunPCA(S3, features = VariableFeatures(object = S3))
S3 = FindNeighbors(S3, dims = 1:10)

mnc1 = as(as.matrix(S3@assays$RNA@counts), 'sparseMatrix')
pdata = new('AnnotatedDataFrame', data = S3@meta.data)
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
cds = reduceDimension(cds, max_components = 3,num_dim = 20,method = 'DDRTree', lambda = 20*ncol(cds))
cds = orderCells(cds)

plot_cell_trajectory(cds, color_by = "sample", cell_size = 0.25) + scale_color_manual(values = palette_84)
plot_cell_trajectory(cds, color_by = "sample", cell_size = 0.1) + facet_wrap(~sample, nrow = 1) + scale_color_manual(values = palette_84)
plot_cell_trajectory(cds, color_by = "State", cell_size = 0.15)+ scale_color_manual(values = palette_84)
plot_cell_trajectory(cds, color_by = "State", cell_size = 0.1) + facet_wrap(~State, nrow = 1) + scale_color_manual(values = palette_84)
plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.3) 

diff_test_res = differentialGeneTest(cds, core = 3, fullModelFormulaStr = "~State")
sig_genes = subset(diff_test_res, qval < 0.01)
sig_genes = sig_genes[order(sig_genes$qval,decreasing=F),]
nrow(sig_genes)

time_test_res = differentialGeneTest(cds, core = 5, fullModelFormulaStr = "~sm.ns(Pseudotime)")
time_genes = subset(time_test_res, qval < 0.01 & use_for_ordering == "TRUE")
nrow(time_genes)

plot_pseudotime_heatmap2=function (cds, cluster_rows = TRUE, hclust_method = "ward.D2", 
                                   num_clusters = 6, hmcols = NULL, add_annotation_row = NULL, 
                                   add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE, 
                                   norm_method = c("log", "vstExprs"), scale_max = 3, scale_min = -3, 
                                   trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = TRUE, 
                                   cores = 1) 
{
  num_clusters <- min(num_clusters, nrow(cds))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime), 
                                         max(pData(cds)$Pseudotime), length.out = 100))
  m <- genSmoothCurves(cds, cores = cores, trend_formula = trend_formula, 
                       relative_expr = T, new_data = newdata)
  m = m[!apply(m, 1, sum) == 0, ]
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs" && is.null(cds@dispFitInfo[["blind"]]$disp_func) == 
      FALSE) {
    m = vstExprs(cds, expr_matrix = m)
  }
  else if (norm_method == "log") {
    m = log10(m + pseudocount)
  }
  m = m[!apply(m, 1, sd) == 0, ]
  m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  m = m[is.na(row.names(m)) == FALSE, ]
  m[is.nan(m)] = 0
  m[m > scale_max] = scale_max
  m[m < scale_min] = scale_min
  heatmap_matrix <- m
  if (return_heatmap) {
    return(heatmap_matrix)
  }
}


exp_matrix = plot_pseudotime_heatmap2(cds[row.names(time_genes),])
site_vec=rep(1,nrow(exp_matrix))
site_vec_1=rep(1,nrow(exp_matrix))
site_vec_2=rep(1,nrow(exp_matrix))
for(i in 1:nrow(exp_matrix)) {
  max_value_1=max(exp_matrix[i,])
  site_1=which(exp_matrix[i,]==max_value_1)
  if(site_1[1]<50){
    site_vec_1[i]=site_1[1]
    max_value_2=exp_matrix[i,site_1[1]:(site_1[1]+5)]
    site_vec_2[i]=mean(max_value_2)
  }
  else{
    site_vec_1[i]=site_1[length(site_1)]
    max_value_2=exp_matrix[i,(site_1[1]-5):site_1[1]]
    site_vec_2[i]=mean(max_value_2)
  }
}

site_vec = as.data.frame(list(site_vec_1,site_vec_2));colnames(site_vec)=c('1','2')
site_vec$order = 1:nrow(exp_matrix)
site_vec[site_vec[,1]>50,] = site_vec[site_vec[,1]>50,][order(site_vec[site_vec[,1]>50,2]),]
site_vec[site_vec[,1]<50,] = site_vec[site_vec[,1]<50,][order(site_vec[site_vec[,1]<50,2],decreasing = T),]
site_vec = site_vec[order(site_vec[,1]),]
site_vec1 = site_vec[,1]
exp_matrix=exp_matrix[site_vec$order,]
bks = seq(-3.1, 3.1, by = 0.1)
hmcols = monocle:::blue2green2red(length(bks) - 1)
split=c(1,1450,2500,nrow(exp_matrix))
pheatmap(exp_matrix,scale = "none", cluster_cols = F, cluster_rows = F, 
         show_rownames = F, show_colnames = F,color = hmcols,
         gaps_row = split[2:(length(split)-1)],border=F)