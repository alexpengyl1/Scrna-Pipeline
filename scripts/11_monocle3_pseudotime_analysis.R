setwd("D:\\BaiduNetdiskDownload\\")

load("scRNA_harmony.Rdata")
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
load("D:\\BaiduNetdiskDownload\\ref_Human_all.RData")

refdata <- ref_Human_all

testdata <- GetAssayData(scRNA_harmony, slot="data")
###把scRNA数据中的seurat_clusters提取出来，注意这里是因子类型的
clusters <- scRNA_harmony@meta.data$seurat_clusters
###开始用singler分析
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)
scRNA_harmony@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

table(scRNA_harmony@meta.data$celltype)

data.sample=scRNA_harmony

DimPlot(data.sample,label=T,group.by = "celltype")
Idents(data.sample) <- data.sample$celltype

anno_sample <- data.sample@meta.data
###monocle3
data <- GetAssayData(data.sample,slot = "counts")
data <- data[rowSums(data>0)>=3,]
dim(data)
cell_metadata <- data.sample@meta.data
gene_annotation <- data.frame(gene_short_name=rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds,num_dim = 50)

cds <- reduce_dimension(cds,preprocess_method = "PCA")
p1 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "celltype")
p1
#########提取seurat和cds的umap信息
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(data.sample,reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "celltype")
p2

cds <- cluster_cells(cds,resolution = 0.01,k=20,random_seed=18,verbose=T)
cds@clusters$UMAP$clusters <- data.sample$celltype_level3_0912
plot_cells(cds,color_cells_by = "partition")

p1 <- plot_cells(cds,group_cells_by = 'cluster')
p1

cds <- learn_graph(cds, verbose =T,
                   use_partition=T,close_loop=F)
p <- plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=T,
                label_leaves=T, label_branch_points=T,cell_size = 0.5,group_label_size=4)

p

cds <- order_cells(cds)

p1 <- plot_cells(cds,color_cells_by = "pseudotime",label_branch_points = FALSE,label_leaves = F)
p1

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
p2 <- plot_cells(cds,color_cells_by = "celltype",
                 label_branch_points = FALSE,label_leaves = F)+
  scale_color_manual(values=my36colors)
p2|p1

ggsave("monocle3_FigS3G.pdf",width = 20,height = 8)

