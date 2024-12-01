###代码 ：  https://github.com/ruoyan-li/RCC-spatial-mapping/blob/main/code/NichNet_Mac_RCC.R
###文章：  https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00548-7#sectitle0030 




##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())


library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)


setwd("D:/shangke/lession17/")
load("scRNA_harmony.Rdata")

library(SingleR)
load("D:/ref.data/ref_Human_all.RData")

refdata <- ref_Human_all
load("scRNA_harmony.Rdata")
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
DimPlot(scRNA_harmony,group.by = "orig.ident")
####################################################################################################################################################
setwd("D:/shangke/lession18/")
##library(devtools)
##install_github("saeyslab/nichenetr")
library(nichenetr)
# weighted_networks列表包含两个数据框，lr_sig是配体-受体权重信号网络，gr是配体-靶基因权重调控网络 
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5]

lr_network = readRDS("lr_network.rds")
head(lr_network)

weighted_networks = readRDS("weighted_networks.rds")
head(weighted_networks$lr_sig) 
scRNA_harmony@meta.data$celltype %>% table()

scRNA_harmony@meta.data$orig.ident %>% table()
DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")


Idents(scRNA_harmony) <- "celltype"


sender_celltypes = c("Chondrocytes", "Fibroblasts", "Monocyte", "T_cells" )
?nichenet_seuratobj_aggregate

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = scRNA_harmony, 
  receiver = "Endothelial_cells", 
  condition_colname = "orig.ident", condition_oi = "sample_11", condition_reference = "sample_3", 
  sender = sender_celltypes, 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
 ##用于预测活性靶基因和构建活性配体-受体网络的配体的数量默认为20个。
# 输出的是一个列表：
nichenet_output %>% names()

#配体活性分析结果。NicheNet做的第一件事，是根据预测的配体活性来确定配体的优先级。使用如下命令查看这些配体的排名:
nichenet_output$ligand_activities
# 查看top20 ligands
nichenet_output$top_ligands

##这些配体由一个或多个输入发送细胞表达。要看哪个细胞群表达了这些配体，你可以运行以下程序:
# 查看top20 ligands在各个细胞亚群中表达情况
DotPlot(scRNA_harmony, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
##RdYlBu也就是Rd红色,Yi黄色,Bu蓝色的过度

##观察这些配体在LCMV感染后是否有差异表达也是很有趣的。
DotPlot(scRNA_harmony, features = nichenet_output$top_ligands %>% rev(), split.by = "orig.ident") + RotatedAxis()
##
VlnPlot(scRNA_harmony, features = c( "CTGF" ,"INHBA","LAMB2","COL4A1" ), split.by = "orig.ident", pt.size = 0, combine = T)
##推断活跃的配体-靶标连接 NicheNet还推断出这些顶级配体的活性靶基因。
##要查看哪个顶级配体被预测调控了哪些差异表达基因的表达，可以运行以下命令来查看热图:

## 查看配体调控靶基因
nichenet_output$ligand_target_heatmap
##
nichenet_output$ligand_target_heatmap + 
  scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.0045,0.009)) + 
  xlab("Endothelial_cells") + 
  ylab("Prioritized  cell ligands")

# 查看top配体调控的靶基因及其评分
ligand.target.matrix <- nichenet_output$ligand_target_matrix

# 查看被配体调控靶基因的表达情况
DotPlot(scRNA_harmony %>% subset(idents = "Endothelial_cells"), 
        features = nichenet_output$top_targets, 
        split.by = "orig.ident") + RotatedAxis()

VlnPlot(scRNA_harmony %>% subset(idents = "Endothelial_cells"), features = nichenet_output$top_targets[1:5], 
        split.by = "orig.ident", pt.size = 0, combine = T, ncol = 8)


## 查看受体情况
# 查看配体-受体互作
nichenet_output$ligand_receptor_heatmap

ligand.receptor.matrix=nichenet_output$ligand_receptor_matrix


DotPlot(scRNA_harmony %>% subset(idents = "Endothelial_cells"), 
        features = nichenet_output$top_receptors, 
        split.by = "orig.ident") + RotatedAxis()

VlnPlot(scRNA_harmony %>% subset(idents = "Endothelial_cells"),  features = nichenet_output$top_receptors[1:5], 
        split.by = "orig.ident",pt.size = 0, combine = T, ncol = 8)


################################################################
