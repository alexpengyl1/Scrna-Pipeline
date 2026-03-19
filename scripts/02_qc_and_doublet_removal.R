remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

dir = c('BC2/')
#names(dir) = c('s1')      
? CreateSeuratObject
counts <- Read10X(data.dir =dir)

## Pre-process Seurat object (standard) 
seu_kidney <- CreateSeuratObject(counts)
seu_kidney <- NormalizeData(seu_kidney)
seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 2000)
seu_kidney <- ScaleData(seu_kidney)
seu_kidney <- RunPCA(seu_kidney)
seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep(seu_kidney, PCs = 1:10, sct = FALSE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)###这个就是pk值

## Assuming 7.5% doublet formation rate - tailor for your dataset 期待的双细胞数量
nExp_poi <- round(0.08*length(seu_kidney@meta.data$orig.ident)*length(seu_kidney@meta.data$orig.ident)/10000)  
#DoubletRate1 = ncol(seu_kidney)*8*1e-6 

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_kidney <- doubletFinder(seu_kidney, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
table(seu_kidney$DF.classifications_0.25_0.09_284) 


seu_kidney$doublet.class <- seu_kidney@meta.data[paste0("DF.classifications_0.25_0.09_",nExp_poi)]
seu_kidney@meta.data[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]] <- NULL
pann <- grep(pattern="pANN", colnames(seu_kidney@meta.data), value=TRUE)
seu_kidney$pANN <- seu_kidney[[pann]]
seu_kidney1 <- subset(seu_kidney, subset = doublet.class == "Singlet")
seu_kidney2 <- subset(seu_kidney, subset = doublet.class != "Doublet")


####这里可精简一下函数，可以减少pca和umap的信息减少数据大小。
seurat_slim <- DietSeurat(
  seu_kidney1,
  counts = TRUE,      # 保留counts
  data = TRUE,        # 保留normalized data
  scale.data = FALSE, # 不保留scaled data
  #dimreducs = c("pca","umap",)# 指定要保留的降维结果
  assays = "RNA"
)


#################多样本去除双细胞
#################多样本去除双细胞
#################多样本去除双细胞
#################多样本去除双细胞
#################多样本去除双细胞
#################多样本去除双细胞
#################多样本去除双细胞
#################多样本去除双细胞
#################多样本去除双细胞
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(DoubletFinder)
set.seed(1)


setwd("/Users/yuepen/Downloads/华哥2024-单细胞高阶班/资料/6")
data_directory=c("BC2/","BC21/")
project_name=c("sample2","sample21")


samples <- project_name

####定义函数make_seurat_object_and_doublet_removal
#################################################
make_seurat_object_and_doublet_removal <- function(data_directory, project_name){
  # function for basic seurat based qc and doubletfinder based doublet removal
  colon.data <- Read10X(data.dir = data_directory)
  currentSample <- CreateSeuratObject(counts = colon.data, project = project_name, min.cells = 3, min.features = 40)
  currentSample[["percent.mt"]] <- PercentageFeatureSet(currentSample, pattern = "^MT-")
  
  # qc plot-pre filtering
  pdf(paste0("./qc_plots_", project_name, "_prefiltered.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05))
  dev.off()
  pdf(paste0("./qc_plots_", project_name, "_prefiltered_no_points.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
  dev.off()
  
  # filter everything to 400 unique genes/cell
  currentSample <- subset(currentSample, subset = nFeature_RNA > 400 & nFeature_RNA < 4000  & percent.mt<20 & nCount_RNA>1000)
  
  # Normalize and make UMAP
  currentSample <- seurat_standard_normalize_and_scale(currentSample, FALSE)
  
  # Run doublet finder
  nExp_poi <- round(0.08*length(currentSample@meta.data$orig.ident)*length(currentSample@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  seu_colon <- doubletFinder(currentSample, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  print(head(seu_colon@meta.data))
  
  # rename columns
  seu_colon$doublet.class <- seu_colon@meta.data[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]]
  seu_colon@meta.data[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]] <- NULL
  pann <- grep(pattern="^pANN", x=names(seu_colon@meta.data), value=TRUE)
  seu_colon$pANN <- seu_colon[[pann]]
  seu_colon[[pann]] <- NULL
  
  
  # plot pre and post doublet finder results
  pdf(paste0("./UMAP_pre_double_removal", project_name, ".pdf"))
  print(DimPlot(seu_colon, reduction = "umap", group.by = "doublet.class", cols = c("#D51F26", "#272E6A")))
  dev.off()
  seu_colon <- subset(seu_colon, subset = doublet.class != "Doublet")
  pdf(paste0("./UMAP_post_double_removal", project_name, ".pdf"))
  print(DimPlot(seu_colon, reduction = "umap", cols = c("#272E6A")))
  dev.off()
  
  # Remove extra stuff and return filtered Seurat object
  seu_colon <- DietSeurat(seu_colon, counts=TRUE, data=TRUE, scale.data=FALSE, assays="RNA")
  return(seu_colon)
}

seurat_standard_normalize_and_scale <- function(colon, cluster, cluster_resolution){
  # colon is seurat object, 
  colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
  colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(colon)
  colon <- ScaleData(colon, features = all.genes)
  colon <- RunPCA(colon, features = VariableFeatures(object = colon))
  if (cluster){
    colon <- FindNeighbors(colon, dims = 1:20)
    colon <- FindClusters(colon, resolution = cluster_resolution)
  }
  colon <- RunUMAP(colon, dims = 1:20)
  return(colon)
}

#################################################
sample1 <- make_seurat_object_and_doublet_removal(data_directory[1], samples[1])


###  多个样本合并 
seu_list <- sample1

####接下来我们会从第二个样本开始
for (i in 2:length(samples)){
  
  
  
  sc.i = make_seurat_object_and_doublet_removal(data_directory[i], samples[i])
  seu_list=merge(seu_list,sc.i)
  
}

table(seu_list$orig.ident)


scRNA_harmony=seu_list
scRNA_harmony  <- NormalizeData(scRNA_harmony ) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
library(harmony)
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")
###问题 harmony之后的数据在哪里？


###一定要指定harmony
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:40) %>% FindClusters(resolution =1)

scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:40)

DimPlot(scRNA_harmony , reduction = "umap",label = T) 
DimPlot(scRNA_harmony, reduction = "umap", split.by ='orig.ident')

DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident')
table(scRNA_harmony$orig.ident)  

#######################铆钉点的方法进行整合
ifnb =seu_list

#ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$orig.ident)

###这是不整合的情况
###这是不整合的情况
###这是不整合的情况
###这是不整合的情况
###这是不整合的情况

# run standard anlaysis workflow
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))


######################这是要整合点情况
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
verbose = FALSE)###这一步是在在layers有多个数据的时候进行的，整合完啦再joinlayers就好了

# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
# Visualization
DimPlot(ifnb, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))

######################


