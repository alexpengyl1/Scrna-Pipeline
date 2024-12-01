setwd("/Users/yuepen/Downloads/华哥2024-单细胞高阶班/资料/6")
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)

x=list.files()
x
#first merhod to load data
dir = c('BC2/', "BC21/")
names(dir) = c('s1',  's2')      
? CreateSeuratObject
counts <- Read10X(data.dir =dir)
scRNA1 = CreateSeuratObject(counts,min.cells = 3, min.features = 200)

#second method to load data

scRNA1[["percent.mt"]] <- PercentageFeatureSet(scRNA1, pattern = "^MT-")

VlnPlot(scRNA1,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt" ), 
        
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 3) 
scRNA2 <- subset(scRNA1, subset = nFeature_RNA > 500 & 
                   percent.mt < 20 &   nCount_RNA > 1000)

ifnb =scRNA2

ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$orig.ident)


# run standard anlaysis workflow
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))



################################
#Perform cca integration

ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
# Visualization
DimPlot(ifnb, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))


#########Perform integration with SCTransform-normalized datasets
#####################
# split datasets and process without integration
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$orig.ident)
ifnb <- SCTransform(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- RunUMAP(ifnb, dims = 1:30)
DimPlot(ifnb, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))

# integrate datasets
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, normalization.method = "SCT", verbose = F)
ifnb <- FindNeighbors(ifnb, reduction = "integrated.dr", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 0.6)


#############################
###harmony整合


scRNA.11=Read10X("BC2/")
scRNA.3=Read10X("BC21/")
scRNA.11 = CreateSeuratObject(scRNA.11 ,project="sample_11",min.cells = 3, min.features = 200)
scRNA.3 = CreateSeuratObject(scRNA.3 ,project="sample_3",min.cells = 3, min.features = 200)

scRNA_harmony <- merge(x=scRNA.11, y=c(scRNA.3 ))
scRNA_harmony[["RNA"]] <- JoinLayers(scRNA_harmony[["RNA"]])######Harmony需要在统一的空间中工作，JoinLayers提供了这个统一的空间，这样Harmony才能有效地进行批次校正
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)


library(harmony)
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")
#system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})

###一定要指定harmony
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 0.5)

scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:15)
?DimPlot
plot1 =DimPlot(scRNA_harmony, reduction = "umap",label = T) 
plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident') 
#combinate
plotc <- plot1+plot2
plotc


#############################
############################
markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)   
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)



scRNA_harmony@meta.data$seurat_clusters

scRNA_harmony <- RenameIdents(scRNA_harmony, "12" = "Macrophage","0"="Macrophage","2"="MSc")
DimPlot(scRNA_harmony,label = T,group.by = "seurat_clusters")
scRNA_harmony@meta.data$celltype=scRNA_harmony@active.ident


