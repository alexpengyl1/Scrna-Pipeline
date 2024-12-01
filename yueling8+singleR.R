setwd("E:\\super.lesson\\lesson7")
###加载所需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)

dir=c("CAFs.1/","CAFs.2/","dapi1/","dapi2/")
dir=c("CAFs.1/","CAFs.2/","dapi1/","dapi2/")
names(dir) = c('CAF1',  'CAF2', 'DapiNeg1',"DapiNeg2")      



scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells = 3, min.features =300)
}



for (i in 1:length(scRNAlist)) {
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst",nfeatures = 3000)
}


library(future)
plan("multiprocess", workers =4)
options(future.globals.maxSize = 25000 * 1024^2)

scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist,anchor.features = 3000)


plan("multiprocess", workers =1)
scRNA1 <- IntegrateData(anchorset = scRNA.anchors)
###################################################################################################################
##如果不整合
scRNA1=merge(scRNAlist[[1]],c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]]))

######################################################################################
DefaultAssay(scRNA1) <- "RNA"
scRNA1 <-  ScaleData(scRNA1,features = rownames(scRNA1))
DefaultAssay(scRNA1) <- "integrated"
scRNA1 <-  ScaleData(scRNA1)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst",nfeatures = 3000)
scRNA1 <- RunPCA(scRNA1, npcs = 20, verbose = T)
# t-SNE and Clustering

scRNA1 <- FindNeighbors(scRNA1, reduction = "pca", dims = 1:7)
scRNA1 <- FindClusters(scRNA1, resolution = 0.4)
scRNA1 <- RunUMAP(scRNA1, reduction = "pca", dims = 1:7)
colnames(scRNA1@meta.data)
scRNA1 <- RunTSNE(scRNA1, dims = 1:7)

DefaultAssay(scRNA1) <- "RNA"
##BiocManager::install("SingleR")
library(SingleR)
load("F:/ref.data/ref_Mouse_all.RData")
refdata <- ref_Mouse
###把rna的转录表达数据提取
testdata <- GetAssayData(scRNA1, slot="data")
###把scRNA数据中的seurat_clusters提取出来，注意这里是因子类型的
clusters <- scRNA1@meta.data$seurat_clusters
###开始用singler分析
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)
scRNA1@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA1@meta.data[which(scRNA1@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}


i=1
scRNA1@meta.data[which(scRNA1@meta.data$seurat_clusters == celltype$ClusterID[1]),'celltype'] <-  "Macrophages"

which(scRNA1@meta.data$seurat_clusters == "0")




DimPlot(scRNA1, reduction = "umap", group.by = "celltype",label = T)
DimPlot(scRNA1, reduction = "umap",label = T)

getwd()
save(scRNA1,file = "scRNA1.rdata")

load("scRNA1.rdata")
#################################################################################################
##scRNA1@active.ident=scRNA1@meta.data$celltype

DimPlot(scRNA1,group.by = "celltype")

DimPlot(scRNA1 )
colnames(scRNA1@meta.data)
Idents(scRNA1)= "celltype"
DimPlot(scRNA1)


degs <- FindAllMarkers(scRNA1, logfc.threshold = 0.5,
                       test.use = "roc", 
                       return.thresh = 0.25, 
                       min.pct = 0.3, only.pos = T) 

degs_sig <- degs %>% 
  filter(pct.1 > 0.3 &
           power > 0.25) %>%
  filter(cluster != "other") %>%
  arrange(cluster, -power)  

# select degs for heatmap
degs_top50 <- degs_sig %>% 
  #  filter(cluster!="other") %>%
  group_by(cluster) %>% 
  top_n(50, power) %>%
  top_n(50, avg_diff) %>%
  arrange(cluster, -power)


avgData <- scRNA1@assays$RNA@data[degs_top50$gene,] %>% 
  apply(1, function(x){
    tapply(x, scRNA1$celltype, mean) # ExpMean
  }) %>% t


phData <- MinMax(scale(avgData), -2, 2) # z-score
rownames(phData) <- 1:nrow(phData)

##install.packages("pheatmap")
library(pheatmap)

phres <- pheatmap(
  phData, 
  color = colorRampPalette(c("darkblue", "white", "red3"))(99), #配色
  scale = "row",
  cluster_rows = F, #不按行聚类
  cluster_cols = T, #按列聚类
  clustering_method = "complete",
  show_rownames = F, #显示cluster名
  annotation_row = data.frame(cluster = degs_top50$cluster), 
)  

phData1= phData[, c("Fibroblasts","Adipocytes","Macrophages","Granulocytes","Monocytes","T cells")] 


phres <- pheatmap(
  phData1, 
  color = colorRampPalette(c("darkblue", "white", "red3"))(99), #配色
  scale = "row",
  cluster_rows = F, #不按行聚类
  cluster_cols = F, #按列聚类
  clustering_method = "complete",
  show_rownames = F, #显示cluster名
  annotation_row = data.frame(cluster = degs_top50$cluster), 
) 

phres
########################################################################################################


load("scRNA1.rdata")

Idents(scRNA1)="celltype"

scRNA1=scRNA2

scRNA1@active.ident=factor(scRNA1@active.ident,levels = c("Adipocytes","Fibroblasts" , "Granulocytes", "Macrophages",   "Monocytes","T cells"))
table(scRNA1@active.ident)
marker=FindAllMarkers(scRNA1,logfc.threshold = 0.5)

write.csv(marker,file = "marker.csv")

top10 <-  marker  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DimPlot(scRNA1, reduction = "tsne", group.by = "celltype", pt.size=1)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
###################################################################################################


DoHeatmap(scRNA1,features = top10$gene,label = F,slot = "scale.data")
table(scRNA1$celltype)
scRNA2=scRNA1

table(scRNA1@active.ident)
scRNA1=subset(scRNA1, downsample = 150)
library(pheatmap)

colanno=scRNA1@meta.data 
colanno$barcode=rownames(colanno)
colanno=colanno%>%arrange(celltype)
rownames(colanno)=colanno$barcode
colanno$barcode=NULL


colanno$celltype=factor(colanno$celltype,levels =)
colanno1=colanno[,6]
colanno1=as.data.frame(colanno1)
rownames(colanno1)=rownames(colanno)
colnames(colanno1)="celltype"

rowanno=top10
colnames(top10)
rowanno=rowanno%>%arrange(cluster)

mat4=scRNA1[["RNA"]]@scale.data[rowanno$gene,rownames(colanno1)]
mat4[mat4>=2.5]=2.5
mat4[mat4 < (-1.5)]= -1  

pheatmap(mat4,cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno1,
         gaps_row=as.numeric(cumsum(table(rowanno$cluster))[-6]),
         gaps_col=as.numeric(cumsum(table(colanno$celltype))[-6]) 
         
)




pheatmap(mat4,cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno1,
         gaps_row=as.numeric(cumsum(table(rowanno$cluster))[-6]),
         gaps_col=as.numeric(cumsum(table(colanno$celltype))[-6]) ,
         filename="marker.heatmap.pdf",width=11,height = 7
)


Idents(scRNA1)="celltype"
table(scRNA1@active.ident)

scRNA1.1=subset(scRNA1, downsample = 100) 



colanno=scRNA1.1@meta.data 
colanno$barcode=rownames(colanno)
colanno=colanno%>%arrange(celltype)
rownames(colanno)=colanno$barcode
colanno$barcode=NULL


colanno$celltype=factor(colanno$celltype,levels = unique(colanno$celltype))
colanno1=colanno[,6]
colanno1=as.data.frame(colanno1)
rownames(colanno1)=rownames(colanno)
colnames(colanno1)="celltype"

rowanno=top10
colnames(top10)
rowanno=rowanno%>%arrange(cluster)

mat4=scRNA1.1[["RNA"]]@scale.data[rowanno$gene,rownames(colanno1)]
mat4[mat4>=2.5]=2.5
mat4[mat4 < (-1.5)]= -1  

pheatmap(mat4,cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno1,
         gaps_row=as.numeric(cumsum(table(rowanno$cluster))[-6]),
         gaps_col=as.numeric(cumsum(table(colanno$celltype))[-6]) 
         
)

################################################################################################


######################################################################################################
VlnPlot(scRNA1, features = top10$gene[1:3], pt.size = 0, ncol = 1)+
  scale_x_discrete("")+
  theme(
    axis.text.x.bottom = element_blank()
  )


top10 <-  marker  %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


library(reshape2)
vln.df=as.data.frame(scRNA1[["RNA"]]@data[top10$gene,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("barcode","exp")

anno=scRNA1@meta.data
anno$barcode=rownames(anno)


vln.df=inner_join(vln.df,anno,by="barcode")
##vln.df$gene=factor(vln.df$gene,levels = top10$gene[1:7]) #为了控制画图的基因顺序

vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_dark()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )


vln.df$celltype_gene=paste(vln.df$celltype,vln.df$gene,sep = "_")
stat.df=as.data.frame(vln.df%>%  dplyr::group_by(celltype,gene)%>%dplyr::summarize(mean=mean(exp)))
colnames(stat.df)


stat.df$celltype_gene=paste(stat.df$celltype,stat.df$gene,sep = "_")
stat.df=stat.df[,c("mean","celltype_gene")]
vln.df=inner_join(vln.df,stat.df,by="celltype_gene")
vln.df$mean=ifelse(vln.df$mean > 3, 3, vln.df$mean)

vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=mean),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_gradient(limits=c(0,3),low = "lightgrey",high = "yellow")+
  scale_x_discrete("")+scale_y_continuous("",expand = c(0.02,0))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1)
  )


###########################################################################################################
## 
DotPlot(scRNA1, features = top10$gene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+theme_bw()
# 

bubble.df=as.matrix(scRNA1[["RNA"]]@data[ top10$gene,])
bubble.df=t(bubble.df)
bubble.df=as.data.frame(scale(bubble.df))
bubble.df$CB=rownames(bubble.df)
meta=scRNA1@meta.data[,c(1,6)]
bubble.df=merge(bubble.df,meta,by.x = "CB",by.y = 0)
bubble.df$CB=NULL

celltype_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(bubble.df$celltype)) {
  bubble.df_small=bubble.df%>%filter(celltype==i)
  for (j in top10$gene) {
    exp_mean=mean(bubble.df_small[,j])
    exp_ratio=sum(bubble.df_small[,j] > min(bubble.df_small[,j])) / length(bubble.df_small[,j])
    celltype_v=append(celltype_v,i)  ##Add elements to a vector.
    gene_v=append(gene_v,j)
    mean_v=append(mean_v,exp_mean)
    ratio_v=append(ratio_v,exp_ratio)
  }
}

plotdf=data.frame(
  celltype=celltype_v,
  gene=gene_v,
  exp=mean_v,
  ratio=ratio_v
)
library(RColorBrewer)
plotdf$celltype=factor(plotdf$celltype,levels = sort(unique(plotdf$celltype)))
plotdf$gene=factor(plotdf$gene,levels = rev(as.character(top10$gene[1:10])))
plotdf$exp=ifelse(plotdf$exp>3,3,plotdf$exp)
plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete(" ")+scale_y_discrete(" ")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  scale_size_continuous(limits = c(0,1))+theme_classic()+
  theme(
    axis.text.x.bottom = element_text(hjust =1, vjust = 1, angle = 90)
  )






################################################################################
scRNA1=scRNA2

FeaturePlot(scRNA1,features = "Rbp1",reduction = "tsne",pt.size = 1,split.by = "orig.ident")


FeaturePlot(scRNA1,features = "Rbp1",reduction = "tsne",pt.size = 1)+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+ 
  theme(  
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),  
    axis.ticks = element_blank(),axis.text = element_blank(),  
    legend.position = "none",  
    plot.title = element_text(hjust = 0.5,size=14)  
  )

#### 

mat1=as.data.frame(scRNA1[["RNA"]]@data["Rbp1",])
colnames(mat1)="exp"
mat2=Embeddings(scRNA1,"tsne")
mat3=merge(mat2,mat1,by="row.names")


mat3%>%ggplot(aes(tSNE_1,tSNE_2))+geom_point(aes(color=exp))+
  scale_color_gradient(low = "grey",high = "purple")+theme_bw()+ 
  theme(  
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),  
    axis.ticks = element_blank(),axis.text = element_blank(),  
    legend.position = "none",  
    plot.title = element_text(hjust = 0.5,size=14)  
  )+
  scale_x_continuous("")+scale_y_continuous("")




library(ggpubr)

VlnPlot(subset(scRNA1,Rbp1 > 0 ), "Rbp1",group.by = "orig.ident")+ 
  theme_bw()+stat_compare_means(method = 't.test')+geom_boxplot()
sc.FURIN=subset(scRNA1,Rbp1 > 0 )
FURIN.exp=sc.FURIN@assays$RNA@data["Rbp1",]
FURIN.exp.m= sc.FURIN@meta.data
table(sc.FURIN$orig.ident)
FURIN.exp.m$exp=FURIN.exp
FURIN.exp.m$Cohort=factor(FURIN.exp.m$orig.ident,levels = c("DapiNeg1","DapiNeg2","CAF1",'CAF2'))
write.csv(FURIN.exp.m,file = "FURIN.exp.m.csv")
colnames(FURIN.exp.m)
library(ggplot2)
library("ggsignif")
my_comparisons=list(c("DapiNeg1","CAF1"),c("DapiNeg2",'CAF2'),
                    c("DapiNeg2" ,"CAF1"))

ggplot(FURIN.exp.m,aes(Cohort,exp,fill=Cohort)) +
  geom_violin()+theme_classic()+
  theme(plot.title=element_text(size = 15,face="bold"),
        axis.text.x=element_text(size=15,face="bold",angle=25,hjust = 1,vjust = 1),
        axis.text.y=element_text(size=18,face="bold"),
        axis.title.x=element_text(size = 0,vjust = 1,face="bold"),
        axis.title.y=element_text(size = 15,face="bold"))+
  labs(y= 'Expression',x="Group")+
  geom_signif(comparisons = my_comparisons,step_increase = 0.1,
              map_signif_level = F,test = t.test,size=1,textsize = 6)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  ggtitle("Rbp1")+
  theme(plot.title = element_text(size =18,hjust = 0.5))+
  guides(fill=guide_legend(title = " "))+ ##更改legend名字
  ##如果是color分组 需要 guides(color=guide_legend(title = "Group")) 
  theme(axis.line.x=element_line(linetype=1,color="black",size=1))+
  theme(axis.line.y=element_line(linetype=1,color="black",size=1))+
  theme(legend.text = element_text(size = 15) ) ##legend字体变大



