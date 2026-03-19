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
DimPlot(scRNA_harmony,group.by = "celltype")

########################################################################################################################
#BiocManager::install("infercnv")
library("infercnv")

setwd("D:/shangke/lession19/") 
pos=read.table("human.gene.positions")
pos1=distinct(pos,V7,.keep_all = TRUE)
rownames(pos1)=pos1$V7
pos2=select(pos1,V7,V2,V3,V4)
View(pos2)
write.table(pos2, 'geneLocate.txt', row.names=F, col.names=F, sep='\t')

scRNA1=scRNA_harmony[,sample(colnames(scRNA_harmony),500)]
exprMatrix <- as.matrix(GetAssayData(scRNA1, slot='counts'))
cellAnnota <- subset(scRNA1@meta.data, select='seurat_clusters')
groupFiles='groupFiles.txt'
dim(exprMatrix)
write.table(cellAnnota,file =" groupFiles.txt",sep = '\t',col.names = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprMatrix,
                                    annotations_file=" groupFiles.txt",
                                    delim="\t",
                                    gene_order_file= "geneLocate.txt",
                                    ref_group_names=NULL)

#ref_group_names参数根据细胞注释文件填写，在示例中，这两种细胞是非恶性细胞，所以作为参照；
#ref_group_names=NULL，则会选用所有细胞的平均表达作为参照，这种模式下，最好能确保细胞内在有足够的差异
###比如如果其中有自己正常细胞t cells，可以写成ref_group_names=t cells

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics 
                             out_dir=  'cnv1/' ,  # dir is auto-created for storing outputs
                             cluster_by_groups=F,   #  cluster_by_groups：先区分细胞来源，再做层次聚类
                             hclust_method="ward.D2", plot_steps=F，
                             write_expr_matrix = T#####这是新版本 需要手动添加的
)


######################################################################################3

cellAnnota <- subset(scRNA1@meta.data, select='celltype')
groupFiles='groupFiles1.txt'
dim(exprMatrix)
write.table(cellAnnota,file =" groupFiles1.txt",sep = '\t',col.names = F)

table(cellAnnota$celltype)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprMatrix,
                                    annotations_file=" groupFiles1.txt",
                                    delim="\t",
                                    gene_order_file= "geneLocate.txt",
                                    ref_group_names="T_cells")

infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics 
                              out_dir=  'cnv.REF/' ,  # dir is auto-created for storing outputs
                              cluster_by_groups=T,   #  cluster_by_groups：先区分细胞来源，再做层次聚类
                              hclust_method="ward.D2", plot_steps=F,denoise=TRUE,
                              HMM=F,  ##特别耗时间
                              num_threads=4，
                              write_expr_matrix = T#####这是新版本 需要手动添加的
)


#################################################################
#接下来我们用infercnv.observations.txt来计算每个细胞的CNV score。
#这里我们根据经验设置一系列的threshold来判断对于每个基因其
#copy+number是不变（Neutral），减少（Copy+loss）或增多（Copy+gain）
#。Copy+neutral+其score为0，Copy+loss和Copy+gain都认为是有event，
#其score不为0。这里简单粗暴的分了五类，具体如何更好的计算score


grp=read.table("cnv.REF/infercnv.observation_groupings.txt",sep = "",header = T)
obs=read.table("cnv.REF/infercnv.observations.txt", header = T,check.names = F)
max(obs)
min(obs)
obs[obs>0.9 & obs<0.93]=2
obs[obs>=0.93 & obs<0.95]=1
obs[obs>=0.95 & obs<1.05]=0
obs[obs>=1.05 & obs<1.07]=1
obs[obs>=1.07 & obs<1.1]=2

scores=as.data.frame(colSums(obs))
scores$cluster=grp$Annotation.Group
colnames(scores)=c("score","cluster")

library(ggpubr)
ggboxplot(scores,"cluster","score",fill = "cluster")


###问题 如果直接计算score可以吗？
grp=read.table("cnv.REF/infercnv.observation_groupings.txt",sep = "",header = T)
obs=read.table("cnv.REF/infercnv.observations.txt", header = T,check.names = F)
scores=as.data.frame(colSums(obs))
scores$cluster=grp$Annotation.Group
colnames(scores)=c("score","cluster")

library(ggpubr)
ggboxplot(scores,"cluster","score",fill = "cluster")


########################################################################################################
##################################这些可视化是旧版的，以下是新版本的
cnv_score_table = data.table::fread("cnv.REF/infercnv.observations.txt",
                                    data.table=F)%<%
  column_to_rownames(var ="v1")

library(scales)
cnvScore <- function(data)
{data <- data %>% as.matrix() %>%
  t()  %>%
  scale()%>%
  rescale(to=c(-1,1))%>%
  t()
cnv_score<- as.data.frame(colSums(data*data))
return(cnv_score)
}
cnv_score <-cnvScore(cnv_score_table)
library(ggpubr)
#将CNV分数与细胞类型信息合井
cnv_score <- cbind(cnv_score, cellAnnota[row.names(cnv_score),])
names(cnv_score)=c("cnv_score","celltype")
#######

color <- ggsci::pal_aaas()(10)
ggplot(cnv_score, aes(reorder(celltype, cnv_score),cnv_score, color = celltype))+
  geomboxplot() +
  # scale color manual(values =color)+
  theme(panel.background = element_blank())+
  theme(axis.line = element_line(colour ="black"))+
  theme(axis.title.x= element_blank())+
  theme(legend.position ="NA")+
  labs(x ="",y="CNV_Scores",title="")+
  theme(axis.title.y=element_text(size=15))+
  theme(axis.text.x =element_text(size=15,angle =45,vjust =1,hjust = 1))+
  stat_compare_means()

###
library(ggplot2)
library(ggridges)
ggplot(cnv_score, aes(x=cnv_score,y= celltype))+ geom_density_ridges()

ggplot(cnv_score, aes(x=cnv_score,y= celltype))+ geom_density_ridges()

ggplot(cnv_score, aes(x=cnv_score,y= celltype, fill =0.5- abs(0.5 - stat(ecdf))))+
  stat_density_ridges(geom ="density_ridges_gradient", calc_ecdf = TRUE)+
  scale_fill_viridis_c(name="Tail probability", direction =-1)+
  theme_bw(base_size-10)+
  theme(axis.text = element_text(colour="black"))

#################################
library(RColorBrewer)
infercnv::plot_cnv(infercnv_obj2, #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = "better_plot",output_format = "pdf", #保存为pdf文件
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2))) #改颜色





