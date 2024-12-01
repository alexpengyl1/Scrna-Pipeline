##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

###加载所需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(ggrepel)
library(reshape2)

setwd("E:/super.lesson/lesson3/")
load("scRNA_harmony.rdata")

table(scRNA_harmony@active.ident)

#######################################################################################
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

t.cells <- subset(scRNA_harmony, idents = "1")
colnames(scRNA_harmony@meta.data)
Idents(t.cells) <- "orig.ident"
table(t.cells@active.ident)
t.cells=subset(t.cells,ident=c("sample_11","sample_3"))

deg <- FindMarkers(t.cells,ident.1 = "sample_11",ident.2 = "sample_3", logfc.threshold =5
                   ,group.by = "orig.ident")
###deg <- FindMarkers(t.cells,ident.1 = c("sample_11","sample_5"),ident.2 = c("sample_3","sample_4"), logfc.threshold =5,group.by = "orig.ident")
###2个分组内有多个样本 使用向量比较就好了
?log1p

avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))###取log+1防止基因表达为0以及表达太高
avg.t.cells$gene <- rownames(avg.t.cells)

colnames(avg.t.cells)



genes.to.label = sample(rownames(deg),50)
colnames(avg.t.cells)=c("sample_11","sample_3","gene")
p1 <- ggplot(avg.t.cells, aes(sample_11,sample_3)) + geom_point() + ggtitle("T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1 

LabelPoints(plot = p1,points = c("RPS4Y1"),repel = TRUE)###查看某个单独的基因

###########################################################

degdf=avg.t.cells
degdf.l=degdf[genes.to.label,]
degdf.p=degdf[-match(genes.to.label,rownames(degdf)),]
degdf.p$Sig="NO_DIFF"
degdf.l$Sig="DIFF"
degdf=rbind(degdf.p,degdf.l)
colnames(degdf)

p=ggplot(degdf, aes(sample_11,sample_3)) +
  geom_point( aes(color=Sig)) +
  
  scale_color_manual(values=c("red","grey"))+
  ggtitle("Stromal Cells")+
  theme(plot.title = element_text(size =18,hjust = 0.5, face = "bold")) +
  
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14)) + 
  
  theme_bw()

p

p+theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))


data_selected <- degdf[genes.to.label,]
p + geom_label_repel(data=data_selected,
                     aes(label=rownames(data_selected)))


?geom_point
p +
  geom_point(size = 6, shape =2, data = data_selected[c("MMP2"),],colour="green") +
  ggrepel::geom_text_repel(
    aes(label =genes.to.label),size =2.5,
    data = data_selected,
    color="blue")+theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  ggtitle("Stromal Cells")+
  theme(plot.title = element_text(size =18,hjust = 0.5))


?geom_label_repel


#######################################################################################
##火山图
table(scRNA_harmony@active.ident)
cd4.naive=subset(scRNA_harmony,idents ="1")

degdf <- FindMarkers(cd4.naive,ident.1 = "sample_11",ident.2 = "sample_3", logfc.threshold = 0.01,group.by = "orig.ident")

degdf$symbol <- rownames(degdf)
logFC_t=0
P.Value_t = 1e-28
degdf$change = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < 0,"down",
                      ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 0,"up","stable"))

p=ggplot(degdf, aes(avg_log2FC,  -log10(p_val_adj))) +
  geom_point(alpha=0.4, size=2.8, aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("green", "grey","red"))+
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()

p=ggplot(degdf, aes(x=avg_log2FC, y= -log10(p_val_adj),size= pct.1)) +
  geom_point(alpha=0.4, aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("green", "grey","red"))+
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()######调整点的大小和pct1相关

p

p+theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))


data_selected <- degdf["TIMP1",]
p + geom_label_repel(data=data_selected,
                     aes(label=rownames(data_selected)))


?geom_point
p +
  geom_point(size = 4, shape =2, data = data_selected,colour="blue") +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = data_selected,
    color="yellow")
`###如何改变点的颜色

#############################################################################################
