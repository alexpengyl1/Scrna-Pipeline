##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)  
library(magrittr)
library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)

setwd("E:/super.lesson/lesson3/")
load("scRNA_harmony.rdata")

library(SingleR)
load("F:/ref.data/ref_Human_all.RData")

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


setwd("E:/super.lesson/lesson9/")
sc.combined=scRNA1

table(sc.combined@meta.data$celltype)
table(sc.combined@meta.data$orig.ident)
type=c("Adipocytes","Fibroblasts","Granulocytes",
       "Monocytes", "T cells","Macrophages")


dir.create("deg/")
setwd("deg/")




library(future)
plan("multiprocess", workers =5)
options(future.globals.maxSize = 2000 * 1024^2)

r.deg=data.frame()
table(sc.combined@meta.data$orig.ident)

Idents(sc.combined)="orig.ident"
sc.combined <-RenameIdents(sc.combined,"CAF1"="CAF","CAF2"="CAF","DapiNeg1"="Dap","DapiNeg2"="Dap")
sc.combined@meta.data$group=sc.combined@active.ident

for (i in 1:length(type)) {
  Idents(sc.combined)="celltype"
  deg=FindMarkers(sc.combined,ident.1 = c("CAF"),ident.2 = c("Dap"),
                  group.by = "group",subset.ident =type[i]   )
  
  write.csv(deg,file = paste0( type[i],'deg.csv') )
  deg$gene=rownames(deg)
  deg$celltype=type[i]
  deg$unm=i-1
  r.deg=rbind(deg,r.deg)
  
}


table(r.deg$unm) 
####明显看到T cells组中CAF只有1个细胞，这导致了统计分析的错误。建议修改代码排除T cells：
type = c("Adipocytes", "Fibroblasts", "Granulocytes", "Monocytes", "Macrophages")
r.deg=data.frame()
for (i in 1:length(type)) {
  Idents(sc.combined) = "celltype"
  deg = FindMarkers(sc.combined, 
                    ident.1 = c("CAF"), 
                    ident.2 = c("Dap"),
                    group.by = "group", 
                    subset.ident = type[i])
  
  write.csv(deg, file = paste0(type[i], 'deg.csv'))
  deg$gene = rownames(deg)
  deg$celltype = type[i]
  deg$unm = i-1
  r.deg = rbind(r.deg,deg)
}
table(r.deg$unm) 
#############################################################################

# 根据自己计算的marker基因数量确定log2FC的阈值，这里先定为1.5
r.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
r.deg$threshold <- as.factor(ifelse(r.deg$avg_log2FC > 0 , 'Up', 'Down'))
dim(r.deg)

r.deg$adj_p_signi <- as.factor(ifelse(r.deg$p_val_adj < 0.01 , 'Highly', 'Lowly'))
r.deg$thr_signi <- paste0(r.deg$threshold, "_", r.deg$adj_p_signi)
r.deg$unm %<>% as.vector(.) %>% as.numeric(.)

##自定义显示想要展示的基因名
##这里挑选log2FC为top5的基因进行展示

top_up_label <- r.deg %>% 
  subset(., threshold%in%"Up") %>% 
  group_by(unm) %>% 
  top_n(n = 5, wt = avg_log2FC) %>% 
  as.data.frame()

top_down_label <- r.deg %>% 
  subset(., threshold %in% "Down") %>% 
  group_by(unm) %>% 
  top_n(n = -5, wt = avg_log2FC) %>% 
  as.data.frame()

top_label <- rbind(top_up_label,top_down_label)
top_label$thr_signi %<>% 
  factor(., levels = c("Up_Highly","Down_Highly","Up_Lowly","Down_Lowly"))

# 保存到文件，便于小白套用格式
write.csv(top_label, "easy_input_label.csv", quote = F)
##也可以基于output_pbmc.markers.csv文件，手动挑选出想要标注名字的基因，
#例如标注参与某一通路的基因，然后将文件命名为easy_input_label.csv


colnames(r.deg)
### 准备绘制暗灰色背景所需数据 
background_position <- r.deg %>%
  dplyr::group_by(unm) %>%
  dplyr::summarise(Min = min(avg_log2FC) - 0.2, Max = max(avg_log2FC) + 0.2) %>%
  as.data.frame()
## `summarise()` ungrouping output (override with `.groups` argument)
background_position$unm %<>% as.vector(.) %>% as.numeric(.)
background_position$start <- background_position$unm - 0.4
background_position$end <- background_position$unm + 0.4

### 准备绘制中间区域cluster彩色bar所需数据
cluster_bar_position <- background_position
cluster_bar_position$start <- cluster_bar_position$unm - 0.5
cluster_bar_position$end <- cluster_bar_position$unm + 0.5
cluster_bar_position$unm %<>% 
  factor(., levels = c(0:max(as.vector(.))))

## 设置填充颜色
cols_thr_signi <- c("Up_Highly" = "#d7301f",
                    "Down_Highly" = "#225ea8",
                    "Up_Lowly" = "black",
                    "Down_Lowly" = "black")
cols_cluster <- c("0" = "#35978f",
                  "1" = "#8dd3c7",
                  "2" = "#ffffb3",
                  "3" = "#bebada",
                  "4" = "#fb8072",
                  "5" = "#80b1d3",
                  "6" = "#fdb462",
                  "7" = "#b3de69","8" = "#b3de89")

p= ggplot() +
  geom_rect(data = background_position, aes(xmin = start, xmax = end, ymin = Min,
                                            ymax = Max),
            fill = "#525252", alpha = 0.1) + ###添加灰色背景色
  geom_jitter(data = r.deg, aes(x =unm, y = avg_log2FC, colour = thr_signi),
              size = 0.5,position = position_jitter(seed = 1)) +
  scale_color_manual(values = cols_thr_signi) +
  scale_x_continuous(limits = c(-0.5, max(r.deg$unm) + 0.5),
                     breaks = seq(0, max(r.deg$unm), 1),
                     label = seq(0, max(r.deg$unm),1)) + #修改坐标轴显示刻度
  
  # 根据top_label标注基因名
  geom_text_repel(data = top_label, aes(x =unm, y = avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), show.legend = F, size = 2.5,
                  box.padding = unit(0, "lines")) +
  
  geom_rect(data = cluster_bar_position, aes(xmin = start, xmax = end, ymin = -0.4,
                                             ymax = 0.4, fill = unm), color = "black", alpha = 1, show.legend = F) +
  scale_fill_manual(values = cols_cluster) +
  labs(x = "Cluster", y = "average log2FC") +
  theme_bw()
p
plot1 <- p + theme(panel.grid.minor = element_blank(), ##去除网格线
                   panel.grid.major = element_blank(),
                   axis.text.y = element_text(colour = 'black', size = 14),
                   axis.text.x = element_text(colour = 'black', size = 14, vjust = 66), #调整x轴坐标,vjust的值按照最终结果稍加调整
                   panel.border = element_blank(), ## 去掉坐标轴
                   axis.ticks.x = element_blank(), ## 去掉的坐标刻度线
                   axis.line.y = element_line(colour = "black")) #添加y轴坐标轴

plot1
ggsave(filename = "deg_pointplot.pdf", plot = plot1, width = 9, height = 6)
