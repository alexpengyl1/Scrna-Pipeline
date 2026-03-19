##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

setwd("E:\\super.lesson\\lesson7")
###加载所需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(ggplot2)
load("scRNA1.rdata")
Idents(scRNA1)="orig.ident"
table(scRNA1$orig.ident)

scRNA1=RenameIdents(scRNA1,"CAF1"="CAF",  "CAF2"="CAF",
                    "DapiNeg1"="DapiNeg", "DapiNeg2"="DapiNeg")


table(scRNA1@active.ident)
scRNA1@meta.data$condition=scRNA1@active.ident

setwd("/Users/yuepen/Downloads/华哥2024单细胞高阶班/资料/12")
##BiocManager::install("AUCell")
library(AUCell)
library(ggplot2)
library(Seurat)
##mn BiocManager::install("clusterProfiler")
library(clusterProfiler)

#sc.id=sample(colnames(scRNA1),1500)
#sc2=scRNA1[,sc.id]
##install.packages("doParallel")
##install.packages("doRNG")
sc2=scRNA1
exp = GetAssayData(sc2,slot = "data")####此处修改了代码，因为seurat5版本多了一个layers如果使用以前的cells_rankings <- AUCell_buildRankings(sc2@assays$RNA@data,  nCores=6, plotStats=TRUE) 会没有基因名和细胞名
cells_rankings <- AUCell_buildRankings(exp,  nCores=6, plotStats=TRUE) 

cells_rankings

c2 <- read.gmt("c2.cp.kegg_medicus.v2023.2.Hs.symbols (2).gmt") 
geneSets <- lapply(unique(c2$term), function(x){print(x);c2$gene[c2$term == x]})
names(geneSets) <- unique(c2$term)
?AUCell_calcAUC
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,nCores =1,
                            aucMaxRank=nrow(cells_rankings)*0.1)
####################################################以下是使用msigdbr构建cells_AUC
##install.packages("msigdbr")
library(msigdbr)
msigdbr_species()

m_df<- msigdbr(species = "Mus musculus",  category = "C2", subcategory = "KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)


cells_AUC <- AUCell_calcAUC(fgsea_sets, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)
#####################################################

grep("OXIDATIVE",rownames(cells_AUC@assays@data$AUC),value = T)

geneSet <- "KEGG_THYROID_CANCER"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
sc2$AUC <- aucs
df<- data.frame(sc2@meta.data, sc2@reductions$umap@cell.embeddings)
colnames(df)

class_avg <- df %>%
  group_by(seurat_clusters) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

ggplot(df, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="D") +
  ggrepel::geom_label_repel(aes(label = seurat_clusters),
                            data = class_avg,
                            size = 3,
                            label.size = 1,
                            segment.color = NA
  )+   theme(legend.position = "none") + theme_bw() + facet_grid(.~condition)

colnames(df)

class_avg <- df %>%
  group_by(celltype) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

ggplot(df, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour = AUC)) + 
  viridis::scale_color_viridis(option="D") +
  ggrepel::geom_label_repel(aes(label = celltype),
                            data = class_avg,
                            size = 3,
                            label.size = 1,
                            segment.color = NA) + 
  theme(legend.position = "none") + 
  theme_bw() + 
  facet_grid(.~condition)


ggplot(data.frame(sc2@meta.data, sc2@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=AUC)
) + geom_point( size=1.5
) + scale_color_viridis(option="A")  + theme_light(base_size = 26) + facet_grid(.~condition)

DimPlot(sc2, reduction = "umap",label = T,split.by = "condition", group.by = "celltype") 
DimPlot(sc2, reduction = "umap", group.by = "celltype", label = FALSE) +
  scale_color_viridis(discrete = TRUE, option = "A") +
  theme_light(base_size = 26)
##############################################以下是寻找比较同一个细胞群在不同condition条件下的差异通路，based on aucell结果

# Subset cells that are Macrophages
macro_cells <- WhichCells(sc2, expression = celltype == "Fibroblasts")
macro_data <- sc2[, macro_cells]

# Get AUCell scores for these cells
macro_auc <- getAUC(cells_AUC)[, macro_cells]

# Split by condition
caf_cells <- WhichCells(macro_data, expression = condition == "CAF")
dapi_cells <- WhichCells(macro_data, expression = condition == "DapiNeg")

# Calculate mean AUC scores for each condition
caf_means <- rowMeans(macro_auc[, caf_cells])
dapi_means <- rowMeans(macro_auc[, dapi_cells])

# Perform Wilcoxon test for each pathway
pathway_pvals <- sapply(rownames(macro_auc), function(pathway) {
  wilcox.test(macro_auc[pathway, caf_cells], 
              macro_auc[pathway, dapi_cells])$p.value
})

# Create results dataframe
results <- data.frame(
  pathway = names(pathway_pvals),
  pvalue = pathway_pvals,
  CAF_mean = caf_means,
  DapiNeg_mean = dapi_means,
  diff = caf_means - dapi_means
)

# Add adjusted p-values
results$padj <- p.adjust(results$pvalue, method = "BH")

# Sort by absolute difference
results <- results[order(abs(results$diff), decreasing = TRUE), ]

# View top differential pathways
head(results)

# Optional: Visualize top pathways
top_n <- 20
top_paths <- head(results, top_n)

ggplot(top_n, aes(x = reorder(pathway, diff), y = diff)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(x = "Pathway", y = "Difference (CAF - DapiNeg)",
       title = "Top Differential Pathways in Macrophages") +
  theme(axis.text.y = element_text(size = 8))


# Get top 20 and bottom 20 pathways by diff
top_20 <- head(results[order(results$diff, decreasing = TRUE), ], 20)
bottom_20 <- head(results[order(results$diff), ], 20)
sig_paths <- rbind(top_20, bottom_20)

# Add significance indicator
sig_paths$significance <- ifelse(sig_paths$padj < 0.05, "FDR < 0.05", "Not Significant")######sig_paths$pvalue

# Visualization 
ggplot(sig_paths, aes(x = reorder(pathway, diff), y = diff)) +
  geom_bar(stat = "identity", aes(fill = significance)) +
  coord_flip() +
  theme_minimal() +
  labs(x = "Pathway", y = "Difference (CAF - DapiNeg)") +
  theme(axis.text.y = element_text(size = 8))

ggplot(sig_paths, aes(x = reorder(pathway, diff), y = diff)) +
  geom_bar(stat = "identity", aes(fill = diff > 0)) +
  scale_fill_manual(values = c("lightblue", "#CD5C5C")) +
  coord_flip() +
  theme_minimal() +
  labs(x = "Pathway", y = "Difference (CAF - DapiNeg)") +
  theme(axis.text.y = element_text(size = 8)) +
  theme(legend.position = "none")

####################################################################################################
#####以下是进行gsva分析
###  gsva
library('GSEABase')
library(GSVA)
library(msigdbr)
m_df<- msigdbr(species ="Mus musculus",  category = "H" )
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
colnames(scRNA1@meta.data)
table(scRNA1$celltype)
Idents(scRNA1)="celltype"
?AverageExpression
exp=AverageExpression(scRNA1) 
exp=exp[["RNA"]]
exp= as.matrix(exp)

?gsva
GSVA_hall <- gsvaParam(expr=as.matrix(exp), 
                       geneSets=geneSets, 
                  kcdf="Gaussian" 
                  #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  ) # 并行线程数目
head(GSVA_hall)
gsva_scores <- gsva(GSVA_hall)
?pheatmap::pheatmap
pheatmap::pheatmap(gsva_scores, #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =F,#列聚类，可以看出样本之间的区分度
                   show_colnames=T,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("#FF7744", "white","#AAAAAA","#0044BB"))(100))






###################这里是进行多组的分析
exp=AverageExpression(scRNA1,group.by = c("condition","celltype"))
exp=exp[["RNA"]]
table(exp@Dimnames[[2]])
counts2=exp[,c(1:3,7:9)]########在这里手动提取需要的分组条件下细胞的种类

GSVA_hall <- gsvaParam(expr=as.matrix(counts2), 
                  geneSets=geneSets, 
                  # 数据为正态分布则T，双峰则F
                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                 ) # 并行线程数目

gsva_scores <- gsva(GSVA_hall)
head(gsva_scores)
pheatmap::pheatmap(gsva_scores, #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =F,#列聚类，可以看出样本之间的区分度
                   
                   show_colnames=T,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("#FF7744", "white","#AAAAAA","#0044BB"))(100))

#################################################################################################




############################################################################################

######以下是对两个分组进行gsva分析，然后完成barplot的差异可视化

##install.packages("msigdbr")



##BiocManager::install("GSVA")
library('GSEABase')
library(GSVA)
library(msigdbr)
m_df<- msigdbr(species ="Homo sapiens",  category = "H" )
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

scRNA1=scRNA_harmony
table(scRNA1$celltype,scRNA1$orig.ident)

Idents(scRNA1)="celltype"
sc.b=subset(scRNA1,ident="Fibroblasts")
table(sc.b$orig.ident)
Idents(sc.b)="orig.ident"
sc.b=subset(sc.b,ident=c("CAF1","DapiNeg1"))
table(sc.b$orig.ident)
exp=sc.b@assays$RNA@data
exp=as.matrix(exp)


GSVA_hall <- gsvaParam(expr=exp, 
                       geneSets=geneSets, 

                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                   ) # 并行线程数目
GSVA_hall <- gsva(GSVA_hall)
head(GSVA_hall)

## limma差异通路分析

#BiocManager::install('limma')
library(limma)
# 设置或导入分组

group <- factor(sc.b@meta.data$orig.ident, levels = c( 'CAF1','DapiNeg1'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall)
design
# Tunor VS Normal
compare <- makeContrasts(CAF1 - DapiNeg1, levels=design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

## 发散条形图绘制
## barplot
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
# 去掉"HALLMARK_"
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
# 新增一列 根据t阈值分类
dat_plot$threshold = factor(ifelse(dat_plot$t  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 排序
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
library(ggplot2)

##install.packages("ggthemes")
library(ggthemes)
# install.packages("ggprism")
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-1,1),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, S2 VS S1') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
# 添加标签

# 小于-1的数量
low1 <- dat_plot %>% filter(t < -21) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于1总数量
high0 <- dat_plot %>% filter(t < 1) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签

p

ggsave("gsva_bar.pdf",p,width = 16,height  = 8)

#############使用pdf函数可以保存所有的图像但是gggave有可能会丢失
pdf("gsva_bar.pdf", width = 10, height = 8, onefile = TRUE)
print(p)
dev.off()

#####################################################################################################


degdf <- FindMarkers(scRNA1,ident.1 = "DapiNeg1",ident.2 = "DapiNeg2", 
                     logfc.threshold = 0.5,group.by = "orig.ident",subset.ident="Fibroblasts")
##BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
library(clusterProfiler)
degs.list=rownames(degdf)
erich.go.BP = enrichGO(gene =degs.list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
dotplot(erich.go.BP,showCategory = 10 )
?dotplot
barplot(erich.go.BP,showCategory = 8)
erich.go.BP=erich.go.BP@result
write.table(erich.go.BP,"4.7.erich.go.BP.deg.con.hf.txt",sep = "\t",col.names = NA)




kegg=read.table("4.7.go.txt",header = T,sep = "\t")
View(kegg)

k = data.frame(kegg)
library(ggplot2)
library(dplyr)
before <- as.numeric(sub("/\\d+$", "", k$GeneRatio))
after <- as.numeric(sub("^\\d+/", "", k$GeneRatio))
k$GeneRatio = before /after
font.size =10

k %>% 
  ## 对进行p值排序
  arrange(p.adjust) %>% 
  ##指定富集的通路数目
  dplyr::slice(1:14) %>% 
  ## 开始ggplot2 作图，其中fct_reorder调整因子level的顺序
  ggplot(aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  ## 画出点图
  geom_point(aes(color=p.adjust, size = Count)) +
  ## 调整颜色，guide_colorbar调整色图的方向
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ## 调整泡泡的大小
  scale_size_continuous(range=c(3, 8))+
  ## 如果用ylab("")或出现左侧空白
  labs(y=NULL) +
  ## 如果没有这一句，上方会到顶
  ggtitle("")+
  ## 设定主题
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))




erich.go.CC = enrichGO(gene =degs.list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.CC,showCategory = 8)
?dotplot
barplot(erich.go.CC,showCategory = 8)

erich.go.MF = enrichGO(gene =degs.list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.MF,showCategory = 8)
?dotplot
barplot(erich.go.MF,showCategory = 8)

keytypes(org.Hs.eg.db)


DEG.entrez_id = mapIds(x = org.Mm.eg.db,
                       keys =  degs.list,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

## install.packages("R.utils")

library(R.utils)

R.utils::setOption("clusterProfiler.download.method","auto")

erich.kegg.res <- enrichKEGG(gene = DEG.entrez_id,
                             organism = "mmu",
                             keyType = "kegg")

dotplot(erich.kegg.res)
kegg.res.df=erich.kegg.res@result
write.table(kegg.res.df,"erich.kegg.res.txt",sep = "\t",col.names = NA)

####
####################################################################################################

##pathway=read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
library(msigdbr)
m_df<- msigdbr(species ="Mus musculus",  category = "H" )
m_df = m_df[,c(3,4)] 
#########enricher可以做很多中富集
?enricher
res=enricher( degs.list,TERM2GENE =m_df )

dotplot(res)


### kegg
m_df.kegg<- msigdbr(species ="Mus musculus",  category = "C2", subcategory = "KEGG" )
m_df.kegg = m_df.kegg[,c(3,4)] #第三列是通路，第四列是基因

res.kegg=enricher( degs.list,TERM2GENE =m_df.kegg )

dotplot( res.kegg)



####################################################################################


##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

setwd("E:/super.lesson/lesson3/")
load("scRNA_harmony.Rdata")
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
#BiocManager::install("SingleR")
library(SingleR)
load("F:/ref.data/ref_Human_all.RData")


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
table(scRNA_harmony@meta.data$seurat_clusters)
sc.t=scRNA_harmony[,rownames(subset(scRNA_harmony@meta.data,celltype=="T_cells"))]  

######################################################################################################

##BiocManager::install("org.Hs.eg.db")
##BiocManager::install("clusterProfiler")

table(scRNA_harmony$orig.ident)
deg =FindMarkers(scRNA_harmony,ident.1 = "sample1",ident.2 = "sample2",group.by = "orig.ident",
                 subset.ident = "0")
degs.list=rownames(deg)
##BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)
library(clusterProfiler)
erich.go.BP = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

dotplot(erich.go.BP,showCategory = 30)

erich.go.BP=erich.go.BP@result

###################################################################################另一个版本的go分析，做圈图
colnames(erich.go.BP)
###install.packages("GOplot")
library(GOplot)
colnames(erich.kegg.res)
setwd("E:/super.lesson/lesson10/")

go1=erich.go.BP[1:3,c(1,2,11,8)]######提取ID，description，geneID，pvalue
rownames(go1)=go1$ID
go=go1
###install.packages("stringr")
library(stringr)
go$geneID=str_replace_all(go$geneID,"/",",")
names(go)=c('ID','Term','Genes','adj_pval')
go$Category="BP"



x1=strsplit(go$Genes[1],split=",",fixed=T)
x2=strsplit(go$Genes[2],split=",",fixed=T)
x3=strsplit(go$Genes[3],split=",",fixed=T)
g1=c(x1[[1]],x2[[1]],x3[[1]])

deg=degdf
genedata1=deg[g1,]   
genedata1$ID=rownames(genedata1)
genedata2=data.frame(ID=genedata1$ID,logFC=genedata1$avg_log2FC)
genedata=data.frame(ID=degs.list,logFC=deg[degs.list,]$avg_log2FC)


circ <- circle_dat(go,genedata2)
##条形图
GOBar(subset(circ, category == 'BP'))
#气泡图
GOBubble(circ, labels = 3)

GOCircle(circ, nsub = 3)


chord<-chord_dat(circ, genedata2)

GOChord(chord, gene.order = 'logFC')

##基因与GO Term的热图(GOHeat)
GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))




#################################################################









