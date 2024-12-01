


##  https://linkinghub.elsevier.com/retrieve/pii/S1535610822005025
## https://github.com/lifei176/Single-cell-multi-omics/blob/main/Cell_cell_communication_CellChat.R

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


setwd("/Users/yuepen/Downloads/华哥2024单细胞高阶班/资料/7")
load("scRNA_harmony.Rdata")

library(SingleR)
load("/Users/yuepen/Downloads/华哥2024单细胞高阶班/资料/7/ref_Human_all.RData")

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
setwd("F:/shangke/lession18/")


install.packages("ggalluvial")

install.packages('NMF')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")

devtools::install_github("sqjin/CellChat")



library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(data.table)
library(ggsci)

#创建cellchat对象
#scRNA_harmony.1= scRNA_harmony[,sample( colnames(scRNA_harmony) ,2000 )] 
##提取表达矩阵和细胞分类信息
data.input <- GetAssayData(scRNA_harmony, slot = "data")
identity <- subset(scRNA_harmony@meta.data, select = "celltype")
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "celltype")


####可选CellChatDB.human, CellChatDB.mouse
CellChatDB <- CellChatDB.human
##下一步不出图的时候运行 dev.new()
showDatabaseCategory(CellChatDB)

##
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)

########在CellChat中，我们还可以先择特定的信息描述细胞间的相互作用，
##可以理解为从特定的侧面来刻画细胞间相互作用，比用一个大的配体库又精细了许多。
##查看可以选择的侧面
unique(CellChatDB$interaction$annotation)
# use Secreted Signaling for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use # set the used database in the object
##################no subset
CellChatDB.use <- subsetDB(CellChatDB)
cellchat@DB <- CellChatDB.use # set the used database in the object

#对表达数据进行预处理

##将信号基因的表达数据进行子集化，以节省计算成本
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)

##相互作用推断
## 1、计算通信概率推断细胞互作的通信网络

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
###如果特定细胞群中只有少数细胞，则过滤掉细胞间的通信
cellchat <- filterCommunication(cellchat, min.cells = 3)

#提取推断出的细胞互作的通信网络数据框，我们提供了一个subsetCommunication 函数，
#可以方便地访问感兴趣的推断的细胞间通信。

##返回一个数据框，包含所有推断的配体/受体级别的细胞-细胞通信。设置slot.name = "netP"以访问信令路径级别的推断通信
df.net <- subsetCommunication(cellchat)

##
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))

##2、在信号通路水平上推断细胞间的通讯
cellchat <- computeCommunProbPathway(cellchat)
##汇总通信概率来计算细胞间的聚合通信网络。
cellchat <- aggregateNet(cellchat)
##3、计算聚合细胞互作通信网络
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")
#左图：外周各种颜色圆圈的大小表示细胞的数量，圈越大，细胞数越多。发出箭头的细胞表达配体，
#箭头指向的细胞表达受体。配体-受体对越多，线越粗。
#右图：互作的概率或者强度值（强度就是概率值相加）


##每个细胞如何跟别的细胞互作（互作的强度或概率图）
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
##每个细胞如何跟别的细胞互作（number+of+interaction图）
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

##可视化每个信号通路
##查看通路

levels(cellchat@idents)            #查看细胞顺序
vertex.receiver = c(1, 2)          #指定靶细胞的索引
cellchat@netP$pathways             #查看富集到的信号通路
pathways.show <- "CCL"            #指定需要展示的通路

?netVisual_aggregate
##层次图
vertex.receiver = seq(1,2) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = "MK",  vertex.receiver = c(1,2),layout="hierarchy")
#在层次图中，实体圆和空心圆分别表示源和目标。圆的大小与每个细胞组的细胞数成比例。线越粗，互作信号越强。
#左图中间的target是我们选定的靶细胞。右图是选中的靶细胞之外的另外一组放在中间看互作。

##圈图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="MK", layout = "circle")

##和弦图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="MK", layout = "chord", vertex.size = groupSize)

##热图
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "MK", color.heatmap = "Reds")
##纵轴是发出信号的细胞，横轴是接收信号的细胞，热图颜色深浅代表信号强度。
##上侧和右侧的柱子是纵轴和横轴强度的累积


#配体-受体层级的可视化（计算各个ligand-receptor+pair对信号通路的贡献）

netAnalysis_contribution(cellchat, signaling = pathways.show)
##也可以看到单个配体-受体对介导的细胞-细胞通信。
#我们提供了一个extractEnrichedLR功能来提取给定信号通路的所有重要相互作用(L-R对)和相关信号基因。
pairLR.MK <- extractEnrichedLR(cellchat, signaling = "MK", geneLR.return = FALSE)
#提取对这个通路贡献最大的配体受体对来展示（也可以选择其他的配体受体对）
LR.show <- pairLR.MK[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,2) # a numeric vector
##层次图
netVisual_individual(cellchat, signaling = "MK",  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout="hierarchy")
##圈图
netVisual_individual(cellchat, signaling = "MK", pairLR.use = LR.show, layout = "circle")
##和弦图
netVisual_individual(cellchat, signaling ="MK", pairLR.use = LR.show, layout = "chord")


##############批量保存
pathway.show.all=cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver=c(1,2,3,4)
dir.create("F:/shangke/lession18/cellchat.res")
setwd("F:/shangke/lession18/cellchat.res/")
for (i in 1:length(pathway.show.all)) {
  
  netVisual(cellchat,signaling = pathway.show.all[i],out.format = c("pdf"),
            vertex.receiver=vertex.receiver,layout="circle")
  plot=netAnalysis_contribution(cellchat,signaling = pathway.show.all[i])
  ggsave(filename = paste0(pathway.show.all[i],".contribution.pdf"),
         plot=plot,width=6,height=4,dpi=300,units="in")
  
}

#####################################################################################
##气泡图
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:5), remove.isolate = FALSE)
##sources.use = 2 是值第二个细胞亚群
netVisual_bubble(cellchat, sources.use =c(1,3), targets.use = c(1:5), remove.isolate = FALSE)
##指定信号通路
cellchat@netP$pathways 
netVisual_bubble(cellchat, sources.use =c(1,3), targets.use =c(1:5),signaling =  c("MK","CCL"), remove.isolate = FALSE)


pairLR  <- extractEnrichedLR(cellchat, signaling =c("MK","CCL"), geneLR.return = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,3), targets.use =c(1:5),pairLR.use =pairLR , remove.isolate = FALSE)

#和弦图               
netVisual_chord_gene(cellchat, sources.use = 2, targets.use = c(1:5), lab.cex = 0.5,legend.pos.y = 30)

##用小提琴图绘制信号基因的表达分布 参与某条信号通路（如TGFb）的所有基因在细胞群中的表达情况展示
plotGeneExpression(cellchat, signaling = "MK")
#默认情况下，plotGeneExpression只显示与推断的重要通信相关的信号基因的表达。
#用户可以通过显示一个信号通路相关的所有信号基因的表达。
plotGeneExpression(cellchat, signaling = "MK", enriched.only = FALSE)
#也可以用气泡图展示
plotGeneExpression(cellchat, signaling = "MK",type = "dot")

###################################################################################################
##可视化配体和受体
## 1、计算网络中心性得分
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
##2、热图  使用热图可视化计算的中心性评分，允许随时识别细胞群的主要信号作用。
netAnalysis_signalingRole_network(cellchat, signaling = "MK", width = 8, height = 2.5, font.size = 10)

##在2D空间中可视化主要的发送者(源)和接收者(目标)。
##我们还提供了另一种直观的方式，使用散点图来可视化2D空间中的主要发送者(源)和接收者(目标)。

##从所有信号通路对聚合的细胞-细胞通信网络的信号作用分析
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
###从所有信号通路对聚合的细胞-细胞通信网络的信号作用分析
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MK", "PARs"))
gg1 + gg2

##识别对某些细胞群的传出或传入信号贡献最大的信号，从所有信号通路对聚合的细胞-细胞通信网络的信号作用分析。

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

########################################################################################################
#细胞通讯模式和信号网络

library(NMF)
library(ggalluvial)
#非负矩阵分解（NMF）识别细胞的通讯模式
##信号输出细胞的模式识别
##计算分解成几个因子(pattern)比较合适（这一步运行比较慢+。在使用NMF对细胞进行亚群细分时，如果不测试的话，最好选择比细胞类型多一点的值）
selectK(cellchat, pattern = "outgoing")


#挑选曲线中第一个出现下降的点（从3就开始下降了）
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
##river plot
library(ggalluvial)
netAnalysis_river(cellchat, pattern = "outgoing")
#气泡图
netAnalysis_dot(cellchat, pattern = "outgoing")
#信号输入细胞的模式识别
selectK(cellchat, pattern = "incoming")

#################################################################################################
#  信号网络聚类
# 1、根据功能相似性来识别信号分组

##reticulate::py_install(packages = 'umap-learn')

## 2、基于结构相似性识别信号分组
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

#############################################################################################
#不同分组之间的配对分析
sc.sp=SplitObject(scRNA_harmony,split.by = "orig.ident")
sc.11=scRNA_harmony[,sample(colnames(sc.sp[["sample_11"]]),1000)]
sc.3=scRNA_harmony[,sample(colnames(sc.sp[["sample_3"]]),1000)]



cellchat.sc11 <- createCellChat(object =sc.11@assays$RNA@data, meta =sc.11@meta.data,  group.by ="celltype")
cellchat.sc3 <- createCellChat(object =sc.3@assays$RNA@data, meta =sc.3@meta.data,  group.by ="celltype")

dir.create("compare")
setwd("compare/")
#########样本1
cellchat=cellchat.sc11 
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc11 = cellchat
#################################样本2
cellchat=cellchat.sc3
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc3 = cellchat
##############################################merge在一起
cc.list=list(SC11=cc.sc11,SC3=cc.sc3)
cellchat=mergeCellChat(cc.list,cell.prefix = T,add.names = names(cc.list))
##可视化
##所有细胞群总体观：通讯数量与强度对比
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "count")
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "weight")
##第一个图展示通讯数量之间的差异，第二个图展示通讯强度之间的差异。 

##数量与强度差异网络图
netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight")
##红色是case相对于control上调的，蓝色是下调的

#数量与强度差异热图
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat,measure = "weight")
#case和control对比，红色是上调，蓝色是下调

#保守和特异性信号通路的识别与可视化
rankNet(cellchat,mode = "comparison",stacked = T,do.stat = T)
rankNet(cellchat,mode = "comparison",stacked =F,do.stat = T)
##左图最下面多个信号通路是case组独有的

##细胞互作数量对比网络图
weight.max=getMaxWeight(cc.list,attribute = c("idents","count"))
netVisual_circle(cc.list[[1]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )

netVisual_circle(cc.list[[2]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )


table(scRNA_harmony@active.ident)
s.cell=c( "Macrophage", "Tissue_stem_cells","Monocyte")
count1=cc.list[[1]]@net$count[s.cell,s.cell]
count2=cc.list[[2]]@net$count[s.cell,s.cell]

netVisual_circle(count1,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )

netVisual_circle(count2,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )



netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)


gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2



###################################

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "LS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "LS",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "NL",ligand.logFC = -0.1, receptor.logFC = -0.1)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2


# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))






