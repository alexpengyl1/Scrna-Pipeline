##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

#install.packages("arrow")
#BiocManager::install(c("AUCell", "RcisTarget"))
#BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost
#BiocManager::install(c("zoo", "mixtools", "rbokeh"))
#BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
#BiocManager::install(c("doMC", "doRNG"))
#devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
#devtools::install_github("aertslab/SCENIC") 
library(ComplexHeatmap)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(harmony)

setwd("D:/shangke/lession16/")

Idents(scRNA_harmony)<- "celltype"

######这里只提取3个细胞群
sc1 = subset(scRNA_harmony,idents=c("Macrophage","Monocyte","Fibroblasts"))
sc2 = subset(sc1,downsample= 150)
View(sc2)

scRNAsub=sc2
#####用于v5
exprMat1 =GetAssayData(scRNAsub,slot = "count")
exprMat1 =as.matrix(exprMat1)
#####用于v4
exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
##设置分析环境

##https://resources.aertslab.org/cistarget/
data(list ="motifAnnotations_hgnc_v9",package='RcisTarget')
motifAnnotations_hgnc=motifAnnotations_hgnc_v9
mydbDIR <- "/Users/yuepen/Downloads"
mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) =c("10kb","500bp")
#初始化 SCENIC 设置,设置分析环境
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=7,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "new")
saveRDS(scenicOptions, "int/scenicOptions.rds")


##==转录调控网络推断==##
##基因过滤
#过滤标准是基因表达量之和>细胞数*3%，且在1%的细胞中表达
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)
##TF-Targets相关性回归分析
exprMat_filtered_log <- log2(exprMat_filtered+1)
#根据表达数据推断潜在的转录因子靶标，使用 GENIE3 或 GRNBoost，
#GENIE3 非常耗时且计算量大（在 3-5k 单元的数据集上需要几个小时或几天的时间）
#GRNboost可在很短的时间内提供与 GENIE3 类似的结果，这儿使用的R，选择GENIC3

##nParts参数，是把表达矩阵分成n份分开计算 
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 50)#####h这一步特别花时间
head("int/1.4_GENIE3_linkList.Rds")
head(readRDS("int/1.4_GENIE3_linkList.Rds") )

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions) #1. 获取共表达模块
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)  #2. 获取regulons
?runSCENIC_2_createRegulons

##==regulon活性评分与可视化==##
##regulons计算AUC值并进行下游分析
#load("scRNA_harmony.Rdata")
#mydbDIR <- "D:/ref.data/"mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
#names(mydbs) <- c("10kb","500bp")
#初始化 SCENIC 设置,设置分析环境
#scenicOptions <- initializeScenic(org="hgnc", 
 #                                 nCores=1,
 #                                 dbDir=mydbDIR, 
 #                                 dbs = mydbs,
 #                                 datasetTitle = "os")
####要是有的时候多核心出错 可以用单核心做不影响结果
########这里是评分的时候我们对所有的细胞进行评分，但是在计算的时候我们会对细胞进行抽样 
library(foreach)
###因为我只提取了3个细胞群为sc1，
#从sc1从我们downsample500来做抽样 减少计算量，现在我比对得分需要给所有的细胞
exprMat_all <- as.matrix(sc1@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
#rm(scRNA_harmony)
###因为我只提取了3个细胞群为sc1

runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all) 






###################################################
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

##准备细胞meta信息，这里就可以使用全部细胞的信息，数量是全部，但是只是基于筛选之前的相同细胞种类
###sc1
cellInfo <- data.frame(sc1@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype_Monaco")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster","celltype")]
saveRDS(cellInfo, file="int/cellInfo.Rds")


###############################################
mydbDIR <- "D:/ref.data/"
mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) <- c("10kb","500bp")
#初始化 SCENIC 设置,设置分析环境
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "os")


library(foreach)
nPcs <- c(5)

fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))

# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="celltype",   cex=.5)




tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
###################
###########################################跑完scenic4之后就到这里
###################################使用loadInt这个函数来读取scenic4生成的文件
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

####现在我们用有了sc1（这是没有sownsample的结果），但是我们想看一下现在的regulon的细胞数和sc1是不是一样的
sctext=sc1
Idents(sctext)="celltype"

x=intersect(colnames(aucell_regulonAUC),row.names(sctext@meta.data))
####下面这一步是可以对aucell_regulonAUC进行subset不过这里我们选取了所有的细胞 ，但是其实是可以subset 的
aucell_regulonAUC_subset= aucell_regulonAUC[,row.names(sctext@meta.data)]

cellInfo.subset=sctext@meta.data

rss <- calcRSS(AUC=getAUC(aucell_regulonAUC_subset), cellAnnotation=cellInfo.subset[colnames(aucell_regulonAUC_subset), "celltype"])
rss <- calcRSS(AUC=getAUC(aucell_regulonAUC_subset), cellAnnotation=cellInfo.subset[colnames(aucell_regulonAUC_subset), "orig.ident"])
#############以上代码的"celltype"还可以替换为分组的信息，也可以吧分组信息和细胞信息 paste0 形成一个新的分组 就可以查看两个分组情况下同一个细胞的转录因子差别
rss = na.omit(rss)
colnames(rss)

rssPlot <- plotRSS(rss)
#####################以上是dotplot，下面我们要做ranking的图，根据分组
# Convert RSS to data frame
B_rss <- as.data.frame(rss)

# Get column names for plot preparation
colnames(B_rss)
celltype <- colnames(B_rss)
rssRanklist <- list()

# Load required library
library(ggrepel)

# Create ranking plot for each cell type
for(i in 1:length(celltype)) {
  # Prepare data for plotting
  data_rank_plot <- cbind(as.data.frame(rownames(B_rss)),
                          as.data.frame(B_rss[,celltype[i]]))
  
  # Set column names
  colnames(data_rank_plot) <- c("TF", "celltype")
  
  # Remove NA values
  data_rank_plot=na.omit(data_rank_plot)
  
  # Order by celltype value in decreasing order
  data_rank_plot <- data_rank_plot[order(data_rank_plot$celltype, decreasing=T),]
  
  # Add rank column
  data_rank_plot$rank <- seq(1, nrow(data_rank_plot))
  
  # Create plot
  p <- ggplot(data_rank_plot, aes(x=rank, y=celltype)) +
    # Add main points
    geom_point(size=3, shape=16, color="#1F77B4", alpha=0.4) +
    # Add highlighted points for top 6
    geom_point(data = data_rank_plot[1:6,],
               size=3, color="#DC050C") +
    # Set theme
    theme_bw() +
    theme(axis.title = element_text(colour = "black", size = 12),
          axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    # Add labels
    labs(x='Regulons Rank', 
         y='Specificity Score',
         title=celltype[i]) +
    # Add text labels for top 6 regulons
    geom_text_repel(data= data_rank_plot[1:6,],
                    aes(label=TF), 
                    color="black",
                    size=3,
                    fontface="italic",
                    box.padding = 0.2,
                    point.padding = 0.3,
                    segment.color = "black",
                    segment.size = 0.3,
                    force = 1,
                    max.iter = 3e3,
                    arrow = arrow(ends="first", length = unit(0.01, "npc")))
  
  # Store the plot in the list
  rssRanklist[[i]] <- p
}

library(cowplot)
plot_grid(rssRanklist[[1]], rssRanklist[[2]], rssRanklist[[3]], ncol=3) ##这是针对3个分组的
plot_grid(rssRanklist[[1]], rssRanklist[[2]],ncol=2)######这是针对2个分组
tf.gene = c(rssRanklist[[1]][["data"]]$TF[1:40], rssRanklist[[2]][["data"]]$TF[1:40])

###导出的时候也要进行修改
write.csv(c(rssRanklist[[1]][["data"]]$TF[1:40], rssRanklist[[2]][["data"]]$TF[1:40], rssRanklist[[3]][["data"]]$TF[1:40]), file = "macro.monocyte。fibroblasts.csv")




#####################################################以上就是不同的组
library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_all, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("JUN","MYC")],], plots="Expression")

#par(bg = "black")
par(mfrow=c(1,2))
regulonNames <- c( "JUN","MYC")
cellCol <- plotEmb_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)


regulonNames <- list( green=c("JUN"),
                      blue=c( "MYC"))
cellCol <- plotEmb_rgb(scenicOptions, regulonNames, aucType="Binary")

regulons <- loadInt(scenicOptions, "regulons")
regulons[c( "JUN","MYC")]

regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="YY1" & highConfAnnot==TRUE]
viewMotifs(tableSubset, options=list(pageLength=5)) 

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="YY1"]
viewMotifs(tableSubset)



regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:20,], name="Regulon activity")


topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$celltype), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
ComplexHeatmap::Heatmap(binaryActPerc_subset[1:20,], name="Regulon activity (%)", col = c("white","pink","red"))


topRegulators <- reshape2::melt(regulonActivity_byCellType_Binarized)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
viewTable(topRegulators)
################################################


rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "celltype"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)


plotRSS_oneSet(rss, setName = "MSC")


library(Seurat)
#scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:16)
dr_coords <- Embeddings(scRNA_harmony, reduction="tsne")
par(mfrow=c(4,4))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, "MYC"), plots = "AUC")


tfs <- c("HDAC2","RAD21","YY1", "SMARCA4")
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")
########################################################################


##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(scRNA_harmony, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNA_harmony, BINmatrix)
scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'scRNAbin.rds')

##利用Seurat可视化AUC
dir.create('scenic_seurat')
#FeaturePlot
colnames(scRNAauc@meta.data)[20:30]
p1 = FeaturePlot(scRNAauc, features="HCFC1_24g", label=T, reduction = 'tsne')
p2 = FeaturePlot(scRNAbin, features="HCFC1_24g", label=T, reduction = 'tsne')
p3 = DimPlot(scRNA_harmony, reduction = 'tsne', group.by = "celltype", label=T)
plotc = p1|p2|p3
plotc



#RidgePlot&VlnPlot
p1 = RidgePlot(scRNAauc, features = "NR3C1_1339g", group.by="celltype") + 
  theme(legend.position='none')
p2 = VlnPlot(scRNAauc, features ="NR3C1_1339g", pt.size = 0, group.by="celltype") + 
  theme(legend.position='none')
plotc = p1 + p2
plotc


library(pheatmap)
cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(cellInfo,select = 'celltype')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
#挑选部分感兴趣的regulons
my.regulons <- rownames(AUCmatrix)[80:100]
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]
#使用regulon原始AUC值绘制热图
pheatmap(myAUCmatrix, show_colnames=F, annotation_col=celltype )
#使用regulon二进制AUC值绘制热图
pheatmap(myBINmatrix, show_colnames=F, annotation_col=celltype,
         color = colorRampPalette(colors = c("white","black"))(100))


