

##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(Seurat)
library(scatterplot3d)

setwd("E:/super.lesson/lesson8/")
load("sc.combined.rdata")


table(sc.combined@meta.data$new.celltype.2)
table(sc.combined@meta.data$type)

#setwd("E:/super.lesson/lesson9/")

##  Fibroblasts

sc.imm=subset(sc.combined,idents = c("Epithelial","Endothelial","Stem_cell"))


Idents(sc.imm)="orig.ident"
sc.imm@active.ident
sc.imm@meta.data$orig.ident=factor(sc.imm@meta.data$orig.ident,levels = c(
  "GF_1","GF_2","GF_3","SPF_1","SPF_2","SPF_3","FMT_1","FMT_2","FMT_3"
))
Idents(sc.imm)="orig.ident"
table(sc.imm@active.ident)

?AverageExpression
x=AverageExpression(sc.imm ,add.ident="new.celltype.2")
x1=AverageExpression(sc.imm  )
x=x$RNA
x1=x1$RNA
pca <- prcomp(t(x))

pca.data=pca$x

rownames(pca.data)

pca.data=as.data.frame(pca.data)
pca.data=pca.data[,1:3]

pca.data$Type = c(rep("GF",9),rep("SPF",9),rep("FMT",9))

s1=strsplit(rownames(pca.data),split = "_",fixed = T)
type=sapply(s1, function(x){x[1]}   )


pca.data$cell.type = c(rep(c("Epithelial","Endothelial","Stem_cell"),9 ))



Type.p=c(rep(11,9),rep(16,9),rep(17,9))
cell.type.p = c(rep(c("blue","red","orange"),9 ))


colors.lib <- c("blue","red","orange")
shapes.lib = c(11,16,17)


################插入部分为claude提供的3d作图函数
library(plotly)

plot_3d_pca <- function(pca_data) {
  cell_colors <- c("Epithelial" = "blue", 
                   "Endothelial" = "red", 
                   "Stem_cell" = "orange")
  
  type_shapes <- c("GF" = "square-open",      # 空心正方形
                   "SPF" = "triangle-up-open", # 空心三角形
                   "FMT" = "circle-open")      # 空心圆形
  
  p <- plot_ly() %>%
    add_trace(
      data = pca_data,
      x = ~PC1, y = ~PC2, z = ~PC3,
      color = ~cell.type,
      symbol = ~Type,
      colors = cell_colors,
      symbols = type_shapes,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 8)
    ) %>%
    layout(
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3")
      ),
      title = "3D PCA Plot"
    )
  
  return(p)
}

p <- plot_3d_pca(pca.data)
p


###############################################以下部分为华哥的
getwd()
setwd("/Users/yuepen/Downloads/华哥2024单细胞高阶班/资料/11/")
setwd("~/Downloads/华哥2024单细胞高阶班/资料/11/")
# 1. Source the function
source( "/Users/yuepen/Downloads/华哥2024单细胞高阶班/资料/11/huage3Dplot.R" )
source("/Users/yuepen/Downloads/华哥2024单细胞高阶班/资料/11/huage3Dplot.R")
# 2. 3D scatter plot


s3d <- scatterplot3d(pca.data[,c("PC1","PC2","PC3")],
                     
                     pch = Type.p,color = cell.type.p,
                     
                     cex.symbols = 1,grid=FALSE, box=FALSE, 
                     
                     main = "3D PCA plot")


legend("topright",
       c("Epithelial","Endothelial","Stem_cell"),
       fill=c('blue',"red","orange"),
       box.col=NA)

legend("topleft",
       c('GF','SPF','FMT'),
       pch = shapes.lib,
       box.col=NA)

# 3. Add grids
addgrids3d(pca.data[,c("PC1","PC2","PC3")], grid = c("xy", "xz", "yz"))

table(sc.combined@meta.data$new.celltype.2)

##  Goblet Paneth cells Enteroendocrine
Idents(sc.combined)="new.celltype.3"
sc.imm=subset(sc.combined,idents = c("Goblet","MSC","T_NK"))


Idents(sc.imm)="orig.ident"
x=AverageExpression(sc.imm ,add.ident="new.celltype.3")
x=x$RNA
pca <- prcomp(t(x))

pca.data=pca$x

rownames(pca.data)

pca.data=as.data.frame(pca.data)
pca.data=pca.data[,1:3]

pca.data$Type = c(rep("GF",9),rep("SPF",9),rep("FMT",9))
pca.data$cell.type = c(rep(c("Goblet","Paneth cells","Enteroendocrine"),9 ))



Type.p=c(rep(11,9),rep(16,9),rep(17,9))
cell.type.p = c(rep(c("blue","red","orange"),9 ))


colors.lib <- c("blue","red","orange")
shapes.lib = c(11,16,17)


getwd()


# 1. Source the function
source('plot.3d.R')
# 2. 3D scatter plot


s3d <- scatterplot3d(pca.data[,c("PC1","PC2","PC3")],
                     
                     pch = Type.p,color = cell.type.p,
                     
                     cex.symbols = 1,grid=FALSE, box=FALSE, 
                     
                     main = "3D PCA plot")




legend("topright",
       c("Goblet","Paneth cells","Enteroendocrine"),
       fill=c('blue',"red","orange"),
       box.col=NA)

legend("topleft",
       c('GF','SPF','FMT'),
       pch = shapes.lib,
       box.col=NA)

# 3. Add grids
addgrids3d(pca.data[,c("PC1","PC2","PC3")], grid = c("xy", "xz", "yz"))

