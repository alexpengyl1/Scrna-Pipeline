# 设置环境变量
Sys.setenv(BOOST_ROOT = "C:/boost_1_87_0")
Sys.setenv(BOOST_LIBRARYDIR = "C:/boost_1_87_0/stage/lib")

# 验证设置
Sys.getenv("BOOST_ROOT")
Sys.getenv("BOOST_LIBRARYDIR")

# 重新尝试安装
devtools::install_github("velocyto-team/velocyto.R")



# 加载必要的包
library(Seurat)
library(velocyto.R)
library(tidyverse)
library(SeuratWrappers)
library(magrittr)
# 读取gene.filtered矩阵
# 10X数据
gene_matrix <- ReadMtx(
  mtx = "gene.filtered/matrix.mtx",
  cells = "gene.filtered/barcodes.tsv",
  features = "velocyto.filtered/features.tsv",  # 使用相同的feature文件
  feature.column = 2  # 使用相同的基因符号列
)


spliced <- ReadMtx(
  mtx = "velocyto.filtered/spliced.mtx",
  cells = "velocyto.filtered/barcodes.tsv",
  features = "velocyto.filtered/features.tsv",
  feature.column = 2  # 使用第2列作为特征名
)

unspliced <- ReadMtx(
  mtx = "velocyto.filtered/unspliced.mtx",
  cells = "velocyto.filtered/barcodes.tsv",
  features = "velocyto.filtered/features.tsv",
  feature.column = 2  # 使用第2列作为特征名
)

# 创建Seurat对象、
velo <- CreateSeuratObject(counts = gene_matrix)
velo[["spliced"]] <- CreateAssayObject(counts = spliced)
velo[["unspliced"]] <- CreateAssayObject(counts = unspliced)


# 在SCTransform之前添加质控
velo <- PercentageFeatureSet(velo, pattern = "^MT-", col.name = "percent.mt")
velo <- subset(velo, subset = percent.mt < 20)  # 过滤线粒体比例高的细胞
# 数据预处理和降维
velo <- velo %>% 
  SCTransform(assay = "RNA") %>% 
  RunPCA(verbose = FALSE)

# 检查PC数量
ElbowPlot(velo, ndims = 50)

# 聚类和降维可视化
nPC <- 1:20
velo <- velo %>% 
  FindNeighbors(dims = nPC) %>% 
  FindClusters() %>% 
  RunUMAP(dims = nPC) %>% 
  RunTSNE(dims = nPC)

# 设置细胞颜色
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = velo)))
names(x = ident.colors) <- levels(x = velo)
cell.colors <- ident.colors[Idents(object = velo)]
names(x = cell.colors) <- colnames(x = velo)
library(SeuratWrappers)

# 检查当前assay名称
DefaultAssay(velo)  # 看看当前的默认assay是什么

# 确认一下所有的assay
Assays(velo)  # 这会显示所有可用的assay

# 设置指定的assay名称
DefaultAssay(velo) <- "RNA"  # 设置RNA作为默认assay


# RNA velocity分析
# 或者直接在运行velocity时指定assay
velo <- RunVelocity(velo, 
                    deltaT = 1, 
                    kCells = 25, 
                    fit.quantile = 0.02, 
                    spliced.average = 0.2, 
                    unspliced.average = 0.05, 
                    ncores = 1,
                    spliced.assay = "spliced",     # 明确指定使用RNA assay作为spliced数据
                    unspliced.assay = "unspliced")

# 可视化
emb <- Embeddings(velo, reduction = "umap")
vel <- Tool(velo, slot = "RunVelocity")

# 全局velocity图
show.velocity.on.embedding.cor(emb = emb, 
                               vel = vel, 
                               n = 200, 
                               scale = "sqrt", 
                               cell.colors = ac(cell.colors, alpha = 0.5), 
                               cex = 0.8, 
                               arrow.scale = 3, 
                               show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, 
                               grid.n = 40, 
                               arrow.lwd = 1, 
                               do.par = FALSE, 
                               cell.border.alpha = 0.1)

p2=DimPlot(velo, reduction = "umap", label = TRUE)
p2
# 使用show.velocity.on.embedding来绘制
show.velocity.on.embedding.cor(
  emb = Embeddings(velo, reduction = "umap"),
  vel = velocyto, 
  n = 200, 
  scale = 'sqrt',
  arrow.scale = 3, 
  show.grid.flow = TRUE
)

# 如果要看特定基因的velocity
gene <- "CAMTA1"  # 替换成您感兴趣的基因名
RunVelocity(velo, 
            deltaT = 1, 
            kCells = 25, 
            fit.quantile = 0.02, 
            old.fit = vel, 
            cell.emb = emb, 
            cell.colors = cell.colors, 
            show.gene = gene, 
            do.par = TRUE)



# 获取velocity分析中使用的基因
vel.genes <- rownames(vel$deltaE)  # 这些是最终用于velocity分析的基因

# 查看基因数量
length(vel.genes)

# 查看前20个基因
head(vel.genes, 20)

# 可以检查这些基因的velocity贡献
deltaE <- vel$deltaE  # 速率变化矩阵
mean.deltaE <- rowMeans(deltaE)  # 计算每个基因的平均速率变化

# 创建一个包含基因名称和其velocity贡献的数据框
gene.velocity <- data.frame(
  gene = vel.genes,
  mean_velocity = mean.deltaE
)

# 按velocity绝对值大小排序
gene.velocity <- gene.velocity[order(abs(gene.velocity$mean_velocity), decreasing = TRUE), ]

# 查看velocity贡献最大的前20个基因
head(gene.velocity, 20)
