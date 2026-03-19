ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
mono.de <- FindMarkers(ifnb, ident.1 = "pDC_STIM", ident.2 = "pDC_CTRL", verbose = FALSE)
head(mono.de, n = 10)

ifnb$donor_id <- "SNG-1488"
ifnb$donor_id <- NULL 

ifnb$donor_id <- sample(1:6, 
                    ncol(ifnb), 
                    replace = TRUE)
# pseudobulk the counts based on donor-condition-celltype
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_ifnb))
pseudo_ifnb$celltype.stim <- paste(pseudo_ifnb$seurat_annotations, pseudo_ifnb$stim, sep = "_")


Idents(pseudo_ifnb) <- "celltype.stim"

bulk.mono.de <- FindMarkers(object = pseudo_ifnb, 
                            ident.1 = "pDC_STIM", 
                            ident.2 = "pDC_CTRL",
                            test.use = "DESeq2")
head(bulk.mono.de, n = 15)

# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.mono.de) <- paste0(names(bulk.mono.de), ".bulk")
bulk.mono.de$gene <- rownames(bulk.mono.de)

names(mono.de) <- paste0(names(mono.de), ".sc")
mono.de$gene <- rownames(mono.de)

merge_dat <- merge(mono.de, bulk.mono.de, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                 merge_dat$p_val.sc < 0.05)]
only_sc <- merge_dat$gene[which(merge_dat$p_val.bulk > 0.05 & 
                                  merge_dat$p_val.sc < 0.05)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                    merge_dat$p_val.sc > 0.05)]
print(paste0('# Common: ',length(common)))
print(paste0('# Only in single-cell: ',length(only_sc)))
print(paste0('# Only in bulk: ',length(only_bulk)))


# create a new column to annotate sample-condition-celltype in the single-cell dataset
ifnb$donor_id.stim <- paste0(ifnb$stim, "-", ifnb$donor_id)

# generate violin plot 
Idents(ifnb) <- "celltype.stim"
print(merge_dat[merge_dat$gene%in%common[1:2],c('gene','p_val.sc.sc','p_val.bulk')])

VlnPlot(ifnb, features = common[1:2], idents = c("pDC_STIM", "pDC_CTRL"), group.by = "stim") 


VlnPlot(ifnb, features = common[1:2], idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "donor_id.stim", ncol = 1) 


###############By contrast, we can examine examples of genes that are only DE under the single-cell analysis.
print(merge_dat[merge_dat$gene%in%c('SRGN','HLA-DRA'),c('gene','p_val.sc','p_val.bulk')])
VlnPlot(ifnb, features <- c('SRGN','HLA-DRA'), idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "stim") 




##############再进行gsea
# 根据log2FC排序
gene_list <- bulk.mono.de$avg_log2FC.bulk
names(gene_list) <- bulk.mono.de$gene
gene_list <- sort(gene_list, decreasing = TRUE)


library(clusterProfiler)
library(msigdbr)
library(enrichplot)

# 获取通路数据库
msig_h <- msigdbr(species = "Homo sapiens", category = "H") # hallmark通路
msig_h$gs_name <- gsub("HALLMARK_", "", msig_h$gs_name)###去除开头的hallmark_
# 或用其他数据库如KEGG
# msig_k <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")

# GSEA分析
gsea_result <- GSEA(
  gene_list,
  TERM2GENE = dplyr::select(msig_h, gs_name, gene_symbol),
  pvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500
)

# 结果可视化
enrichplot::dotplot(gsea_result)



# 2. 查看所有通路
head(gsea_result)

# 3. 绘制单个通路的富集图
# gseaplot2可以同时展示富集分数曲线、hits位置和排序基因表达量
gseaplot2(gsea_result, 
          geneSetID = c(1,2),    # 可以是通路名称或序号
          title = "Pathway Name",
          pvalue_table = TRUE,  # 显示p值等统计信息
          ES_geom = "line")    # 富集分数曲线类型

# 4. 保存图片
ggsave("gsea_pathway.pdf", width = 8, height = 6)
