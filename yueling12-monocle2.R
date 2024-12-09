
### https://github.com/Swarbricklab-code/BrCa_cell_atlas/tree/main/monocle_analysis_stromal_cells
### https://www.nature.com/articles/s41588-021-00911-1


setwd("D:\\BaiduNetdiskDownload\\")

load("scRNA_harmony.Rdata")
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
load("D:\\BaiduNetdiskDownload\\ref_Human_all.RData")

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


#BiocManager::install("monocle")
library(monocle)


#scRNA.Osteoclastic=subset(scRNA_harmony,ident=c(0,1))
scRNA.Osteoclastic=sc.t
##data <- as(as.matrix(scRNA.Osteoclastic@assays$RNA@counts), 'sparseMatrix') 这是v4版本的 以下是 v5
data =GetAssayData(scRNA.Osteoclastic,slot = "counts")
?new
pd <- new('AnnotatedDataFrame', data = scRNA.Osteoclastic@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())



monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))

HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)


HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')


HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")


plot_cell_trajectory(HSMM, color_by = "State")

plot_cell_trajectory(HSMM, color_by = "Pseudotime")

plot_cell_trajectory(HSMM, color_by = "celltype")
plot_cell_trajectory(HSMM, color_by = "orig.ident")
plot_cell_trajectory(HSMM, color_by = "State") +
  facet_wrap(~State, nrow = 3)


###################################################################
blast_genes <- row.names(subset(fData(HSMM),
                                gene_short_name %in% c("GAPDH", "RORA")))
plot_genes_jitter(HSMM[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)



HSMM_expressed_genes <-  row.names(subset(fData(HSMM),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("YWHAB", "GAPDH", "TNNC1")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")



plot_genes_in_pseudotime(cds_subset, color_by =  "State")
plot_genes_in_pseudotime(cds_subset, color_by =  "celltype")


#################这是查看 目的基因在 不同的state条件下的表达差异
genes <- c("TNNT2", "TNNC1", "CDK1")
p1 <- plot_genes_jitter(HSMM[genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(HSMM[genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(HSMM[genes,], color_by = "State")
plotc <- p1|p2|p3
plotc

to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("MYH3", "MEF2C", "CCNB2", "TNNT1")))
cds_subset <- HSMM[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_in_pseudotime(cds_subset, color_by ="State")

###########################################


marker_genes <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5",
                                                        "ANPEP", "PDGFRA","MYOG",
                                                        "TPM1",  "TPM2",  "MYH2",
                                                        "MYH3",  "NCAM1", "TNNT1",
                                                        "TNNT2", "TNNC1", "CDK1",
                                                        "CDK2",  "CCNB1", "CCNB2",
                                                        "CCND1", "CCNA1", "ID1")))

diff_test_res <- differentialGeneTest(HSMM[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 1))#####可以调整q的大小来控制基因数量
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                        num_clusters = 5,#######可以设置cluster的数量
                        cores = 1,
                        show_rownames = T)

################也可以使用全部的基因来做差异分析
diff_test_res <- differentialGeneTest(HSMM,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
gene_to_cluster = diff_test_res %>%    # 从diff_test_res开始
  arrange(qval) %>%                  # 按qval值排序
  head(50) %>%                       # 取前50行
  pull(gene_short_name)              # 提取gene_short_name列的值

##sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))#####可以调整q的大小来控制基因数量，但是有时候基因太多了只能这样用上面的代码来做
plot_pseudotime_heatmap(HSMM[gene_to_cluster,],
                        num_clusters = 5,#######可以设置cluster的数量
                        cores = 20,
                        show_rownames = T)




##################################################################
####beam分支点 为什么细胞群的转化 过程中的基因是什么
####发育起点是可以通过HSMM <- orderCells(HSMM,root_state = n)来调整的
####发育起点是可以通过HSMM <- orderCells(HSMM,root_state = n)来调整的
####发育起点是可以通过HSMM <- orderCells(HSMM,root_state = n)来调整的
####发育起点是可以通过HSMM <- orderCells(HSMM,root_state = n)来调整的
plot_cell_trajectory(HSMM, color_by = "State")
###BEAM_res <- BEAM(cds[ordergene,], branch_point = 2, cores = 2, progenitor_method = 'duplicate')
BEAM_res <- BEAM(HSMM, branch_point = 1, cores =10, progenitor_method = 'duplicate')
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,
                                                  qval < 1e-10)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)


genes <- row.names(subset(fData(HSMM),
                          gene_short_name %in% c( "MEF2C", "CCNB2", "TNNT1")))

plot_genes_branched_pseudotime(HSMM[genes,],
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1)
HSMM <- orderCells(HSMM,root_state = 3)
