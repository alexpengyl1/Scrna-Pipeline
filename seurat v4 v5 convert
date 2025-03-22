# https://github.com/satijalab/seurat/issues/7409
# convert a v5 assay to a v3 ,v4 assay
pbmc3k[["RNA3"]] <- as(object = pbmc3k[["RNA"]], Class = "Assay")
 
​pbmc3k[["RNA4"]] <- as(object = pbmc3k[["RNA"]], Class = "Assay4")
 
# convert a v3 assay to a v5 assay
pbmc3k[["RNA5"]] <- as(object = pbmc3k[["RNA3"]], Class = "Assay5")
