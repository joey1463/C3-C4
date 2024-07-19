# This script will subcluster mesophyll 10X nuclei from rice, using the genes found differentially expressed as variable features. 

library(Seurat)

#setwd("~/Desktop/Salk/Project - Light/R Files/")
load("L3_rice_mesophyll_10X.RData")
load("L3_rice_mesophyll_DE_Significant.RData")

# Redo Clusters with Variable Features
VariableFeatures(subcluster_rice_combined_10X)<-rownames(DE_Significant)
subcluster_rice_combined_10X <- SCTransform(subcluster_rice_combined_10X, vars.to.regress = "nCount_RNA", verbose = TRUE)
subcluster_rice_combined_10X <- RunPCA(object = subcluster_rice_combined_10X, npcs = 30, verbose = FALSE)
subcluster_rice_combined_10X <- RunUMAP(object = subcluster_rice_combined_10X, reduction = "pca", dims = 1:15) 
subcluster_rice_combined_10X <- FindNeighbors(object = subcluster_rice_combined_10X, reduction = "pca", dims = 1:15)
subcluster_rice_combined_10X <- FindClusters(subcluster_rice_combined_10X, resolution = 0) 
Idents(subcluster_rice_combined_10X)<-subcluster_rice_combined_10X$light_time 
DefaultAssay(object = subcluster_rice_combined_10X) <- "RNA" 

#load("L3_rice_mesophyll_reclustered.RData")
DimPlot(object = subcluster_rice_combined_10X, reduction = "umap", label = TRUE, raster = FALSE, pt.size = 0.5, cols =plasma(8, direction = 1)) + NoLegend()

save(subcluster_rice_combined_10X, file = "L3_rice_mesophyll_reclustered.RData")
