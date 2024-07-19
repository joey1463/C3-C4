# This script will subcluster bundle sheath 10X nuclei from sorghum, using the genes found differentially expressed as variable features. 

library(Seurat)

load("L3_sorghum_bundle_sheath_10X.RData")
load("L3_sorghum_bundle_sheath_DE_Significant.RData")

# Redo Clusters with Variable Features
VariableFeatures(subcluster_sorghum_combined_10X)<-rownames(DE_Significant)
subcluster_sorghum_combined_10X <- SCTransform(subcluster_sorghum_combined_10X, vars.to.regress = "nCount_RNA", verbose = TRUE)
subcluster_sorghum_combined_10X <- RunPCA(object = subcluster_sorghum_combined_10X, npcs = 30, verbose = FALSE)
subcluster_sorghum_combined_10X <- RunUMAP(object = subcluster_sorghum_combined_10X, reduction = "pca", dims = 1:15) 
subcluster_sorghum_combined_10X <- FindNeighbors(object = subcluster_sorghum_combined_10X, reduction = "pca", dims = 1:15)
subcluster_sorghum_combined_10X <- FindClusters(subcluster_sorghum_combined_10X, resolution = 0) 
Idents(subcluster_sorghum_combined_10X)<-subcluster_sorghum_combined_10X$light_time 
DimPlot(object = subcluster_sorghum_combined_10X, reduction = "umap", label = FALSE, raster = FALSE, pt.size = 0.5, cols =plasma(8, direction = 1)) + NoLegend()
FeaturePlot(object = subcluster_sorghum_combined_10X, features = "Sobic.003G036200", cols = c("grey80","hotpink4"),  pt.size = 0.5, raster=FALSE)

save(subcluster_sorghum_combined_10X, file = "L3_sorghum_bundle_sheath_reclustered.RData")
#load("L3_sorghum_bundle_sheath_reclustered.RData")
