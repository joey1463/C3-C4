# This script will take mesophyll, epidermal (guard, epidermis), and vasculature (bundle sheath, phloem, xylem) clusters and re-cluster them. 

library(Seurat)
load("L1_rice_combined.RData")

# Mesophyll
subcluster_rice_combined<-rice_combined[,WhichCells(rice_combined, idents = c(0,5))] 
DefaultAssay(object = subcluster_rice_combined) <- "integrated"
subcluster_rice_combined <- ScaleData(object = subcluster_rice_combined, verbose = TRUE)
subcluster_rice_combined <- RunPCA(object = subcluster_rice_combined, npcs = 30, verbose = FALSE)
subcluster_rice_combined <- RunUMAP(object = subcluster_rice_combined, reduction = "pca", dims = 1:30)
subcluster_rice_combined <- FindNeighbors(object = subcluster_rice_combined, reduction = "pca", dims = 1:30)
subcluster_rice_combined <- FindClusters(subcluster_rice_combined, resolution = 0.4)
save(subcluster_rice_combined, file = "L2_rice_mesophyll_combined.RData")

# Vasculature
subcluster_rice_combined<-rice_combined[,WhichCells(rice_combined, idents = c(9,6,10,11,12,16))]
DefaultAssay(object = subcluster_rice_combined) <- "integrated"
subcluster_rice_combined <- ScaleData(object = subcluster_rice_combined, verbose = TRUE)
subcluster_rice_combined <- RunPCA(object = subcluster_rice_combined, npcs = 30, verbose = FALSE)
subcluster_rice_combined <- RunUMAP(object = subcluster_rice_combined, reduction = "pca", dims = 1:30)
subcluster_rice_combined <- FindNeighbors(object = subcluster_rice_combined, reduction = "pca", dims = 1:30)
subcluster_rice_combined <- FindClusters(subcluster_rice_combined, resolution = 0.4)
save(subcluster_rice_combined, file = "L2_rice_vasculature_combined.RData")

# Epidermis
subcluster_rice_combined<-rice_combined[,WhichCells(rice_combined, idents = c(3,15,14,4,7))]
DefaultAssay(object = subcluster_rice_combined) <- "integrated"
subcluster_rice_combined <- ScaleData(object = subcluster_rice_combined, verbose = TRUE)
subcluster_rice_combined <- RunPCA(object = subcluster_rice_combined, npcs = 30, verbose = FALSE)
subcluster_rice_combined <- RunUMAP(object = subcluster_rice_combined, reduction = "pca", dims = 1:30)
subcluster_rice_combined <- FindNeighbors(object = subcluster_rice_combined, reduction = "pca", dims = 1:30)
subcluster_rice_combined <- FindClusters(subcluster_rice_combined, resolution = 0.4)
save(subcluster_rice_combined, file = "L2_rice_epidermis_combined.RData")
