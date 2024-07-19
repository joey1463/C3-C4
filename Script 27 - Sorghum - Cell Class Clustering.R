# This script will take mesophyll, epidermal (guard, epidermis), and vasculature (bundle sheath, phloem, xylem) clusters and re-cluster them. 


library(Seurat)
load("L1_sorghum_combined.RData")

# EPIDERMIS
subcluster_sorghum_combined<-sorghum_combined[,WhichCells(sorghum_combined, idents = c(6,14,17,8))]
DefaultAssay(object = subcluster_sorghum_combined) <- "integrated"
subcluster_sorghum_combined <- ScaleData(object = subcluster_sorghum_combined, verbose = TRUE)
subcluster_sorghum_combined <- RunPCA(object = subcluster_sorghum_combined, npcs = 30, verbose = FALSE)
subcluster_sorghum_combined <- RunUMAP(object = subcluster_sorghum_combined, reduction = "pca", dims = 1:30)
subcluster_sorghum_combined <- FindNeighbors(object = subcluster_sorghum_combined, reduction = "pca", dims = 1:30)
subcluster_sorghum_combined <- FindClusters(subcluster_sorghum_combined, resolution = 0.4)
save(subcluster_sorghum_combined, file = "L2_sorghum_epidermis_combined.RData")

# VASCULATURE
subcluster_sorghum_combined<-sorghum_combined[,WhichCells(sorghum_combined, idents = c(15,4,13,12,9,10))]
DefaultAssay(object = subcluster_sorghum_combined) <- "integrated"
subcluster_sorghum_combined <- ScaleData(object = subcluster_sorghum_combined, verbose = TRUE)
subcluster_sorghum_combined <- RunPCA(object = subcluster_sorghum_combined, npcs = 30, verbose = FALSE)
subcluster_sorghum_combined <- RunUMAP(object = subcluster_sorghum_combined, reduction = "pca", dims = 1:30)
subcluster_sorghum_combined <- FindNeighbors(object = subcluster_sorghum_combined, reduction = "pca", dims = 1:30)
subcluster_sorghum_combined <- FindClusters(subcluster_sorghum_combined, resolution = 0.4)
save(subcluster_sorghum_combined, file = "L2_sorghum_vasculature_combined.RData")

# MESOPHYLL
subcluster_sorghum_combined<-sorghum_combined[,WhichCells(sorghum_combined, idents = c(0,1))]
DefaultAssay(object = subcluster_sorghum_combined) <- "integrated"
subcluster_sorghum_combined <- ScaleData(object = subcluster_sorghum_combined, verbose = TRUE)
subcluster_sorghum_combined <- RunPCA(object = subcluster_sorghum_combined, npcs = 30, verbose = FALSE)
subcluster_sorghum_combined <- RunUMAP(object = subcluster_sorghum_combined, reduction = "pca", dims = 1:30)
subcluster_sorghum_combined <- FindNeighbors(object = subcluster_sorghum_combined, reduction = "pca", dims = 1:30)
subcluster_sorghum_combined <- FindClusters(subcluster_sorghum_combined, resolution = 0.4)
save(subcluster_sorghum_combined, file = "L2_sorghum_mesophyll_combined.RData")
