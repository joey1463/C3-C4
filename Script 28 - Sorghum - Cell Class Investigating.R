# This script will visualize cell class clustertring 


# VASCULATURE
setwd("~/Desktop/Salk/Project - Light/R Files")
load("L2_sorghum_vasculature_combined.RData")
DimPlot(object = subcluster_sorghum_combined, reduction = "umap", label=TRUE, raster=FALSE) + NoLegend()
DimPlot(object = subcluster_sorghum_combined, reduction = "umap", label=TRUE, split.by = 'assay_type')
DimPlot(object = subcluster_sorghum_combined, reduction = "umap", label=TRUE, split.by = 'light_time')
DefaultAssay(object = subcluster_sorghum_combined) <- "RNA" 
FeaturePlot(object = subcluster_sorghum_combined, features = "nCount_RNA", cols = c("grey", "red"), raster=TRUE) 
FeaturePlot(object = subcluster_sorghum_combined, features = "Sobic.001G488700", max.cutoff = 20, min.cutoff = 0, cols = c("grey", "red"), raster=FALSE)

# MESOPHYLL
setwd("~/Desktop/Salk/Project - Light/R Files")
load("L2_sorghum_mesophyll_combined.RData")
DimPlot(object = subcluster_sorghum_combined, reduction = "umap", label=TRUE, raster=FALSE) + NoLegend()
DimPlot(object = subcluster_sorghum_combined, reduction = "umap", label=TRUE, split.by = 'assay_type')
DimPlot(object = subcluster_sorghum_combined, reduction = "umap", label=TRUE, split.by = 'light_time')
DefaultAssay(object = subcluster_sorghum_combined) <- "RNA" 
FeaturePlot(object = subcluster_sorghum_combined, features = "nCount_RNA", cols = c("grey", "red"), raster=FALSE) 

# EPIDERMIS
setwd("~/Desktop/Salk/Project - Light/R Files")
load("L2_sorghum_epidermis_combined.RData")
DimPlot(object = subcluster_sorghum_combined, reduction = "umap", label=TRUE, raster=FALSE) + NoLegend()
DefaultAssay(object = subcluster_sorghum_combined) <- "RNA" 
FeaturePlot(object = subcluster_sorghum_combined, features = "nCount_RNA", cols = c("grey", "red"), raster=TRUE) 
FeaturePlot(object = subcluster_sorghum_combined, features = "Sobic.010G164500", raster=FALSE) # PDF1
