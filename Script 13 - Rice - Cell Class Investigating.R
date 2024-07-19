# This script will visualize cell class clustertring 

# VASCULATURE
setwd("~/Desktop/Salk/Project - Light/R Files")
load("L2_rice_vasculature_combined.RData")
DimPlot(object = subcluster_rice_combined, reduction = "umap", label=TRUE, raster=TRUE) + NoLegend()
DimPlot(object = subcluster_rice_combined, reduction = "umap", label=TRUE, split.by = 'assay_type')
DimPlot(object = subcluster_rice_combined, reduction = "umap", label=TRUE, split.by = 'light_time')
DefaultAssay(object = subcluster_rice_combined) <- "RNA" 
FeaturePlot(object = subcluster_rice_combined, features = Validated_Markers$Gene_ID[c(28,34,35)], max.cutoff = 20, min.cutoff = 0, cols = c("grey", "red"), raster=FALSE) 

# MESOPHYLL
setwd("~/Desktop/Salk/Project - Light/R Files")
load("L2_rice_mesophyll_combined.RData")
DimPlot(object = subcluster_rice_combined, reduction = "umap", label=TRUE, raster=FALSE) + NoLegend()
DimPlot(object = subcluster_rice_combined, reduction = "umap", label=TRUE, split.by = 'assay_type')
DimPlot(object = subcluster_rice_combined, reduction = "umap", label=TRUE, split.by = 'light_time')
DefaultAssay(object = subcluster_rice_combined) <- "RNA" 
FeaturePlot(object = subcluster_rice_combined, features = 'xxx', max.cutoff = 20, min.cutoff = 0, cols = c("grey", "red"), raster=FALSE) 

# EPIDERMIS
setwd("~/Desktop/Salk/Project - Light/R Files")
load("L2_rice_epidermis_combined.RData")
DimPlot(object = subcluster_rice_combined, reduction = "umap", label=TRUE, raster=FALSE) + NoLegend()
DimPlot(object = subcluster_rice_combined, reduction = "umap", label=TRUE, split.by = 'assay_type')
DimPlot(object = subcluster_rice_combined, reduction = "umap", label=TRUE, split.by = 'light_time')
DefaultAssay(object = subcluster_rice_combined) <- "RNA" 
FeaturePlot(object = subcluster_rice_combined, features = "LOC-Os12g22284", max.cutoff = 20, min.cutoff = 0, cols = c("grey", "red"), raster=FALSE) 
