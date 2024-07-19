# This script will read in the clustered Sorghum Atlas (RNA data only). The ful dataset is quite large, thus the "subsample" data may be used instead on local computer. 

setwd("~/Desktop/Salk/Project - Light/R files")

load("L1_sorghum_combined_subsample.RData")
sorghum_combined<-sorghum_combined_subsample
remove(sorghum_combined_subsample)
DefaultAssay(object = sorghum_combined) <- "RNA"

# Metrics
length(sorghum_combined$nFeature_RNA) 
median(sorghum_combined$nFeature_RNA) 
median(sorghum_combined$nCount_RNA)
DimPlot(object = sorghum_combined, reduction = "umap", label=TRUE, raster = TRUE) + NoLegend()
DimPlot(object = sorghum_combined, reduction = "umap", label=TRUE, raster = TRUE, split.by = 'Time') + NoLegend()
DimPlot(object = sorghum_combined, reduction = "umap", label=TRUE, raster = TRUE, split.by = 'assay_type') + NoLegend()

# Feature Plots 
FeaturePlot(object = sorghum_combined, features = "Sobic.003G036200", raster=FALSE, reduction = "umap", cols = c("grey90","red3")) #NADPME
FeaturePlot(object = sorghum_combined, features = "Sobic.003G234200", raster=FALSE, reduction = "umap", cols = c("grey90","red3")) #CA 
