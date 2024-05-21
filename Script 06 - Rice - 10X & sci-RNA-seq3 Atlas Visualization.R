# This script will read in the clustered Rice Atlas (RNA data only). The ful dataset is quite large, thus the "subsample" data may be used instead on local computer. 

setwd("~/Desktop/Salk/Project - Light/R files")

load("L1_rice_combined_subsample.RData") # subsetted data for lower-memory environments
rice_combined<-rice_combined_subsample
remove(rice_combined_subsample)

#load("L1_rice_combined.RData") # full data set

DefaultAssay(object = rice_combined) <- "RNA"

# Metrics
length(rice_combined$nFeature_RNA) 
median(rice_combined$nFeature_RNA) 
median(rice_combined$nCount_RNA) 
DimPlot(object = rice_combined, reduction = "umap", label=TRUE, raster = TRUE) + NoLegend()
DimPlot(object = rice_combined, reduction = "umap", label=TRUE, raster = TRUE, split.by = 'assay_type')
DimPlot(object = rice_combined, reduction = "umap", label=TRUE, raster = TRUE, split.by = 'light_time')
DimPlot(object = rice_combined, reduction = "umap", label=TRUE, raster = TRUE, split.by = 'Time')
