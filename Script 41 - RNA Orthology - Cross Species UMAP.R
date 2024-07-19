# This script will read in the cross-species clustered data, and then label nuclei based off their annotation from either the rice or sorghum RNA atlas. 

setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_species_combined.RData")
DefaultAssay(object = species_combined) <- "RNA"

# Plot Unlabelled
length(species_combined$nFeature_RNA) 
median(species_combined$nFeature_RNA) 
median(species_combined$nCount_RNA) 
DimPlot(object = species_combined, reduction = "umap", label=TRUE, raster = TRUE) + NoLegend()
DimPlot(object = species_combined, reduction = "umap", label=TRUE, raster = FALSE, split.by = 'Time')

# Labeling with Rice Cell Types
load("L1_rice_combined_CI_labelled.RData")
rice_combined_CI_labelled<-rice_combined_CI_labelled[,substr(colnames(rice_combined_CI_labelled),6,7) == '48']
rice_combined_CI_labelled<-RenameCells(rice_combined_CI_labelled, new.names = str_sub(colnames(rice_combined_CI_labelled), end=-4))
species_combined<-species_combined[,colnames(rice_combined_CI_labelled)]
species_combined$cell_type<-Idents(rice_combined_CI_labelled)
Idents(species_combined)<-species_combined$cell_type
colors<-as.data.frame(table(Idents(species_combined)))
colors[,3]<-c("springgreen4","grey90","navajowhite3","grey90","deepskyblue4","steelblue1","turquoise3","tan3","grey90")
DimPlot(object = species_combined, reduction = "umap", label=FALSE, raster = FALSE, cols = colors[,3]) + theme(legend.position = "none")
FeaturePlot(object = species_combined, features = "LOC-Os01g09320", raster=FALSE,cols = c("grey90","steelblue3")) #NADPME
remove(species_combined)
remove(rice_combined_CI_labelled)

# Labeling with Sorghum Cell Types
load("L1_species_combined.RData")
load("L1_sorghum_combined_CI_labelled.RData")
sorghum_combined_CI_labelled<-sorghum_combined_CI_labelled[,substr(colnames(sorghum_combined_CI_labelled),6,7) == '48']
sorghum_combined_CI_labelled<-RenameCells(sorghum_combined_CI_labelled, new.names = str_sub(colnames(sorghum_combined_CI_labelled), end=-4))
species_combined<-species_combined[,colnames(sorghum_combined_CI_labelled)]
#sorghum_combined_CI_labelled$cell_type<-ifelse(substr(sorghum_combined_CI_labelled$cell_type,1,4) == "Sorg", "unknown", sorghum_combined_CI_labelled$cell_type)
species_combined$cell_type<-Idents(sorghum_combined_CI_labelled)
Idents(species_combined)<-species_combined$cell_type
colors<-as.data.frame(table(Idents(species_combined)))
colors[,3]<-c("springgreen4","grey90","steelblue1","navajowhite3","deepskyblue4","turquoise3","tan3")
DimPlot(object = species_combined, reduction = "umap", label=FALSE, raster = FALSE, cols = colors[,3])  + theme(legend.position = "none")
FeaturePlot(object = species_combined, features = "LOC-Os01g09320", max.cutoff = 4, raster=FALSE,cols = c("grey90","steelblue3")) #NADPME
remove(species_combined)
remove(sorghum_combined_CI_labelled)
