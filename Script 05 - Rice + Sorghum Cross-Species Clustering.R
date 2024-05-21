# This script will take the OrthoMCL Rice-Sorghum ortholog matches and select the first orthologs listed. This file is required for cross-species clustering. 
# This script will then convert sorghum to rice orthologs and cluster both species' 48h time point in the same Seurat object. 

library(Seurat)
library(Matrix)
library(stringr)

orthology<-read.table('Rice_Sorghum_Orthology.txt',header=TRUE, sep="\t") # can be found in Mendeley folder
orthology_refined<-as.data.frame(orthology)
orthology_refined$Sorghum<-substr(orthology_refined$Sorghum, 1,16)
orthology_refined$Rice<-substr(orthology_refined$Rice, 1,14)
colnames(orthology_refined)<-c("Orthogroup","sorghum_group","rice_group")
save(orthology_refined, file = "Orthology_1_to_1.RData")

#setwd("~/Desktop/Salk/Project - Light/R files")
load("Orthology_1_to_1.RData")

load("CI_full_matrix_rice_48hr.RData")
New_Names<-rep(NA,length(rownames(full_matrix_rice_48hr)))
for (i in 1:length(rownames(full_matrix_rice_48hr))) {if ((rownames(full_matrix_rice_48hr)[i] %in% orthology_refined$rice_group) == TRUE) {New_Names[i]<- subset(orthology_refined, orthology_refined$rice_group == rownames(full_matrix_rice_48hr)[i])$rice_group }}
New_Names[is.na(New_Names)] <- "none"
Rice_48_CI<-full_matrix_rice_48hr[substr(New_Names,1,1) =="L",]
remove(full_matrix_rice_48hr)

load("CI_full_matrix_sorghum_48hr.RData")
New_Names<-rep(NA,length(rownames(full_matrix_sorghum_48hr)))
for (i in 1:length(rownames(full_matrix_sorghum_48hr))) {if ((rownames(full_matrix_sorghum_48hr)[i] %in% orthology_refined$sorghum_group) == TRUE) {New_Names[i]<- subset(orthology_refined, orthology_refined$sorghum_group == rownames(full_matrix_sorghum_48hr)[i])$rice_group }}
New_Names[is.na(New_Names)] <- "none"
Sorghum_48_CI<-full_matrix_sorghum_48hr[substr(New_Names,1,1) =="L",]
rownames(Sorghum_48_CI)<-New_Names[substr(New_Names,1,1) =="L"]
remove(full_matrix_sorghum_48hr)

Read_CI<-function(path, species, minimum){
  DayX_Light <- CreateSeuratObject(counts = path, min.cells = 10, min.features = 200)
  DayX_Light$Time <- species
  DayX_Light <- subset(x = DayX_Light, nFeature_RNA > minimum & nFeature_RNA < 3000)
  DayX_Light <- ScaleData(object = DayX_Light ,verbose = TRUE) 
  DayX_Light <- NormalizeData(object = DayX_Light, verbose = FALSE)
  DayX_Light <- FindVariableFeatures(object = DayX_Light, selection.method = "vst", nfeatures = 2000)
  DayX_Light}

Rice_48_CI<-Read_CI(Rice_48_CI, minimum = 450, species ="Rice_CI")
Sorghum_48_CI<-Read_CI(Sorghum_48_CI, minimum = 400, species ="Sorghum_CI")

# Integration
plant_anchors <- FindIntegrationAnchors(object.list = list(Rice_48_CI, Sorghum_48_CI), dims = 1:30)
species_combined_integrated <- IntegrateData(anchorset = plant_anchors, dims = 1:30)

remove(plant_anchors)
remove(Rice_48_CI, Sorghum_48_CI)

DefaultAssay(object = species_combined_integrated) <- "integrated"
species_combined_integrated <- ScaleData(object = species_combined_integrated, vars.to.regress = c('nCount_RNA'),verbose = TRUE)
species_combined_integrated <- RunPCA(object = species_combined_integrated, npcs = 50, verbose = FALSE)
species_combined_integrated <- RunUMAP(object = species_combined_integrated, reduction = "pca", dims = 1:50)
species_combined_integrated <- FindNeighbors(object = species_combined_integrated, reduction = "pca", dims = 1:50)
species_combined_integrated <- FindClusters(species_combined_integrated, resolution = 0.4)  

species_combined<-species_combined_integrated
length(species_combined$nFeature_RNA) # no. of nuceli
median(species_combined$nFeature_RNA) # median no. of genes.
median(species_combined$nCount_RNA) # no. of reads 

save(species_combined, file = "L1_species_combined.RData")
