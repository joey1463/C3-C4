# This script will read in and cluster the Sorghum 10X-Multiome RNA data. 

Read_in_RNA<-function(path, treatment){
  TimeX_Light <- Read10X_h5(path)
  TimeX_Light <- CreateSeuratObject(counts = TimeX_Light$`Gene Expression`, assay = "RNA", min.cells = 20, min.features = 100)
  TimeX_Light  <- subset(x = TimeX_Light, nFeature_RNA > 300 & nFeature_RNA < 5000) 
  TimeX_Light <- NormalizeData(object = TimeX_Light, verbose = FALSE)
  TimeX_Light <- FindVariableFeatures(object = TimeX_Light, selection.method = "vst", nfeatures = 6000)
  DefaultAssay(object = TimeX_Light) <- "RNA"
  TimeX_Light$Treatment<- treatment
  TimeX_Light}
  
# Rep 1
Sorghum_Light_1_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/1_Sorghum_Light/filtered_feature_bc_matrix.h5","Sorghum_Light_1_RNA")
Sorghum_Light_5_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/5_Sorghum_Light/filtered_feature_bc_matrix.h5","Sorghum_Light_5_RNA")
Sorghum_Dark_2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/2_Sorghum_Dark/filtered_feature_bc_matrix.h5","Sorghum_Dark_2_RNA")
Sorghum_Dark_6_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/6_Sorghum_Dark/filtered_feature_bc_matrix.h5","Sorghum_Dark_6_RNA")

# Rep 2
Sorghum_Dark_1_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Dark_1/filtered_feature_bc_matrix.h5", 'Sorghum_Dark_1_Rep2_RNA' )
Sorghum_Dark_2_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Dark_2/filtered_feature_bc_matrix.h5", 'Sorghum_Dark_2_Rep2_RNA' )
Sorghum_Dark_3_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Dark_3/filtered_feature_bc_matrix.h5", 'Sorghum_Dark_3_Rep2_RNA' )
Sorghum_Dark_4_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Dark_4/filtered_feature_bc_matrix.h5", 'Sorghum_Dark_4_Rep2_RNA' )

Sorghum_Light_1_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Light_1/filtered_feature_bc_matrix.h5", 'Sorghum_Light_1_Rep2_RNA' )
Sorghum_Light_3_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Light_3/filtered_feature_bc_matrix.h5", 'Sorghum_Light_3_Rep2_RNA' )
Sorghum_Light_4_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Light_4/filtered_feature_bc_matrix.h5", 'Sorghum_Light_4_Rep2_RNA' )

setwd("~/Desktop/Salk/Project - Light/R files")
save(Sorghum_Light_1_RNA,file="Sorghum_Light_1_RNA.RData")
save(Sorghum_Light_5_RNA,file="Sorghum_Light_5_RNA.RData")
save(Sorghum_Dark_2_RNA,file="Sorghum_Dark_2_RNA.RData")
save(Sorghum_Dark_6_RNA,file="Sorghum_Dark_6_RNA.RData")
save(Sorghum_Dark_1_Rep2_RNA,file="Sorghum_Dark_1_Rep2_RNA.RData")
save(Sorghum_Dark_2_Rep2_RNA,file="Sorghum_Dark_2_Rep2_RNA.RData")
save(Sorghum_Dark_3_Rep2_RNA,file="Sorghum_Dark_3_Rep2_RNA.RData")
save(Sorghum_Dark_4_Rep2_RNA,file="Sorghum_Dark_4_Rep2_RNA.RData")
save(Sorghum_Light_1_Rep2_RNA,file="Sorghum_Light_1_Rep2_RNA.RData")
save(Sorghum_Light_3_Rep2_RNA,file="Sorghum_Light_3_Rep2_RNA.RData")
save(Sorghum_Light_4_Rep2_RNA,file="Sorghum_Light_4_Rep2_RNA.RData")

plant_anchors <- FindIntegrationAnchors(object.list = list(Sorghum_Light_1_RNA, Sorghum_Light_5_RNA, Sorghum_Dark_2_RNA, Sorghum_Dark_6_RNA,
                                                           Sorghum_Dark_1_Rep2_RNA,Sorghum_Dark_2_Rep2_RNA,Sorghum_Dark_3_Rep2_RNA,Sorghum_Dark_4_Rep2_RNA,
                                                           Sorghum_Light_1_Rep2_RNA,Sorghum_Light_3_Rep2_RNA,Sorghum_Light_4_Rep2_RNA), dims = 1:30)
sorghum_combined_integrated <- IntegrateData(anchorset = plant_anchors, dims = 1:30)


DefaultAssay(object = sorghum_combined_integrated) <- "integrated"
sorghum_combined_integrated <- ScaleData(object = sorghum_combined_integrated, vars.to.regress = c('nCount_RNA'),verbose = TRUE)
sorghum_combined_integrated <- RunPCA(object = sorghum_combined_integrated, npcs = 50, verbose = FALSE)
sorghum_combined_integrated <- RunUMAP(object = sorghum_combined_integrated, reduction = "pca", dims = 1:50)
sorghum_combined_integrated <- FindNeighbors(object = sorghum_combined_integrated, reduction = "pca", dims = 1:50)
sorghum_combined_integrated <- FindClusters(sorghum_combined_integrated, resolution = 0.4) 
sorghum_multiome_RNA_combined <- sorghum_combined_integrated

# Metrics
length(sorghum_combined_integrated$nFeature_RNA) 
median(sorghum_combined_integrated$nFeature_RNA) 
median(sorghum_combined_integrated$nCount_RNA) 
DimPlot(object = sorghum_combined_integrated, reduction = "umap", label=TRUE, raster = FALSE) + NoLegend()
DimPlot(object = sorghum_combined_integrated, reduction = "umap", label=TRUE, raster = FALSE, split.by='Treatment') + NoLegend()

setwd("~/Desktop/Salk/Project - Light/R files")
save(sorghum_multiome_RNA_combined,file="L1_sorghum_multiome_RNA_combined.RData")
