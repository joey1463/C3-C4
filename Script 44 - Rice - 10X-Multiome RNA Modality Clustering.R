# This script will read in and cluster the Rice 10X-Multiome RNA data. 

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
Rice_Light_3_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/3_Rice_Light/filtered_feature_bc_matrix.h5", 'Rice_Light_3_RNA')
Rice_Light_7_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/7_Rice_Light/filtered_feature_bc_matrix.h5", 'Rice_Light_7_RNA')
Rice_Dark_4_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/4_Rice_Dark/filtered_feature_bc_matrix.h5", 'Rice_Dark_4_RNA')
Rice_Dark_8_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/8_Rice_Dark/filtered_feature_bc_matrix.h5", 'Rice_Dark_8_RNA')

# Rep 2
Rice_Dark_1_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Dark_1/filtered_feature_bc_matrix.h5", 'Rice_Dark_1_Rep2_RNA' )
Rice_Dark_2_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Dark_2/filtered_feature_bc_matrix.h5", 'Rice_Dark_2_Rep2_RNA' )
Rice_Dark_3_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Dark_3/filtered_feature_bc_matrix.h5", 'Rice_Dark_3_Rep2_RNA' )

Rice_Light_1_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Light_1/filtered_feature_bc_matrix.h5", 'Rice_Light_1_Rep2_RNA' )
Rice_Light_2_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Light_2/filtered_feature_bc_matrix.h5", 'Rice_Light_2_Rep2_RNA' )
Rice_Light_3_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Light_3/filtered_feature_bc_matrix.h5", 'Rice_Light_3_Rep2_RNA' )
Rice_Light_4_Rep2_RNA <-Read_in_RNA("~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Light_4/filtered_feature_bc_matrix.h5", 'Rice_Light_4_Rep2_RNA' )

setwd("~/Desktop/Salk/Project - Light/R files/")
save(Rice_Light_3_RNA,file="Rice_Light_3_RNA.RData")
save(Rice_Light_7_RNA,file="Rice_Light_7_RNA.RData")
save(Rice_Dark_4_RNA,file="Rice_Dark_4_RNA.RData")
save(Rice_Dark_8_RNA,file="Rice_Dark_8_RNA.RData")
save(Rice_Dark_1_Rep2_RNA,file="Rice_Dark_1_Rep2_RNA.RData")
save(Rice_Dark_2_Rep2_RNA,file="Rice_Dark_2_Rep2_RNA.RData")
save(Rice_Dark_3_Rep2_RNA,file="Rice_Dark_3_Rep2_RNA.RData")
save(Rice_Light_1_Rep2_RNA,file="Rice_Light_1_Rep2_RNA.RData")
save(Rice_Light_2_Rep2_RNA,file="Rice_Light_2_Rep2_RNA.RData")
save(Rice_Light_3_Rep2_RNA,file="Rice_Light_3_Rep2_RNA.RData")
save(Rice_Light_4_Rep2_RNA,file="Rice_Light_4_Rep2_RNA.RData")

plant_anchors <- FindIntegrationAnchors(object.list = list(Rice_Light_3_RNA,
                                                           Rice_Light_7_RNA,
                                                           Rice_Dark_4_RNA,
                                                           Rice_Dark_8_RNA,
                                                           # Rep 2
                                                           Rice_Dark_1_Rep2_RNA,
                                                           Rice_Dark_2_Rep2_RNA,
                                                           Rice_Dark_3_Rep2_RNA,
                                                           Rice_Light_1_Rep2_RNA,
                                                           Rice_Light_2_Rep2_RNA,
                                                           Rice_Light_3_Rep2_RNA,
                                                           Rice_Light_4_Rep2_RNA), dims = 1:30)

rice_combined_integrated <- IntegrateData(anchorset = plant_anchors, dims = 1:30)

DefaultAssay(object = rice_combined_integrated) <- "integrated"
rice_combined_integrated <- ScaleData(object = rice_combined_integrated, vars.to.regress = c('nCount_RNA'),verbose = TRUE)
rice_combined_integrated <- RunPCA(object = rice_combined_integrated, npcs = 50, verbose = FALSE)
rice_combined_integrated <- RunUMAP(object = rice_combined_integrated, reduction = "pca", dims = 1:50)
rice_combined_integrated <- FindNeighbors(object = rice_combined_integrated, reduction = "pca", dims = 1:50)
rice_combined_integrated <- FindClusters(rice_combined_integrated, resolution = 0.5)  

# Metrics
length(rice_combined_integrated$nFeature_RNA) 
median(rice_combined_integrated$nFeature_RNA) 
median(rice_combined_integrated$nCount_RNA) 
DimPlot(object = rice_combined_integrated, reduction = "umap", label=TRUE, raster = FALSE) + NoLegend()
DimPlot(object = rice_combined_integrated, reduction = "umap", label=TRUE, raster = FALSE, split.by='Treatment') + NoLegend()
rice_multiome_RNA_combined<-rice_combined_integrated

setwd("~/Desktop/Salk/Project - Light/R files")
save(rice_multiome_RNA_combined,file="L1_rice_multiome_RNA_combined.RData")
