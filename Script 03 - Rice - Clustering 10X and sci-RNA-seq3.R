# This script will take sci-RNA-seq3 and 10X-RNA Rice experiments and integrate them into one Seurat object. 

library(Seurat)
library(Matrix)
library(stringr)
library(DoubletFinder)
load("CI_full_matrix_rice.RData")
load("CI_full_matrix_rice_48hr.RData")

#nFeature_RNA > 300
Read_10X<-function(path, name){
  TimeX_Light <- Read10X(data.dir = path)
  rownames(TimeX_Light)<-substring(rownames(TimeX_Light),1,16)
  TimeX_Light <- CreateSeuratObject(counts = TimeX_Light, min.cells = 200, min.features = 200) 
  TimeX_Light$Time <- name
  TimeX_Light  <- subset(x = TimeX_Light, nFeature_RNA > 300 & nFeature_RNA < 5000) 
  TimeX_Light <- NormalizeData(object = TimeX_Light, verbose = FALSE)
  TimeX_Light <- FindVariableFeatures(object = TimeX_Light, selection.method = "vst", nfeatures = 6000)
  DefaultAssay(object = TimeX_Light) <- "RNA"
  TimeX_Light  <- ScaleData(object = TimeX_Light , vars.to.regress = c('nCount_RNA'),verbose = TRUE)
  TimeX_Light  <- RunPCA(object = TimeX_Light , npcs = 30, verbose = FALSE)
  TimeX_Light <- RunUMAP(object = TimeX_Light, reduction = "pca", dims = 1:30)
  TimeX_Light<- FindNeighbors(object = TimeX_Light, reduction = "pca", dims = 1:30)
  TimeX_Light<- FindClusters(TimeX_Light, resolution = 0.4)
  TimeX_Light}

Remove_Doublets<-function(sample_input, collision_rate){
  # Calculate pK
  sample_input_sweep <- paramSweep_v3(sample_input , PCs = 1:30, sct = FALSE)
  sample_input_sweep_stats <- summarizeSweep(sample_input_sweep, GT = FALSE)
  sample_input_sweep_stats_pK <- find.pK(sample_input_sweep_stats)
  pK_pull<-as.data.frame(find.pK(sample_input_sweep_stats))
  pK_pull<-as.numeric(as.vector((subset(pK_pull, pK_pull$BCmetric == max(pK_pull$BCmetric))$pK)))
    
  ## Homotypic Doublet Proportion Estimate 
  annotations <- sample_input@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)      
  nExp_poi <- round(collision_rate*nrow(sample_input@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 
  
  ## Run DoubletFinder with varying classification stringencies 
  sample_input_doublet_finder <- doubletFinder_v3(sample_input, PCs = 1:30, pN = 0.25, pK = pK_pull, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(sample_input_doublet_finder @meta.data)[8]<-'doublet_call'
  
  doublet_cluster<-DimPlot(object = sample_input_doublet_finder , reduction = "umap", group.by = colnames(sample_input_doublet_finder @meta.data)[8], raster= TRUE)
  jpeg('doublet_cluster.jpg')
  plot(doublet_cluster)
  new_sample_doublets_removed<-subset(x = sample_input_doublet_finder, doublet_call == 'Singlet')
  new_sample_doublets_removed}

Rice_0<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_0_Rice/10X_0_Rice/outs/filtered_feature_bc_matrix","Rice_0")
Rice_0.5<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_0.5_Rice/10X_0point5_Rice/outs/filtered_feature_bc_matrix","Rice_0.5")
Rice_1<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_1_Rice/10X_1_Rice/outs/filtered_feature_bc_matrix","Rice_01")
Rice_2<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_2_Rice/10X_2_Rice/outs/filtered_feature_bc_matrix","Rice_02")
Rice_4<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_4_Rice/10X_4_Rice/outs/filtered_feature_bc_matrix","Rice_04")
Rice_6<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_6_Rice/10X_6_Rice/outs/filtered_feature_bc_matrix","Rice_06")
Rice_8<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_8_Rice/10X_8_Rice/outs/filtered_feature_bc_matrix","Rice_08")
Rice_12<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_12_Rice/10X_12_Rice/outs/filtered_feature_bc_matrix","Rice_12")

Rice_0<-Remove_Doublets(Rice_0,0.1576)
Rice_0.5<-Remove_Doublets(Rice_0.5,0.1030)
Rice_1<-Remove_Doublets(Rice_1,0.1950)
Rice_2<-Remove_Doublets(Rice_2,0.0491)
Rice_4<-Remove_Doublets(Rice_4,0.0990)
Rice_6<-Remove_Doublets(Rice_6,0.0779)
Rice_8<-Remove_Doublets(Rice_8,0.1521)
Rice_12<-Remove_Doublets(Rice_12,0.0527)

Read_CI<-function(path, minimum){
  DayX_Light <- CreateSeuratObject(counts = path, min.cells = 10, min.features = 10)
  DayX_Light$Time <- "CI"
  DayX_Light <- subset(x = DayX_Light, nFeature_RNA > minimum & nFeature_RNA < 3000)
  DayX_Light <- ScaleData(object = DayX_Light ,verbose = TRUE) 
  DayX_Light <- NormalizeData(object = DayX_Light, verbose = FALSE)
  DayX_Light <- FindVariableFeatures(object = DayX_Light, selection.method = "vst", nfeatures = 2000)
  DayX_Light}

### CI
Rice_0_CI<-Read_CI(full_matrix_rice[,(substr(colnames(full_matrix_rice),1,8) == "Rice_000")], minimum = 250)
Rice_0.5_CI<-Read_CI(full_matrix_rice[,(substr(colnames(full_matrix_rice),1,8) == "Rice_0.5")], minimum = 250)
Rice_1_CI<-Read_CI(full_matrix_rice[,(substr(colnames(full_matrix_rice),1,8) == "Rice_001")], minimum = 200)
Rice_2_CI<-Read_CI(full_matrix_rice[,(substr(colnames(full_matrix_rice),1,8) == "Rice_002")], minimum = 200)
Rice_4_CI<-Read_CI(full_matrix_rice[,(substr(colnames(full_matrix_rice),1,8) == "Rice_004")], minimum = 200)
Rice_6_CI<-Read_CI(full_matrix_rice[,(substr(colnames(full_matrix_rice),1,8) == "Rice_006")], minimum = 250)
Rice_8_CI<-Read_CI(full_matrix_rice[,(substr(colnames(full_matrix_rice),1,8) == "Rice_008")], minimum = 250)
Rice_12_CI<-Read_CI(full_matrix_rice[,(substr(colnames(full_matrix_rice),1,8) == "Rice_012")], minimum = 250)
Rice_48_CI<-Read_CI(full_matrix_rice_48hr, minimum = 450)

### Merge
rice_combined_merged <- merge(Rice_0, y = c(Rice_0.5, Rice_1, Rice_2, Rice_4, Rice_6, Rice_8, Rice_12,
                                                  Rice_0_CI, Rice_0.5_CI, Rice_1_CI, Rice_2_CI, Rice_4_CI, Rice_6_CI, Rice_8_CI, Rice_12_CI,
                                                  Rice_48_CI))
rice_combined_merged <- NormalizeData(object = rice_combined_merged , verbose = FALSE)
rice_combined_merged <- FindVariableFeatures(object = rice_combined_merged , selection.method = "vst", nfeatures = 2000) 
save(rice_combined_merged, file = "rice_combined_merged.RData")

### Subset by Variable Features
Variable_Features<-VariableFeatures(rice_combined_merged)
length(Variable_Features)
Rice_0<-Rice_0[Variable_Features,]
Rice_0.5<-Rice_0.5[Variable_Features,]
Rice_1<-Rice_1[Variable_Features,]
Rice_2<-Rice_2[Variable_Features,]
Rice_4<-Rice_4[Variable_Features,]
Rice_6<-Rice_6[Variable_Features,]
Rice_8<-Rice_8[Variable_Features,]
Rice_12<-Rice_12[Variable_Features,]
Rice_0_CI<-Rice_0_CI[Variable_Features,]
Rice_0.5_CI<-Rice_0.5_CI[Variable_Features,]
Rice_1_CI<-Rice_1_CI[Variable_Features,]
Rice_2_CI<-Rice_2_CI[Variable_Features,]
Rice_4_CI<-Rice_4_CI[Variable_Features,]
Rice_6_CI<-Rice_6_CI[Variable_Features,]
Rice_8_CI<-Rice_8_CI[Variable_Features,]
Rice_12_CI<-Rice_12_CI[Variable_Features,]
Rice_48_CI<-Rice_48_CI[Variable_Features,]

# Integration
plant_anchors <- FindIntegrationAnchors(object.list = list(Rice_0, Rice_0.5, Rice_1, Rice_2, Rice_4, Rice_6, Rice_8, Rice_12,
                                                           Rice_0_CI, Rice_0.5_CI, Rice_1_CI, Rice_2_CI, Rice_4_CI, Rice_6_CI, Rice_8_CI, Rice_12_CI, Rice_48_CI), reference= c(1:8), dims = 1:30)

rice_combined_integrated <- IntegrateData(anchorset = plant_anchors, dims = 1:30)

remove(plant_anchors)
remove(Rice_0, Rice_0.5, Rice_1, Rice_2,Rice_4, Rice_6, Rice_8, Rice_12, 
       Rice_0_CI, Rice_0.5_CI, Rice_1_CI, Rice_2_CI, Rice_4_CI, Rice_6_CI, Rice_8_CI, Rice_12_CI, Rice_48_CI)

DefaultAssay(object = rice_combined_integrated) <- "integrated"
rice_combined_integrated <- ScaleData(object = rice_combined_integrated, vars.to.regress = c('nCount_RNA'),verbose = TRUE)
rice_combined_integrated <- RunPCA(object = rice_combined_integrated, npcs = 50, verbose = FALSE)
rice_combined_integrated <- RunUMAP(object = rice_combined_integrated, reduction = "pca", dims = 1:50)
rice_combined_integrated <- FindNeighbors(object = rice_combined_integrated, reduction = "pca", dims = 1:50)
rice_combined_integrated <- FindClusters(rice_combined_integrated, resolution = 0.3)  # if change resolution here, make sure you change it down where you transfer identities!
save(rice_combined_integrated, file = "rice_combined_integrated.RData")

# Transfering UMAP and Meta Data
rice_combined_merged[['umap']] <- rice_combined_integrated[['umap']]
rice_combined_merged[['integrated']] <- rice_combined_integrated[['integrated']]
rice_combined_merged@meta.data$seurat_clusters <- rice_combined_integrated@meta.data$seurat_clusters
rice_combined_merged@meta.data$integrated_snn_res.0.3 <-rice_combined_integrated@meta.data$integrated_snn_res.0.3
rice_combined<-rice_combined_merged # save over

# Factors
light_time<-substring(colnames(rice_combined),6,8)
for (i in 1:length(light_time)) { if (substring(light_time[i],1,1) %in% c('G','C','A','T')) {light_time[i] <- rice_combined$Time[i] }}
light_time<-gsub("Rice_", "", light_time)
light_time<-as.numeric(gsub("_", "",light_time))
rice_combined$light_time <- light_time

assay_type<-rice_combined$Time
for (i in 1:length(assay_type)) { if (substr(assay_type[i],1,3) == 'Ric')  { assay_type[i]<- '10X' }}
rice_combined$assay_type <- assay_type

Idents(rice_combined)<-rice_combined@meta.data$integrated_snn_res.0.3
table(Idents(rice_combined))
rice_combined<-rice_combined[,WhichCells(rice_combined, idents = c(0:18))] # remove tiny clusters

length(rice_combined$nFeature_RNA) # no. of nuceli
median(rice_combined$nFeature_RNA) # median no. of genes.
median(rice_combined$nCount_RNA) # no. of reads 

save(rice_combined, file = "L1_rice_combined.RData")
rice_combined_subsample <- rice_combined[, sample(colnames(rice_combined), size = 80000, replace=F)]
save(rice_combined_subsample, file = "L1_rice_combined_subsample.RData")
jpeg('rice_combined.jpg')
plot(DimPlot(object = rice_combined, reduction = "umap", label = TRUE))
