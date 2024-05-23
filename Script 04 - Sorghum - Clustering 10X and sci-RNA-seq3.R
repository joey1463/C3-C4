# This script will take sci-RNA-seq3 and 10X-RNA Sorghum experiments and integrate them into one Seurat object. 


library(Seurat)
library(Matrix)
library(stringr)
library(DoubletFinder)
load("CI_full_matrix_sorghum.RData")
load("CI_full_matrix_sorghum_48hr.RData")

### 10X 
Read_10X<-function(path, name){
  TimeX_Light <- Read10X(data.dir = path)
  rownames(TimeX_Light)<-substring(rownames(TimeX_Light),1,16)
  TimeX_Light <- CreateSeuratObject(counts = TimeX_Light, min.cells = 200, min.features = 200) 
  TimeX_Light$Time <- name
  TimeX_Light <- subset(x = TimeX_Light, nFeature_RNA > 200 & nFeature_RNA < 5000) 
  TimeX_Light <- NormalizeData(object = TimeX_Light, verbose = FALSE)
  TimeX_Light <- FindVariableFeatures(object = TimeX_Light, selection.method = "vst", nfeatures = 6000)
  DefaultAssay(object = TimeX_Light) <- "RNA"
  TimeX_Light <- ScaleData(object = TimeX_Light , vars.to.regress = c('nCount_RNA'),verbose = TRUE)
  TimeX_Light <- RunPCA(object = TimeX_Light , npcs = 30, verbose = FALSE)
  TimeX_Light <- RunUMAP(object = TimeX_Light, reduction = "pca", dims = 1:30)
  TimeX_Light <- FindNeighbors(object = TimeX_Light, reduction = "pca", dims = 1:30)
  TimeX_Light <- FindClusters(TimeX_Light, resolution = 0.4)
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

Sorghum_0<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_0_Sorghum/10X_0_Sorghum/outs/filtered_feature_bc_matrix","Sorghum_0")
Sorghum_0.5<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_0.5_Sorghum/10X_0point5_Sorghum/outs/filtered_feature_bc_matrix","Sorghum_0.5")
Sorghum_1<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_1_Sorghum/10X_1_Sorghum/outs/filtered_feature_bc_matrix","Sorghum_01")
Sorghum_2<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_2_Sorghum/10X_2_Sorghum/outs/filtered_feature_bc_matrix","Sorghum_02")
Sorghum_4<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_4_Sorghum/10X_4_Sorghum/outs/filtered_feature_bc_matrix","Sorghum_04")
Sorghum_6<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_6_Sorghum/10X_6_Sorghum/outs/filtered_feature_bc_matrix","Sorghum_06")
Sorghum_8<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_8_Sorghum/10X_8_Sorghum/outs/filtered_feature_bc_matrix","Sorghum_08")
Sorghum_12<-Read_10X("/gale/netapp/home/jswift/analysis/Light_10X/10X_12_Sorghum/10X_12_Sorghum/outs/filtered_feature_bc_matrix","Sorghum_12")

Sorghum_0<-Remove_Doublets(Sorghum_0,0.0475)
Sorghum_0.5<-Remove_Doublets(Sorghum_0.5,0.0893)
Sorghum_1<-Remove_Doublets(Sorghum_1,0.1030)
Sorghum_2<-Remove_Doublets(Sorghum_2,0.1802)
Sorghum_4<-Remove_Doublets(Sorghum_4,0.1384)
Sorghum_6<-Remove_Doublets(Sorghum_6,0.0241)
Sorghum_8<-Remove_Doublets(Sorghum_8,0.1521)
Sorghum_12<-Remove_Doublets(Sorghum_12,0.0527)

### CI
Read_CI<-function(path, minimum){
  DayX_Light <- CreateSeuratObject(counts = path, min.cells = 10, min.features = 10)
  DayX_Light$Time <- "CI"
  DayX_Light <- subset(x = DayX_Light, nFeature_RNA > minimum & nFeature_RNA < 3000)
  DayX_Light <- ScaleData(object = DayX_Light ,verbose = TRUE) #  vars.to.regress = c('nCount_RNA')
  DayX_Light <- NormalizeData(object = DayX_Light, verbose = FALSE)
  DayX_Light <- FindVariableFeatures(object = DayX_Light, selection.method = "vst", nfeatures = 2000)
  DayX_Light}

### CI
Sorghum_0_CI<-Read_CI(full_matrix_sorghum[,(substr(colnames(full_matrix_sorghum),1,8) == "Sorg_000")], minimum = 200)
Sorghum_0.5_CI<-Read_CI(full_matrix_sorghum[,(substr(colnames(full_matrix_sorghum),1,8) == "Sorg_0.5")], minimum = 200)
Sorghum_1_CI<-Read_CI(full_matrix_sorghum[,(substr(colnames(full_matrix_sorghum),1,8) == "Sorg_001")], minimum = 200)
Sorghum_2_CI<-Read_CI(full_matrix_sorghum[,(substr(colnames(full_matrix_sorghum),1,8) == "Sorg_002")], minimum = 200)
Sorghum_4_CI<-Read_CI(full_matrix_sorghum[,(substr(colnames(full_matrix_sorghum),1,8) == "Sorg_004")], minimum = 200)
Sorghum_6_CI<-Read_CI(full_matrix_sorghum[,(substr(colnames(full_matrix_sorghum),1,8) == "Sorg_006")], minimum = 200)
Sorghum_8_CI<-Read_CI(full_matrix_sorghum[,(substr(colnames(full_matrix_sorghum),1,8) == "Sorg_008")], minimum = 200)
Sorghum_12_CI<-Read_CI(full_matrix_sorghum[,(substr(colnames(full_matrix_sorghum),1,8) == "Sorg_012")], minimum = 200)
Sorghum_48_CI<-Read_CI(full_matrix_sorghum_48hr, minimum = 400)

### Merge
sorghum_combined_merged <- merge(Sorghum_0, y = c(Sorghum_0.5, Sorghum_1, Sorghum_2, Sorghum_4, Sorghum_6, Sorghum_8, Sorghum_12,
                                                  Sorghum_0_CI, Sorghum_0.5_CI, Sorghum_1_CI, Sorghum_2_CI, Sorghum_4_CI, Sorghum_6_CI, Sorghum_8_CI, Sorghum_12_CI,
                                                  Sorghum_48_CI))
sorghum_combined_merged <- NormalizeData(object = sorghum_combined_merged , verbose = FALSE)
sorghum_combined_merged <- FindVariableFeatures(object = sorghum_combined_merged , selection.method = "vst", nfeatures = 2000) 
save(sorghum_combined_merged, file = "sorghum_combined_merged.RData")

### Subset by Variable Features
Variable_Features<-VariableFeatures(sorghum_combined_merged)
length(Variable_Features)
Sorghum_0<-Sorghum_0[Variable_Features,]
Sorghum_0.5<-Sorghum_0.5[Variable_Features,]
Sorghum_1<-Sorghum_1[Variable_Features,]
Sorghum_2<-Sorghum_2[Variable_Features,]
Sorghum_4<-Sorghum_4[Variable_Features,]
Sorghum_6<-Sorghum_6[Variable_Features,]
Sorghum_8<-Sorghum_8[Variable_Features,]
Sorghum_12<-Sorghum_12[Variable_Features,]
Sorghum_0_CI<-Sorghum_0_CI[Variable_Features,]
Sorghum_0.5_CI<-Sorghum_0.5_CI[Variable_Features,]
Sorghum_1_CI<-Sorghum_1_CI[Variable_Features,]
Sorghum_2_CI<-Sorghum_2_CI[Variable_Features,]
Sorghum_4_CI<-Sorghum_4_CI[Variable_Features,]
Sorghum_6_CI<-Sorghum_6_CI[Variable_Features,]
Sorghum_8_CI<-Sorghum_8_CI[Variable_Features,]
Sorghum_12_CI<-Sorghum_12_CI[Variable_Features,]
Sorghum_48_CI<-Sorghum_48_CI[Variable_Features,]

### Integrate
plant_anchors <- FindIntegrationAnchors(object.list = list(Sorghum_0, Sorghum_0.5, Sorghum_1, Sorghum_2, Sorghum_4, Sorghum_6, Sorghum_8, Sorghum_12,
                                                           Sorghum_0_CI, Sorghum_0.5_CI, Sorghum_1_CI, Sorghum_2_CI, Sorghum_4_CI, Sorghum_6_CI, Sorghum_8_CI, Sorghum_12_CI, Sorghum_48_CI), reference= c(1:8), dims = 1:30)

sorghum_combined_integrated <- IntegrateData(anchorset = plant_anchors, dims = 1:30)

remove(plant_anchors)
remove(Sorghum_0, Sorghum_0.5, Sorghum_1, Sorghum_2,Sorghum_4, Sorghum_6, Sorghum_8, Sorghum_12,
       Sorghum_0_CI, Sorghum_0.5_CI, Sorghum_1_CI, Sorghum_2_CI, Sorghum_4_CI, Sorghum_6_CI, Sorghum_8_CI, Sorghum_12_CI, Sorghum_48_CI)

DefaultAssay(object = sorghum_combined_integrated) <- "integrated"
sorghum_combined_integrated <- ScaleData(object = sorghum_combined_integrated, vars.to.regress = c('nCount_RNA'),verbose = TRUE)
sorghum_combined_integrated <- RunPCA(object = sorghum_combined_integrated, npcs = 50, verbose = FALSE)
sorghum_combined_integrated <- RunUMAP(object = sorghum_combined_integrated, reduction = "pca", dims = 1:50)
sorghum_combined_integrated <- FindNeighbors(object = sorghum_combined_integrated, reduction = "pca", dims = 1:50)
sorghum_combined_integrated <- FindClusters(sorghum_combined_integrated, resolution = 0.5) # 0.5 is important, if change resolution here, make sure you change it for transferring across down below!
save(sorghum_combined_integrated, file = "sorghum_combined_integrated.RData")

# Transfering UMAP and Meta Data
sorghum_combined_merged[['umap']] <- sorghum_combined_integrated[['umap']]
sorghum_combined_merged[['integrated']] <- sorghum_combined_integrated[['integrated']]
sorghum_combined_merged@meta.data$seurat_clusters <- sorghum_combined_integrated@meta.data$seurat_clusters
sorghum_combined_merged@meta.data$integrated_snn_res.0.5 <- sorghum_combined_integrated@meta.data$integrated_snn_res.0.5
sorghum_combined<- sorghum_combined_merged  # save over

# Factors
light_time<-substring(colnames(sorghum_combined),6,8)
for (i in 1:length(light_time)) { if (substring(light_time[i],1,1) %in% c('G','C','A','T')) {light_time[i] <- sorghum_combined$Time[i] }}
sorghum_combined$light_time <- light_time

assay_type<-sorghum_combined$Time
for (i in 1:length(assay_type)) { if (substr(assay_type[i],1,3) == 'Sor')  { assay_type[i]<- '10X' }}
sorghum_combined$assay_type <- assay_type

Idents(sorghum_combined)<-sorghum_combined@meta.data$integrated_snn_res.0.5
table(Idents(sorghum_combined))
sorghum_combined<-sorghum_combined[,WhichCells(sorghum_combined, idents = c(0:18))] # remove tiny clusters

length(sorghum_combined$nFeature_RNA) # no. of nuceli
median(sorghum_combined$nFeature_RNA) # median no. of genes.
median(sorghum_combined$nCount_RNA) # no. of reads 

save(sorghum_combined, file = "L1_sorghum_combined.RData")
sorghum_combined_subsample <- sorghum_combined[, sample(colnames(sorghum_combined), size = 80000, replace=F)]
save(sorghum_combined_subsample, file = "L1_sorghum_combined_subsample.RData")
jpeg('sorghum_combined.jpg')
plot(DimPlot(object = sorghum_combined, reduction = "umap", label = TRUE))
