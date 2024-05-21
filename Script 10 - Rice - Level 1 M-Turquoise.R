# This script will read in the 10X-RNA nuclei data sourced from the rice mTurqouise line. The script combines 2 replicates and clusters these data. 

Rice_mT_1 <- Read10X(data.dir = "~/Desktop/Salk/Project - Light/Nuceli 10X/Sequencing/Rice_mTurquoise_trial2_rep1") # files can be found in Mendeley folder
Rice_mT_1 <- CreateSeuratObject(counts = Rice_mT_1, min.cells = 3, min.features = 200) 
Rice_mT_1 <- ScaleData(object = Rice_mT_1,verbose = TRUE)
Rice_mT_1 <- NormalizeData(object = Rice_mT_1 , verbose = FALSE)
Rice_mT_1  <- FindVariableFeatures(object = Rice_mT_1 , selection.method = "vst", nfeatures = 2000)
Rice_mT_1$Time <- 'Rice_mT1'

Rice_mT_2 <- Read10X(data.dir = "~/Desktop/Salk/Project - Light/Nuceli 10X/Sequencing/Rice_mTurquoise_trial2_rep2") # files can be found in Mendeley folder
Rice_mT_2 <- CreateSeuratObject(counts = Rice_mT_2, min.cells = 3, min.features = 200) 
Rice_mT_2 <- ScaleData(object = Rice_mT_2,verbose = TRUE)
Rice_mT_2 <- NormalizeData(object = Rice_mT_2 , verbose = FALSE)
Rice_mT_2  <- FindVariableFeatures(object = Rice_mT_2 , selection.method = "vst", nfeatures = 2000)
Rice_mT_2$Time <- 'Rice_mT2'

plant_anchors <- FindIntegrationAnchors(object.list = list(Rice_mT_1, Rice_mT_2), dims = 1:30) 
rice_mT_combined <- IntegrateData(anchorset = plant_anchors, dims = 1:30)

DefaultAssay(object = rice_mT_combined) <- "integrated"
rice_mT_combined <- ScaleData(object = rice_mT_combined, vars.to.regress = c('nCount_RNA'),verbose = TRUE)
rice_mT_combined <- RunPCA(object = rice_mT_combined, npcs = 30, verbose = FALSE)
rice_mT_combined <- RunUMAP(object = rice_mT_combined, reduction = "pca", dims = 1:30)
rice_mT_combined <- FindNeighbors(object = rice_mT_combined, reduction = "pca", dims = 1:30)
rice_mT_combined <- FindClusters(rice_mT_combined, resolution = 0.3) 

DefaultAssay(object = rice_mT_combined) <- "RNA"
DimPlot(object = rice_mT_combined, reduction = "umap", label = TRUE)
DimPlot(object = rice_mT_combined, reduction = "umap", label = TRUE, split.by = 'Time')
FeaturePlot(object = rice_mT_combined, features = "INTACT-eGFP", cols = c("grey90","green4"),max.cutoff = 2)
FeaturePlot(object = rice_mT_combined, features = "H2B-mTurquoise2", cols = c("grey90","blue4"),max.cutoff = 2)
VlnPlot(object = rice_mT_combined, features = "H2B-mTurquoise2")
VlnPlot(object = rice_mT_combined, features = "INTACT-eGFP")
FeaturePlot(object = rice_mT_combined, features = Validated_Markers$Gene_ID[c(28,34,35)], cols = c("grey80","blue4"),max.cutoff = 1)

# Find All Markers
setwd("~/Desktop/")
mT_Markers<-list(NA)
for (i in 1:13){
  markers<-FindMarkers(rice_mT_combined, ident.1 = as.character(i-1), verbose = FALSE, only.pos = TRUE)
  markers<-subset(markers, markers$p_val_adj<0.01)
  markers[,6]<-markers$pct.1/markers$pct.2
  mT_Markers[[i]]<-markers
  write.table(markers, file = paste("mTurq_clusters_",i-1,"_markers.txt",sep=''),quote=F,sep="\t")}

# Make Nice Plot 
FeaturePlot(object = rice_mT_combined, features = "H2B-mTurquoise2", cols = c("grey85","dodgerblue3"),max.cutoff = 1, pt.size = 0.5)
FeaturePlot(object = rice_mT_combined, features = "INTACT-eGFP", cols = c("grey80","green4"),max.cutoff = 1, pt.size = 0.5)

setwd("~/Desktop/Salk/Project - Light/R files")
save(rice_mT_combined, file = "L1_rice_mT_combined.RData")
save(mT_Markers, file = "L1_all_mT_Markers.RData")
