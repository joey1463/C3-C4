# This script will calculate the top 200 cluster specific markers for each cluster. 

library(Seurat)
load("L1_sorghum_combined.RData") 

# Find Markers 
for (i in 0:18){
  markers<-FindMarkers(sorghum_combined, ident.1 = as.character(i), verbose = FALSE, only.pos = TRUE)
  markers[,6]<-markers$pct.1/markers$pct.2
  markers<-markers[order(markers$p_val_adj,decreasing = FALSE),]
  markers<-markers[1:200,]
  write.table(markers, file = paste("L1_sorghum_cluster_",i,"_markers.txt",sep=''),quote=F,sep="\t")}
