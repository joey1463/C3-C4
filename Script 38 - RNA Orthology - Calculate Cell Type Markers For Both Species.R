# This script will calculate cell type specific marker genes for each of the 6 major cell types in each species.

library(Seurat)
library(stringr)
load("L1_rice_combined_labelled.RData")
load("L1_sorghum_combined_labelled.RData")

cell_types<-c("mesophyll","bundle_sheath","guard","epidermis","phloem","xylem")

rice_markers<-list(NA)
for (i in 1:6){
  print(i)
  markers<-FindMarkers(rice_combined_labelled, ident.1 = cell_types[i], verbose = FALSE, only.pos = TRUE, min.pct = 0.1)
  markers[,6]<-markers$pct.1/markers$pct.2
  markers<-markers[order(markers$p_val_adj,decreasing = FALSE),]
  markers<-subset(markers, markers$V6>1) 
  rice_markers[[i]]<-markers}

sorghum_markers<-list(NA)
for (i in 1:6){
  print(i)
  markers<-FindMarkers(sorghum_combined_labelled, ident.1 = cell_types[i], verbose = FALSE, only.pos = TRUE, min.pct = 0.1)
  markers[,6]<-markers$pct.1/markers$pct.2
  markers<-markers[order(markers$p_val_adj,decreasing = FALSE),]
  markers<-subset(markers, markers$V6>1) 
  sorghum_markers[[i]]<-markers}

save(rice_markers, file = "L3_Orthology_rice_markers.RData")
save(sorghum_markers, file = "L3_Orthology_sorghum_markers.RData")

# Export
setwd("~/Desktop/")
for (i in 1:6){
rice_markers_all<-rice_markers[[i]]
rice_markers_all<-subset(rice_markers_all, rice_markers_all$p_val_adj<0.01)
write.table(rice_markers_all, file = paste("rice_markers_all",i,"_markers.txt",sep=''),quote=F,sep="\t")}
for (i in 1:6){
sorghum_markers_all<-sorghum_markers[[i]]
sorghum_markers_all<-subset(sorghum_markers_all, sorghum_markers_all$p_val_adj<0.01)
write.table(sorghum_markers_all, file = paste("sorghum_markers_all",i,"_markers.txt",sep=''),quote=F,sep="\t")}
