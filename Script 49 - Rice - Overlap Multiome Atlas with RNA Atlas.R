# This script will compute cell type specific gene expression markers within the Rice Mutiome data, and compare resulting markers with those from the RNA Atlas data

setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_rice_multiome_RNA_combined_with_ATAC.RData")
DefaultAssay(rice_multiome_RNA_combined) <- "RNA"

cell_types<-c("mesophyll", "bundle_sheath","guard","epidermis","phloem","xylem")

# Find Markers of Multiome RNA 
Rice_Multiome_Markers<-list(NA)
for (i in 1:6){
  print(i)
  markers<-FindMarkers(rice_multiome_RNA_combined, ident.1 =cell_types[i], verbose = FALSE, only.pos = TRUE)
  markers[,6]<-markers$pct.1/markers$pct.2
  markers<-markers[order(markers$p_val_adj,decreasing = FALSE),]
  markers<-subset(markers, markers$p_val_adj<0.01)
  Rice_Multiome_Markers[[i]]<-markers}

# Load Markers of RNA Atlas 
load("L3_Orthology_rice_markers.RData")
Rice_RNA_Markers<-rice_markers


# Call Overlap
call_overlap_fisher<-function(rice_RNA, rice_Multiome){
  rice_RNA_list<-subset(rice_RNA, rice_RNA$p_val_adj<0.01 & rice_RNA$V6 > 1.3)
  rice_Multiome_list<-subset(rice_Multiome, rice_Multiome$p_val_adj<0.01 & rice_Multiome$V6 > 1.3)

  intersect<-sum(as.numeric((rownames(rice_RNA_list) %in% rownames(rice_Multiome_list))))
  rice_RNA_length<-length(rice_RNA_list[,1])
  rice_Multiome_length<-length(rice_Multiome_list[,1])
  matrix_to_test<-matrix(c(intersect,(rice_RNA_length-intersect),(rice_Multiome_length-intersect),30467),2,2) # using all expressed genes in RNA object
  p_value<-as.numeric(fisher.test(matrix_to_test, alternative='greater')[1])
  #intersect
  -log(p_value,10)}

call_marker_overlap<-matrix(NA, nrow=6, ncol=6)
for (a in 1:6) { 
  for (b in 1:6) {
     call_marker_overlap[a,b]<- call_overlap_fisher(Rice_RNA_Markers[[a]], Rice_Multiome_Markers[[b]])  }}
call_marker_overlap

# Make Heatmap
call_marker_overlap <- ifelse(call_marker_overlap > 50 , 50, call_marker_overlap) # Set upper limit identity
rownames(call_marker_overlap)<-cell_types
colnames(call_marker_overlap)<-cell_types
plot_matrix<-melt(call_marker_overlap)
colnames(plot_matrix)<-c("RNA","Multiome","log_10_p")
g1<- qplot(RNA, Multiome, fill=log_10_p, data=plot_matrix, geom='tile') + theme(axis.text = element_text(angle = 90, vjust = 0.5, hjust=1))
g1 + scale_fill_gradient(low="grey100", high="red4") + labs(x = 'RNA', y = 'Multiome') + theme_bw()
