# This script will compute cell type specific gene expression markers within the Sorghum Mutiome data, and compare resulting markers with those from the RNA Atlas data

setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_sorghum_multiome_RNA_combined_with_ATAC.RData")
DefaultAssay(sorghum_multiome_RNA_combined) <- "RNA"

cell_types<-c("mesophyll", "bundle_sheath","guard","epidermis","phloem","xylem")

# Find Markers of Multiome
Sorghum_Multiome_Markers<-list(NA)
for (i in 1:6){
  print(i)
  markers<-FindMarkers(sorghum_multiome_RNA_combined, ident.1 =cell_types[i], verbose = FALSE, only.pos = TRUE)
  markers[,6]<-markers$pct.1/markers$pct.2
  markers<-markers[order(markers$p_val_adj,decreasing = FALSE),]
  markers<-subset(markers, markers$p_val_adj<0.01)
  Sorghum_Multiome_Markers[[i]]<-markers}

# Load Markers of RNA Atlas 
load("L3_Orthology_sorghum_markers.RData")
Sorghum_RNA_Markers<-sorghum_markers

# Call Overlap
call_overlap_fisher<-function(sorghum_RNA, sorghum_Multiome){
  sorghum_RNA_list<-subset(sorghum_RNA, sorghum_RNA$p_val_adj<0.01 & sorghum_RNA$V6 > 1.8)
  sorghum_Multiome_list<-subset(sorghum_Multiome, sorghum_Multiome$p_val_adj<0.01 & sorghum_Multiome$V6 > 1.8)

  intersect<-sum(as.numeric((rownames(sorghum_RNA_list) %in% substr(rownames(sorghum_Multiome_list),1,16))))
  sorghum_RNA_length<-length(sorghum_RNA_list[,1])
  sorghum_Multiome_length<-length(sorghum_Multiome_list[,1])
  matrix_to_test<-matrix(c(intersect,(sorghum_RNA_length-intersect),(sorghum_Multiome_length-intersect),30467),2,2) # using all expressed genes in RNA object
  p_value<-as.numeric(fisher.test(matrix_to_test, alternative='greater')[1])
  -log(p_value,10)}

call_marker_overlap<-matrix(NA, nrow=6, ncol=6)
for (a in 1:6) { 
  for (b in 1:6) {
     call_marker_overlap[a,b]<- call_overlap_fisher(Sorghum_RNA_Markers[[a]], Sorghum_Multiome_Markers[[b]])  }}
call_marker_overlap
call_marker_overlap[is.infinite(call_marker_overlap)] <- 300

# Make Heatmap
#call_marker_overlap <- ifelse(call_marker_overlap > 100 , 100, call_marker_overlap) # Set upper limit identity
rownames(call_marker_overlap)<-cell_types
colnames(call_marker_overlap)<-cell_types
plot_matrix<-melt(call_marker_overlap)
colnames(plot_matrix)<-c("RNA","Multiome","log_10_p")
g1<- qplot(RNA, Multiome, fill=log_10_p, data=plot_matrix, geom='tile') + theme(axis.text = element_text(angle = 90, vjust = 0.5, hjust=1))
g1 + scale_fill_gradient(low="grey100", high="red4") + labs(x = 'Multiome', y = 'RNA') + theme_bw()
