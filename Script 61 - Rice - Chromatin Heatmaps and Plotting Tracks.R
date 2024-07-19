# This script will call cell type specific accessible chromatin from the multiome dataset. 
# It will also create a heatmap that visualizes these patterns of accessibility. 

setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_rice_multiome_RNA_combined_with_ATAC.RData")

# Calling Cell Type Peaks
Cell_Type_Peaks<-list(NA)
cell_types<-c("mesophyll","bundle_sheath","epidermis","guard","phloem","xylem")
for (i in 1:6){ 
  print(i)
  response_cluster<-FindMarkers(rice_multiome_RNA_combined, ident.1 = cell_types[i], min.pct = 0.05, test.use = 'LR', only.pos = TRUE) 
  response_cluster<-subset(response_cluster, response_cluster$p_val_adj<0.05)
  response_cluster[,6]<-response_cluster$pct.1/response_cluster$pct.2
  Cell_Type_Peaks[[i]]<-response_cluster}
save(Cell_Type_Peaks, file = "L1_rice_cell_type_peaks.RData")
load("L1_rice_cell_type_peaks.RData")

# Visualize Rubisco Subunit
rice_multiome_RNA_combined<-rice_multiome_RNA_combined[,WhichCells(rice_multiome_RNA_combined, idents = c("mesophyll","bundle_sheath","phloem","xylem","epidermis","guard"))]
Idents(rice_multiome_RNA_combined) <- factor(rice_multiome_RNA_combined@active.ident, c("mesophyll","bundle_sheath","phloem","xylem","epidermis","guard"))
p<- CoveragePlot(object = rice_multiome_RNA_combined, region = "LOC-Os12g19470", extend.upstream = 800, extend.downstream = 0)  # PS.calvin cycle.rubisco small subunit 
p & scale_fill_manual(values = c("springgreen4","steelblue1", "deepskyblue4","turquoise3","navajowhite3","tan3")) 


# Calling All Peaks
DefaultAssay(rice_multiome_RNA_combined) <- "ATAC"
load("L1_rice_multiome_Average_Peaks.RData")
Average_Peaks<-Average_Peaks[,order(colnames(Average_Peaks))]
Peaks_Close_Feature<-ClosestFeature(rice_multiome_RNA_combined, rownames(rice_multiome_RNA_combined))
Peaks_Close_Feature<-subset(Peaks_Close_Feature, Peaks_Close_Feature$distance<2000)
DE_Peaks_Average<-merge(Average_Peaks,Peaks_Close_Feature, by.x='row.names', by.y='query_region')
rownames(DE_Peaks_Average)<-DE_Peaks_Average[,1]
DE_Peaks_Average<-DE_Peaks_Average[,-1]
DE_Peaks_Average_Rice<-DE_Peaks_Average
save(DE_Peaks_Average_Rice, file = "DE_Peaks_Average_Rice_peaks.RData")
load("DE_Peaks_Average_Rice_peaks.RData")

# Looking at Peak Heatmap
all_peaks<-unique(c(rownames(Cell_Type_Peaks[[1]],Cell_Type_Peaks[[2]]),rownames(Cell_Type_Peaks[[3]]),rownames(Cell_Type_Peaks[[4]]),rownames(Cell_Type_Peaks[[5]]),rownames(Cell_Type_Peaks[[6]])))
DE_Peaks_Average_Rice<-subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% all_peaks)
to_heatmap<-DE_Peaks_Average_Rice[,c("mesophyll","bundle_sheath","phloem","xylem","epidermis","guard")]
to_heatmap$max<-apply(to_heatmap, 1, max, na.rm=TRUE)
mesophyll_heatmap<-subset(to_heatmap, to_heatmap$mesophyll == to_heatmap[,7])[,1:6]
bundle_sheath_heatmap<-subset(to_heatmap, to_heatmap$bundle_sheath == to_heatmap[,7])[,1:6]
phloem_heatmap<-subset(to_heatmap, to_heatmap$phloem == to_heatmap[,7])[,1:6]
xylem_heatmap<-subset(to_heatmap, to_heatmap$xylem == to_heatmap[,7])[,1:6]
epidermis_heatmap<-subset(to_heatmap, to_heatmap$epidermis == to_heatmap[,7])[,1:6]
guard_heatmap<-subset(to_heatmap, to_heatmap$guard == to_heatmap[,7])[,1:6]
to_heatmap_rearragned<-as.matrix(rbind(mesophyll_heatmap, bundle_sheath_heatmap, phloem_heatmap, xylem_heatmap, epidermis_heatmap, guard_heatmap))
heatmap.2(to_heatmap_rearragned, Rowv=NA, Colv=NA,dendrogram="row", scale="row",col=viridis(25, direction = 1), tracecol=NA,margins = c(8, 16))

# Output Gene Names For GO Terms
Gene_Output<-matrix(NA, nrow=length(DE_Peaks_Average_Rice[,1]), ncol=6)
colnames(Gene_Output)<-c("mesophyll","bundle_sheath","epidermis","guard","phloem","xylem")
Gene_Output[1:length(subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% rownames(Cell_Type_Peaks[[1]]))$transcript_id),1]<-subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% rownames(Cell_Type_Peaks[[1]]))$transcript_id
Gene_Output[1:length(subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% rownames(Cell_Type_Peaks[[2]]))$transcript_id),2]<-subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% rownames(Cell_Type_Peaks[[2]]))$transcript_id
Gene_Output[1:length(subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% rownames(Cell_Type_Peaks[[3]]))$transcript_id),3]<-subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% rownames(Cell_Type_Peaks[[3]]))$transcript_id
Gene_Output[1:length(subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% rownames(Cell_Type_Peaks[[4]]))$transcript_id),4]<-subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% rownames(Cell_Type_Peaks[[4]]))$transcript_id
Gene_Output[1:length(subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% rownames(Cell_Type_Peaks[[5]]))$transcript_id),5]<-subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% rownames(Cell_Type_Peaks[[5]]))$transcript_id
Gene_Output[1:length(subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% rownames(Cell_Type_Peaks[[6]]))$transcript_id),6]<-subset(DE_Peaks_Average_Rice, rownames(DE_Peaks_Average_Rice) %in% rownames(Cell_Type_Peaks[[6]]))$transcript_id

setwd("~/Desktop/")
write.table(Gene_Output, file = "Rice_Chromatin_Unique_Peaks.txt",quote=F,sep="\t")
