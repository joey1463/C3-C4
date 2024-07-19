# This script will calculate which rice bundle sheath specific cell type markers are lost in the sorghum bundle sheath, and vice versa

setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_species_combined.RData")
DefaultAssay(object = species_combined) <- "RNA"

# Labeling One UMAP for Both Species To Generate Normalized Heatmap 
load("L1_rice_combined_CI_labelled.RData")

rice_combined_CI_labelled<-rice_combined_CI_labelled[,substr(colnames(rice_combined_CI_labelled),6,7) == '48']
rice_combined_CI_labelled<-RenameCells(rice_combined_CI_labelled, new.names = str_sub(colnames(rice_combined_CI_labelled), end=-4))
rice_labels<-Idents(rice_combined_CI_labelled)
remove(rice_combined_CI_labelled) # remove, otherwise memory issues

load("L1_sorghum_combined_CI_labelled.RData")
sorghum_combined_CI_labelled<-sorghum_combined_CI_labelled[,substr(colnames(sorghum_combined_CI_labelled),6,7) == '48']
sorghum_combined_CI_labelled<-RenameCells(sorghum_combined_CI_labelled, new.names = str_sub(colnames(sorghum_combined_CI_labelled), end=-4))
sorghum_labels<-Idents(sorghum_combined_CI_labelled)
remove(sorghum_combined_CI_labelled) # remove, otherwise memory issues

# Apply Labels to Existing Data
current_names<-colnames(species_combined)
new_names<-rep("unknown",length(colnames(species_combined)))
for (i in 1:length(current_names)) {if ( current_names[i] %in% names(rice_labels) == TRUE) { new_names[i]<- as.character(subset(rice_labels, names(rice_labels) %in% current_names[i])) }}
for (i in 1:length(current_names)) {if ( current_names[i] %in% names(sorghum_labels) == TRUE) { new_names[i]<- as.character(subset(sorghum_labels, names(sorghum_labels) %in% current_names[i])) }}
species_combined$designation<-new_names
Idents(species_combined)<-paste(str_sub(colnames(species_combined), start=1, end=4),species_combined$designation,sep="_")
DimPlot(object = species_combined, reduction = "umap", label=FALSE, raster = TRUE, split.by = 'Time')

DefaultAssay(object = species_combined) <- "RNA"
species_combined<-species_combined[,WhichCells(species_combined, idents = as.character(unique(Idents(species_combined)))[c(3:7,9,11:16)])]
cross_species_cell_type_average<-as.matrix((x = AverageExpression(object = species_combined, verbose = FALSE)$RNA))

# Save Data
#setwd("~/Desktop/Salk/Project - Light/R files")
#save(cross_species_cell_type_average, file = "L1_cross_species_cell_type_average.RData")
load("L1_cross_species_cell_type_average.RData")
cross_species_cell_type_average<-as.data.frame(cross_species_cell_type_average)

# Finding Gain and Loss Sorghum BS Markers
DefaultAssay(object = species_combined) <- "integrated" # if doesn't work, reload original data.
species_combined <- FindClusters(species_combined, resolution = 0.8)
DimPlot(object = species_combined, reduction = "umap", label=TRUE, raster = FALSE, split.by = 'Time')

# Barplot of Species Representation
raw_table<-table(as.character(species_combined@active.ident), species_combined$Time)
normalized_table<-matrix(NA,nrow=nrow(raw_table),ncol=ncol(raw_table))
rownames(normalized_table)<-rownames(raw_table)
colnames(normalized_table)<-colnames(raw_table)
for (i in 1:ncol(raw_table)) {normalized_table[,i]<-as.numeric(raw_table[,i]/sum(raw_table[,i]))}
markers_overlap<-as.data.frame(normalized_table[,2]/normalized_table[,1])
markers_overlap[,2]<-rownames(markers_overlap)
colnames(markers_overlap)<-c("percentage","cluster")
markers_overlap$cluster<-as.numeric(markers_overlap$cluster)
p1<- ggplot(data=markers_overlap, aes( x= cluster, y=percentage)) + geom_bar(stat="identity", fill="grey70") + theme_bw()
p2<- p1 + labs(x = "") + labs(y = "") + theme(axis.text=element_text(size=12),axis.title=element_text(size=12)) + scale_x_continuous(markers_overlap, breaks = markers_overlap$cluster)

# Find Markers of Sorghum BS (Cluster 10)
DefaultAssay(object = species_combined) <- "RNA"
species_combined_sorghum<-species_combined[,species_combined$Time == "Sorghum_CI"] # works better when consider sorghum alone
BS_markers_sorghum<-FindMarkers(species_combined_sorghum, ident.1 = 10, verbose = FALSE, only.pos = TRUE)
BS_markers_sorghum<-subset(BS_markers_sorghum, BS_markers_sorghum$p_val_adj<0.05)

# Find Markers of Rice BS (Cluster 17)
species_combined_rice<-species_combined[,species_combined$Time == "Rice_CI"] # works better when consider rice alone
BS_markers_rice<-FindMarkers(species_combined_rice, ident.1 = 18, verbose = FALSE, only.pos = TRUE) 
BS_markers_rice<-subset(BS_markers_rice, BS_markers_rice$p_val_adj<0.05)

# Genes 'Gained' by Sorghum BS
cell_type_average_subset<-subset(cross_species_cell_type_average, rownames(cross_species_cell_type_average) %in% rownames(BS_markers_sorghum))
cell_type_average_subset$maximum <- apply(cell_type_average_subset, 1, max, na.rm=TRUE)
cell_type_average_subset<-as.matrix(subset(cell_type_average_subset, cell_type_average_subset$maximum == cell_type_average_subset$Sorg_bundle_sheath)[1:12]) # important to include
cell_type_average_subset<-cell_type_average_subset[,order(colnames(cell_type_average_subset))]
heatmap.2(cell_type_average_subset, Rowv=as.dendrogram(hclust(as.dist(1-cor(t(cell_type_average_subset), method="pearson")))), Colv=NA, col=plasma(50, direction = 1), trace="none", scale="row",margins = c(8, 16))

# Plot Top 50 'Gained' Markers
cell_type_average_subset_top<-subset(cell_type_average_subset, rownames(cell_type_average_subset) %in% rownames(BS_markers_sorghum)[1:80]) #  TOP 50 FOR FIGURE
heatmap.2(cell_type_average_subset_top, Rowv=as.dendrogram(hclust(as.dist(1-cor(t(cell_type_average_subset_top), method="pearson")))), Colv=NA, col=colorRampPalette(c("blue", "black", "yellow"))(n = 100), trace="none", scale="row",margins = c(8, 16))
#setwd("~/Desktop/")
#write.table(cell_type_average_subset, file = "BS_markers_sorghum_gain.txt",quote=F,sep="\t")


# Genes 'Lost' by Sorghum BS
cell_type_average_subset<-subset(cross_species_cell_type_average, rownames(cross_species_cell_type_average) %in% rownames(BS_markers_rice))
cell_type_average_subset$maximum <- apply(cell_type_average_subset, 1, max, na.rm=TRUE)
cell_type_average_subset<-as.matrix(subset(cell_type_average_subset, cell_type_average_subset$maximum == cell_type_average_subset$Rice_bundle_sheath)[1:12])
cell_type_average_subset<-cell_type_average_subset[,order(colnames(cell_type_average_subset))]
heatmap.2(cell_type_average_subset, Rowv=as.dendrogram(hclust(as.dist(1-cor(t(cell_type_average_subset), method="pearson")))), Colv=NA, col=plasma(50, direction = 1), trace="none", scale="row",margins = c(8, 16))

# Plot Top 50 'Lost' Markers
cell_type_average_subset_top<-subset(cell_type_average_subset, rownames(cell_type_average_subset) %in% rownames(BS_markers_rice)[1:80]) #  TOP 50 FOR FIGURE
heatmap.2(cell_type_average_subset_top, Rowv=as.dendrogram(hclust(as.dist(1-cor(t(cell_type_average_subset_top), method="pearson")))), Colv=NA, col=colorRampPalette(c("blue", "black", "yellow"))(n = 100), trace="none", scale="row",margins = c(8, 16))
#setwd("~/Desktop/")
#write.table(cell_type_average_subset, file = "BS_markers_sorghum_loss.txt",quote=F,sep="\t")
