# This script will read in the significant light-responsive differentially expressed genes found in the mesophyll or bundle sheath cell type of rice. 
# The script will then subset these genes by those involved in photosynthesis, and then make a heatmap. 

# All DE Genes
setwd("~/Desktop/Salk/Project - Light/R Files/")
my_palette<-colorRampPalette(c("blue", "black", "yellow"))(n = 100)

chosen_cluster<-c("L3_rice_mesophyll","L3_rice_bundle_sheath")
DE_genes_rice<-list(NA)
for (i in 1:length(chosen_cluster)){
load(paste(chosen_cluster[i],"_DE_Significant.RData",sep=''))
DE_Significant<-subset(DE_Significant, DE_Significant$padj<0.01) 
DE_genes_rice[[i]]<-rownames(DE_Significant)}
all_genes<-unique(c(DE_genes_rice[[1]],DE_genes_rice[[2]]))

cell_type_subset_sig_rice<-subset(cell_type_subset, rownames(cell_type_subset) %in% all_genes)
cell_type_subset_sig_photo<-merge(rice_lightdark_extensive,cell_type_subset_sig_rice,by='row.names')
rownames(cell_type_subset_sig_photo)<-cell_type_subset_sig_photo[,1]
heatmap_cluster<-as.matrix(cell_type_subset_sig_photo[,5:52])
heatmap<-heatmap.2(heatmap_cluster, Rowv=as.dendrogram(hclust(as.dist(1-cor(t(heatmap_cluster), method="pearson")))), Colv=NA, col=my_palette, trace="none", scale="row",margins = c(8, 16))
heatmap_cluster_ordered<-as.data.frame(heatmap_cluster[rownames(heatmap_cluster)[heatmap$rowInd],])

rice_lightdark_extensive_calvin<-subset(rice_lightdark_extensive, rice_lightdark_extensive$category %in% c('calvin_cycle','photorespiration'))
rice_lightdark_extensive_light<-subset(rice_lightdark_extensive, rice_lightdark_extensive$category %in% c('light_reaction'))
heatmap_cluster_ordered[,49]<-ifelse(rownames(heatmap_cluster_ordered) %in% rice_lightdark_extensive_calvin$gene.ID, 2, 0)
heatmap_cluster_ordered[,50]<-ifelse(rownames(heatmap_cluster_ordered) %in% rice_lightdark_extensive_light$gene.ID, 1, 0)
heatmap_cluster_ordered[,51]<-heatmap_cluster_ordered[,49]+heatmap_cluster_ordered[,50]
heatmap_cluster_ordered<-subset(heatmap_cluster_ordered, heatmap_cluster_ordered$V51 >0)
heatmap_cluster_ordered<-heatmap_cluster_ordered[order(heatmap_cluster_ordered$V51),]
heatmap_cluster_ordered_rice<-heatmap_cluster_ordered
heatmap.2(as.matrix(heatmap_cluster_ordered[,1:48]), Rowv=NA, Colv=NA, col=my_palette, trace="none", scale="row", margins = c(8, 16))

# Plot Category
heatmap_category<-as.matrix(cbind(heatmap_cluster_ordered$V51,heatmap_cluster_ordered$V51))
heatmap.2(heatmap_category, Rowv=NA, Colv=NA, tracecol=NA,col=c("green","blue"),margins = c(8, 16))

# For Supplementary
heatmap_sup<-merge(rice_lightdark_extensive, heatmap_cluster_ordered_rice, by='row.names')
#setwd("~/Desktop")
#write.table(heatmap_sup,file="heatmap_sup.txt",quote=F,sep="\t")
