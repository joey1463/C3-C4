# This script does 2 things. It will:
# 1. Compute the mesophyll specific and bundle sheath specific genes within the multiome data. 
# 2. Visualize the DOF family gene expression patterns within these 2 cell types. 

# Compute M and BS specific gene expression patterns in Rice

setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_rice_multiome_RNA_combined_with_ATAC.RData")
DefaultAssay(rice_multiome_RNA_combined) <- "RNA"
#DimPlot(object = rice_multiome_RNA_combined, reduction = "umap", label=TRUE, raster = FALSE) + NoLegend()

# Use Find Markers()
markers<-FindMarkers(rice_multiome_RNA_combined, ident.1 = "mesophyll", ident.2 = "bundle_sheath", verbose = FALSE, only.pos = FALSE)
markers[,6]<-markers$pct.1/markers$pct.2
markers<-markers[order(markers$p_val_adj,decreasing = FALSE),]

markers<-subset(markers, markers$p_val_adj<0.01)
rice_M_specific_genes<-subset(markers, markers$V6>1.25)
rice_BS_specific_genes<-subset(markers, markers$V6<(1/1.25))

setwd("~/Desktop/Salk/Project - Light/R files")
save(rice_M_specific_genes, file="rice_M_specific_genes.RData")
save(rice_BS_specific_genes, file="rice_BS_specific_genes.RData")


# Compute DOF Gene Expression Patterns
setwd("~/Desktop/Salk/Project - Light/R files")
Rice_TFs<-read.table("Osj_TF_list.txt",header=TRUE)
Rice_DOFs<-unique(subset(Rice_TFs, Rice_TFs[,2] == "Dof")[,1])
load("L1_rice_multiome_Average_Expression.RData")
DOF_names<-read.table("Rice DOF Names.txt",header=TRUE, sep="\t")
rownames(DOF_names)<-DOF_names[,2]

# Plot Average Expression
Average_Expression<-Average_Expression[,c(1,7)] # need this from first principles
plot_matrix<-subset(Average_Expression, rownames(Average_Expression) %in% Rice_DOFs)
plot_matrix<-merge(plot_matrix, DOF_names, by='row.names')
rownames(plot_matrix)<-plot_matrix$Name
plot_matrix<-plot_matrix[,2:3]
plot_matrix[15,]<-c(0,100) # fake to get scale right 
plot_matrix<-as.matrix(plot_matrix)
plot_matrix<-t(apply(plot_matrix, 1, function(x)(x/(max(x)))))
plot_matrix<-plot_matrix[order(plot_matrix[,1]/plot_matrix[,2]),]
heatmap.2(plot_matrix, Rowv=NA, Colv=NA,dendrogram="row", scale=NULL, col=colorRampPalette(c("blue","black", "yellow"))(n = 100), tracecol=NA,margins = c(8, 16),cexRow=2) # ! Make Sure You Remove Fake Last Row !
