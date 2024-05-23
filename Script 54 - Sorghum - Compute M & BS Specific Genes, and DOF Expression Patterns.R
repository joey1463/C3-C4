# This script does 2 things. It will:
# 1. Compute the mesophyll specific and bundle sheath specific genes within the multiome data. 
# 2. Visualize the DOF family gene expression patterns within these 2 cell types. 

# Compute BS-specific gene expression patterns in Sorghum

setwd("~/Desktop/Salk/Project - Light/R files/Multiome")
load("L1_sorghum_multiome_RNA_combined_with_ATAC.RData")
DefaultAssay(sorghum_multiome_RNA_combined) <- "RNA"

# Use FindMarkers()
markers<-FindMarkers(sorghum_multiome_RNA_combined, ident.1 = "bundle_sheath", ident.2 = "mesophyll", verbose = FALSE, only.pos = FALSE)
markers[,6]<-markers$pct.1/markers$pct.2
markers<-markers[order(markers$p_val_adj,decreasing = FALSE),]
markers<-subset(markers, markers$p_val_adj<0.001)
sorghum_BS_specific_genes<-subset(markers, markers$V6>1.5)

setwd("~/Desktop/Salk/Project - Light/R files")
save(sorghum_BS_specific_genes, file="sorghum_BS_specific_genes.RData")

# Compute DOF Gene Expression Patterns

setwd("~/Desktop/Salk/Project - Light/R files")
Sorghum_TFs<-read.table("Sbi_TF_list.txt",header=TRUE) # TFs
Sorghum_DOFs<-unique(subset(Sorghum_TFs, Sorghum_TFs$Family == "Dof")[,1])
load("L1_sorghum_multiome_Average_Expression.RData")

# Plot Average Expression
Average_Expression<-Average_Expression[,c(2,5)] # need this from first principles
plot_matrix<-subset(Average_Expression, substr(rownames(Average_Expression),1,16) %in% Sorghum_DOFs)
plot_matrix<-as.matrix(plot_matrix)
plot_matrix<-t(apply(plot_matrix, 1, function(x)(x/(max(x)))))
heatmap.2(plot_matrix, Rowv=NA, Colv=NA,dendrogram="row", scale=NULL, col=colorRampPalette(c("blue","black", "yellow"))(n = 100), tracecol=NA,margins = c(8, 16))


# Computing DOF Rice Orthologs
setwd("~/Desktop/Salk/Project - Light/R files")
orthology<-read.table('Rice_Sorghum_Orthology.txt',header=TRUE, sep="\t")
Rice_DOFs<-str_replace(Rice_DOFs,"-","_")
Rice_Sorghum_DOFs<-matrix(NA, nrow=length(Rice_DOFs),ncol=3)
for (i in 1:length(Rice_DOFs)) {Rice_Sorghum_DOFs[i,]<-as.character(orthology[grep(Rice_DOFs[i], orthology$Rice),])}
rownames(Rice_Sorghum_DOFs)<-Rice_DOFs
colnames(Rice_Sorghum_DOFs)[1:3]<-colnames(orthology)
Rice_Sorghum_DOFs<-as.data.frame(Rice_Sorghum_DOFs)
Rice_Sorghum_DOFs<-subset(Rice_Sorghum_DOFs, substr(Rice_Sorghum_DOFs$Orthogroup,1,1) =="O")
#setwd("~/Desktop/")
#write.table(Rice_Sorghum_DOFs, file = 'Rice_Sorghum_DOFs_orthologs.txt',quote=F,sep="\t")
