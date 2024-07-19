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

# Computing DOF Rice Orthologs
setwd("~/Desktop/Salk/Project - Light/R files")
DOF_names<-read.table("Rice DOF Names.txt",header=TRUE, sep="\t")
rownames(DOF_names)<-DOF_names[,2]
orthology<-read.table('Rice_Sorghum_Orthology.txt',header=TRUE, sep="\t")
Rice_TFs<-read.table("Osj_TF_list.txt",header=TRUE)
Rice_DOFs<-unique(subset(Rice_TFs, Rice_TFs[,2] == "Dof")[,1])
Rice_DOFs<-str_replace(Rice_DOFs,"-","_")
Rice_Sorghum_DOFs<-matrix(NA, nrow=length(Rice_DOFs),ncol=3)
for (i in 1:length(Rice_DOFs)) {Rice_Sorghum_DOFs[i,]<-as.character(orthology[grep(Rice_DOFs[i], orthology$Rice),])}
rownames(Rice_Sorghum_DOFs)<-Rice_DOFs
colnames(Rice_Sorghum_DOFs)[1:3]<-colnames(orthology)
Rice_Sorghum_DOFs<-as.data.frame(Rice_Sorghum_DOFs)
Rice_Sorghum_DOFs<-subset(Rice_Sorghum_DOFs, substr(Rice_Sorghum_DOFs$Orthogroup,1,1) =="O")
rownames(Rice_Sorghum_DOFs)<-str_replace(rownames(Rice_Sorghum_DOFs),"_","-")
Rice_Sorghum_DOFs<-merge(Rice_Sorghum_DOFs, DOF_names, by='row.names')

# Compute DOF Gene Expression Patterns
setwd("~/Desktop/Salk/Project - Light/R files")
Sorghum_TFs<-read.table("Sbi_TF_list.txt",header=TRUE) # TFs
Sorghum_DOFs<-unique(subset(Sorghum_TFs, Sorghum_TFs$Family == "Dof")[,1])
load("L1_sorghum_multiome_Average_Expression.RData")

# Plot Average Expression
Average_Expression<-Average_Expression[,c(2,5)] # need this from first principles
plot_matrix<-subset(Average_Expression, substr(rownames(Average_Expression),1,16) %in% Sorghum_DOFs)
for(i in c(2:14,16:21)) { plot_matrix[i,3:9]<-(Rice_Sorghum_DOFs[grep(substr(rownames(plot_matrix[i,]),1,16), Rice_Sorghum_DOFs$Sorghum),])} # hard coded
plot_matrix$Name[1]<-"no ortholog"
plot_matrix$Name[15]<-"no ortholog "
plot_matrix$Name[21]<-"OsDof5 "
rownames(plot_matrix)<-plot_matrix$Name
plot_matrix<-as.matrix(plot_matrix[,1:2])
plot_matrix<-t(apply(plot_matrix, 1, function(x)(x/(max(x)))))
plot_matrix<-plot_matrix[order(plot_matrix[,1]/plot_matrix[,2]),]
heatmap.2(plot_matrix, Rowv=NA, Colv=NA,dendrogram="row", scale=NULL, col=colorRampPalette(c("blue","black", "yellow"))(n = 100), tracecol=NA,margins = c(8, 16),cexRow=2)

setwd("~/Desktop/")
write.table(plot_matrix,file="DOF_orthologs.txt",quote=F,sep="\t")

