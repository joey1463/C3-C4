# This script will analyze the mesophyll and bundle sheath gene expression patterns within Multiome datsets, and detect
# which genes are differentially partititoned, and which are consistently partitioned. 

# Load Orthology
setwd("~/Desktop/Salk/Project - Light/R files")
orthology<-read.table('Rice_Sorghum_Orthology.txt',header=TRUE, sep="\t")

### Differentially Partitioned BS Genes

# Load Rice-M Specific Genes
load("rice_M_specific_genes.RData")
rownames(rice_M_specific_genes)<-str_replace(rownames(rice_M_specific_genes),"-","_")
load("L1_rice_multiome_Average_Expression.RData")
rownames(Average_Expression)<-str_replace(rownames(Average_Expression),"-","_")
Average_Expression_Rice<-Average_Expression[,c(1,7)]
rice_M_specific_genes<-merge(rice_M_specific_genes,  Average_Expression_Rice, by='row.names')
rownames(rice_M_specific_genes)<-rice_M_specific_genes[,1]
rice_M_specific_genes<-rice_M_specific_genes[,c(3,6,8,9)]

# Load Sorghum-BS Specific Genes
load("sorghum_BS_specific_genes.RData")
rownames(sorghum_BS_specific_genes)<-substr(rownames(sorghum_BS_specific_genes),1,16)
load("L1_sorghum_multiome_Average_Expression.RData")
rownames(Average_Expression)<-substr(rownames(Average_Expression), 1,16)
Average_Expression_Sorghum<-Average_Expression[,c(2,5)]
sorghum_BS_specific_genes<-merge(sorghum_BS_specific_genes,  Average_Expression_Sorghum, by='row.names')
rownames(sorghum_BS_specific_genes)<-sorghum_BS_specific_genes[,1]
sorghum_BS_specific_genes<-sorghum_BS_specific_genes[,c(3,6,8,9)]

# Compute Rice Orthogroups
rice_orthology_M<-matrix(NA,nrow=length(rice_M_specific_genes[,1]),ncol=3)
for (i in 1:length(rice_M_specific_genes[,1])) {rice_orthology_M[i,]<-as.character(orthology[grep(rownames(rice_M_specific_genes)[i], orthology$Rice),])}
rice_orthology_M<-cbind(rice_orthology_M, rice_M_specific_genes)
colnames(rice_orthology_M)[1:3]<-colnames(orthology)
rice_orthology_M<-subset(rice_orthology_M, substr(rice_orthology_M$Orthogroup,1,1) =="O" )
rice_orthology_M$Rice<-rownames(rice_orthology_M)

# Compute Sorghum Orthogroups
sorghum_orthology_BS<-matrix(NA,nrow=length(sorghum_BS_specific_genes[,1]),ncol=3)
for (i in 1:length(sorghum_BS_specific_genes[,1])) {sorghum_orthology_BS[i,]<-as.character(orthology[grep(rownames(sorghum_BS_specific_genes)[i], orthology$Sorghum),])}
sorghum_orthology_BS<-cbind(sorghum_orthology_BS, sorghum_BS_specific_genes)
colnames(sorghum_orthology_BS)[1:3]<-colnames(orthology)
sorghum_orthology_BS<-subset(sorghum_orthology_BS, substr(sorghum_orthology_BS$Orthogroup,1,1) =="O" )
sorghum_orthology_BS$Sorghum<-rownames(sorghum_orthology_BS)

# Merge By Orthogroups
Swap_Genes<-merge(x = rice_orthology_M, y = sorghum_orthology_BS, by='Orthogroup') 
Swap_Genes$Sorghum.x<-Swap_Genes$Sorghum.y
Swap_Genes$Rice.y<-Swap_Genes$Rice.x

setwd("~/Desktop/Salk/Project - Light/R files")
save(Swap_Genes, file="Swap_Genes.RData")
#write.table(Swap_Genes,file="Swap_Genes.txt",quote=F,sep="\t")


### Consistently BS Partitioned Genes 

# Load Rice-BS Specific Genes
setwd("~/Desktop/Salk/Project - Light/R files")
load("rice_BS_specific_genes.RData")
rownames(rice_BS_specific_genes)<-str_replace(rownames(rice_BS_specific_genes),"-","_")

# Compute Rice Orthogroups
rice_orthology_BS<-matrix(NA,nrow=length(rice_BS_specific_genes[,1]),ncol=3)
for (i in 1:length(rice_BS_specific_genes[,1])) {rice_orthology_BS[i,]<-as.character(orthology[grep(rownames(rice_BS_specific_genes)[i], orthology$Rice),])}
rice_orthology_BS<-cbind(rice_orthology_BS, rice_BS_specific_genes)
colnames(rice_orthology_BS)[1:3]<-colnames(orthology)
rice_orthology_BS<-subset(rice_orthology_BS, substr(rice_orthology_BS$Orthogroup,1,1) =="O" )
rice_orthology_BS$Rice<-rownames(rice_orthology_BS)

# Merge By Orthogroups
Consistent_Genes<-merge(x = rice_orthology_BS, y = sorghum_orthology_BS, by='Orthogroup') 
Consistent_Genes$Sorghum.x<-Consistent_Genes$Sorghum.y
save(Consistent_Genes, file="Consistent_Genes.RData")

# Save Data
Swap_Genes_1<-cbind(Swap_Genes[,1:3],"differential")
colnames(Swap_Genes_1)<-c("Orthogroup","Sorghum","Rice","Pattern")
Consistent_Genes_1<-cbind(Consistent_Genes[,1:3],"consistent")
colnames(Consistent_Genes_1)<-c("Orthogroup","Sorghum","Rice","Pattern")
Target_Genes<-rbind(Swap_Genes_1, Consistent_Genes_1)
save(Target_Genes, file="Target_Genes.RData")

# Plotting
load("Swap_Genes.RData")
Swap_Genes_Rice<-Swap_Genes[!duplicated(Swap_Genes$Rice.x), ]
plot_matrix<-t(apply(as.matrix(Swap_Genes_Rice[,c(6,7)]), 1, function(x)(x/(max(x)))))
heatmap.2(plot_matrix, Rowv=NA, Colv=NA,dendrogram=NULL, scale=NULL, col=colorRampPalette(c("blue","black", "yellow"))(n = 100), tracecol=NA,margins = c(8, 16))

Swap_Genes_Sorghum<-Swap_Genes[!duplicated(Swap_Genes$Sorghum.x), ]
plot_matrix<-t(apply(as.matrix(Swap_Genes_Sorghum[,c(12,13)]), 1, function(x)(x/(max(x)))))
heatmap.2(plot_matrix, Rowv=NA, Colv=NA,dendrogram=NULL, scale=NULL, col=colorRampPalette(c("blue","black", "yellow"))(n = 100), tracecol=NA,margins = c(8, 16))

# Check Photosynthesis Genes
# Gene List Rice
setwd("~/Desktop/Salk/Project - Light/R files")
rice_lightdark_extensive<-read.table('extensive_photosynthesis_Rice.txt',header=TRUE, sep="\t")
rownames(rice_lightdark_extensive)<-str_replace(rice_lightdark_extensive[,1],"-","_")

# Rice Function
Swap_Genes_Rice[,14]<-NA
for (i in 1:length(Swap_Genes_Rice[,1])) {if (Swap_Genes_Rice$Rice.x[i] %in% rownames(rice_lightdark_extensive))
  {Swap_Genes_Rice[i,14]<- subset(rice_lightdark_extensive, rownames(rice_lightdark_extensive) == Swap_Genes_Rice$Rice.x[i])$process}}

# Gene List Sorghum
sorghum_lightdark_extensive<-read.table('extensive_photosynthesis_Sorghum.txt',header=TRUE, sep="\t")
sorghum_lightdark_extensive<-sorghum_lightdark_extensive[!duplicated(sorghum_lightdark_extensive[,1]),]
rownames(sorghum_lightdark_extensive)<-sorghum_lightdark_extensive[,1]

# Sorghum Function
Swap_Genes_Sorghum[,14]<-NA
for (i in 1:length(Swap_Genes_Sorghum[,1])) {if (Swap_Genes_Sorghum$Sorghum.x[i] %in% rownames(sorghum_lightdark_extensive))
  {Swap_Genes_Sorghum[i,14]<- subset(sorghum_lightdark_extensive, rownames(sorghum_lightdark_extensive) ==Swap_Genes_Sorghum$Sorghum.x[i])$process}}
