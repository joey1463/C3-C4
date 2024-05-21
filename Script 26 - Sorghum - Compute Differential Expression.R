# This script will calcualte differentially expressed genes in response to light for the 6 major cell types. 
# For each cell type, it will save both the list of differentially expressed genes, as well as 10X nuclei object. 

library(Seurat)
library(DESeq2)
library(EDASeq)

subcluster_and_DE<-function(subcluster_sorghum_combined, target_subclusters, chosen_cluster){
  
  subcluster_sorghum_combined<-subcluster_sorghum_combined[,WhichCells(subcluster_sorghum_combined, idents = target_subclusters)] 
  subcluster_sorghum_combined_CI<-subcluster_sorghum_combined[,subcluster_sorghum_combined$assay_type %in% c("CI")]
  subcluster_sorghum_combined_10X<-subcluster_sorghum_combined[,subcluster_sorghum_combined$assay_type %in% c("10X")]
  
  Idents(subcluster_sorghum_combined_CI)<-subcluster_sorghum_combined_CI$light_time
  Idents(subcluster_sorghum_combined_10X)<-subcluster_sorghum_combined_10X$light_time

  DefaultAssay(object = subcluster_sorghum_combined_CI) <- "RNA"
  DefaultAssay(object = subcluster_sorghum_combined_10X) <- "RNA"
  
  # Pairwise DE
  response_cluster<-FindMarkers(subcluster_sorghum_combined_10X, ident.1 = 'Sorghum_12' , ident.2 = 'Sorghum_0', verbose = FALSE)
  Pairwise_Significant<-subset(response_cluster, response_cluster$p_val_adj<0.01)
  
  # Linear Model on Pseudo-Bulk
  counts_10X<-as.matrix(GetAssayData(object = subcluster_sorghum_combined_10X, slot = "counts"))
  counts_CI<-as.matrix(GetAssayData(object = subcluster_sorghum_combined_CI, slot = "counts"))

  time_points<-c(0,0.5,1,2,4,6,8,12)
  
  pseudobulk_10X<-AggregateExpression(object = subcluster_sorghum_combined_10X, slot = "counts", group.by = "light_time")
  pseudobulk_10X<-as.data.frame(pseudobulk_10X[[1]])
  colnames(pseudobulk_10X)<-paste(time_points,"10X", sep="_")
  pseudobulk_10X<-pseudobulk_10X[rowSums(pseudobulk_10X == 0) < 5, ] # restricting drop out, setting on the 10X data
  
  pseudobulk_CI<-AggregateExpression(object = subcluster_sorghum_combined_CI, slot = "counts", group.by = "light_time")
  pseudobulk_CI<-as.data.frame(pseudobulk_CI[[1]])
  pseudobulk_CI<-pseudobulk_CI[,c(1:8)]
  colnames(pseudobulk_CI)<-paste(time_points,"CI", sep="_")
  pseudobulk_CI<-subset(pseudobulk_CI, rownames(pseudobulk_CI) %in% rownames(pseudobulk_10X))

  pseudo_bulk<-cbind(pseudobulk_10X, pseudobulk_CI)
  maxes<-rowSums(pseudo_bulk) 
  #hist(log(maxes,2),breaks=1000, xlim=c(1,15), ylim=c(0,100))
  pseudo_bulk<-pseudo_bulk[maxes>128,]

  pseudo_bulk_pairwise<-subset(pseudo_bulk, rownames(pseudo_bulk) %in% rownames(Pairwise_Significant))
  PC_call<-betweenLaneNormalization(as.matrix(pseudo_bulk_pairwise+1), which="full",offset=FALSE, round=TRUE)
  PCs<-((prcomp(t(log(PC_call,2)))$x)[,1:3]) # First PC is sequencing depth between techniques
  #plot(PCs[,1],PCs[,2])
  #text(PCs[,1],PCs[,2], labels=rownames(PCs), cex= 0.7)

  # Linear
  Conditions<-data.frame(time= c(time_points, time_points), technique=c(rep("10X",8),rep("CI",8)),row.names = colnames(pseudo_bulk))
  DEmat<-DESeq(DESeqDataSetFromMatrix(countData = pseudo_bulk,colData = Conditions, design = ~ time + technique)) 
  pseudo_bulk_normalized <- counts(estimateSizeFactors(DEmat), normalized=TRUE)
  Significance_Values<-results(DESeq(DEmat, test='LRT', reduced = ~ technique))[,c(2,6)]
  Linear_Significant_Linear<-subset(Significance_Values, Significance_Values$padj<0.01)
  
  # First PC
  Conditions<-data.frame(time= PCs[,2], technique=c(rep("10X",8),rep("CI",8)),row.names = colnames(pseudo_bulk))
  DEmat<-DESeq(DESeqDataSetFromMatrix(countData = pseudo_bulk,colData = Conditions, design = ~ time + technique)) 
  Significance_Values<-results(DESeq(DEmat, test='LRT', reduced = ~ technique))[,c(2,6)]
  Linear_Significant_PC1<-subset(Significance_Values, Significance_Values$padj<0.01)

  # Second PC
  Conditions<-data.frame(time= PCs[,3], technique=c(rep("10X",8),rep("CI",8)),row.names = colnames(pseudo_bulk))
  DEmat<-DESeq(DESeqDataSetFromMatrix(countData = pseudo_bulk,colData = Conditions, design = ~ time + technique)) 
  Significance_Values<-results(DESeq(DEmat, test='LRT', reduced = ~ technique))[,c(2,6)]
  Linear_Significant_PC2<-subset(Significance_Values, Significance_Values$padj<0.01)
  
  Linear_Significant<-rbind(Linear_Significant_Linear,Linear_Significant_PC1,Linear_Significant_PC2)
  Linear_Significant<-Linear_Significant[unique(rownames(Linear_Significant)),]
  Linear_Significant<-Linear_Significant[order(Linear_Significant$padj, decreasing = FALSE),]
  #
  Pairwise_Excluded<-Pairwise_Significant[!(rownames(Pairwise_Significant) %in% rownames(Linear_Significant)),]
  Pairwise_To_Add<-subset(Significance_Values, rownames(Significance_Values) %in% rownames(Pairwise_Excluded))
  DE_Significant<-rbind(Linear_Significant,Pairwise_To_Add)
  
  # Generating FC values
  Conditions<-data.frame(time= time_points, technique=c(rep("10X",8)),row.names = colnames(pseudo_bulk)[1:8])
  DEmat<-DESeq(DESeqDataSetFromMatrix(countData = pseudo_bulk[,1:8],colData = Conditions, design = ~ time)) 
  FC_Values<-results(DESeq(DEmat, test='LRT', reduced = ~ 1))[,c(2,6)]
  DE_Significant_FC<-merge(as.data.frame(DE_Significant), as.data.frame(FC_Values), by='row.names')
  rownames(DE_Significant_FC)<-DE_Significant_FC[,1]
  DE_Significant_FC<-DE_Significant_FC[,c(4,3)]
  colnames(DE_Significant_FC)<-c("log2FoldChange","padj")
  DE_Significant<-DE_Significant_FC
  
  # Saving Output
  save(subcluster_sorghum_combined_10X, file = paste(paste(chosen_cluster, sep="", collapse=""),"10X.RData",sep='_'))
  save(DE_Significant, file = paste(paste(chosen_cluster, sep="", collapse=""),"DE_Significant.RData",sep='_'))}

load("L1_sorghum_combined.RData")

chosen_cluster<-'L3_sorghum_guard'
target_subclusters<-c(14) 
subcluster_and_DE(sorghum_combined, target_subclusters, chosen_cluster)

chosen_cluster<-'L3_sorghum_epidermis'
target_subclusters<-c(6,8,17) 
subcluster_and_DE(sorghum_combined, target_subclusters, chosen_cluster)

chosen_cluster<-'L3_sorghum_mesophyll'
target_subclusters<-c(0,1) 
subcluster_and_DE(sorghum_combined, target_subclusters, chosen_cluster)

chosen_cluster<-'L3_sorghum_bundle_sheath'
target_subclusters<-c(4,15)
subcluster_and_DE(sorghum_combined, target_subclusters, chosen_cluster)

chosen_cluster<-'L3_sorghum_xylem'
target_subclusters<-c(10) 
subcluster_and_DE(sorghum_combined, target_subclusters, chosen_cluster)

chosen_cluster<-'L3_sorghum_phloem'
target_subclusters<-c(12,13,9)
subcluster_and_DE(sorghum_combined, target_subclusters, chosen_cluster)
