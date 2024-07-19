# This script will specifically insepct the differential parititoning patterns of mesophyll vs. bundle sheath, 
# and group gene pairs based on whether they are consistently or differentially partitioned across species. 


# Load Orthology
setwd("~/Desktop/Salk/Project - Light/R files/")
orthology<-read.table('Rice_Sorghum_Orthology.txt',header=TRUE, sep="\t")

# Read in just M and BS data
setwd("~/Desktop/Salk/Project - Light/R files")
load("L3_Sorghum_Pair_Permutations.RData")
load("L3_Rice_Pair_Permutations.RData")

partition_1<-Sorghum_Pair_Permutations[[1]]
partition_2<-Rice_Pair_Permutations[[1]]

# Sorghum
partitioning_1<-cell_type_subset_sorghum[,str_detect(colnames(cell_type_subset_sorghum), "mesophyll")]
partitioning_2<-cell_type_subset_sorghum[,str_detect(colnames(cell_type_subset_sorghum), "bundle_sheath")]
partitioning_sorghum<-cbind(partitioning_1,partitioning_2)
cell_type_subset_sig_sorghum<-subset(partitioning_sorghum, rownames(partitioning_sorghum) %in% partition_1)

# Rice
rownames(cell_type_subset_rice)<-str_replace(rownames(cell_type_subset_rice),"-","_")
partitioning_1<-cell_type_subset_rice[,str_detect(colnames(cell_type_subset_rice), "mesophyll")]
partitioning_2<-cell_type_subset_rice[,str_detect(colnames(cell_type_subset_rice), "bundle_sheath")]
partitioning_rice<-cbind(partitioning_1,partitioning_2)
cell_type_subset_sig_rice<-subset(partitioning_rice, rownames(partitioning_rice) %in% str_replace(partition_2,"-","_"))

# Rice orthology compute
rice_orthology<-matrix(NA,nrow=length(cell_type_subset_sig_rice[,1]),ncol=3)
for (i in 1:length(cell_type_subset_sig_rice[,1])) {rice_orthology[i,]<-as.character(orthology[grep(rownames(cell_type_subset_sig_rice)[i], orthology$Rice),])}
rice_orthology<-cbind(rice_orthology, cell_type_subset_sig_rice)
colnames(rice_orthology)[1:3]<-colnames(orthology)
rice_orthology<-subset(rice_orthology, substr(rice_orthology$Orthogroup,1,1) =="O" )
rice_orthology$Rice<-rownames(rice_orthology)

# Sorghum orthology compute
sorghum_orthology<-matrix(NA,nrow=length(cell_type_subset_sig_sorghum[,1]),ncol=3)
for (i in 1:length(cell_type_subset_sig_sorghum[,1])) {sorghum_orthology[i,]<-as.character(orthology[grep(rownames(cell_type_subset_sig_sorghum)[i], orthology$Sorghum),])}
sorghum_orthology<-cbind(sorghum_orthology, cell_type_subset_sig_sorghum)
colnames(sorghum_orthology)[1:3]<-colnames(orthology)
sorghum_orthology<-subset(sorghum_orthology, substr(sorghum_orthology$Orthogroup,1,1) =="O" )
sorghum_orthology$Sorghum<-rownames(sorghum_orthology)

# Merge by orthogroups
category_genes_1<-merge(x = rice_orthology, y = sorghum_orthology, by='Orthogroup') 
category_genes_1$Sorghum.x<-category_genes_1$Sorghum.y
category_genes_1<-category_genes_1[,-c(20:21)]

# Normalizing gene expression across species
category_genes_rice<-category_genes_1[,c(4:19)]
category_genes_rice_normalized<-t(apply(category_genes_rice, 1, function(x)(x-min(x))/(max(x)-min(x))))
category_genes_sorghum<-category_genes_1[,c(20:35)]
category_genes_sorghum_normalized<-t(apply(category_genes_sorghum, 1, function(x)(x-min(x))/(max(x)-min(x))))
cluster_output<-as.data.frame(cbind(category_genes_rice_normalized, category_genes_sorghum_normalized))
rownames(cluster_output)<-paste(category_genes_1$Rice.x, category_genes_1$Sorghum.x)
# heatmap of everything together deleted
cluster_output[,33]<-category_genes_1$Orthogroup

# Heatmaps of each class
cluster_output_simple<-matrix(NA, nrow=length(cluster_output[,1]), ncol=4)
rownames(cluster_output_simple)<-rownames(cluster_output)
for (i in 1:length(cluster_output[,1])) {cluster_output_simple[i,1]<- mean(as.numeric(cluster_output[i,1:8]), na.rm = TRUE)
                                        cluster_output_simple[i,2]<-  mean(as.numeric(cluster_output[i,9:16]), na.rm = TRUE) 
                                        cluster_output_simple[i,3]<-  mean(as.numeric(cluster_output[i,17:24]), na.rm = TRUE) 
                                        cluster_output_simple[i,4]<-  mean(as.numeric(cluster_output[i,25:32]), na.rm = TRUE) }
cluster_output_simple<- cluster_output_simple[complete.cases(cluster_output_simple),]
colnames(cluster_output_simple)<-c("rice_cell_type_1","rice_cell_type_2","sorghum_cell_type_1","sorghum_cell_type_2")

focus_partition_MtoBS<-subset(cluster_output_simple, cluster_output_simple[,1]>cluster_output_simple[,2] & cluster_output_simple[,3]<cluster_output_simple[,4] )
focus_partition_MtoBS<-subset(cluster_output, rownames(cluster_output) %in% rownames(focus_partition_MtoBS))
heatmap.2(as.matrix(focus_partition_MtoBS[,1:32]), Rowv=as.dendrogram(hclust(as.dist(1-cor(t(as.matrix(focus_partition_MtoBS[,1:32])), method="pearson")))), Colv=NA, col=colorRampPalette(c("blue", "black", "yellow"))(n = 100), trace="none", scale="row",margins = c(8, 16))

focus_partition_BStoM<-subset(cluster_output_simple, cluster_output_simple[,1]<cluster_output_simple[,2] & cluster_output_simple[,3]>cluster_output_simple[,4] )
focus_partition_BStoM<-subset(cluster_output, rownames(cluster_output) %in% rownames(focus_partition_BStoM))
heatmap.2(as.matrix(focus_partition_BStoM[,1:32]), Rowv=as.dendrogram(hclust(as.dist(1-cor(t(as.matrix(focus_partition_BStoM[,1:32])), method="pearson")))), Colv=NA, col=colorRampPalette(c("blue", "black", "yellow"))(n = 100), trace="none", scale="row",margins = c(8, 16))

focus_partition_MtoM<-subset(cluster_output_simple, cluster_output_simple[,1]>cluster_output_simple[,2] & cluster_output_simple[,3]>cluster_output_simple[,4] )
focus_partition_MtoM<-subset(cluster_output, rownames(cluster_output) %in% rownames(focus_partition_MtoM))
heatmap.2(as.matrix(focus_partition_MtoM[,1:32]), Rowv=as.dendrogram(hclust(as.dist(1-cor(t(as.matrix(focus_partition_MtoM[,1:32])), method="pearson")))), Colv=NA, col=colorRampPalette(c("blue", "black", "yellow"))(n = 100), trace="none", scale="row",margins = c(8, 16))

focus_partition_BStoBS<-subset(cluster_output_simple, cluster_output_simple[,1]<cluster_output_simple[,2] & cluster_output_simple[,3]<cluster_output_simple[,4] )
focus_partition_BStoBS<-subset(cluster_output, rownames(cluster_output) %in% rownames(focus_partition_BStoBS))
heatmap.2(as.matrix(focus_partition_BStoBS[,1:32]), Rowv=as.dendrogram(hclust(as.dist(1-cor(t(as.matrix(focus_partition_BStoBS[,1:32])), method="pearson")))), Colv=NA, col=colorRampPalette(c("blue", "black", "yellow"))(n = 100), trace="none", scale="row",margins = c(8, 16))

# Consistent partitioning - overlap with photosynthesis gene annotation
focus_partition_MtoM$category<-ifelse(substr(rownames(focus_partition_MtoM),1,14) %in% str_replace(rownames(rice_lightdark_extensive),"-","_"), rice_lightdark_extensive$process, "NA")
focus_partition_BStoBS$category<-ifelse(substr(rownames(focus_partition_BStoBS),1,14) %in% str_replace(rownames(rice_lightdark_extensive),"-","_"), rice_lightdark_extensive$process, "NA")

# Differential partitioning - overlap with photosynthesis gene annotation
focus_partition_MtoBS$category<-ifelse(substr(rownames(focus_partition_MtoBS),1,14) %in% str_replace(rownames(rice_lightdark_extensive),"-","_"), rice_lightdark_extensive$process, "NA")
focus_partition_BStoM$category<-ifelse(substr(rownames(focus_partition_BStoM),1,14) %in% str_replace(rownames(rice_lightdark_extensive),"-","_"), rice_lightdark_extensive$process, "NA")



## Supplementary Data on Orthology Complexity 
represented_orthogroups<-unique(cluster_output[,33])
pairs_count<-NA
pair_combinations<-matrix(NA,ncol=4,nrow=length(represented_orthogroups))

for (i in 1:length(represented_orthogroups)) {
pairs_count[i]<-length(rownames(subset(cluster_output, cluster_output[,33] == represented_orthogroups[i])))
pairs_themselves<-rownames(subset(cluster_output, cluster_output[,33] == represented_orthogroups[i]))
pair_combinations[i,1]<-sum(as.numeric(pairs_themselves %in% rownames(focus_partition_MtoM)))
pair_combinations[i,2]<-sum(as.numeric(pairs_themselves %in% rownames(focus_partition_BStoBS)))
pair_combinations[i,3]<-sum(as.numeric(pairs_themselves %in% rownames(focus_partition_BStoM)))
pair_combinations[i,4]<-sum(as.numeric(pairs_themselves %in% rownames(focus_partition_MtoBS)))}
complexity_orthology<-as.data.frame(cbind(pairs_count,pair_combinations))
rownames(complexity_orthology)<-represented_orthogroups
