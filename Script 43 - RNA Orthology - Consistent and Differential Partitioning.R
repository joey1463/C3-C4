# This script will compare differentially partitioned genes acros species. 
# For example, it will compare how mesophyll vs. bundle sheath differentially partitioned genes in rice compare to those in sorghum.
# This is iterated through all possible 15 paris. 

# Load Orthology
setwd("~/Desktop/Salk/Project - Light/R files/")
orthology<-read.table('Rice_Sorghum_Orthology.txt',header=TRUE, sep="\t")

setwd("~/Desktop/Salk/Project - Light/R files/")
load("L3_Sorghum_Pair_Permutations.RData")
load("L3_Rice_Pair_Permutations.RData")

cell_types<-c("mesophyll","bundle_sheath","guard","epidermis","phloem","xylem")
possible_pairs<-combn(cell_types,2)
Diff_Partitioning_Outcome<-function(partition_1, partition_2, pair_column, function_call){

#partition_1<-Sorghum_Pair_Permutations[[1]]
#partition_2<-Rice_Pair_Permutations[[1]]
#pair_column<-1
  
# Sorghum
partitioning_1<-cell_type_subset_sorghum[,str_detect(colnames(cell_type_subset_sorghum), possible_pairs[1,pair_column])]
partitioning_2<-cell_type_subset_sorghum[,str_detect(colnames(cell_type_subset_sorghum), possible_pairs[2,pair_column])]
partitioning_sorghum<-cbind(partitioning_1,partitioning_2)
cell_type_subset_sig_sorghum<-subset(partitioning_sorghum, rownames(partitioning_sorghum) %in% partition_1)

# Rice
rownames(cell_type_subset_rice)<-str_replace(rownames(cell_type_subset_rice),"-","_")
partitioning_1<-cell_type_subset_rice[,str_detect(colnames(cell_type_subset_rice), possible_pairs[1,pair_column])]
partitioning_2<-cell_type_subset_rice[,str_detect(colnames(cell_type_subset_rice), possible_pairs[2,pair_column])]
partitioning_rice<-cbind(partitioning_1,partitioning_2)
cell_type_subset_sig_rice<-subset(partitioning_rice, rownames(partitioning_rice) %in% str_replace(partition_2,"-","_"))

# Compute Rice Orthogroups
rice_orthology<-matrix(NA,nrow=length(cell_type_subset_sig_rice[,1]),ncol=3)
for (i in 1:length(cell_type_subset_sig_rice[,1])) {rice_orthology[i,]<-as.character(orthology[grep(rownames(cell_type_subset_sig_rice)[i], orthology$Rice),])}
rice_orthology<-cbind(rice_orthology, cell_type_subset_sig_rice)
colnames(rice_orthology)[1:3]<-colnames(orthology)
rice_orthology<-subset(rice_orthology, substr(rice_orthology$Orthogroup,1,1) =="O" )
rice_orthology$Rice<-rownames(rice_orthology)

# Compute Sorghum Orthogroups
sorghum_orthology<-matrix(NA,nrow=length(cell_type_subset_sig_sorghum[,1]),ncol=3)
for (i in 1:length(cell_type_subset_sig_sorghum[,1])) {sorghum_orthology[i,]<-as.character(orthology[grep(rownames(cell_type_subset_sig_sorghum)[i], orthology$Sorghum),])}
sorghum_orthology<-cbind(sorghum_orthology, cell_type_subset_sig_sorghum)
colnames(sorghum_orthology)[1:3]<-colnames(orthology)
sorghum_orthology<-subset(sorghum_orthology, substr(sorghum_orthology$Orthogroup,1,1) =="O" )
sorghum_orthology$Sorghum<-rownames(sorghum_orthology)

# Merge By Orthogroups
category_genes_1<-merge(x = rice_orthology, y = sorghum_orthology, by='Orthogroup') 
category_genes_1$Sorghum.x<-category_genes_1$Sorghum.y
category_genes_1<-category_genes_1[,-c(20:21)]

# Calculate Fischer
rice_length<-length(rice_orthology[,1])
sorghum_length<-length(sorghum_orthology[,1])
intersect<-length(category_genes_1[,1])
#intersect<-length(unique(category_genes_1$Sorghum.x)) # older way of doing it when forcing unique names
matrix_to_test<-matrix(c(intersect,(sorghum_length-intersect),(rice_length-intersect),(14444)),2,2) # number of total unique orthogroups
p_value<-as.numeric(fisher.test(matrix_to_test, alternative='greater')[1])
overlap_and_pvalue<-c(intersect,p_value)

# Normalizing Read Counts Across Species
category_genes_rice<-category_genes_1[,c(4:19)]
category_genes_rice_normalized<-t(apply(category_genes_rice, 1, function(x)(x-min(x))/(max(x)-min(x))))
category_genes_sorghum<-category_genes_1[,c(20:35)]
category_genes_sorghum_normalized<-t(apply(category_genes_sorghum, 1, function(x)(x-min(x))/(max(x)-min(x))))
cluster_output<-as.data.frame(cbind(category_genes_rice_normalized, category_genes_sorghum_normalized))
rownames(cluster_output)<-paste(category_genes_1$Rice.x, category_genes_1$Sorghum.x)
cluster_output<-as.matrix(cluster_output)

#quartz()
#heatmap.2(cluster_output, Rowv=as.dendrogram(hclust(as.dist(1-cor(t(cluster_output), method="pearson")))), Colv=NA, col=colorRampPalette(c("blue", "black", "yellow"))(n = 100), trace="none", scale="row",margins = c(8, 16))

# Calulate Consistent and Differential Partitioning
cluster_output_simple<-matrix(NA, nrow=length(cluster_output[,1]), ncol=4)
rownames(cluster_output_simple)<-rownames(cluster_output)
for (i in 1:length(cluster_output[,1])) {cluster_output_simple[i,1]<-  mean(cluster_output[i,1:8], na.rm = TRUE)  # doesnt really matter if mean or sum
                                         cluster_output_simple[i,2]<-  mean(cluster_output[i,9:16], na.rm = TRUE) 
                                         cluster_output_simple[i,3]<-  mean(cluster_output[i,17:24], na.rm = TRUE) 
                                         cluster_output_simple[i,4]<-  mean(cluster_output[i,25:32], na.rm = TRUE) }
cluster_output_simple<- cluster_output_simple[complete.cases(cluster_output_simple),]
colnames(cluster_output_simple)<-c("rice_cell_type_1","rice_cell_type_2","sorghum_cell_type_1","sorghum_cell_type_2")

same_partition<-subset(cluster_output_simple, cluster_output_simple[,1]>cluster_output_simple[,2] & cluster_output_simple[,3]>cluster_output_simple[,4] | 
                                              cluster_output_simple[,1]<cluster_output_simple[,2] & cluster_output_simple[,3]<cluster_output_simple[,4] )
no_same_partiton<-length(same_partition[,1])

dif_partition<-subset(cluster_output_simple, cluster_output_simple[,1]>cluster_output_simple[,2] & cluster_output_simple[,3]<cluster_output_simple[,4] |
                                             cluster_output_simple[,1]<cluster_output_simple[,2] & cluster_output_simple[,3]>cluster_output_simple[,4])
no_dif_partition<-length(dif_partition[,1])

partitioning_numbers<-c(no_same_partiton, no_dif_partition)
output<-c(overlap_and_pvalue, partitioning_numbers)
output
}

partitioning_outcome<-matrix(NA,ncol=4, nrow=15)
for (a in 1:15) {print(a)
                 partitioning_outcome[a,]<-Diff_Partitioning_Outcome(Sorghum_Pair_Permutations[[a]], Rice_Pair_Permutations[[a]], pair_column = a)}

#setwd("~/Desktop/Salk/Project - Light/R files/")
#save(partitioning_outcome, file = "L3_orthology_partitioning_outcome.RData")

setwd("~/Desktop/Salk/Project - Light/R files/")
load("L3_orthology_partitioning_outcome.RData")
partitioning_outcome<-as.data.frame(partitioning_outcome)

# Barplot of Overlap
cell_types_short<-c("M","BS","G","E","P","X")
possible_pairs_short<-combn(cell_types_short,2)
names<-NA
for (b in 1:15){names[b]<-paste(possible_pairs_short[1,b],possible_pairs_short[2,b],sep=" vs ")}
partitioning_outcome[,5]<-names
partitioning_outcome[,6]<- -log(partitioning_outcome[,2],2)
partitioning_outcome[!is.finite(-log(partitioning_outcome[,2],2)),6] <- 1000
partitioning_outcome[,7]<-partitioning_outcome[,4]/partitioning_outcome[,1]
partitioning_outcome[,8]<-partitioning_outcome[,3]/partitioning_outcome[,1]
order<-c("M vs BS","BS vs P","BS vs X","BS vs E","BS vs G","G vs E","G vs X","G vs P","P vs X","E vs P","E vs X","M vs G","M vs E","M vs P","M vs X")

p1<- ggplot(data=partitioning_outcome, aes( x= factor(partitioning_outcome$V5, level=order), y=partitioning_outcome$V1)) + geom_bar(stat="identity", fill="grey50") 
p2<- p1 + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_y_continuous(position = "right") 
p3<- p2 + labs(x = "") + labs(y = "") + theme(axis.text=element_text(size=12),axis.title=element_text(size=12)) #+ scale_x_discrete(limits = positions)
p3
