# This script will take all genes found differentially expressed in response to light in each cell type, and cluster their expression patterns
# to find dominant gene expression trends.
# These trends are then visualized using a Z-score plot. 

setwd("~/Desktop/Salk/Project - Light/R Files")
chosen_cluster<-c("L3_sorghum_mesophyll","L3_sorghum_bundle_sheath","L3_sorghum_guard","L3_sorghum_phloem","L3_sorghum_epidermis","L3_sorghum_xylem")

DE_genes<-list(NA)
for (i in 1:length(chosen_cluster)){
load(paste(chosen_cluster[i],"_DE_Significant.RData",sep=''))
DE_Significant<-subset(DE_Significant, DE_Significant$padj<2 & DE_Significant$log2FoldChange>0) 
DE_genes[[i]]<-rownames(DE_Significant)}


all_genes<-unique(c(DE_genes[[1]],DE_genes[[2]],DE_genes[[3]],DE_genes[[4]],DE_genes[[5]],DE_genes[[6]]))

test_genes<-subset(time_average, rownames(time_average) %in% all_genes)

plot(hclust(as.dist(1-cor(t(test_genes))), "ave"))
clustering_1<-hclust(as.dist(1-cor(t(test_genes))), "ave")
table(cutree(clustering_1, h=0.7)) #0.6 # 0.7 for Figure 5 &  0.6 for Sup figure
cut_groups<-cutree(clustering_1, h=0.7) #0.6 # 0.7 for Figure 5 &  0.6 for Sup figure
table(cut_groups)
order<-order(table(cut_groups), decreasing=TRUE)

cluster_groups<-list(NA)
for (i in 1:length(unique(cut_groups))) {groupings<-subset(test_genes, rownames(test_genes) %in% names(cut_groups[cut_groups==i]))
                                         cluster_groups[[i]]<-as.data.frame(groupings)
                                         cluster_groups[[i]][,length(test_genes[1,])+1]<-c(rep(i,length(cluster_groups[[i]][,1])))}

plot_list<-list(NA)
for (i in 1:length(order)) {
cluster_group<-cluster_groups[[i]][,c(1:length(test_genes[1,]))]
Z <- as.data.frame(t(scale(t(cluster_group))))
Z_average<-colSums(Z)/length(Z[,1])
Z<-rbind(Z, Z_average)
rownames(Z)[length(Z[,1])]<-"average"
plot_matrix<-melt(t(Z))
factored<-as.character(plot_matrix[,1])
time<-substr(factored,9,11)
time<-str_remove(time, "L") 
time<-as.numeric(str_remove(time, "_"))
cell_type<-substr(factored,22,36)
cell_type<-str_remove(cell_type, "m_") 
cell_type<-str_remove(cell_type, "_") 
cell_type<-str_remove(cell_type, "_")
plot_matrix[,4]<-time
plot_matrix[,5]<-paste(cell_type,plot_matrix[,2])
plot_matrix[,1]<-cell_type
plot_matrix[,6]<-as.numeric(plot_matrix[,2] == 'average')
colnames(plot_matrix)<-c("cell_type","gene","expression","time","gene_tissue","color")
plot_matrix<-subset(plot_matrix, plot_matrix$gene %in% c('average'))
title<-paste(dim(cluster_group)[1],"genes",sep=" ")
plot_list[[i]]<-ggplot(data=plot_matrix, aes(x=time, y=expression, group=gene_tissue)) + geom_smooth(method = "loess", se = FALSE, aes(color=cell_type)) +  theme_bw() +   scale_x_continuous(breaks=c(0,4,8,12)) + ylab("Z-score") + xlab("time (h)") + labs(title=title)}# } # + labs(title=title)


# Figure 3
grid.arrange(plot_list[[order[1]]],plot_list[[order[6]]],plot_list[[order[8]]],plot_list[[order[7]]],plot_list[[order[4]]],plot_list[[order[5]]],nrow=2) 

# Supplementary Fig for Downregulated
grid.arrange(plot_list[[order[1]]],plot_list[[order[2]]],plot_list[[order[4]]],plot_list[[order[6]]],plot_list[[order[7]]],plot_list[[order[9]]],nrow=2)
