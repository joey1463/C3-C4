# This script will compute plot light-responsive changes in chromatin accessibility for select genes.
# It will also assess the changes in accessible chromatin among canonical photosynthesis genes within different cell types. 


setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_sorghum_multiome_RNA_combined_with_ATAC.RData")
Idents(sorghum_multiome_RNA_combined)<-paste(Idents(sorghum_multiome_RNA_combined), str_remove(substr(sorghum_multiome_RNA_combined$dataset,1,5), "_"), sep="_")

# Testing Canonicals
small_plot_sorghum<-sorghum_multiome_RNA_combined[,WhichCells(sorghum_multiome_RNA_combined, idents = c("mesophyll_Light","mesophyll_Dark","bundle_sheath_Light","bundle_sheath_Dark"))]
Idents(small_plot_sorghum) <- factor(small_plot_sorghum@active.ident, c("mesophyll_Dark","bundle_sheath_Dark","mesophyll_Light","bundle_sheath_Light"))
#save(small_plot_sorghum, file="small_plot_sorghum.RData")
#load("small_plot_sorghum.RData")

p<- CoveragePlot(object = small_plot_sorghum, region = "Sobic.005G042000.v3.2", extend.upstream = 0, extend.downstream = 500) # Rubisco subunit
p<- CoveragePlot(object = small_plot_sorghum, region = "Sobic.002G036000.v3.2", extend.upstream = 500, extend.downstream = -1000) # glycolate oxydase
p & scale_fill_manual(values = c("green4","blue4", "green3","blue2")) 

# Calling All Light-Dark Peaks
DefaultAssay(sorghum_multiome_RNA_combined) <- "ATAC" 
Average_Peaks<-as.data.frame(AverageExpression(object = sorghum_multiome_RNA_combined)[[3]]) 
Average_Peaks<-Average_Peaks[,order(colnames(Average_Peaks))]
Peaks_Close_Feature<-ClosestFeature(sorghum_multiome_RNA_combined, rownames(sorghum_multiome_RNA_combined))
Peaks_Close_Feature<-subset(Peaks_Close_Feature, Peaks_Close_Feature$distance<2000)
Light_Peaks_Average<-merge(Average_Peaks,Peaks_Close_Feature, by.x='row.names', by.y='query_region')
rownames(Light_Peaks_Average)<-Light_Peaks_Average[,1]
Light_Peaks_Average<-Light_Peaks_Average[,-1]
#save(Light_Peaks_Average, file="Light_Peaks_Average_Sorghum.RData")
#load("Light_Peaks_Average_Sorghum.RData")

# Canonical Photosynthesis Genes in Violin Plot
heatmap_cluster_ordered_sorghum # heatmap 
heatmap_ATAC<-subset(Light_Peaks_Average, Light_Peaks_Average$transcript_id %in% rownames(heatmap_cluster_ordered_sorghum))
heatmap_ATAC<-subset(heatmap_ATAC, heatmap_ATAC$distance<2000)
heatmap_ATAC<-heatmap_ATAC[!duplicated(heatmap_ATAC$transcript_id),]
rownames(heatmap_ATAC)<-heatmap_ATAC$transcript_id
heatmap_ATAC<-as.matrix(heatmap_ATAC[,c(1:10,13,14)])
heatmap_cluster_ordered_subset<-subset(heatmap_cluster_ordered_sorghum, rownames(heatmap_cluster_ordered_sorghum) %in% rownames(heatmap_ATAC))
heatmap_ATAC<-heatmap_ATAC[rownames(heatmap_cluster_ordered_subset),]
heatmap_ATAC<-heatmap_ATAC[,c(7,1,8,2,5,6,9,10,3,4,11,12)] 
heatmap_ATAC<-heatmap_ATAC[,c(1,3,2,4,5:12)] # reorder
heatmap_ATAC<-as.data.frame(heatmap_ATAC)
plot_matrix<-melt(heatmap_ATAC)
plot_matrix[,3]<-rep(rownames(heatmap_ATAC),12)
colnames(plot_matrix)<-c("condition","value")
plot_matrix<-subset(plot_matrix, plot_matrix$value<2)
plot_matrix$value<-(plot_matrix$value-min(plot_matrix$value))/(max(plot_matrix$value)-min(plot_matrix$value))

ggplot(plot_matrix, aes(x=as.factor(condition), y=value,  color = as.factor(condition))) +  theme_bw() + geom_violin(trim=TRUE) + stat_summary(fun = "median",geom = "crossbar", width = 0.2, colour = "black") + geom_jitter(shape=1, position=position_jitter(0.2), size = 0.2, alpha=1) + ylim(c(0,1)) +  ylab("") +  xlab("") + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + scale_color_manual(values=rep(c("blue4", "yellow4"),6)) + scale_y_continuous(breaks=c(0,1))

# reporting the p value for claim that difference between mesophyll light and dark as well as BS light and dark is significant. 
t.test(subset(plot_matrix, plot_matrix$condition == 'mesophyll_Light')$value, 
       subset(plot_matrix, plot_matrix$condition == 'mesophyll_Dark')$value, alternative= 'greater') # 9.753e-16
t.test(subset(plot_matrix, plot_matrix$condition == 'bundle_sheath_Light')$value, 
       subset(plot_matrix, plot_matrix$condition == 'bundle_sheath_Dark')$value, alternative= 'greater') #  4.69e-10

# Pairwise t-test
t.test(subset(plot_matrix, plot_matrix$condition == 'guard_Light')$value, 
       subset(plot_matrix, plot_matrix$condition == 'guard_Dark')$value, alternative= 'greater') 
t.test(subset(plot_matrix, plot_matrix$condition == 'phloem_Light')$value, 
       subset(plot_matrix, plot_matrix$condition == 'phloem_Dark')$value, alternative= 'greater') 
t.test(subset(plot_matrix, plot_matrix$condition == 'epidermis_Light')$value, 
       subset(plot_matrix, plot_matrix$condition == 'epidermis_Dark')$value, alternative= 'greater')
t.test(subset(plot_matrix, plot_matrix$condition == 'xylem_Light')$value, 
       subset(plot_matrix, plot_matrix$condition == 'xylem_Dark')$value, alternative= 'greater') 
