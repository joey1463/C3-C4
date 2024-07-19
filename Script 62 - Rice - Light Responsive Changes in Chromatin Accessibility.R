# This script will compute plot light-responsive changes in chromatin accessibility for select genes.
# It will also assess the changes in accessible chromatin among canonical photosynthesis genes within different cell types. 

setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_rice_multiome_RNA_combined_with_ATAC.RData")
Idents(rice_multiome_RNA_combined)<-paste(Idents(rice_multiome_RNA_combined), str_remove(substr(rice_multiome_RNA_combined$dataset,1,5), "_"), sep="_")

# Plotting Canonicals (will need to assemble fresh)
small_plot_rice<-rice_multiome_RNA_combined[,WhichCells(rice_multiome_RNA_combined, idents = c("mesophyll_Light","mesophyll_Dark","bundle_sheath_Light","bundle_sheath_Dark"))]
Idents(small_plot_rice) <- factor(small_plot_rice@active.ident, c("mesophyll_Dark","bundle_sheath_Dark","mesophyll_Light","bundle_sheath_Light"))
save(small_plot_rice, file="small_plot_rice.RData")
load("small_plot_rice.RData")

p<- CoveragePlot(object = small_plot_rice ,region = "LOC-Os12g17600",expression.assay = "RNA",extend.upstream = 0,extend.downstream = 500) # RBCS2
p<- CoveragePlot(object = small_plot_rice ,region = "LOC-Os03g57220",expression.assay = "RNA",extend.upstream = 500,extend.downstream = -1000) # glycolate oxydase 
p & scale_fill_manual(values = c("green4","blue4", "green3","blue2")) 


# Calling All Light-Dark Peaks (must be done on re-assembled data)
DefaultAssay(rice_multiome_RNA_combined) <- "ATAC"
Average_Peaks<-as.data.frame(AverageExpression(object = rice_multiome_RNA_combined)[[3]]) 
Average_Peaks<-Average_Peaks[,order(colnames(Average_Peaks))]
Peaks_Close_Feature<-ClosestFeature(rice_multiome_RNA_combined, rownames(rice_multiome_RNA_combined))
Peaks_Close_Feature<-subset(Peaks_Close_Feature, Peaks_Close_Feature$distance<2000)
Light_Peaks_Average<-merge(Average_Peaks,Peaks_Close_Feature, by.x='row.names', by.y='query_region')
rownames(Light_Peaks_Average)<-Light_Peaks_Average[,1]
Light_Peaks_Average<-Light_Peaks_Average[,-1]
#save(Light_Peaks_Average, file="Light_Peaks_Average_Rice.RData")
#load("Light_Peaks_Average_Rice.RData")

# Canonicals in Violin Plot
heatmap_cluster_ordered_rice # Computed from RNA-atlas
heatmap_ATAC<-subset(Light_Peaks_Average, Light_Peaks_Average$transcript_id %in% rownames(heatmap_cluster_ordered_rice))
heatmap_ATAC<-subset(heatmap_ATAC, heatmap_ATAC$distance<2000)
heatmap_ATAC<-heatmap_ATAC[!duplicated(heatmap_ATAC$transcript_id),]
rownames(heatmap_ATAC)<-heatmap_ATAC$transcript_id
heatmap_ATAC<-heatmap_ATAC[,c(9,10,1,2,7,8,11,12,3,4,15,16)]
plot_matrix<-melt(as.data.frame(heatmap_ATAC))
plot_matrix[,3]<-rep(rownames(heatmap_ATAC),12)
colnames(plot_matrix)<-c("condition","value")
plot_matrix<-subset(plot_matrix, plot_matrix$value<15)
plot_matrix$value<-(plot_matrix$value-min(plot_matrix$value))/(max(plot_matrix$value)-min(plot_matrix$value))
ggplot(plot_matrix, aes(x=as.factor(condition), y=value,  color = as.factor(condition))) +  theme_bw() + geom_violin(trim=TRUE) + stat_summary(fun = "median",geom = "crossbar", width = 0.2, colour = "black") + geom_jitter(shape=1, position=position_jitter(0.2), size = 0.2, alpha=1) + ylim(c(0,1)) +  ylab("") +  xlab("") + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + scale_color_manual(values=rep(c("blue4", "yellow4"),6))+ scale_y_continuous(breaks=c(0,1))


# reporting the p value for claim that difference between mesophyll light and dark is marginal
t.test(subset(plot_matrix, plot_matrix$condition == 'mesophyll_Light')$value, 
       subset(plot_matrix, plot_matrix$condition == 'mesophyll_Dark')$value, alternative= 'greater') # 0.05292

# reporting the highest p-value for claim that mesophyll is more acccessible compared to other cell types
t.test(subset(plot_matrix, plot_matrix$condition %in% c('mesophyll_Light','mesophyll_Dark'))$value, 
       subset(plot_matrix, plot_matrix$condition %in% c('bundle_sheath_Light','bundle_sheath_Dark'))$value,alternative= 'greater') # 4.01e-16
t.test(subset(plot_matrix, plot_matrix$condition %in% c('mesophyll_Light','mesophyll_Dark'))$value, 
       subset(plot_matrix, plot_matrix$condition %in% c('guard_Light','guard_Dark'))$value,alternative= 'greater') #  7.467e-16
t.test(subset(plot_matrix, plot_matrix$condition %in% c('mesophyll_Light','mesophyll_Dark'))$value, 
       subset(plot_matrix, plot_matrix$condition %in% c('epidermis_Light','epidermis_Dark'))$value,alternative= 'greater') # 7.62e-16
t.test(subset(plot_matrix, plot_matrix$condition %in% c('mesophyll_Light','mesophyll_Dark'))$value, 
       subset(plot_matrix, plot_matrix$condition %in% c('xylem_Light','xylem_Dark'))$value,alternative= 'greater') # 7.115e-16
t.test(subset(plot_matrix, plot_matrix$condition %in% c('mesophyll_Light','mesophyll_Dark'))$value, 
       subset(plot_matrix, plot_matrix$condition %in% c('phloem_Light','phloem_Dark'))$value,alternative= 'greater') # < 2.2e-16
