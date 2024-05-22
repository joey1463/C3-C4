# This script will take the BS cluster specific markers from the mTurquoise line, and overlap this list with cluster specific markers from the Rice atlas.   

setwd("~/Desktop/Salk/Project - Light/R files/")
load("L1_rice_mT_combined.RData")
load("L1_all_mT_Markers.RData")
Rice_BS_Markers<-read.table("mTurq_clusters_0_markers.txt",header=TRUE, sep="\t")  #new
DimPlot(object = rice_mT_combined, reduction = "umap", label=TRUE, raster=FALSE) + NoLegend()
DefaultAssay(object = rice_mT_combined) <- "RNA"
Rice_BS_Markers[,6]<-Rice_BS_Markers$pct.1/Rice_BS_Markers$pct.2

# To Show that rice atlas cluster 10 is best match for mT Cluster 0, look at overlap
setwd("~/Desktop/Salk/Project - Light/R files/cluster_markers")
markers_overlap<-NA
for (i in 0:18){
L1_Markers<-read.table(paste("L1_rice_cluster_",i,"_markers.txt",sep=''))
L1_Markers<-subset(L1_Markers, L1_Markers$V6>2)
mT_Markers_subset<-subset(mT_Markers[[1]], mT_Markers[[1]]$V6>2)
outcome<-as.numeric(rownames(L1_Markers) %in% rownames(mT_Markers_subset))
markers_overlap[i+1]<-(sum(outcome)/length(outcome))}

plot_matrix<-as.data.frame(cbind(markers_overlap*100, c(0:18)))
p1<- ggplot(data=plot_matrix, aes( x= plot_matrix$V2, y=plot_matrix$V1)) + geom_bar(stat="identity", fill="grey70") + theme_bw()
p2<- p1 + labs(x = "") + labs(y = "") + theme(axis.text=element_text(size=12),axis.title=element_text(size=12)) #+ scale_x_discrete(limits = positions)
p2 + scale_x_continuous(plot_matrix$V2, breaks = plot_matrix$V2) + xlab("")

# Module Score Plot
rice_combined <- AddModuleScore(object = rice_combined, features = list(c(rownames(test_set))), name = "BS_Markers")
FeaturePlot(object = rice_combined , features = 'BS_Markers1', cols = c("grey85","dodgerblue3"), min.cutoff = 0, max.cutoff = 1, raster=FALSE)
