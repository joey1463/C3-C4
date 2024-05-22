# This script will take the labeled Sorghum atlas data, and create pseduobulked expression profiles for each cell type at each time point. 
#Additionally, this script will isolate the 48h time point for use in downstream orthology analysses

library(Seurat)

# Create read count for downstream pseudobulk analyses (time points 0 - 12 h)

load("L1_sorghum_combined_labelled.RData")
sorghum_combined_labelled<-sorghum_combined_labelled[,sorghum_combined_labelled$assay_type %in% c("10X")]
sorghum_combined_labelled$cell_types<-Idents(sorghum_combined_labelled)
sorghum_combined_labelled<-sorghum_combined_labelled[,sorghum_combined_labelled$cell_types %in% c("mesophyll","bundle_sheath","xylem","phloem","guard","epidermis")]

Idents(sorghum_combined_labelled)<-paste(sorghum_combined_labelled$light_time,sorghum_combined_labelled$cell_types,sep='_')
sorghum_time_average<-as.data.frame(x = AverageExpression(object = sorghum_combined_labelled, verbose = FALSE)$RNA)
sorghum_time_aggregate<-AggregateExpression(object = sorghum_combined_labelled, slot = "counts")
sorghum_time_aggregate<-as.data.frame(sorghum_time_aggregate[[1]])

save(sorghum_time_average, file="L3_sorghum_time_average_10X.RData")
save(sorghum_time_aggregate, file="L3_sorghum_time_aggregate_10X.RData")


# Create read count for downstream othology analyses (time points 48 h)

load("L1_sorghum_combined_labelled.RData")
sorghum_combined_labelled<-sorghum_combined_labelled[,sorghum_combined_labelled$assay_type %in% c("CI")]
sorghum_combined_CI_labelled<-sorghum_combined_labelled[,substr(colnames(sorghum_combined_labelled),6,7) == '48'] 
save(sorghum_combined_CI_labelled, file="L1_sorghum_combined_CI_labelled.RData")
