# This script will take the labeled Rice atlas data, and create pseduobulked expression profiles for each cell type at each time point. 
#Additionally, this script will isolate the 48h time point for use in downstream orthology analysses

library(Seurat)

# Create read count for downstream pseudobulk analyses (time points 0 - 12 h)

load("L1_rice_combined_labelled.RData")
rice_combined_labelled<-rice_combined_labelled[,rice_combined_labelled$assay_type %in% c("10X")]
rice_combined_labelled$cell_types<-Idents(rice_combined_labelled)
rice_combined_labelled<-rice_combined_labelled[,rice_combined_labelled$cell_types %in% c("mesophyll","bundle_sheath","xylem","phloem","guard","epidermis")]

Idents(rice_combined_labelled)<-paste(rice_combined_labelled$light_time,rice_combined_labelled$cell_types,sep='_')
rice_time_average<-as.data.frame(x = AverageExpression(object = rice_combined_labelled, verbose = FALSE)$RNA)
rice_time_aggregate<-AggregateExpression(object = rice_combined_labelled, slot = "counts")
rice_time_aggregate<-as.data.frame(rice_time_aggregate[[1]])

save(rice_time_average, file="L3_rice_time_average_10X.RData")
save(rice_time_aggregate, file="L3_rice_time_aggregate_10X.RData")


# Create read count for downstream othology analyses (time points 48 h)

load("L1_rice_combined_labelled.RData")
rice_combined_labelled<-rice_combined_labelled[,rice_combined_labelled$assay_type %in% c("CI")]
rice_combined_CI_labelled<-rice_combined_labelled[,substr(colnames(rice_combined_labelled),6,7) == '48'] 
save(rice_combined_CI_labelled, file="L1_rice_combined_CI_labelled.RData")
