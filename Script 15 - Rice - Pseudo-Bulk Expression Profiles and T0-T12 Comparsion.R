# This script will compute pseudo-bulked transcriptional profiles of each cell type
# This script will also detect which genes are differentially partitioned between the mesophyll and bundle sheath cell types under etiolated conditions. 

library(Seurat)

# Generate pseudo-bulked expression profiles from 10X Data
load("L1_rice_combined.RData")
rice_combined<-rice_combined[,rice_combined$assay_type %in% c("10X")]
cell_type_identification<-colnames(rice_combined)

chosen_cluster<-c("L3_rice_mesophyll","L3_rice_bundle_sheath","L3_rice_xylem","L3_rice_phloem","L3_rice_guard","L3_rice_epidermis")

for (i in 1:length(chosen_cluster)){
  print(i)
  load(paste(chosen_cluster[i],"_10X.RData",sep=''))
  hits<-colnames(rice_combined) %in% colnames(subcluster_rice_combined_10X)
  co_ords<-which(hits == TRUE)
  cell_type_identification[co_ords] <- chosen_cluster[i]}

# Calculate Average and Aggregate Expression on 10X Data
rice_combined$cell_type<-cell_type_identification
rice_combined<-rice_combined[,rice_combined$cell_type %in% chosen_cluster]
Idents(rice_combined)<-paste(rice_combined$light_time,rice_combined$cell_type,sep='_')
rice_time_average<-as.data.frame(x = AverageExpression(object = rice_combined, verbose = FALSE)$RNA)
rice_time_aggregate<-AggregateExpression(object = rice_combined, slot = "counts")
rice_time_aggregate<-as.data.frame(rice_time_aggregate[[1]])

save(rice_time_average, file="L3_rice_time_average_10X.RData")
save(rice_time_aggregate, file="L3_rice_time_aggregate_10X.RData")


# Compute T0 (etiolated) and T12 Response in the mesophyll and bundle sheath
Idents(rice_combined)<-rice_combined$cell_type

rice_combined_zero<-rice_combined[,rice_combined$Time %in% c("Rice_0")]
rice_combined_twelve<-rice_combined[,rice_combined$Time %in% c("Rice_12")]

T0_response <- FindMarkers(rice_combined_zero, ident.1 = 'L3_rice_bundle_sheath' , ident.2 = 'L3_rice_mesophyll', verbose = FALSE, logfc.threshold = 0.05, min.pct=0.05)  
T12_response <- FindMarkers(rice_combined_twelve, ident.1 = 'L3_rice_bundle_sheath' , ident.2 = 'L3_rice_mesophyll', verbose = FALSE, logfc.threshold = 0.05, min.pct=0.05)  

save(T0_response,file="L3_rice_T0_response_10X.RData")
save(T12_response,file="L3_rice_T12_response_10X.RData")
