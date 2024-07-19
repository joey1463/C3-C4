# This script will compute pseudo-bulked transcriptional profiles of each cell type
# This script will also detect which genes are differentially partitioned between the mesophyll and bundle sheath cell types under etiolated conditions. 

library(Seurat)

### 10X Data
load("L1_sorghum_combined.RData")
sorghum_combined<-sorghum_combined[,sorghum_combined$assay_type %in% c("10X")]
cell_type_identification<-colnames(sorghum_combined)

chosen_cluster<-c("L3_sorghum_mesophyll","L3_sorghum_bundle_sheath","L3_sorghum_xylem","L3_sorghum_phloem","L3_sorghum_guard","L3_sorghum_epidermis")

for (i in 1:length(chosen_cluster)){
  print(i)
  load(paste(chosen_cluster[i],"_10X.RData",sep=''))
  hits<-colnames(sorghum_combined) %in% colnames(subcluster_sorghum_combined_10X)
  co_ords<-which(hits == TRUE)
  cell_type_identification[co_ords] <- chosen_cluster[i]}

# Get Average and Aggregate Expression on 10X Data
sorghum_combined$cell_type<-cell_type_identification
sorghum_combined<-sorghum_combined[,sorghum_combined$cell_type %in% chosen_cluster]
Idents(sorghum_combined)<-paste(sorghum_combined$light_time,sorghum_combined$cell_type,sep='_')
sorghum_time_average<-as.data.frame(x = AverageExpression(object = sorghum_combined, verbose = FALSE)$RNA)
sorghum_time_aggregate<-AggregateExpression(object = sorghum_combined, slot = "counts")
sorghum_time_aggregate<-as.data.frame(sorghum_time_aggregate[[1]])

save(sorghum_time_average, file="L3_sorghum_time_average_10X.RData")
save(sorghum_time_aggregate, file="L3_sorghum_time_aggregate_10X.RData")

# T0 and T12 Response
Idents(sorghum_combined)<-sorghum_combined$cell_type

sorghum_combined_zero<-sorghum_combined[,sorghum_combined$Time %in% c("Sorghum_0")]
sorghum_combined_twelve<-sorghum_combined[,sorghum_combined$Time %in% c("Sorghum_12")]

T0_response <- FindMarkers(sorghum_combined_zero, ident.1 = 'L3_sorghum_bundle_sheath' , ident.2 = 'L3_sorghum_mesophyll', verbose = FALSE, logfc.threshold = 0.05, min.pct=0.05)  
T12_response <- FindMarkers(sorghum_combined_twelve, ident.1 = 'L3_sorghum_bundle_sheath' , ident.2 = 'L3_sorghum_mesophyll', verbose = FALSE, logfc.threshold = 0.05, min.pct=0.05)  

save(T0_response,file="L3_sorghum_T0_response_10X.RData")
save(T12_response,file="L3_sorghum_T12_response_10X.RData")
