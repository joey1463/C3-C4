# This script will calculate differential expression patterns that exist between mesophyll and bundle sheath cell types in rice at either the 0 hour or 12 hour time point

library(Seurat)

load("L1_rice_combined_labelled.RData")
rice_combined_labelled<-rice_combined_labelled[,rice_combined_labelled$assay_type %in% c("10X")]

rice_combined_zero<-rice_combined_labelled[,rice_combined_labelled$Time %in% c("Rice_0")]
rice_combined_twelve<-rice_combined_labelled[,rice_combined_labelled$Time %in% c("Rice_12")]

# Perform Pairwise Tests
T0_response <- FindMarkers(rice_combined_zero, ident.1 = 'bundle_sheath' , ident.2 = 'mesophyll', verbose = FALSE, logfc.threshold = 0.05, min.pct=0.05)  
T12_response <- FindMarkers(rice_combined_twelve, ident.1 = 'bundle_sheath' , ident.2 = 'mesophyll', verbose = FALSE, logfc.threshold = 0.05, min.pct=0.05)  

save(T0_response,file="L3_rice_T0_response_10X.RData")
save(T12_response,file="L3_rice_T12_response_10X.RData")
