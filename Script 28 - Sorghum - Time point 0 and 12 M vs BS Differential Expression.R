# This script will calculate differential expression patterns that exist between mesophyll and bundle sheath cell types in sorghum at either the 0 hour or 12 hour time point

library(Seurat)

load("L1_sorghum_combined_labelled.RData")
sorghum_combined_labelled<-sorghum_combined_labelled[,sorghum_combined_labelled$assay_type %in% c("10X")]

sorghum_combined_zero<-sorghum_combined_labelled[,sorghum_combined_labelled$Time %in% c("Sorghum_0")]
sorghum_combined_twelve<-sorghum_combined_labelled[,sorghum_combined_labelled$Time %in% c("Sorghum_12")]

# Perform Pairwise Tests
T0_response <- FindMarkers(sorghum_combined_zero, ident.1 = 'bundle_sheath' , ident.2 = 'mesophyll', verbose = FALSE, logfc.threshold = 0.05, min.pct=0.05)  
T12_response <- FindMarkers(sorghum_combined_twelve, ident.1 = 'bundle_sheath' , ident.2 = 'mesophyll', verbose = FALSE, logfc.threshold = 0.05, min.pct=0.05)  

save(T0_response,file="L3_sorghum_T0_response_10X.RData")
save(T12_response,file="L3_sorghum_T12_response_10X.RData")
