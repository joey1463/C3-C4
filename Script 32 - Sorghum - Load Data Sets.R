# This script will read in the pseudo-bulked expression profiles, as well as read in the names of photosynthesis related genes and transcription factors.

setwd("~/Desktop/Salk/Project - Light/R files")

# Time Average
load('L3_sorghum_time_average_10X.RData')
time_average<-sorghum_time_average[,order(colnames(sorghum_time_average))]
cell_type_mesophyll<-time_average[,grepl("mesophyll", colnames(time_average), fixed = TRUE)]
cell_type_bundle_sheath<-time_average[,grepl("bundle_sheath", colnames(time_average), fixed = TRUE)]
cell_type_guard<-time_average[,grepl("guard", colnames(time_average), fixed = TRUE)]
cell_type_phloem<-time_average[,grepl("phloem", colnames(time_average), fixed = TRUE)]
cell_type_epidermis<-time_average[,grepl("epidermis", colnames(time_average), fixed = TRUE)]
cell_type_xylem<-time_average[,grepl("xylem", colnames(time_average), fixed = TRUE)]
cell_type_subset<-cbind(cell_type_mesophyll,cell_type_bundle_sheath,cell_type_guard,cell_type_phloem,cell_type_epidermis,cell_type_xylem)
cell_type_subset_sorghum<-cell_type_subset

# Time Aggregate 
load('L3_sorghum_time_aggregate_10X.RData')
time_aggregate<-sorghum_time_aggregate[,order(colnames(sorghum_time_aggregate))]
time_aggregate<-subset(time_aggregate, rownames(time_aggregate) %in% rownames(time_average))

# Photosynthesis Genes
sorghum_lightdark_extensive<-read.table('extensive_photosynthesis_Sorghum.txt',header=TRUE, sep="\t")
sorghum_lightdark_extensive<-sorghum_lightdark_extensive[!duplicated(sorghum_lightdark_extensive[,1]),]
rownames(sorghum_lightdark_extensive)<-sorghum_lightdark_extensive[,1]
Sorghum_TFs<-read.table("Sbi_TF_list.txt",header=TRUE) # TFs
Sorghum_DOFs<-unique(subset(Sorghum_TFs, Sorghum_TFs$Family == "Dof")[,1])

# Export Gene Lists to Supplementary
chosen_cluster<-c("L3_sorghum_mesophyll","L3_sorghum_bundle_sheath","L3_sorghum_guard","L3_sorghum_phloem","L3_sorghum_epidermis","L3_sorghum_xylem")
DE_genes<-list(NA)
for (i in 1:length(chosen_cluster)){
setwd("~/Desktop/Salk/Project - Light/R Files")
load(paste(chosen_cluster[i],"_DE_Significant.RData",sep=''))
DE_Significant<-subset(DE_Significant, DE_Significant$padj<1)  
setwd("~/Desktop/")
write.table(DE_Significant,file=paste(chosen_cluster[i],".txt",sep=''),quote=F,sep="\t")}
