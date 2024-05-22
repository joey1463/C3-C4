# This script will read in the pseudo-bulked expression profiles, as well as read in the names of photosynthesis related genes and transcription factors.

setwd("~/Desktop/Salk/Project - Light/R Files")

# Time Average
load('L3_rice_time_average_10X.RData')
time_average<-rice_time_average[,order(colnames(rice_time_average))]
cell_type_mesophyll<-time_average[,grepl("mesophyll", colnames(time_average), fixed = TRUE)]
cell_type_bundle_sheath<-time_average[,grepl("bundle_sheath", colnames(time_average), fixed = TRUE)]
cell_type_guard<-time_average[,grepl("guard", colnames(time_average), fixed = TRUE)]
cell_type_phloem<-time_average[,grepl("phloem", colnames(time_average), fixed = TRUE)]
cell_type_epidermis<-time_average[,grepl("epidermis", colnames(time_average), fixed = TRUE)]
cell_type_xylem<-time_average[,grepl("xylem", colnames(time_average), fixed = TRUE)]
cell_type_subset<-cbind(cell_type_mesophyll,cell_type_bundle_sheath,cell_type_guard,cell_type_phloem,cell_type_epidermis,cell_type_xylem)
cell_type_subset<-cell_type_subset[,c(c(1:3,5:8,4),c(1:3,5:8,4)+8,c(1:3,5:8,4)+16,c(1:3,5:8,4)+24,c(1:3,5:8,4)+32,c(1:3,5:8,4)+40)]
cell_type_subset_rice<-cell_type_subset

# Time Aggregate
load('L3_rice_time_aggregate_10X.RData')
time_aggregate<-rice_time_aggregate[,order(colnames(rice_time_aggregate))]
time_aggregate<-subset(time_aggregate, rownames(time_aggregate) %in% rownames(time_average))

# Gene Lists
setwd("~/Desktop/Salk/Project - Light/R files/")
rice_lightdark_extensive<-read.table('extensive_photosynthesis_Rice.txt',header=TRUE, sep="\t")
rownames(rice_lightdark_extensive)<-rice_lightdark_extensive[,1]
Rice_TFs<-read.table("Osj_TF_list.txt",header=TRUE)
Rice_DOFs<-unique(subset(Rice_TFs, Rice_TFs[,2] == "Dof")[,1])
Rice_TFs<-unique(Rice_TFs[,1])

# Export Supplementary Gene List Files
chosen_cluster<-c("L3_rice_mesophyll","L3_rice_bundle_sheath","L3_rice_guard","L3_rice_phloem","L3_rice_epidermis","L3_rice_xylem")
DE_genes<-list(NA)
for (i in 1:length(chosen_cluster)){
setwd("~/Desktop/Salk/Project - Light/R Files")
load(paste(chosen_cluster[i],"_DE_Significant.RData",sep=''))
DE_Significant<-subset(DE_Significant, DE_Significant$padj<1) 
setwd("~/Desktop/")
write.table(DE_Significant,file=paste(chosen_cluster[i],".txt",sep=''),quote=F,sep="\t")}
