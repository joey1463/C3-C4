# This script will use motif scanning to count how many Dof.2 motifs are present within differentially partitioned genes in Rice. 

# Read in required datasets
setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_rice_multiome_Average_Peaks.RData")
load("L1_rice_multiome_RNA_combined_with_ATAC_with_ChromVAR.RData")
load("Target_Genes.RData")
load("L1_cis_sig_list_cross_species.RData")
Dof2_motif <- MotifDb::MotifDb %>%  MotifDb::query(cis_sig_list_identity[[5]][22]) %>%  universalmotif::convert_motifs() %>%  .[[5]]  # Dof2 motif

# Calling All Peaks
Average_Peaks<-Average_Peaks[,order(colnames(Average_Peaks))]
Peaks_Close_Feature<-ClosestFeature(rice_multiome_RNA_combined, rownames(rice_multiome_RNA_combined))
Peaks_Close_Feature<-subset(Peaks_Close_Feature, Peaks_Close_Feature$distance<1500) 
DE_Peaks_Average_Rice<-merge(Average_Peaks,Peaks_Close_Feature, by.x='row.names', by.y='query_region')
DE_Peaks_Average_Rice$transcript_id<-str_replace(DE_Peaks_Average_Rice$transcript_id,"_","-")
Target_Peaks<-subset(DE_Peaks_Average_Rice, DE_Peaks_Average_Rice$transcript_id %in% str_replace(Target_Genes$Rice,"_","-"))

# Count how many Dof2 sites 
# Counting is done using FIMO with a sigificance threshold of 0.005. Count is compared to a manual count of 'AAAG' occurances as a sanity check
Counts_Outcome<-as.data.frame(matrix(NA, nrow=length(Target_Peaks[,1]), ncol=2))
colnames(Counts_Outcome)<-c("manual_count","fimo_count")

for (i in 1:length(Target_Peaks[,1])){
target_peak<-as.character(sub("-", ":", Target_Peaks$Row.names[i], fixed = TRUE))
sequences <- target_peak %>% get_sequence(BSgenome.Osativa.MSU.1) 
manual_count<-str_count(as.character(sequences[[1]]), pattern = 'AAAG') + str_count(as.character(sequences[[1]]), pattern = 'CTTT')
fimo_count<-rbind(as.data.frame(runFimo(sequences, Dof2_motif,thresh=0.005)))
fimo_count<-fimo_count[!duplicated(fimo_count$start), ]
fimo_count<- dim(fimo_count)[1]
Counts_Outcome[i,]<-c(manual_count, fimo_count)}

Target_Peaks<-cbind(Target_Peaks, Counts_Outcome)
Counts_Outcome$fimo_count == Counts_Outcome$manual_count 

# Create DOF Count Table 
gene_list<-unique(Target_Peaks$transcript_id)
DOF_count_by_gene<-as.data.frame(matrix(NA, nrow=length(gene_list),ncol=2))
rownames(DOF_count_by_gene)<-gene_list
for (i in 1:length(gene_list)){
DOF_count_by_gene[i,1]<-sum(subset(Target_Peaks, Target_Peaks$transcript_id == gene_list[i])[,26])
DOF_count_by_gene[i,2]<-sum(subset(Target_Peaks, Target_Peaks$transcript_id == gene_list[i])[,27])}
colnames(DOF_count_by_gene)<-c("manual_count_sorg","fimo_count_sorg")
DOF_count_by_gene_Rice<-DOF_count_by_gene
save(DOF_count_by_gene_Rice, file="DOF_count_by_gene_Rice.RData")
