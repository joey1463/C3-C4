# This script will compute enriched motifs found within the accessbile chromatin of differentially partitioned genes in Sorghum

# Read In Data
setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_sorghum_multiome_RNA_combined_with_ATAC_with_ChromVAR.RData")
DefaultAssay(sorghum_multiome_RNA_combined) <- "ATAC"
load("Swap_Genes.RData")

# Run Motifs
target_list<-paste(Swap_Genes$Sorghum.x,".v3.2",sep="")
sorghum_multiome_RNA_combined <- LinkPeaks(object = sorghum_multiome_RNA_combined, peak.assay = "ATAC", expression.assay = "RNA", genes.use = target_list)

# Pull out Peaks
DefaultAssay(sorghum_multiome_RNA_combined) <- "ATAC"
Links(sorghum_multiome_RNA_combined)
my_hits<-as.data.frame(Links(sorghum_multiome_RNA_combined))
my_hits[,11]<-p.adjust(my_hits$pvalue, method = "fdr")
my_hits<-subset(my_hits, my_hits[,11]<0.001) 
my_hits<-subset(my_hits, my_hits$score>0)
my_hits<-my_hits[ !duplicated(my_hits$peak),]
rownames(my_hits)<-my_hits$peak
my_hits$gene<-substr(my_hits$gene,1,16)
#save(my_hits, file="my_hits_sorghum_swap.RData")
#load("my_hits_sorghum_swap.RData")

# Call Motifs
enriched_motifs <- FindMotifs(object = sorghum_multiome_RNA_combined, features = rownames(my_hits))
#save(enriched_motifs, file="enriched_motifs_sorghum_swap.RData")
enriched_motifs <-subset(enriched_motifs, enriched_motifs$pvalue<0.01)
MotifPlot(object =sorghum_multiome_RNA_combined,motifs = rownames(enriched_motifs)[1:12])
#load("enriched_motifs_sorghum_swap.RData")

# Iterate to Find Stability
iterations<-matrix(nrow=10,ncol=100)
for (i in 1:100){
enriched_motifs <- FindMotifs(object = sorghum_multiome_RNA_combined, features = rownames(my_hits))
enriched_motifs <-subset(enriched_motifs, enriched_motifs$pvalue<0.01)[1:10,]
iterations[,i]<- rownames(enriched_motifs)/100}
#save(iterations, file="iterations_enriched_motifs_sorghum_swap.RData")
#load(file="iterations_enriched_motifs_sorghum_swap.RData")

# Rank Motifs
all_motifs<-unique(as.vector(iterations))
all_ranks<-NA
for (b in 1:length(all_motifs)){
rank<-NA
for (a in 1:100){
rank[a]<-match(all_motifs[b],iterations[,a])}
rank[is.na(rank)] <- 11 # giving rank 11 to those not found
all_ranks[b]<-sum(rank)}

motif_outcome<-as.data.frame(cbind(all_motifs, all_ranks))
motif_outcome[,2]<-as.numeric(motif_outcome[,2])
motif_outcome<-motif_outcome[order(motif_outcome[,2]),]
motif_outcome_top<-motif_outcome
MotifPlot(object =sorghum_multiome_RNA_combined,motifs = motif_outcome_top[,1])

# Export with JASPAR description
JASPAR_Plants<-read.table("JASPAR_Plant_TFs.txt",header=TRUE, sep="\t")
rownames(JASPAR_Plants)<-JASPAR_Plants[,1]
rownames(motif_outcome)<-motif_outcome[,1]
motif_outcome <- merge(JASPAR_Plants,motif_outcome,by='row.names',all.y = TRUE)
