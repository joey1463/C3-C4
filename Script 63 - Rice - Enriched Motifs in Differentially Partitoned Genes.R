# This script will compute enriched motifs found within the accessbile chromatin of differentially partitioned genes in Rice. 

setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_rice_multiome_RNA_combined_with_ATAC_with_ChromVAR.RData")
DefaultAssay(rice_multiome_RNA_combined) <- "ATAC"
load("Swap_Genes.RData")

# Run Motifs
target_list<-str_replace(Swap_Genes$Rice.x,"_","-")
rice_multiome_RNA_combined <- LinkPeaks(object = rice_multiome_RNA_combined, peak.assay = "ATAC", expression.assay = "RNA", genes.use = target_list)

# Pull out Peaks
DefaultAssay(rice_multiome_RNA_combined) <- "ATAC"
Links(rice_multiome_RNA_combined)
my_hits<-as.data.frame(Links(rice_multiome_RNA_combined))
my_hits[,11]<-p.adjust(my_hits$pvalue, method = "fdr")
my_hits<-subset(my_hits, my_hits[,11]<0.001) 
my_hits<-subset(my_hits, my_hits$score>0)
my_hits<-my_hits[ !duplicated(my_hits$peak),]
rownames(my_hits)<-my_hits$peak
#save(my_hits, file="my_hits_rice_swap.RData")
#load("my_hits_rice_swap.RData")

# Call Motifs
enriched_motifs <- FindMotifs(object = rice_multiome_RNA_combined, features = rownames(my_hits))
#save(enriched_motifs, file="enriched_motifs_rice_swap.RData")
#load("enriched_motifs_rice_swap.RData")
enriched_motifs <-subset(enriched_motifs, enriched_motifs$pvalue<0.01)
MotifPlot(object =rice_multiome_RNA_combined,motifs = rownames(enriched_motifs)[1:12])

# Iterate to Ensure Stability
iterations<-matrix(nrow=10,ncol=100)
for (i in 1:100){
enriched_motifs <- FindMotifs(object = rice_multiome_RNA_combined, features = rownames(my_hits))
enriched_motifs <-subset(enriched_motifs, enriched_motifs$pvalue<0.01)[1:10,]
iterations[,i]<- rownames(enriched_motifs)}
#save(iterations, file="iterations_enriched_motifs_rice_swap.RData")
#load(file="iterations_enriched_motifs_rice_swap.RData")

# Rank Motifs
all_motifs<-unique(as.vector(iterations))
all_ranks<-NA
for (b in 1:length(all_motifs)){
rank<-NA
for (a in 1:100){ rank[a]<-match(all_motifs[b],iterations[,a])}
                  rank[is.na(rank)] <- 11 # giving rank 11 to those not found
                  all_ranks[b]<-sum(rank)/100}

motif_outcome<-as.data.frame(cbind(all_motifs, all_ranks))
motif_outcome[,2]<-as.numeric(motif_outcome[,2])
motif_outcome<-motif_outcome[order(motif_outcome[,2]),]
motif_outcome_top<-motif_outcome[1:10,]
MotifPlot(object =rice_multiome_RNA_combined,motifs = motif_outcome_top[,1])

# Export with JASPAR description
JASPAR_Plants<-read.table("JASPAR_Plant_TFs.txt",header=TRUE, sep="\t")
rownames(JASPAR_Plants)<-JASPAR_Plants[,1]
rownames(motif_outcome)<-motif_outcome[,1]
motif_outcome <- merge(JASPAR_Plants,motif_outcome,by='row.names',all.y = TRUE)
