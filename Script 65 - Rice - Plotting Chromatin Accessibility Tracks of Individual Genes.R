# This script will plot accessbile chromatin tracks within mesophyll and bundle sheath cell types for candidate genes in Rice

# Read in Data
setwd("~/Desktop/Salk/Project - Light/R files")
load("small_plot_rice.RData")
Idents(small_plot_rice)<-substr(Idents(small_plot_rice),1,9)
Idents(small_plot_rice) <- factor(small_plot_rice@active.ident, c("mesophyll","bundle_sh"))
load("L1_cis_sig_list_cross_species.RData")
Dof2_motif <- MotifDb::MotifDb %>%  MotifDb::query(cis_sig_list_identity[[5]][22]) %>%  universalmotif::convert_motifs() %>%  .[[5]]  

# Calling All Peaks
load("L1_rice_multiome_Average_Peaks.RData")
load("L1_rice_multiome_RNA_combined_with_ATAC_with_ChromVAR.RData")
Average_Peaks<-Average_Peaks[,order(colnames(Average_Peaks))]
Peaks_Close_Feature<-ClosestFeature(rice_multiome_RNA_combined, rownames(rice_multiome_RNA_combined))
Peaks_Close_Feature<-subset(Peaks_Close_Feature, Peaks_Close_Feature$distance<1500)
DE_Peaks_Average_Rice<-merge(Average_Peaks,Peaks_Close_Feature, by.x='row.names', by.y='query_region')
DE_Peaks_Average_Rice$transcript_id<-str_replace(DE_Peaks_Average_Rice$transcript_id,"_","-")
remove(rice_multiome_RNA_combined)


# GAPDH 
p<- CoveragePlot(object = small_plot_rice,region = "LOC-Os04g38600",expression.assay = "RNA",extend.upstream = 1500,extend.downstream = 400) 
p & scale_fill_manual(values = c("green4","blue4")) 
subset(DE_Peaks_Average_Rice, DE_Peaks_Average_Rice$gene_id == "LOC-Os04g38600")
peak_1<-"Chr4:22937829-22938032" # in the gene
sequences <- peak_1 %>%  get_sequence(BSgenome.Osativa.MSU.1) 
str_count(as.character(sequences[[1]]), pattern = 'AAAG') + str_count(as.character(sequences[[1]]), pattern = 'CTTT')
motif_bundle<-rbind(as.data.frame(runFimo(sequences, Dof2_motif,thresh=0.005)))
dim(motif_bundle)[1]
peak_2<-"Chr4:22936400-22937186" # promoter
sequences <- peak_2 %>%  get_sequence(BSgenome.Osativa.MSU.1) 
str_count(as.character(sequences[[1]]), pattern = 'AAAG') + str_count(as.character(sequences[[1]]), pattern = 'CTTT')
motif_bundle<-rbind(as.data.frame(runFimo(sequences, Dof2_motif,thresh=0.005)))
dim(motif_bundle)[1]

# NADP-ME 
p<- CoveragePlot(object = small_plot_rice,region = "LOC-Os01g09320",expression.assay = "RNA",extend.upstream = -4500,extend.downstream = 400) 
p & scale_fill_manual(values = c("green4","blue4")) 
subset(DE_Peaks_Average_Rice, DE_Peaks_Average_Rice$gene_id == "LOC-Os01g09320")
peak_1<-"Chr1:4743653-4744777" 
sequences <- peak_1 %>%  get_sequence(BSgenome.Osativa.MSU.1) 
str_count(as.character(sequences[[1]]), pattern = 'AAAG') + str_count(as.character(sequences[[1]]), pattern = 'CTTT')
motif_bundle<-rbind(as.data.frame(runFimo(sequences, Dof2_motif,thresh=0.005)))
dim(motif_bundle)[1]
