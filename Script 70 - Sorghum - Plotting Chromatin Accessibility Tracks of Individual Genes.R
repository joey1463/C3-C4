# This script will plot accessbile chromatin tracks within mesophyll and bundle sheath cell types for candidate genes in Sorghu,

# Read in Data
setwd("~/Desktop/Salk/Project - Light/R files")
load("small_plot_sorghum.RData")
Idents(small_plot_sorghum)<-substr(Idents(small_plot_sorghum),1,9)
Idents(small_plot_sorghum) <- factor(small_plot_sorghum@active.ident, c("mesophyll","bundle_sh"))
load("L1_cis_sig_list_cross_species.RData")
Dof2_motif <- MotifDb::MotifDb %>%  MotifDb::query(cis_sig_list_identity[[5]][22]) %>%  universalmotif::convert_motifs() %>%  .[[5]]  

# Calling All Peaks
load("L1_sorghum_multiome_Average_Peaks.RData")
load("L1_sorghum_multiome_RNA_combined_with_ATAC.RData")
Average_Peaks<-Average_Peaks[,order(colnames(Average_Peaks))]
Peaks_Close_Feature<-ClosestFeature(sorghum_multiome_RNA_combined, rownames(sorghum_multiome_RNA_combined))
Peaks_Close_Feature<-subset(Peaks_Close_Feature, Peaks_Close_Feature$distance<1500)
DE_Peaks_Average_sorghum<-merge(Average_Peaks,Peaks_Close_Feature, by.x='row.names', by.y='query_region')
DE_Peaks_Average_sorghum$transcript_id<-str_replace(DE_Peaks_Average_sorghum$transcript_id,"_","-")



# GAPDH 
p<- CoveragePlot(object = small_plot_sorghum,region = "Sobic.006G105900.v3.2",expression.assay = "RNA",extend.upstream = 2000,extend.downstream = -2000) 
p<- CoveragePlot(object = small_plot_sorghum,region = "Sobic.006G105900.v3.2",expression.assay = "RNA",extend.upstream = 2000,extend.downstream = -500) 
p & scale_fill_manual(values = c("green4","blue4"))
subset(DE_Peaks_Average_sorghum, DE_Peaks_Average_sorghum$gene_id  == "Sobic.006G105900.v3.2")
peak_1<-"Chr06:47594873-47595180" # 1000 bp upstream
sequences <- peak_1 %>%  get_sequence(BSgenome.SBicolor.454.2) 
str_count(as.character(sequences[[1]]), pattern = 'AAAG') + str_count(as.character(sequences[[1]]), pattern = 'CTTT')
motif_bundle<-rbind(as.data.frame(runFimo(sequences,Dof2_motif,thresh=0.005)))
dim(motif_bundle)[1]
peak_2<-"Chr06:47595429-47596593" # in gene
sequences <- peak_2 %>%  get_sequence(BSgenome.SBicolor.454.2) 
str_count(as.character(sequences[[1]]), pattern = 'AAAG') + str_count(as.character(sequences[[1]]), pattern = 'CTTT')
motif_bundle<-rbind(as.data.frame(runFimo(sequences, Dof2_motif,thresh=0.005)))
dim(motif_bundle)[1]
peak_3<-"Chr06:47596835-47597499" # in gene
sequences <- peak_3 %>%  get_sequence(BSgenome.SBicolor.454.2) 
str_count(as.character(sequences[[1]]), pattern = 'AAAG') + str_count(as.character(sequences[[1]]), pattern = 'CTTT')
motif_bundle<-rbind(as.data.frame(runFimo(sequences, Dof2_motif,thresh=0.005)))
dim(motif_bundle)[1]

# NADPME 
p<- CoveragePlot(object = small_plot_sorghum,region = "Sobic.003G036200.v3.2",expression.assay = "RNA",extend.upstream = 500,extend.downstream = -3500) # NADPME
p & scale_fill_manual(values = c("green4","blue4"))
subset(DE_Peaks_Average_sorghum, DE_Peaks_Average_sorghum$gene_id  == "Sobic.003G036200.v3.2")
peak_1<-"Chr03:3316880-3317834 " # at promoter and into gene
sequences <- peak_1 %>%  get_sequence(BSgenome.SBicolor.454.2) 
str_count(as.character(sequences[[1]]), pattern = 'AAAG') + str_count(as.character(sequences[[1]]), pattern = 'CTTT')
motif_bundle<-rbind(as.data.frame(runFimo(sequences, Dof2_motif,thresh=0.005)))
dim(motif_bundle)[1]
