# We computed the genes that were bundle sheath partitioned in both rice and sorghum (consistently partitioned). 
# This script reads in the orthologs of these genes from other C3 grasses -  Barley, Brachy and Claxum orthologs - and uses
# AME from the MEME suite to assess cis-regulatory enrichment within their promoters. 

# Read In Required Data
setwd("~/Desktop/Salk/Project - Light/R files")

Claxum_Targets<-read.table("Claxum_Targets.txt",header=TRUE)
Hvulgare_Targets<-read.table("Hvulgare_Targets.txt",header=TRUE)
Bdistachyon_Targets<-read.table("Bdistachyon_Targets.txt",header=TRUE)

Claxum_gff3<-read.table("~/Desktop/Salk/Project - Light/R files/Claxum_690_v1.1.gene.txt")
Osativa_gff3 <- read.table('~/Desktop/Salk/Project - Light/R files/MSU_rice_onlygene_edit.gtf')
Sbicolor_gff3 <- read.table('~/Desktop/Salk/Project - Light/R files/Sbicolor_454_v3.1.1.gene_exons_genes_only.gtf')
Hvulgare_gff3 <- read.table('~/Desktop/Salk/Project - Light/R files/Hvulgare_462_r1.gene_exons.gff3')
Bdistachyon_gff3 <- read.table('~/Desktop/Salk/Project - Light/R files/Bdistachyon_314_v3.1.gene_exons.gff3')

load("L1_cis_sig_list_cross_species.RData")
motif_file <- "~/Desktop/Salk/Project - Light/R files/JASPAR_plant_core.meme"
breadth<-999

# Claxum 
Claxum_Hits<-subset(Claxum_gff3, Claxum_gff3[,9] %in% Claxum_Targets$Ortholog_ID)
Claxum_Hits<-subset(Claxum_Hits, Claxum_Hits[,3] == "mRNA")
Claxum_Hits<-Claxum_Hits[!duplicated(Claxum_Hits[,9]),]
rownames(Claxum_Hits)<-Claxum_Hits[,9]

Claxum_Hits_positive<-subset(Claxum_Hits, Claxum_Hits[,7] == "+")
ranges<- paste(Claxum_Hits_positive[,4]-breadth,Claxum_Hits_positive[,4],sep="-")
ranges<- paste(Claxum_Hits_positive[,1],ranges,sep=":")
Claxum_Hits_positive[,11]<-ranges

Claxum_Hits_negative<-subset(Claxum_Hits, Claxum_Hits[,7] == "-")
ranges<- paste(Claxum_Hits_negative[,5],Claxum_Hits_negative[,5]+breadth,sep="-") 
ranges<- paste(Claxum_Hits_negative[,1],ranges,sep=":")
Claxum_Hits_negative[,11]<-ranges

Claxum_Hits <- rbind(Claxum_Hits_positive, Claxum_Hits_negative)
sequences <- Claxum_Hits[,11] %>% get_sequence(Claxum.version.1.0) 
Claxum_outcome <- as.data.frame(runAme(sequences, control = "shuffle", database = motif_file))
Claxum_outcome$Dof<-ifelse(Claxum_outcome$motif_id %in% cis_sig_list_identity[[5]], "2", "1")
Claxum_outcome<-subset(Claxum_outcome, Claxum_outcome$adj.pvalue<0.01)
Claxum_outcome$adj.pvalue_log<- -log(Claxum_outcome$adj.pvalue, 10)
Claxum_plot<-Claxum_outcome[,c(3,18,19)]
Claxum_plot$species<-"C_Claxum"


# Rice
consistent_genes<-subset(DOF_comparison, DOF_comparison$pattern == 'consistent') # read in from previous script
Sbicolor_Hits<-subset(Sbicolor_gff3, substr(Sbicolor_gff3[,13],1,16) %in% consistent_genes$sorghum_gene)
Osativa_Hits<-subset(Osativa_gff3, Osativa_gff3[,13] %in% str_replace(consistent_genes$rice_gene,"_","-"))
Osativa_Hits<-subset(Osativa_Hits, Osativa_Hits[,3] == "transcript")
Osativa_Hits<-Osativa_Hits[!duplicated(Osativa_Hits[,13]),]
rownames(Osativa_Hits)<-Osativa_Hits[,13]

Osativa_Hits_positive<-subset(Osativa_Hits, Osativa_Hits[,7] == "+")
ranges<- paste(Osativa_Hits_positive[,4]-breadth,Osativa_Hits_positive[,4],sep="-") 
ranges<- paste(Osativa_Hits_positive[,1],ranges,sep=":")
Osativa_Hits_positive[,11]<-ranges

Osativa_Hits_negative<-subset(Osativa_Hits, Osativa_Hits[,7] == "-")
ranges<- paste(Osativa_Hits_negative[,5],Osativa_Hits_negative[,5]+breadth,sep="-") 
ranges<- paste(Osativa_Hits_negative[,1],ranges,sep=":")
Osativa_Hits_negative[,11]<-ranges

Osativa_Hits <- rbind(Osativa_Hits_positive, Osativa_Hits_negative)
sequences <- Osativa_Hits[,11] %>% get_sequence(BSgenome.Osativa.MSU.1) 
Osativa_outcome <- as.data.frame(runAme(sequences, control = "shuffle", database = motif_file))
Osativa_outcome$Dof<-ifelse(Osativa_outcome$motif_id %in% cis_sig_list_identity[[5]], "2", "1")
Osativa_outcome<-subset(Osativa_outcome, Osativa_outcome$adj.pvalue<0.01)
Osativa_outcome$adj.pvalue_log<- -log(Osativa_outcome$adj.pvalue, 10)
Osativa_plot<-Osativa_outcome[,c(3,18,19)]
Osativa_plot$species<-"A_Osativa"




# Sorghum
consistent_genes<-subset(DOF_comparison, DOF_comparison$pattern == 'consistent') # read in from previous script
Sbicolor_Hits<-subset(Sbicolor_gff3, substr(Sbicolor_gff3[,13],1,16) %in% consistent_genes$sorghum_gene)
Sbicolor_Hits<-subset(Sbicolor_Hits, Sbicolor_Hits[,3] == "transcript")
Sbicolor_Hits<-Sbicolor_Hits[!duplicated(Sbicolor_Hits[,13]),]
rownames(Sbicolor_Hits)<-Sbicolor_Hits[,13]

Sbicolor_Hits_positive<-subset(Sbicolor_Hits, Sbicolor_Hits[,7] == "+")
ranges<- paste(Sbicolor_Hits_positive[,4]-breadth,Sbicolor_Hits_positive[,4],sep="-")
ranges<- paste(Sbicolor_Hits_positive[,1],ranges,sep=":")
Sbicolor_Hits_positive[,11]<-ranges

Sbicolor_Hits_negative<-subset(Sbicolor_Hits, Sbicolor_Hits[,7] == "-")
ranges<- paste(Sbicolor_Hits_negative[,5],Sbicolor_Hits_negative[,5]+breadth,sep="-")
ranges<- paste(Sbicolor_Hits_negative[,1],ranges,sep=":")
Sbicolor_Hits_negative[,11]<-ranges

Sbicolor_Hits <- rbind(Sbicolor_Hits_positive, Sbicolor_Hits_negative)
sequences <- Sbicolor_Hits[,11] %>% get_sequence(BSgenome.SBicolor.454.2) 
Sbicolor_outcome <- as.data.frame(runAme(sequences, control = "shuffle", database = motif_file))
Sbicolor_outcome$Dof<-ifelse(Sbicolor_outcome$motif_id %in% cis_sig_list_identity[[5]], "2", "1")
Sbicolor_outcome<-subset(Sbicolor_outcome, Sbicolor_outcome$adj.pvalue<0.01)
Sbicolor_outcome$adj.pvalue_log<- -log(Sbicolor_outcome$adj.pvalue, 10)
Sbicolor_plot<-Sbicolor_outcome[,c(3,18,19)]
Sbicolor_plot$species<-"B_Sbicolor"

# Barley
Hvulgare_Hits<-subset(Hvulgare_gff3, substr(Hvulgare_gff3[,9],4,19)  %in% Hvulgare_Targets$Ortholog_ID)
Hvulgare_Hits<-subset(Hvulgare_Hits, Hvulgare_Hits[,3] == "mRNA")
Hvulgare_Hits<- Hvulgare_Hits %>% filter(grepl("longest=1", V9))
Hvulgare_Hits<-Hvulgare_Hits[!duplicated(Hvulgare_Hits[,9]),]
rownames(Hvulgare_Hits)<-substr(Hvulgare_Hits[,9],4,19) 

Hvulgare_Hits_positive<-subset(Hvulgare_Hits, Hvulgare_Hits[,7] == "+")
ranges<- paste(Hvulgare_Hits_positive[,4]-breadth,Hvulgare_Hits_positive[,4],sep="-") 
ranges<- paste(Hvulgare_Hits_positive[,1],ranges,sep=":")
Hvulgare_Hits_positive[,10]<-ranges

Hvulgare_Hits_negative<-subset(Hvulgare_Hits, Hvulgare_Hits[,7] == "-")
ranges<- paste(Hvulgare_Hits_negative[,5],Hvulgare_Hits_negative[,5]+breadth,sep="-") 
ranges<- paste(Hvulgare_Hits_negative[,1],ranges,sep=":")
Hvulgare_Hits_negative[,10]<-ranges

Hvulgare_Hits <- rbind(Hvulgare_Hits_positive, Hvulgare_Hits_negative)
sequences <- Hvulgare_Hits[,10] %>% get_sequence(Hvulgare.version.1.0) 
Hvulgare_outcome <- as.data.frame(runAme(sequences, control = "shuffle", database = motif_file))
Hvulgare_outcome$Dof<-ifelse(Hvulgare_outcome$motif_id %in% cis_sig_list_identity[[5]], "2", "1")
Hvulgare_outcome<-subset(Hvulgare_outcome, Hvulgare_outcome$adj.pvalue<0.01)
Hvulgare_outcome$adj.pvalue_log<- -log(Hvulgare_outcome$adj.pvalue, 10)
Hvulgare_plot<-Hvulgare_outcome[,c(3,18,19)]
Hvulgare_plot$species<-"D_Hvulgare"


# Brachy
Bdistachyon_Hits<-subset(Bdistachyon_gff3, substr(Bdistachyon_gff3[,9],4,15)  %in% Bdistachyon_Targets$Ortholog_ID)
Bdistachyon_Hits<-subset(Bdistachyon_Hits, Bdistachyon_Hits[,3] == "mRNA")
Bdistachyon_Hits<- Bdistachyon_Hits %>% filter(grepl("longest=1", V9))
Bdistachyon_Hits<-Bdistachyon_Hits[!duplicated(Bdistachyon_Hits[,9]),]
rownames(Bdistachyon_Hits)<-substr(Bdistachyon_Hits[,9],4,15) 

Bdistachyon_Hits_positive<-subset(Bdistachyon_Hits, Bdistachyon_Hits[,7] == "+")
ranges<- paste(Bdistachyon_Hits_positive[,4]-breadth,Bdistachyon_Hits_positive[,4],sep="-")
ranges<- paste(Bdistachyon_Hits_positive[,1],ranges,sep=":")
Bdistachyon_Hits_positive[,10]<-ranges

Bdistachyon_Hits_negative<-subset(Bdistachyon_Hits, Bdistachyon_Hits[,7] == "-")
ranges<- paste(Bdistachyon_Hits_negative[,5],Bdistachyon_Hits_negative[,5]+breadth,sep="-")
ranges<- paste(Bdistachyon_Hits_negative[,1],ranges,sep=":")
Bdistachyon_Hits_negative[,10]<-ranges

Bdistachyon_Hits <- rbind(Bdistachyon_Hits_positive, Bdistachyon_Hits_negative)
sequences <- Bdistachyon_Hits[,10] %>% get_sequence(BDistachyon.version.1.0) 
Bdistachyon_outcome <- as.data.frame(runAme(sequences, control = "shuffle", database = motif_file))
Bdistachyon_outcome$Dof<-ifelse(Bdistachyon_outcome$motif_id %in% cis_sig_list_identity[[5]], "2", "1")
Bdistachyon_outcome<-subset(Bdistachyon_outcome, Bdistachyon_outcome$adj.pvalue<0.01)
Bdistachyon_outcome$adj.pvalue_log<- -log(Bdistachyon_outcome$adj.pvalue, 10)
Bdistachyon_plot<-Bdistachyon_outcome[,c(3,18,19)]
Bdistachyon_plot$species<-"E_Bdistachyon"


# Plot All Together
plot_matrix<-rbind(Claxum_plot,Osativa_plot,Sbicolor_plot,Hvulgare_plot,Bdistachyon_plot)
ggplot(plot_matrix, aes(x=species, y=adj.pvalue_log, color = Dof)) + geom_point(aes(size=Dof)) + theme_bw()+ scale_color_manual(values = c("grey80","red3")) +theme(legend.position = 'none') + scale_size_manual(values = c(2,4))+ theme(axis.text=element_text(size=12),axis.title=element_text(size=12))  + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Export
#setwd("~/Desktop/")
#write.table(Claxum_outcome,file="Claxum_outcome_elements.txt",quote=F,sep="\t")
#write.table(Osativa_outcome,file="Osativa_outcome_elements.txt",quote=F,sep="\t")
#write.table(Sbicolor_outcome,file="Sbicolor_outcome_elements.txt",quote=F,sep="\t")
#write.table(Hvulgare_outcome,file="Hvulgare_outcome_elements.txt",quote=F,sep="\t")
#write.table(Bdistachyon_outcome,file="Bdistachyon_outcome_elements.txt",quote=F,sep="\t")
