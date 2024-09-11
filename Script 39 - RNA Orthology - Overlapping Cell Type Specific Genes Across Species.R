# This script will read in the cell-type specific marker genes for each species, and overlap these datasets across species. 
# The script will calculate the signficance of these overlaps using a Fisher test, and create a heatmap of the reuslt. 

# Load Orthology and Marker Sets
setwd("~/Desktop/Salk/Project - Light/R files")
orthology<-read.table('Rice_Sorghum_Orthology.txt',header=TRUE, sep="\t")
load("L3_Orthology_rice_markers.RData")
load("L3_Orthology_sorghum_markers.RData")

# Make 4-by-4 Matrix and Save Orthogroups in Overlap
output_positive<-matrix(NA,nrow=6, ncol=6)
output_fisher<-matrix(NA,nrow=6, ncol=6)
rice_orthogroups<-list(NA)
sorghum_orthogroups<-list(NA)

for (b in 1:6){
print(b)
sorghum_markers_all<-sorghum_markers[[b]]
sorghum_markers_sig<-subset(sorghum_markers_all, sorghum_markers_all$p_val_adj<0.01)
sorghum_cell_type<-rownames(sorghum_markers_sig) # put in cell type of interest
sorghum_orthology<-matrix(NA,nrow=length(sorghum_cell_type),ncol=3)
for (i in 1:length(sorghum_cell_type)) {sorghum_orthology[i,]<-as.character(orthology[grep(sorghum_cell_type[i], orthology$Sorghum),])}
rownames(sorghum_orthology)<-sorghum_cell_type
colnames(sorghum_orthology)[1:3]<-colnames(orthology)
sorghum_orthology<-as.data.frame(sorghum_orthology)
sorghum_orthology<-subset(sorghum_orthology, substr(sorghum_orthology$Orthogroup,1,1) =="O")
sorghum_orthology$Sorghum<-rownames(sorghum_orthology)

for (a in 1:6){
rice_markers_all<-rice_markers[[a]]
rice_markers_sig<-subset(rice_markers_all, rice_markers_all$p_val_adj<0.01)
rice_cell_type<-str_replace(rownames(rice_markers_sig) ,"-","_") # rice cell type
rice_orthology<-matrix(NA,nrow=length(rice_cell_type),ncol=3)
for (i in 1:length(rice_cell_type)) {rice_orthology[i,]<-as.character(orthology[grep(rice_cell_type[i], orthology$Rice),])}
rownames(rice_orthology)<-rice_cell_type
colnames(rice_orthology)[1:3]<-colnames(orthology)
rice_orthology<-as.data.frame(rice_orthology)
rice_orthology<-subset(rice_orthology, substr(rice_orthology$Orthogroup,1,1) =="O")

rice_orthogroups[[length(rice_orthogroups)+1]]<-subset(rice_orthology, rice_orthology$Orthogroup %in% sorghum_orthology$Orthogroup)
sorghum_orthogroups[[length(sorghum_orthogroups)+1]]<-subset(sorghum_orthology, sorghum_orthology$Orthogroup %in% rice_orthology$Orthogroup)

intersect<-sum(as.numeric(sorghum_orthology$Orthogroup %in% rice_orthology$Orthogroup))
sorghum_length<-length(sorghum_orthology$Orthogroup)
rice_length<-length(rice_orthology$Orthogroup)

matrix_to_test<-matrix(c(intersect,(sorghum_length-intersect),(rice_length-intersect),(14444)),2,2) # no. of orthogroups total is background
p_value<-as.numeric(fisher.test(matrix_to_test, alternative='greater')[1])
total<-length(as.numeric(sorghum_orthology$Orthogroup %in% rice_orthology$Orthogroup))

output_positive[a,b]<-intersect
output_fisher[a,b]<-p_value}}

# Heatmap on log2 p-value
plot_matrix <- output_fisher
cell_types<-c("mesophyll","bundle_sheath","guard","epidermis","phloem","xylem")
rownames(plot_matrix)<-cell_types
colnames(plot_matrix)<-cell_types
plot_matrix<-melt(plot_matrix)
plot_matrix$value<-p.adjust(plot_matrix$value, method = "fdr") # correct p-values
plot_matrix$value <- -log(plot_matrix$value,10)
plot_matrix$value <- ifelse(plot_matrix$value > 100 , 100, plot_matrix$value) # Set upper limit
plot_matrix[,4]<-c(rep(2,6),rep(3,6),rep(1,6),rep(1,6),rep(1,6),rep(1,6))
colnames(plot_matrix)<-c("rice","sorghum","log_10_p","color")
plot_matrix <- plot_matrix %>% mutate(value_rounded = sprintf("%.2f", log_10_p))
g1<- qplot(as.factor(sorghum), as.factor(rice), fill=log_10_p, data=plot_matrix, geom='tile') + theme(axis.text = element_text(angle = 90, vjust = 0.5, hjust=1))
g1 + scale_fill_gradient(low="grey100", high="hotpink4") + labs(x = 'sorghum', y = 'rice') + theme_bw() #+ geom_text(aes(label = value_rounded), size = 3, color = "black", vjust = 0.5, hjust = 0.5)


# Heatmap on number
plot_matrix <- output_positive
cell_types<-c("mesophyll","bundle_sheath","guard","epidermis","phloem","xylem")
rownames(plot_matrix)<-cell_types
colnames(plot_matrix)<-cell_types
plot_matrix<-melt(plot_matrix)
plot_matrix[,4]<-c(rep(2,6),rep(3,6),rep(1,6),rep(1,6),rep(1,6),rep(1,6))
colnames(plot_matrix)<-c("rice","sorghum","overlap","color")
g1<- qplot(as.factor(sorghum), as.factor(rice), fill=overlap, data=plot_matrix, geom='tile') + theme(axis.text = element_text(angle = 90, vjust = 0.5, hjust=1))
g1 + scale_fill_gradient(low="grey100", high="hotpink4") + labs(x = 'sorghum', y = 'rice') + theme_bw() + geom_text(aes(label = overlap), size = 3, color = "black", vjust = 0.5, hjust = 0.5)

save(output_positive, file = "L3_Orthology_output_positive.RData")
save(output_fisher, file = "L3_Orthology_output_fisher.RData")
save(rice_orthogroups, file = "L3_Orthology_rice_orthogroups.RData")
save(sorghum_orthogroups, file = "L3_Orthology_sorghum_orthogroups.RData")
