# This script will make a Sankey plot displaying how cell-type marker genes compare across species. 

# Plot Table
setwd("~/Desktop/Salk/Project - Light/R files")
load("L3_Orthology_output_fisher.RData")
load("L3_Orthology_output_positive.RData")
load("L3_Orthology_rice_orthogroups.RData")
load("L3_Orthology_sorghum_orthogroups.RData")

# Sankey Plot
cell_types<-c("mesophyll","bundle_sheath","guard","epidermis","phloem","xylem")
rownames(output_positive)<-cell_types
colnames(output_positive)<-cell_types
plot_matrix<-melt(output_positive)
plot_matrix[,4]<-c(rep(2,6),rep(3,6),rep(1,6),rep(1,6),rep(1,6),rep(1,6))
colnames(plot_matrix)<-c("rice","sorghum","number","color")
ggplot(as.data.frame(plot_matrix),
       aes(y = number, axis1 = factor(rice, level=unique(plot_matrix$rice)[c(4,3,1,2,5,6)]), axis2 = factor(sorghum, level=unique(plot_matrix$sorghum)[c(4,3,1,2,5,6)]))) +
  geom_alluvium(aes(fill =  as.factor(color)), width = 1/12) +  # fill = Admit - TO COLOR LABELS as.factor(color)
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Rice", "Sorghum"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") + theme_minimal()
ggplot(as.data.frame(plot_matrix),
       aes(y = number, axis1 = factor(rice, level=unique(plot_matrix$rice)[c(4,3,1,2,5,6)]), axis2 = factor(sorghum, level=unique(plot_matrix$sorghum)[c(4,3,1,2,5,6)]))) +
  geom_alluvium(aes(fill =  as.factor(color)), width = 1/12)  + geom_stratum(width = 1/12, fill = "grey50", color = "grey") + scale_fill_manual(values=c("grey79", "green4", "blue4")) + theme_bw()

# Exporting Specific Sankey Lists
new_sorghum_orthogroups<-list(NA)
for (i in 1:36) {new_sorghum_orthogroups[[i]]<-sorghum_orthogroups[[i+1]]}

new_rice_orthogroups<-list(NA)
for (i in 1:36) {new_rice_orthogroups[[i]]<-rice_orthogroups[[i+1]]}

cell_types<-c("mesophyll","bundle sheath","guard","epidermis","phloem","xylem")
list_table<-as.data.frame(cbind(c(rep(cell_types[1],6),rep(cell_types[2],6),rep(cell_types[3],6),rep(cell_types[4],6),rep(cell_types[5],6),rep(cell_types[6],6)),cell_types))
colnames(list_table)<-c("sorghum","rice")

# OK so now you want to pull out a specific set
subset(list_table, list_table$rice == 'guard' &list_table$sorghum == 'guard') # get list number

setwd("~/Desktop")
for (i in 1:36){
name<-paste(list_table[i,2],list_table[i,1],sep='_')
output_table<-subset(new_sorghum_orthogroups[[i]], new_sorghum_orthogroups[[i]]$Orthogroup %in% new_rice_orthogroups[[i]]$Orthogroup)
write.table(output_table, file = paste(name,".txt",sep=''),quote=F,sep="\t")}
