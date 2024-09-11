# This script will read in the cis-regulatory elements that are enriched in a cell-type specific or light-specific manner, 
# and overlap these cis-elements across species. The statistical significance of this overlap is determined using a Fisher test. 

setwd("~/Desktop/Salk/Project - Light/R files")
load("rice_cis_elements.RData")
load("sorghum_cis_elements.RData")
load("rice_cis_elements_light.RData")
load("sorghum_cis_elements_light.RData")
cell_types<-rev(c("mesophyll","bundle_sheath","phloem","xylem","epidermis","guard"))


call_cis_fisher<-function(rice_list, sorghum_list){
  rice_list<-subset(rice_list, rice_list$p_val_adj<0.05)
  sorghum_list<-subset(sorghum_list, sorghum_list$p_val_adj<0.05)
  rice_list<-rice_list[1:25,] # take top 25
  sorghum_list<-sorghum_list[1:25,] # take top 25
  intersect<-sum(as.numeric((rownames(rice_list) %in% rownames(sorghum_list))))
  rice_length<-length(rice_list[,1])
  sorghum_length<-length(sorghum_list[,1])
  matrix_to_test<-matrix(c(intersect,(sorghum_length-intersect),(rice_length-intersect),(530)),2,2) # you have 530 cis elements
  p_value<-as.numeric(fisher.test(matrix_to_test, alternative='greater')[1])}

# Identity Cis-Element Grid Plot
cis_sig_table_identity<-matrix(NA, nrow=6, ncol=6)
for (a in 1:6) {  for (b in 1:6) { cis_sig_table_identity[a,b]<- call_cis_fisher(rice_cis_elements[[a]],sorghum_cis_elements[[b]])  }}
output_fisher<-cis_sig_table_identity
rownames(output_fisher)<-cell_types
colnames(output_fisher)<-cell_types
output_fisher<-output_fisher[,rev(colnames(output_fisher))]
plot_matrix<-melt(output_fisher)
plot_matrix$value<-p.adjust(plot_matrix$value, method = 'fdr')
plot_matrix$value<- -log(plot_matrix$value,10)
plot_matrix$value<- ifelse(plot_matrix$value > 10 , 10, plot_matrix$value) # Set upper limit identity
plot_matrix[,4]<-c(rep(2,6),rep(3,6),rep(1,6),rep(1,6),rep(1,6),rep(1,6))
colnames(plot_matrix)<-c("rice","sorghum","log_10_p","color")
g1<- qplot(sorghum, rice, fill=log_10_p, data=plot_matrix, geom='tile') + theme(axis.text = element_text(angle = 90, vjust = 0.5, hjust=1))
g1 + scale_fill_gradient(low="grey100", high="hotpink4") + labs(x = 'sorghum', y = 'rice') + theme_bw()



# Light Cis-Element Grid Plot
cis_sig_table_light<-matrix(NA, nrow=6, ncol=6)
for (a in 1:6) { for (b in 1:6) {cis_sig_table_light[a,b]<- call_cis_fisher(rice_cis_elements_light[[a]],sorghum_cis_elements_light[[b]])  }}
output_fisher<-cis_sig_table_light
rownames(output_fisher)<-cell_types
colnames(output_fisher)<-cell_types
output_fisher<-output_fisher[,rev(colnames(output_fisher))]
plot_matrix<-melt(output_fisher)
plot_matrix$value<-p.adjust(plot_matrix$value, method = 'fdr')
plot_matrix$value<- -log(plot_matrix$value,10)
plot_matrix$value<- ifelse(plot_matrix$value > 20 , 20, plot_matrix$value)  # Set upper limit light
plot_matrix[,4]<-c(rep(2,6),rep(3,6),rep(1,6),rep(1,6),rep(1,6),rep(1,6))
colnames(plot_matrix)<-c("rice","sorghum","log_10_p","color")
g1<- qplot(sorghum, rice, fill=log_10_p, data=plot_matrix, geom='tile') + theme(axis.text = element_text(angle = 90, vjust = 0.5, hjust=1))
g1 + scale_fill_gradient(low="grey100", high="hotpink4") + labs(x = 'sorghum', y = 'rice') + theme_bw()



## Print Cis-Elements themselves
call_cis_elements<-function(rice_list, sorghum_list){
  rice_list<-subset(rice_list, rice_list$p_val_adj<0.05)
  sorghum_list<-subset(sorghum_list, sorghum_list$p_val_adj<0.05)
  rice_list<-rice_list[1:25,]
  sorghum_list<-sorghum_list[1:25,]
  outcome<-subset(rice_list, rownames(rice_list) %in% rownames(sorghum_list))
  outcome<-rownames(outcome)
  outcome}

# Identity (Dont need nested forloop as simply going across diagonal)
cis_sig_list_identity<-list(NA)
for (a in 1:6) { cis_sig_list_identity[[a]]<- call_cis_elements(rice_cis_elements[[a]],sorghum_cis_elements[[a]]) }
#save(cis_sig_list_identity, file="L1_cis_sig_list_cross_species.RData")


# Light (Since there is no apparent cross species conservation, we print all cis elements together)
cis_sig_list_light<-list(NA)
for (a in 1:6) { for (b in 1:6) { cis_sig_list_light[[length(cis_sig_list_light)+1]]<- call_cis_elements(rice_cis_elements_light[[a]],sorghum_cis_elements_light[[b]])  }}
cis_sig_list_light<- cis_sig_list_light[-c(1, 23)]  # careful, need to remove empty list sets
cis_sig_list_all<-NA 
for (i in c(1:35)) {cis_sig_list_all[length(cis_sig_list_all):(length(cis_sig_list_all)+length(cis_sig_list_light[[i]])-1)]<- cis_sig_list_light[[i]]}
cis_sig_list_light<-unique(cis_sig_list_all)
#save(cis_sig_list_light, file="L1_cis_sig_list_light_dominant_cross_species.RData")
