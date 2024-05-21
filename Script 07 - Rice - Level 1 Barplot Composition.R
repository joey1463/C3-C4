# This script will make barplots that compare cluster representation with assay type and time point. 

raw_table<-table(as.character(rice_combined@active.ident), paste(rice_combined$assay_type, rice_combined$light_time, sep="_"))
normalized_table<-matrix(NA,nrow=nrow(raw_table),ncol=ncol(raw_table))
rownames(normalized_table)<-rownames(raw_table)
colnames(normalized_table)<-colnames(raw_table)
for (i in 1:ncol(raw_table)) {normalized_table[,i]<-as.numeric(raw_table[,i]/sum(raw_table[,i]))}

normalized_table_1<-melt(normalized_table)
colnames(normalized_table_1)<-c("cluster","day","value")
normalized_table_1[,1]<-as.factor(normalized_table_1[,1])

b1 <- ggplot(normalized_table_1, aes(fill=cluster, y=value, x=as.factor(day))) + geom_bar(position="stack", stat="identity")
b2 <- b1 + labs(x = "") + labs(y = "") +  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
b3 <- b2 + theme(axis.text=element_text(size=12),axis.title=element_text(size=12)) 
b3
