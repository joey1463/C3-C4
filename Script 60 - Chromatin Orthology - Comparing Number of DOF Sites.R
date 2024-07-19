# This script will read in the number of DOF sites found associated with each differentially partitioned gene, and then compare 
# these counts across species. Enrichment in one species vs. another is computed using a binomial test. 

setwd("~/Desktop/Salk/Project - Light/R files")
load("Target_Genes.RData")
load("DOF_count_by_gene_Sorghum.RData")
load("DOF_count_by_gene_Rice.RData")
rownames(DOF_count_by_gene_Rice)<-str_replace(rownames(DOF_count_by_gene_Rice),"-","_")

DOF_comparison<-as.data.frame(matrix(NA, nrow=length(Target_Genes[,1]),ncol=5))

for (i in 1:length(Target_Genes[,1])){
DOF_comparison[i,1]<-DOF_count_by_gene_Sorghum[ Target_Genes[i,'Sorghum'],2]
DOF_comparison[i,2]<-DOF_count_by_gene_Rice[ Target_Genes[i,'Rice'],2]}
DOF_comparison[,3]<-Target_Genes$Sorghum
DOF_comparison[,4]<-Target_Genes$Rice
DOF_comparison[,5]<-Target_Genes$Pattern
colnames(DOF_comparison)<-c("sorghum_DOF_count","rice_DOF_count","sorghum_gene","rice_gene","pattern")
DOF_comparison$rice_DOF_count[is.na(DOF_comparison$rice_DOF_count)] <- 0
DOF_comparison$sorghum_DOF_count[is.na(DOF_comparison$sorghum_DOF_count)] <- 0

subset(DOF_comparison, DOF_comparison$sorghum_gene == "Sobic.009G189400") # SiR, "LOC-Os05g42350" S - 13, R - 19
subset(DOF_comparison, DOF_comparison$sorghum_gene == "Sobic.003G036200") # NADPME, "LOC-Os01g09320" S - 3, R - 3
subset(DOF_comparison, DOF_comparison$sorghum_gene == "Sobic.006G105900") # GAPDH,  "LOC_Os04g38600" & "LOC_Os03g03720" S - 9, R - 9 & 1
subset(DOF_comparison, DOF_comparison$sorghum_gene == "Sobic.005G056400") # Aldolase,  "LOC_Os11g07020" S - 32, R - 21

summary(subset(DOF_comparison, DOF_comparison$pattern =="differential")$sorghum_DOF_count > subset(DOF_comparison, DOF_comparison$pattern =="differential")$rice_DOF_count)
binom.test(34, (34+14), p = 0.5, alternative = c("two.sided"), conf.level = 0.99)

summary(subset(DOF_comparison, DOF_comparison$pattern =="consistent")$sorghum_DOF_count > subset(DOF_comparison, DOF_comparison$pattern =="consistent")$rice_DOF_count)
binom.test(44, (44+52), p = 0.5, alternative = c("two.sided"), conf.level = 0.99)

# Plot Diffrential scatterplot
DOF_comparison_diff<-subset(DOF_comparison, DOF_comparison$pattern %in% c("differential"))
p1<- ggplot(data=DOF_comparison_diff, aes(x=rice_DOF_count, y=sorghum_DOF_count, group=pattern))  + geom_point(size=2,aes(color=as.factor(pattern))) + theme_bw() + scale_color_manual(values = c("red3","red3")) + geom_abline(slope=1,  intercept=0) + xlim(c(-0,40)) + ylim(c(0,40)) + xlab("") + ylab("") + theme(legend.position = 'none') + theme(axis.text=element_text(size=15),axis.title=element_text(size=15)) + theme(axis.text=element_text(color="black"))
p1

# Plot Consistent scatterplot
DOF_comparison_consis<-subset(DOF_comparison, DOF_comparison$pattern %in% c("consistent"))
p2<- ggplot(data=DOF_comparison_consis, aes(x=rice_DOF_count, y=sorghum_DOF_count, group=pattern))  + geom_point(size=2,aes(color=as.factor(pattern))) + theme_bw() + scale_color_manual(values = c("blue3","blue3")) + geom_abline(slope=1,  intercept=0) + xlim(c(-0,40)) + ylim(c(0,40)) +  xlab("") + ylab("") + theme(legend.position = 'none') + theme(axis.text=element_text(size=15),axis.title=element_text(size=15))

grid.arrange(p1, p2, nrow=1)

# Export
setwd("~/Desktop/")
write.table(DOF_comparison,file="DOF_comparison.txt",quote=F,sep="\t")
