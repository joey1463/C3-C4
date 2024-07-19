# This script will take unnormalized pseduobulked transcriptional profiles from each cell type, and perform 
# an ANCOVA analysis using cell type and light response as factorial and quantitative variables respectively.
# The ANCOVA analysis is performed on each possible pair of cell types. 

Iterate_ANCOVA<-function(cell_type_1, cell_type_2){

time_aggregate_cell_type_1<-time_aggregate[,str_detect(colnames(time_aggregate), cell_type_1)]
time_aggregate_cell_type_2<-time_aggregate[,str_detect(colnames(time_aggregate), cell_type_2)]

time_aggregate_subset<-cbind(time_aggregate_cell_type_1,time_aggregate_cell_type_2)
time_aggregate_subset<-subset(time_aggregate_subset, rownames(time_aggregate_subset) %in% rownames(time_average))
maxes<-rowSums(time_aggregate_subset) 
time_aggregate_subset<-time_aggregate_subset[maxes>100,] 
time_aggregate_subset<-time_aggregate_subset[rowSums(time_aggregate_subset == 0) < 6, ] # restrict drop out

time_points<-c(0,0.5,1,2,4,6,8,12)
Conditions<-data.frame(time= rep(time_points,2), tissue=c(rep("cell_type_1",8),rep("cell_type_2",8)),row.names = colnames(time_aggregate_subset))
print(Conditions)
cutoff<-0.05 

# Determine if signficant in response to both Tissue and Time variables
DEmat<-DESeq(DESeqDataSetFromMatrix(countData = time_aggregate_subset,colData = Conditions, design = ~ time + tissue)) 
Significance_Values<-as.data.frame(cbind(results(DESeq(DEmat, test='LRT', reduced = ~ time))[,c(2,6)],
                                         results(DESeq(DEmat, test='LRT', reduced = ~ tissue))[,c(2,6)]))
colnames(Significance_Values)<-c('Tissue_FC','Tissue_pvalue','Time_FC','Time_pvalue')
rownames(Significance_Values)<-rownames(time_aggregate_subset)
Significance_Values[is.na(Significance_Values)] <- 1
Tissue_Time_Significant<-subset(Significance_Values, Significance_Values$Tissue_pvalue<cutoff & Significance_Values$Time_pvalue<cutoff)

# Tissue only
DEmat<-DESeq(DESeqDataSetFromMatrix(countData = time_aggregate_subset,colData = Conditions, design = ~ tissue))
Tissue_Significance_Values<-as.data.frame(results(DESeq(DEmat, test='LRT', reduced = ~ 1)))[,c(2,6)]
Tissue_Significance_Values[is.na(Tissue_Significance_Values)] <- 1
colnames(Tissue_Significance_Values)<-c('Tissue_FC','Tissue_pvalue')
rownames(Tissue_Significance_Values)<-rownames(time_aggregate_subset)

# Time only
DEmat<-DESeq(DESeqDataSetFromMatrix(countData = time_aggregate_subset,colData = Conditions, design = ~ time))
Time_Significance_Values<-as.data.frame(results(DESeq(DEmat, test='LRT', reduced = ~ 1)))[,c(2,6)]
Time_Significance_Values[is.na(Time_Significance_Values)] <- 1
colnames(Time_Significance_Values)<-c('Time_FC','Time_pvalue')
rownames(Time_Significance_Values)<-rownames(time_aggregate_subset)

# Place Tissue or Time FC and p-values onto those that are significant for both
Tissue_or_Time_1<-cbind(Tissue_Significance_Values,Time_Significance_Values)
Tissue_or_Time_2<-merge(Tissue_or_Time_1, Tissue_Time_Significant,by="row.names")[,-1]
rownames(Tissue_or_Time_2)<-merge(Tissue_or_Time_1, Tissue_Time_Significant,by="row.names")[,1]
Tissue_Time_Significant_updated<-Tissue_or_Time_2[,c(1,6,3,8)]
colnames(Tissue_Time_Significant_updated)<-c('Tissue_FC','Tissue_pvalue','Time_FC','time_pvalue')

# Add in those genes that are not dynamic with respct to time 
Tissue_Significant_updated<-subset(Tissue_Significance_Values, Tissue_Significance_Values$Tissue_pvalue < cutoff)
Tissue_Significant_updated<-Tissue_Significant_updated[!(rownames(Tissue_Significant_updated) %in% rownames(Tissue_Time_Significant_updated)), ]

# Make final export
Sorghum_ANCOVA_all_Tissue<-unique(c(rownames(Tissue_Significant_updated),rownames(Tissue_Time_Significant_updated)))
Sorghum_ANCOVA_all_Tissue}

cell_types<-c("mesophyll","bundle_sheath","guard","epidermis","phloem","xylem")

possible_pairs<-combn(cell_types,2)

Sorghum_Pair_Permutations<-list(NA)
for (i in 1:15){ Sorghum_Pair_Permutations[[i]]<-Iterate_ANCOVA(possible_pairs[1,i],possible_pairs[2,i])}

setwd("~/Desktop/Salk/Project - Light/R files")
save(Sorghum_Pair_Permutations, file = "L3_Sorghum_Pair_Permutations.RData")
