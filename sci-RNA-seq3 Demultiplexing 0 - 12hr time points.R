# This script will take read counts for 0-12h time point sci-RNA-seq3 experiments and assemble them into a counts matrix for Seurat. This does it for each species seperately. 

library(Matrix)
library(stringr)

### Read in Barcodes
# RT Barcodes - all barcodes are of length 10
RT_Barcodes<-read.table("/gale/raidix/rdx-7/jswift/R_CI/text_files/sc_RT_plate.txt",header=TRUE, sep="\t")
RT_Barcode<-NA
for (i in 1:384) {RT_Barcode[i]<-substring(as.character(RT_Barcodes[i,4]),22,31)}
RT_Barcodes[,5]<-RT_Barcode 

# LIG Barcodes - some LIG barcodes are 9 bp, others are 10bp
LIG_Barcodes<-read.table("/gale/raidix/rdx-7/jswift/R_CI/text_files/sc_ligation_plate.txt",header=TRUE, sep="\t")
LIG_Barcode<-NA
for (i in 1:384) {if (nchar(as.character(LIG_Barcodes[i,4])) == 53) {LIG_Barcode[i]<-substring(as.character(LIG_Barcodes[i,4]),44,53)} 
                  if (nchar(as.character(LIG_Barcodes[i,4])) == 51) {LIG_Barcode[i]<-substring(as.character(LIG_Barcodes[i,4]),43,51)}}
LIG_Barcodes[,5]<-LIG_Barcode

### Generate Matrix
generate_matrix<-function(gene_count_folder, name_start, name_stop){
  gene_matrix <- paste(gene_count_folder, "/count.MM", sep = "")
  df_gene <- paste(gene_count_folder, "/gene_name_annotate.txt", sep = "")
  df_cell <- paste(gene_count_folder, "/cell_annotate.txt", sep = "")

  # this is your entire gene list    
  df_gene <- read.csv(df_gene, header = F) 
  rownames(df_gene) <- df_gene$gene_id
  colnames(df_gene) <- c("gene_id", "gene_type", "exon_intron", "gene_name", "index")
  
  # this is your entire detected index list
  df_cell <- read.csv(df_cell, header = F)
  colnames(df_cell) <- c("sample", "index")
  rownames(df_cell) <- df_cell$sample 
  
  # this is your condensed matrix of counts
  gene_matrix <- read.csv(gene_matrix, header = F)
  colnames(gene_matrix)<- c("gene_name","index","count")
  remove<-c("nothing","nothing_intron")
  gene_matrix<-gene_matrix[!gene_matrix[,1] %in% remove, ]
  
  full_names<-as.character(gene_matrix$gene_name)
  full_genes<-NA
  for(i in 1:length(full_names)){full_genes[i] <- substring(full_names[i],name_start,name_stop)}
  gene_matrix[,1]<-as.factor(full_genes)
  gene_matrix<-gene_matrix[order(gene_matrix[,1], decreasing = FALSE),]
  gene_matrix[,4]<-gene_matrix[,1]
  gene_matrix[,1]<-as.numeric(gene_matrix$gene_name)
  colnames(gene_matrix)<-c("gene_identifier","index","count","gene_name")
  gene_matrix<-gene_matrix[complete.cases(gene_matrix), ]
  
  # generate your sparse matrix
  gene_count <- sparseMatrix(i = gene_matrix$gene_identifier, j = gene_matrix$index, x = gene_matrix$count)
  colnames(gene_count) <- rownames(df_cell) # which infact, these align with your df_cell
  rownames(gene_count) <- unique(gene_matrix$gene_name) # each row is in order of the gene 'identifier' which is in gene alphabetical order 
  gene_count}

rice_matrix_lists<-list(NA)
sorghum_matrix_lists<-list(NA)

for (i in 1:29) {rice_matrix_lists[[i]] <- generate_matrix(gene_count_folder=paste("/gale/netapp/home/jswift/analysis/CI_Light/CI_Light_rice_run1/split_",i,"/Arabidopsis_gene_count",sep=''), name_start=1,name_stop=14)}
for (i in 1:31) {sorghum_matrix_lists[[i]] <- generate_matrix(gene_count_folder=paste("/gale/netapp/home/jswift/analysis/CI_Light/CI_Light_sorghum_run1/split_",i,"/Arabidopsis_gene_count",sep=''), name_start=1,name_stop=16)}

print('no. nuclei/combinations in each rice split (this is raw count; before low count removal)')
for (i in 1:29) {print(length(rice_matrix_lists[[i]][1,]))}
print('no. nuclei/combinations in each sorghum split (this is raw count; before low count removal)')
for (i in 1:31) {print(length(sorghum_matrix_lists[[i]][1,]))}

# Make all rows the same
all_rice_gene_names <- NA
for (i in 1:29) {all_rice_gene_names <- c(all_rice_gene_names, rownames(rice_matrix_lists[[i]]))}
all_rice_gene_names <- unique(all_rice_gene_names[-1])

all_sorghum_gene_names <- NA
for (i in 1:31) {all_sorghum_gene_names <- c(all_sorghum_gene_names, rownames(sorghum_matrix_lists[[i]]))}
all_sorghum_gene_names <- unique(all_sorghum_gene_names[-1])

adding_missing_genes<-function(gene_names,input_matrix){
  missing_genes<-setdiff(gene_names,rownames(input_matrix))
  filler_matrix<-matrix(0,nrow=length(missing_genes),ncol=length(colnames(input_matrix)))
  rownames(filler_matrix)<-missing_genes
  colnames(filler_matrix)<-colnames(input_matrix)
  output<-rbind(input_matrix,filler_matrix)
  output<-output[gene_names,]
  output}

rice_matrix_lists_sized<-list(NA)
for (i in 1:29) {rice_matrix_lists_sized[[i]] <- adding_missing_genes(all_rice_gene_names, rice_matrix_lists[[i]])}

sorghum_matrix_lists_sized<-list(NA)
for (i in 1:31) {sorghum_matrix_lists_sized[[i]] <- adding_missing_genes(all_sorghum_gene_names, sorghum_matrix_lists[[i]])}

# Export for Barcode Report
save(rice_matrix_lists_sized, file = "rice_matrix_lists_sized.RData")
save(sorghum_matrix_lists_sized, file = "sorghum_matrix_lists_sized.RData")

RT_Identifier<-read.table("/gale/ddn/ddn_neomorph/jswift/analysis/R_Light/text_files/Sorghum_Rice_Identifier.txt",header=TRUE, sep="\t") 

Identity<-paste(RT_Identifier$Species, RT_Identifier$Time, RT_Identifier$Replicate, sep='_')
RT_Identifier[,8]<-Identity

get_new_identifiers<-function(input_matrix){

  # get PCR well and barcode name
  PCR_list<-""
  barcode_list<-""
  LIG_list<-""
  RT_list<-""
  fullnames<-colnames(input_matrix)
  
  # Call Barcodes
  for(i in 1:length(fullnames))
    {fullname <- fullnames[i]
     
     # PCR Barcode
     PCRwell <- substring(fullname, str_locate(fullname, "_combo_index")[1]-3 , str_locate(fullname, "_combo_index")[1]-1 )
     PCR_list <- c(PCR_list, PCRwell)
     
     # Barcode
     barcode <- substring(fullname,nchar(fullname)-19,nchar(fullname)) 
     barcode <- str_remove(barcode, "[.]") # some barcodes are 19 and some are 20
     barcode_list <- c(barcode_list, barcode)
     
     # LIG Barcode
     if (nchar(barcode) == 20) Ligation_name<-as.character(LIG_Barcodes[LIG_Barcodes[,5] %in% substr(barcode,1,10),1])
     if (nchar(barcode) == 19) Ligation_name<-as.character(LIG_Barcodes[LIG_Barcodes[,5] %in% substr(barcode,1,9),1])
     LIG_list<- c(LIG_list,Ligation_name)
     
     # RT Barcode
     if (nchar(barcode) == 20) RT_name<-as.character(RT_Barcodes[RT_Barcodes[,5] %in% substr(barcode,11,20),1])
     if (nchar(barcode) == 19) RT_name<-as.character(RT_Barcodes[RT_Barcodes[,5] %in% substr(barcode,10,19),1])
     RT_list<- c(RT_list,RT_name)
     }
  
  barcode_list<-barcode_list[-1]
  PCR_list<-PCR_list[-1]
  LIG_list<-LIG_list[-1]
  RT_list<-RT_list[-1]
  RT_numeric<-as.numeric(str_remove_all(RT_list, "[sc_ligation_RT_]"))
  LIG_numeric<-as.numeric(str_remove_all(LIG_list, "[sc_ligation_]"))

  # Rename RT with sample identity
  RT_name<-NA
  for (i in 1:length(RT_list)) {RT_name[i] <- subset(RT_Identifier, RT_Identifier$Index %in% RT_list[i])[,8] }
  
  # Rename PCR, LIG and RT 
  PCR_name<-paste('PCR',PCR_list,sep='_')
  LIG_name<-paste('LIG',LIG_numeric,sep='_')
  RT_name_2<-paste('RT',formatC(RT_numeric, width=3, flag="0"),sep='_')
  cell_identifier<-paste(RT_name, PCR_name, LIG_name, RT_name_2,sep='_')
  cell_identifier}

for (i in 1:29) {colnames(rice_matrix_lists_sized[[i]]) <- get_new_identifiers(rice_matrix_lists_sized[[i]])}
for (i in 1:31) {colnames(sorghum_matrix_lists_sized[[i]]) <- get_new_identifiers(sorghum_matrix_lists_sized[[i]])}

# Remove Control Wells in Sorghum, and removing RT barcodes designated to Rice plate
sorghum_matrix_lists_sized_clean<-list(NA)
for (i in 1:31) {column_names<-colnames(sorghum_matrix_lists_sized[[i]])
                control_columns<-colnames(sorghum_matrix_lists_sized[[i]][, substring(column_names,6,8) == 'CCC'])
                opposite_species_columns<-as.numeric(substring(column_names, nchar(column_names)-2, nchar(column_names))) %in% c(1:96)
                opposite_species_columns<-colnames(sorghum_matrix_lists_sized[[i]][,opposite_species_columns])
                sorghum_matrix_lists_sized_clean[[i]] <- sorghum_matrix_lists_sized[[i]][ , ! colnames(sorghum_matrix_lists_sized[[i]]) %in% control_columns ] 
                sorghum_matrix_lists_sized_clean[[i]] <- sorghum_matrix_lists_sized_clean[[i]][ , ! colnames(sorghum_matrix_lists_sized_clean[[i]]) %in% opposite_species_columns ] }

# Remove Control Wells in Rice, and removing RT barcodes designated to Sorghum plate
rice_matrix_lists_sized_clean<-list(NA)
for (i in 1:29) {column_names<-colnames(rice_matrix_lists_sized[[i]])
                control_columns<-colnames(rice_matrix_lists_sized[[i]][, substring(column_names,6,8) == 'CCC'])
                opposite_species_columns<-as.numeric(substring(column_names, nchar(column_names)-2, nchar(column_names))) %in% c(97:192)
                opposite_species_columns<-colnames(rice_matrix_lists_sized[[i]][,opposite_species_columns])
                rice_matrix_lists_sized_clean[[i]] <- rice_matrix_lists_sized[[i]][ , ! colnames(rice_matrix_lists_sized[[i]]) %in% control_columns ] 
                rice_matrix_lists_sized_clean[[i]] <- rice_matrix_lists_sized_clean[[i]][ , ! colnames(rice_matrix_lists_sized_clean[[i]]) %in% opposite_species_columns ] }

sorghum_matrix_lists_sized<-sorghum_matrix_lists_sized_clean
rice_matrix_lists_sized<-rice_matrix_lists_sized_clean

# Checking No Duplicate Names Sorghum
print('checking sorghum no duplicates (last number should be zero)')
all_well_names_sorghum <- NA
for (i in 1:31) {all_well_names_sorghum <- c(all_well_names_sorghum, colnames(sorghum_matrix_lists_sized[[i]]))}
length(all_well_names_sorghum)
length(unique(all_well_names_sorghum))
all_well_names_sorghum[duplicated(all_well_names_sorghum)]

# Checking No Duplicate Names Rice
print('checking rice no duplicates (last number should be zero)')
all_well_names_rice <- NA
for (i in 1:29) {all_well_names_rice <- c(all_well_names_rice, colnames(rice_matrix_lists_sized[[i]]))}
length(all_well_names_rice)
length(unique(all_well_names_rice))
all_well_names_rice[duplicated(all_well_names_rice)]

# Resubsetting
full_matrix_rice <- matrix(NA,nrow=length(all_rice_gene_names),ncol=1)
for (i in 1:29) {full_matrix_rice <- cbind(full_matrix_rice, rice_matrix_lists_sized[[i]]) }
full_matrix_rice<-full_matrix_rice[,-1]

full_matrix_sorghum <- matrix(NA,nrow=length(all_sorghum_gene_names),ncol=1)
for (i in 1:31) {full_matrix_sorghum <- cbind(full_matrix_sorghum, sorghum_matrix_lists_sized[[i]]) }
full_matrix_sorghum <-full_matrix_sorghum [,-1]

save(full_matrix_rice, file = "CI_full_matrix_rice.RData")
save(full_matrix_sorghum, file = "CI_full_matrix_sorghum.RData")
