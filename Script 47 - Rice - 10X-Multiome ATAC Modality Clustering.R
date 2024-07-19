# This script will read in and cluster the Rice 10X-Multiome ATAC data. 

# MACS2 called peaks
Rice_Light_3_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Rice_Light_3_bed.txt", header = TRUE)
Rice_Light_7_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Rice_Light_7_bed.txt", header = TRUE)
Rice_Dark_4_bed <-  read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Rice_Dark_4_bed.txt", header = TRUE)
Rice_Dark_8_bed <-  read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Rice_Dark_8_bed.txt", header = TRUE)

Rice_Dark_1_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Rice_Dark_1_Rep2_bed.txt", header = TRUE)
Rice_Dark_2_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Rice_Dark_2_Rep2_bed.txt", header = TRUE)
Rice_Dark_3_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Rice_Dark_3_Rep2_bed.txt", header = TRUE)

Rice_Light_1_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Rice_Light_1_Rep2_bed.txt", header = TRUE)
Rice_Light_2_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Rice_Light_2_Rep2_bed.txt", header = TRUE)
Rice_Light_3_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Rice_Light_3_Rep2_bed.txt", header = TRUE)
Rice_Light_4_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Rice_Light_4_Rep2_bed.txt", header = TRUE)

Rice_Light_3_ranges <- makeGRangesFromDataFrame(Rice_Light_3_bed)
Rice_Light_7_ranges <- makeGRangesFromDataFrame(Rice_Light_7_bed)
Rice_Dark_4_ranges <- makeGRangesFromDataFrame(Rice_Dark_4_bed)
Rice_Dark_8_ranges <- makeGRangesFromDataFrame(Rice_Dark_8_bed)

Rice_Dark_1_Rep2_ranges <- makeGRangesFromDataFrame(Rice_Dark_1_Rep2_bed)
Rice_Dark_2_Rep2_ranges <- makeGRangesFromDataFrame(Rice_Dark_2_Rep2_bed)
Rice_Dark_3_Rep2_ranges <- makeGRangesFromDataFrame(Rice_Dark_3_Rep2_bed)

Rice_Light_1_Rep2_ranges <- makeGRangesFromDataFrame(Rice_Light_1_Rep2_bed)
Rice_Light_2_Rep2_ranges <- makeGRangesFromDataFrame(Rice_Light_2_Rep2_bed)
Rice_Light_3_Rep2_ranges <- makeGRangesFromDataFrame(Rice_Light_3_Rep2_bed)
Rice_Light_4_Rep2_ranges <- makeGRangesFromDataFrame(Rice_Light_4_Rep2_bed)

Rice_peaks_combined <- reduce(x = c(Rice_Light_3_ranges, Rice_Light_7_ranges, Rice_Dark_4_ranges, Rice_Dark_8_ranges, 
                                    Rice_Dark_1_Rep2_ranges, Rice_Dark_2_Rep2_ranges, Rice_Dark_3_Rep2_ranges,
                                    Rice_Light_1_Rep2_ranges, Rice_Light_2_Rep2_ranges, Rice_Light_3_Rep2_ranges, Rice_Light_4_Rep2_ranges))
peakwidths <- width(Rice_peaks_combined)
Rice_peaks_combined <- Rice_peaks_combined[peakwidths  < 10000 & peakwidths > 20]

# Get Selected Nuclei
setwd("~/Desktop/Salk/Project - Light/R files")
load("Rice_Light_3_RNA.RData")
load("Rice_Light_7_RNA.RData")
load("Rice_Dark_4_RNA.RData")
load("Rice_Dark_8_RNA.RData")

load("Rice_Dark_1_Rep2_RNA.RData")
load("Rice_Dark_2_Rep2_RNA.RData")
load("Rice_Dark_3_Rep2_RNA.RData")
load("Rice_Light_1_Rep2_RNA.RData")
load("Rice_Light_2_Rep2_RNA.RData")
load("Rice_Light_3_Rep2_RNA.RData")
load("Rice_Light_4_Rep2_RNA.RData")

# Get Fragments
Rice_Light_3_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/3_Rice_Light/atac_fragments.tsv.gz", cells = colnames(Rice_Light_3_RNA))   
Rice_Light_7_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/7_Rice_Light/atac_fragments.tsv.gz", cells = colnames(Rice_Light_7_RNA)) 
Rice_Dark_4_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/4_Rice_Dark/atac_fragments.tsv.gz", cells = colnames(Rice_Dark_4_RNA)) 
Rice_Dark_8_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/8_Rice_Dark/atac_fragments.tsv.gz", cells = colnames(Rice_Dark_8_RNA)) 

Rice_Dark_1_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Dark_1/atac_fragments.tsv.gz", cells = colnames(Rice_Dark_1_Rep2_RNA))   
Rice_Dark_2_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Dark_2/atac_fragments.tsv.gz", cells = colnames(Rice_Dark_2_Rep2_RNA))   
Rice_Dark_3_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Dark_3/atac_fragments.tsv.gz", cells = colnames(Rice_Dark_3_Rep2_RNA))   

Rice_Light_1_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Light_1/atac_fragments.tsv.gz", cells = colnames(Rice_Light_1_Rep2_RNA))   
Rice_Light_2_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Light_2/atac_fragments.tsv.gz", cells = colnames(Rice_Light_2_Rep2_RNA))   
Rice_Light_3_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Light_3/atac_fragments.tsv.gz", cells = colnames(Rice_Light_3_Rep2_RNA))   
Rice_Light_4_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Rice_Light_4/atac_fragments.tsv.gz", cells = colnames(Rice_Light_4_Rep2_RNA))   

# Create Counts
Rice_Light_3_counts <- FeatureMatrix(fragments = Rice_Light_3_fragments, features = Rice_peaks_combined, cells = colnames(Rice_Light_3_RNA)) 
Rice_Light_7_counts <- FeatureMatrix(fragments = Rice_Light_7_fragments, features = Rice_peaks_combined, cells = colnames(Rice_Light_7_RNA)) 
Rice_Dark_4_counts <- FeatureMatrix(fragments = Rice_Dark_4_fragments, features = Rice_peaks_combined, cells = colnames(Rice_Dark_4_RNA)) 
Rice_Dark_8_counts <- FeatureMatrix(fragments = Rice_Dark_8_fragments, features = Rice_peaks_combined, cells = colnames(Rice_Dark_8_RNA)) 

Rice_Dark_1_Rep2_counts <- FeatureMatrix(fragments = Rice_Dark_1_Rep2_fragments, features = Rice_peaks_combined, cells = colnames(Rice_Dark_1_Rep2_RNA)) 
Rice_Dark_2_Rep2_counts <- FeatureMatrix(fragments = Rice_Dark_2_Rep2_fragments, features = Rice_peaks_combined, cells = colnames(Rice_Dark_2_Rep2_RNA)) 
Rice_Dark_3_Rep2_counts <- FeatureMatrix(fragments = Rice_Dark_3_Rep2_fragments, features = Rice_peaks_combined, cells = colnames(Rice_Dark_3_Rep2_RNA)) 

Rice_Light_1_Rep2_counts <- FeatureMatrix(fragments = Rice_Light_1_Rep2_fragments, features = Rice_peaks_combined, cells = colnames(Rice_Light_1_Rep2_RNA)) 
Rice_Light_2_Rep2_counts <- FeatureMatrix(fragments = Rice_Light_2_Rep2_fragments, features = Rice_peaks_combined, cells = colnames(Rice_Light_2_Rep2_RNA)) 
Rice_Light_3_Rep2_counts <- FeatureMatrix(fragments = Rice_Light_3_Rep2_fragments, features = Rice_peaks_combined, cells = colnames(Rice_Light_3_Rep2_RNA)) 
Rice_Light_4_Rep2_counts <- FeatureMatrix(fragments = Rice_Light_4_Rep2_fragments, features = Rice_peaks_combined, cells = colnames(Rice_Light_4_Rep2_RNA)) 

# Create Objects
Rice_Light_3_assay <- CreateChromatinAssay(Rice_Light_3_counts, fragments = Rice_Light_3_fragments)
Rice_Light_3 <- CreateSeuratObject(Rice_Light_3_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata

Rice_Light_7_assay <- CreateChromatinAssay(Rice_Light_7_counts, fragments = Rice_Light_7_fragments)
Rice_Light_7 <- CreateSeuratObject(Rice_Light_7_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata

Rice_Dark_4_assay <- CreateChromatinAssay(Rice_Dark_4_counts, fragments = Rice_Dark_4_fragments)
Rice_Dark_4 <- CreateSeuratObject(Rice_Dark_4_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata

Rice_Dark_8_assay <- CreateChromatinAssay(Rice_Dark_8_counts, fragments = Rice_Dark_8_fragments)
Rice_Dark_8 <- CreateSeuratObject(Rice_Dark_8_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata


Rice_Dark_1_Rep2_assay <- CreateChromatinAssay(Rice_Dark_1_Rep2_counts, fragments = Rice_Dark_1_Rep2_fragments)
Rice_Dark_1_Rep2 <- CreateSeuratObject(Rice_Dark_1_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata

Rice_Dark_2_Rep2_assay <- CreateChromatinAssay(Rice_Dark_2_Rep2_counts, fragments = Rice_Dark_2_Rep2_fragments)
Rice_Dark_2_Rep2 <- CreateSeuratObject(Rice_Dark_2_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata

Rice_Dark_3_Rep2_assay <- CreateChromatinAssay(Rice_Dark_3_Rep2_counts, fragments = Rice_Dark_3_Rep2_fragments)
Rice_Dark_3_Rep2 <- CreateSeuratObject(Rice_Dark_3_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata

Rice_Light_1_Rep2_assay <- CreateChromatinAssay(Rice_Light_1_Rep2_counts, fragments = Rice_Light_1_Rep2_fragments)
Rice_Light_1_Rep2 <- CreateSeuratObject(Rice_Light_1_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Light_metadata

Rice_Light_2_Rep2_assay <- CreateChromatinAssay(Rice_Light_2_Rep2_counts, fragments = Rice_Light_2_Rep2_fragments)
Rice_Light_2_Rep2 <- CreateSeuratObject(Rice_Light_2_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Light_metadata

Rice_Light_3_Rep2_assay <- CreateChromatinAssay(Rice_Light_3_Rep2_counts, fragments = Rice_Light_3_Rep2_fragments)
Rice_Light_3_Rep2 <- CreateSeuratObject(Rice_Light_3_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Light_metadata

Rice_Light_4_Rep2_assay <- CreateChromatinAssay(Rice_Light_4_Rep2_counts, fragments = Rice_Light_4_Rep2_fragments)
Rice_Light_4_Rep2 <- CreateSeuratObject(Rice_Light_4_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Light_metadata


# Merging
Rice_Light_3$dataset <- 'Light_3'
Rice_Light_7$dataset <- 'Light_7'
Rice_Dark_4$dataset <- 'Dark_4'
Rice_Dark_8$dataset <- 'Dark_8'

Rice_Dark_1_Rep2$dataset <- 'Dark_1_rep2'
Rice_Dark_2_Rep2$dataset <- 'Dark_2_rep2'
Rice_Dark_3_Rep2$dataset <- 'Dark_3_rep2'

Rice_Light_1_Rep2$dataset <- 'Light_1_rep2'
Rice_Light_2_Rep2$dataset <- 'Light_2_rep2'
Rice_Light_3_Rep2$dataset <- 'Light_3_rep2'
Rice_Light_4_Rep2$dataset <- 'Light_4_rep2'

# merge all datasets, adding a cell ID to make sure cell names are unique
rice_multiome_atac_combined <- merge(x = Rice_Light_3, y = list(Rice_Light_7, Rice_Dark_4, Rice_Dark_8,
                                                                Rice_Dark_1_Rep2, Rice_Dark_2_Rep2, Rice_Dark_3_Rep2,
                                                                Rice_Light_1_Rep2, Rice_Light_2_Rep2, Rice_Light_3_Rep2, Rice_Light_4_Rep2)) 
rice_multiome_atac_combined <- RunTFIDF(rice_multiome_atac_combined, method = 3)
rice_multiome_atac_combined <- FindTopFeatures(rice_multiome_atac_combined, min.cutoff = 20)
rice_multiome_atac_combined <- RunSVD(rice_multiome_atac_combined) 
rice_multiome_atac_combined <- RunUMAP(rice_multiome_atac_combined, dims = 2:30, reduction = 'lsi')
rice_multiome_atac_combined <- FindNeighbors(object = rice_multiome_atac_combined, reduction = 'lsi', dims = 2:30)
rice_multiome_atac_combined <- FindClusters(object = rice_multiome_atac_combined, verbose = FALSE, algorithm = 3, resolution = 1.2)
#DimPlot(rice_multiome_atac_combined, pt.size = 1, raster=TRUE)

setwd("~/Desktop/Salk/Project - Light/R files/")
save(rice_multiome_atac_combined,file="L1_rice_multiome_atac_combined.RData")
