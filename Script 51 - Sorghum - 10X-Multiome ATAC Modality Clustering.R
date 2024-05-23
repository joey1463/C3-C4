# This script will read in and cluster the Sorghum 10X-Multiome ATAC data. 

# MACS2 called peaks
Sorghum_Light_1_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Sorghum_Light_1_bed.txt",header = TRUE)
Sorghum_Light_5_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Sorghum_Light_5_bed.txt",header = TRUE)
Sorghum_Dark_2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Sorghum_Dark_2_bed.txt",header = TRUE)
Sorghum_Dark_6_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Sorghum_Dark_6_bed.txt",header = TRUE)

Sorghum_Dark_1_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Sorghum_Dark_1_Rep2_bed.txt", header = TRUE)
Sorghum_Dark_2_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Sorghum_Dark_2_Rep2_bed.txt", header = TRUE)
Sorghum_Dark_3_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Sorghum_Dark_3_Rep2_bed.txt", header = TRUE)
Sorghum_Dark_4_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Sorghum_Dark_4_Rep2_bed.txt", header = TRUE)

Sorghum_Light_1_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Sorghum_Light_1_Rep2_bed.txt", header = TRUE)
Sorghum_Light_3_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Sorghum_Light_3_Rep2_bed.txt", header = TRUE)
Sorghum_Light_4_Rep2_bed <- read.table(file = "~/Desktop/Salk/Project - Light/R Files/MACS peaks/L1_Sorghum_Light_4_Rep2_bed.txt", header = TRUE)

Sorghum_Light_1_ranges <- makeGRangesFromDataFrame(Sorghum_Light_1_bed)
Sorghum_Light_5_ranges <- makeGRangesFromDataFrame(Sorghum_Light_5_bed)
Sorghum_Dark_2_ranges <- makeGRangesFromDataFrame(Sorghum_Dark_2_bed)
Sorghum_Dark_6_ranges <- makeGRangesFromDataFrame(Sorghum_Dark_6_bed)

Sorghum_Dark_1_Rep2_ranges <- makeGRangesFromDataFrame(Sorghum_Dark_1_Rep2_bed)
Sorghum_Dark_2_Rep2_ranges <- makeGRangesFromDataFrame(Sorghum_Dark_2_Rep2_bed)
Sorghum_Dark_3_Rep2_ranges <- makeGRangesFromDataFrame(Sorghum_Dark_3_Rep2_bed)
Sorghum_Dark_4_Rep2_ranges <- makeGRangesFromDataFrame(Sorghum_Dark_4_Rep2_bed)

Sorghum_Light_1_Rep2_ranges <- makeGRangesFromDataFrame(Sorghum_Light_1_Rep2_bed)
Sorghum_Light_3_Rep2_ranges <- makeGRangesFromDataFrame(Sorghum_Light_3_Rep2_bed)
Sorghum_Light_4_Rep2_ranges <- makeGRangesFromDataFrame(Sorghum_Light_4_Rep2_bed)

Sorghum_peaks_combined <- reduce(x = c(Sorghum_Light_1_ranges, Sorghum_Light_5_ranges, Sorghum_Dark_2_ranges, Sorghum_Dark_6_ranges,
                                       Sorghum_Dark_1_Rep2_ranges, Sorghum_Dark_2_Rep2_ranges, Sorghum_Dark_3_Rep2_ranges, Sorghum_Dark_4_Rep2_ranges,
                                       Sorghum_Light_1_Rep2_ranges, Sorghum_Light_3_Rep2_ranges, Sorghum_Light_4_Rep2_ranges))

peakwidths <- width(Sorghum_peaks_combined)
Sorghum_peaks_combined <- Sorghum_peaks_combined[peakwidths  < 10000 & peakwidths > 20]

# Get Selected Nuclei
setwd("~/Desktop/Salk/Project - Light/R files")
load("Sorghum_Light_1_RNA.RData")
load("Sorghum_Light_5_RNA.RData")
load("Sorghum_Dark_2_RNA.RData")
load("Sorghum_Dark_6_RNA.RData")

load("Sorghum_Dark_1_Rep2_RNA.RData")
load("Sorghum_Dark_2_Rep2_RNA.RData")
load("Sorghum_Dark_3_Rep2_RNA.RData")
load("Sorghum_Dark_4_Rep2_RNA.RData")
load("Sorghum_Light_1_Rep2_RNA.RData")
load("Sorghum_Light_3_Rep2_RNA.RData")
load("Sorghum_Light_4_Rep2_RNA.RData")

# Get Fragments
Sorghum_Light_1_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/1_Sorghum_Light/atac_fragments.tsv.gz", cells = colnames(Sorghum_Light_1_RNA)) 
Sorghum_Light_5_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/5_Sorghum_Light/atac_fragments.tsv.gz", cells = colnames(Sorghum_Light_5_RNA)) 
Sorghum_Dark_2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/2_Sorghum_Dark/atac_fragments.tsv.gz", cells = colnames(Sorghum_Dark_2_RNA)) 
Sorghum_Dark_6_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/RNA_ATAC Alignment/6_Sorghum_Dark/atac_fragments.tsv.gz", cells = colnames(Sorghum_Dark_6_RNA)) 

Sorghum_Dark_1_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Dark_1/atac_fragments.tsv.gz", cells = colnames(Sorghum_Dark_1_Rep2_RNA))   
Sorghum_Dark_2_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Dark_2/atac_fragments.tsv.gz", cells = colnames(Sorghum_Dark_2_Rep2_RNA))   
Sorghum_Dark_3_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Dark_3/atac_fragments.tsv.gz", cells = colnames(Sorghum_Dark_3_Rep2_RNA))   
Sorghum_Dark_4_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Dark_4/atac_fragments.tsv.gz", cells = colnames(Sorghum_Dark_4_Rep2_RNA))   

Sorghum_Light_1_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Light_1/atac_fragments.tsv.gz", cells = colnames(Sorghum_Light_1_Rep2_RNA))   
Sorghum_Light_3_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Light_3/atac_fragments.tsv.gz", cells = colnames(Sorghum_Light_3_Rep2_RNA))   
Sorghum_Light_4_Rep2_fragments <- CreateFragmentObject(path = "~/Desktop/Salk/Project - Light/Multiome/Rep 2/Sorghum_Light_4/atac_fragments.tsv.gz", cells = colnames(Sorghum_Light_4_Rep2_RNA))   


# Create Counts
Sorghum_Light_1_counts <- FeatureMatrix(fragments = Sorghum_Light_1_fragments, features = Sorghum_peaks_combined, cells = colnames(Sorghum_Light_1_RNA)) 
Sorghum_Light_5_counts <- FeatureMatrix(fragments = Sorghum_Light_5_fragments, features = Sorghum_peaks_combined, cells = colnames(Sorghum_Light_5_RNA)) 
Sorghum_Dark_2_counts <- FeatureMatrix(fragments = Sorghum_Dark_2_fragments, features = Sorghum_peaks_combined, cells = colnames(Sorghum_Dark_2_RNA)) 
Sorghum_Dark_6_counts <- FeatureMatrix(fragments = Sorghum_Dark_6_fragments, features = Sorghum_peaks_combined, cells = colnames(Sorghum_Dark_6_RNA)) 

Sorghum_Dark_1_Rep2_counts <- FeatureMatrix(fragments = Sorghum_Dark_1_Rep2_fragments, features = Sorghum_peaks_combined, cells = colnames(Sorghum_Dark_1_Rep2_RNA)) 
Sorghum_Dark_2_Rep2_counts <- FeatureMatrix(fragments = Sorghum_Dark_2_Rep2_fragments, features = Sorghum_peaks_combined, cells = colnames(Sorghum_Dark_2_Rep2_RNA)) 
Sorghum_Dark_3_Rep2_counts <- FeatureMatrix(fragments = Sorghum_Dark_3_Rep2_fragments, features = Sorghum_peaks_combined, cells = colnames(Sorghum_Dark_3_Rep2_RNA)) 
Sorghum_Dark_4_Rep2_counts <- FeatureMatrix(fragments = Sorghum_Dark_4_Rep2_fragments, features = Sorghum_peaks_combined, cells = colnames(Sorghum_Dark_4_Rep2_RNA)) 

Sorghum_Light_1_Rep2_counts <- FeatureMatrix(fragments = Sorghum_Light_1_Rep2_fragments, features = Sorghum_peaks_combined, cells = colnames(Sorghum_Light_1_Rep2_RNA)) 
Sorghum_Light_3_Rep2_counts <- FeatureMatrix(fragments = Sorghum_Light_3_Rep2_fragments, features = Sorghum_peaks_combined, cells = colnames(Sorghum_Light_3_Rep2_RNA)) 
Sorghum_Light_4_Rep2_counts <- FeatureMatrix(fragments = Sorghum_Light_4_Rep2_fragments, features = Sorghum_peaks_combined, cells = colnames(Sorghum_Light_4_Rep2_RNA)) 

# Create Objects
Sorghum_Light_1_assay <- CreateChromatinAssay(Sorghum_Light_1_counts, fragments = Sorghum_Light_1_fragments)
Sorghum_Light_1 <- CreateSeuratObject(Sorghum_Light_1_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata

Sorghum_Light_5_assay <- CreateChromatinAssay(Sorghum_Light_5_counts, fragments = Sorghum_Light_5_fragments)
Sorghum_Light_5 <- CreateSeuratObject(Sorghum_Light_5_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata

Sorghum_Dark_2_assay <- CreateChromatinAssay(Sorghum_Dark_2_counts, fragments = Sorghum_Dark_2_fragments)
Sorghum_Dark_2 <- CreateSeuratObject(Sorghum_Dark_2_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata

Sorghum_Dark_6_assay <- CreateChromatinAssay(Sorghum_Dark_6_counts, fragments = Sorghum_Dark_6_fragments)
Sorghum_Dark_6 <- CreateSeuratObject(Sorghum_Dark_6_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata


Sorghum_Dark_1_Rep2_assay <- CreateChromatinAssay(Sorghum_Dark_1_Rep2_counts, fragments = Sorghum_Dark_1_Rep2_fragments)
Sorghum_Dark_1_Rep2 <- CreateSeuratObject(Sorghum_Dark_1_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata

Sorghum_Dark_2_Rep2_assay <- CreateChromatinAssay(Sorghum_Dark_2_Rep2_counts, fragments = Sorghum_Dark_2_Rep2_fragments)
Sorghum_Dark_2_Rep2 <- CreateSeuratObject(Sorghum_Dark_2_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata

Sorghum_Dark_3_Rep2_assay <- CreateChromatinAssay(Sorghum_Dark_3_Rep2_counts, fragments = Sorghum_Dark_3_Rep2_fragments)
Sorghum_Dark_3_Rep2 <- CreateSeuratObject(Sorghum_Dark_3_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata

Sorghum_Dark_4_Rep2_assay <- CreateChromatinAssay(Sorghum_Dark_4_Rep2_counts, fragments = Sorghum_Dark_4_Rep2_fragments)
Sorghum_Dark_4_Rep2 <- CreateSeuratObject(Sorghum_Dark_4_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Dark_metadata


Sorghum_Light_1_Rep2_assay <- CreateChromatinAssay(Sorghum_Light_1_Rep2_counts, fragments = Sorghum_Light_1_Rep2_fragments)
Sorghum_Light_1_Rep2 <- CreateSeuratObject(Sorghum_Light_1_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Light_metadata

Sorghum_Light_3_Rep2_assay <- CreateChromatinAssay(Sorghum_Light_3_Rep2_counts, fragments = Sorghum_Light_3_Rep2_fragments)
Sorghum_Light_3_Rep2 <- CreateSeuratObject(Sorghum_Light_3_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Light_metadata

Sorghum_Light_4_Rep2_assay <- CreateChromatinAssay(Sorghum_Light_4_Rep2_counts, fragments = Sorghum_Light_4_Rep2_fragments)
Sorghum_Light_4_Rep2 <- CreateSeuratObject(Sorghum_Light_4_Rep2_assay, assay = "ATAC") # meta.data = Sorghum_Light_metadata



# Merging
Sorghum_Light_1$dataset <- 'Light_1'
Sorghum_Light_5$dataset <- 'Light_5'
Sorghum_Dark_2$dataset <- 'Dark_2'
Sorghum_Dark_6$dataset <- 'Dark_6'

Sorghum_Dark_1_Rep2$dataset <- 'Dark_1_rep2'
Sorghum_Dark_2_Rep2$dataset <- 'Dark_2_rep2'
Sorghum_Dark_3_Rep2$dataset <- 'Dark_3_rep2'
Sorghum_Dark_4_Rep2$dataset <- 'Dark_4_rep2'

Sorghum_Light_1_Rep2$dataset <- 'Light_1_rep2'
Sorghum_Light_3_Rep2$dataset <- 'Light_3_rep2'
Sorghum_Light_4_Rep2$dataset <- 'Light_4_rep2'


# merge all datasets, adding a cell ID to make sure cell names are unique
sorghum_multiome_atac_combined <- merge(x = Sorghum_Light_1, y = list(Sorghum_Light_5, Sorghum_Dark_2, Sorghum_Dark_6,
                                                                      Sorghum_Dark_1_Rep2, Sorghum_Dark_2_Rep2, Sorghum_Dark_3_Rep2,
                                                                      Sorghum_Dark_4_Rep2,
                                                                      Sorghum_Light_1_Rep2, Sorghum_Light_3_Rep2, Sorghum_Light_4_Rep2)) 

sorghum_multiome_atac_combined <- RunTFIDF(sorghum_multiome_atac_combined, method = 3)
sorghum_multiome_atac_combined <- FindTopFeatures(sorghum_multiome_atac_combined, min.cutoff = 20)
sorghum_multiome_atac_combined <- RunSVD(sorghum_multiome_atac_combined) #DepthCor(sorghum_multiome_atac_combined)
sorghum_multiome_atac_combined <- RunUMAP(sorghum_multiome_atac_combined, dims = 2:30, reduction = 'lsi')
sorghum_multiome_atac_combined <- FindNeighbors(object = sorghum_multiome_atac_combined, reduction = 'lsi', dims = 2:30)
sorghum_multiome_atac_combined <- FindClusters(object = sorghum_multiome_atac_combined, verbose = FALSE, algorithm = 3, resolution = 1.2)
#DimPlot(sorghum_multiome_atac_combined, pt.size = 1, raster=TRUE)

setwd("~/Desktop/Salk/Project - Light/R files")
save(sorghum_multiome_atac_combined,file="L1_sorghum_multiome_atac_combined.RData")
