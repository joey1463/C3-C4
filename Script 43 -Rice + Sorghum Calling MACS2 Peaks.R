# This script will call peaks using MACS2 on each of the different 10X-multiome libraries. 

library(hdf5r)
library(Seurat)
library(Signac)

MACS_multiome <- function(counts_path, frag_path, name){
  counts <- Read10X_h5(counts_path)
  frag_path <- frag_path
  
  # create a Seurat object containing the RNA adata
  multiome <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA")
  
  # create ATAC assay and add it to the object
  multiome[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = frag_path)
  
  multiome <- subset(x = multiome, subset = nCount_ATAC < 50000 & nCount_RNA < 8000 & nCount_ATAC > 1000 & nCount_RNA > 400) 
  
  # call peaks using MACS2
  DefaultAssay(multiome) <- "ATAC"
  peaks <- CallPeaks(multiome)
  write.table(peaks,file=paste("L1_",name,".txt",sep=""),quote=F,sep="\t")}

# Rep 1 Rice
MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Rice_Dark_4/Rice_Dark_4/outs/filtered_feature_bc_matrix.h5",
                                      frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Rice_Dark_4/Rice_Dark_4/outs/atac_fragments.tsv.gz",
                                      name = "Rice_Dark_4_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Rice_Dark_8/Rice_Dark_8/outs/filtered_feature_bc_matrix.h5",
                                      frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Rice_Dark_8/Rice_Dark_8/outs/atac_fragments.tsv.gz",
                                      name = "Rice_Dark_8_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Rice_Light_3/Rice_Light_3/outs/filtered_feature_bc_matrix.h5",
                                      frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Rice_Light_3/Rice_Light_3/outs/atac_fragments.tsv.gz",
                                      name = "Rice_Light_3_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Rice_Light_7/Rice_Light_7/outs/filtered_feature_bc_matrix.h5",
                                      frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Rice_Light_7/Rice_Light_7/outs/atac_fragments.tsv.gz",
                                      name = "Rice_Light_7_bed")


# Rep 2 Rice
MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Dark_1/Rice_Dark_1/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Dark_1/Rice_Dark_1/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Rice_Dark_1_Rep2_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Dark_2/Rice_Dark_2/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Dark_2/Rice_Dark_2/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Rice_Dark_2_Rep2_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Dark_3/Rice_Dark_3/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Dark_3/Rice_Dark_3/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Rice_Dark_3_Rep2_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Light_1/Rice_Light_1/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Light_1/Rice_Light_1/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Rice_Light_1_Rep2_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Light_2/Rice_Light_2/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Light_2/Rice_Light_2/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Rice_Light_2_Rep2_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Light_3/Rice_Light_3/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Light_3/Rice_Light_3/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Rice_Light_3_Rep2_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Light_4/Rice_Light_4/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Rice_Light_4/Rice_Light_4/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Rice_Light_4_Rep2_bed")

# Rep 1 Sorghum
Sorghum_Dark_2_Multiome <- MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Sorghum_Dark_2/Sorghum_Dark_2/outs/filtered_feature_bc_matrix.h5",
                                      frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Sorghum_Dark_2/Sorghum_Dark_2/outs/atac_fragments.tsv.gz",
                                      name = "Sorghum_Dark_2_bed")

Sorghum_Dark_6_Multiome <- MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Sorghum_Dark_6/Sorghum_Dark_6/outs/filtered_feature_bc_matrix.h5",
                                      frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Sorghum_Dark_6/Sorghum_Dark_6/outs/atac_fragments.tsv.gz",
                                      name = "Sorghum_Dark_6_bed")
Sorghum_Light_1_Multiome <- MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Sorghum_Light_1/Sorghum_Light_1/outs/filtered_feature_bc_matrix.h5",
                                      frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Sorghum_Light_1/Sorghum_Light_1/outs/atac_fragments.tsv.gz",
                                      name = "Sorghum_Light_1_bed")

Sorghum_Light_5_Multiome <- MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Sorghum_Light_5/Sorghum_Light_5/outs/filtered_feature_bc_matrix.h5",
                                      frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/RNA_ATAC/Sorghum_Light_5/Sorghum_Light_5/outs/atac_fragments.tsv.gz",
                                      name = "Sorghum_Light_5_bed")


# Rep 2 Sorghum
MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Dark_1/Sorghum_Dark_1/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Dark_1/Sorghum_Dark_1/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Sorghum_Dark_1_Rep2_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Dark_2/Sorghum_Dark_2/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Dark_2/Sorghum_Dark_2/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Sorghum_Dark_2_Rep2_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Dark_3/Sorghum_Dark_3/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Dark_3/Sorghum_Dark_3/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Sorghum_Dark_3_Rep2_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Dark_4/Sorghum_Dark_4/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Dark_4/Sorghum_Dark_4/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Sorghum_Dark_4_Rep2_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Light_1/Sorghum_Light_1/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Light_1/Sorghum_Light_1/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Sorghum_Light_1_Rep2_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Light_3/Sorghum_Light_3/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Light_3/Sorghum_Light_3/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Sorghum_Light_3_Rep2_bed")

MACS_multiome(counts_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Light_4/Sorghum_Light_4/outs/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5",
                                          frag_path = "/gale/netapp/home/jswift/analysis/Light_Multiome/Rep_2/Sorghum_Light_4/Sorghum_Light_4/outs/filtered_feature_bc_matrix/atac_fragments.tsv.gz",
                                           name = "Sorghum_Light_4_Rep2_bed")
