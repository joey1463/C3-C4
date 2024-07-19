# This script will combine assembled RNA and ATAC multiome data into one object, and identify cell types. 

setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_sorghum_multiome_RNA_combined.RData")
# Pause
load("L1_sorghum_multiome_atac_combined.RData")

summary(colnames(sorghum_multiome_atac_combined) == colnames(sorghum_multiome_RNA_combined))
sorghum_multiome_RNA_combined[['ATAC']] <- sorghum_multiome_atac_combined[['ATAC']] # Transfer Information
sorghum_multiome_RNA_combined[['ATAC_umap']] <- sorghum_multiome_atac_combined[['umap']] # Transfer Information
sorghum_multiome_RNA_combined[['dataset']] <- sorghum_multiome_atac_combined[['dataset']] # Transfer Information
DefaultAssay(sorghum_multiome_RNA_combined) <- "ATAC"
remove(sorghum_multiome_atac_combined)

DimPlot(object = sorghum_multiome_RNA_combined, reduction = "ATAC_umap", label=TRUE, raster = FALSE) + NoLegend()
DimPlot(object = sorghum_multiome_RNA_combined, reduction = "umap", label=TRUE, raster = FALSE) + NoLegend()
DimPlot(object = sorghum_multiome_RNA_combined, reduction = "umap", split.by = 'dataset', label=TRUE, raster = FALSE) + NoLegend()

# Make your GTF file 
gene_coords <- rtracklayer::import('~/Desktop/Salk/Project - Light/R files/Sbicolor_454_v3.1.1.gene_exons_genes_only.gtf')
gene_coords@elementMetadata[,7]<-gene_coords@elementMetadata$type
gene_coords@elementMetadata[,8]<-gene_coords@elementMetadata$gene_id
colnames(gene_coords@elementMetadata)[7] <- "gene_biotype"
colnames(gene_coords@elementMetadata)[8] <- "gene_name"
gene_coords@elementMetadata[,9]<-gene_coords@elementMetadata$gene_id
colnames(gene_coords@elementMetadata)[9] <- "tx_id"
gene_coords@elementMetadata$transcript_id<-substr(gene_coords@elementMetadata$transcript_id,1,16)
gene.ranges <- gene_coords[gene_coords$gene_biotype %in% c('transcript','exon'), ]
tss.ranges <- GRanges(seqnames = seqnames(gene.ranges), ranges = IRanges(start = start(gene.ranges), width = 2), strand = strand(gene.ranges))
gene_coords@elementMetadata[,10]<-tss.ranges
colnames(gene_coords@elementMetadata)[10] <- "tss.ranges"
Annotation(sorghum_multiome_RNA_combined) <- gene_coords

# Call Cell Types
sorghum_multiome_RNA_combined<-RenameIdents(sorghum_multiome_RNA_combined,                     `0`   = "unknown",
                                                                                               `1`   = "mesophyll",
                                                                                               `2`   = "mesophyll",
                                                                                               `3`   = "mesophyll",
                                                                                               `4`   = "guard", 
                                                                                               `5`   = "epidermis",
                                                                                               `6`   = "bundle_sheath",
                                                                                               `7`   = "phloem",
                                                                                               `8`   = "unknown", 
                                                                                               `9`   = "unknown",   
                                                                                              `10`   = "unknown",
                                                                                              `11`   = "unknown",
                                                                                              `12`   = "xylem",
                                                                                              `13`   = "epidermis", 
                                                                                              `14`   = "epidermis", 
                                                                                              `15`   = "unknown", 
                                                                                              `16`   = "unknown",
                                                                                              `17`   = "bundle_sheath",
                                                                                              `18`   = "unknown", 
                                                                                              `19`   = "unknown") 
# L1 UMAP 
colors<-as.data.frame(table(Idents(sorghum_multiome_RNA_combined)))
colors[,3]<-c("grey90", "springgreen4","tan3","navajowhite3","steelblue1","deepskyblue4","turquoise3")
DimPlot(object = sorghum_multiome_RNA_combined, reduction = "umap", label=FALSE, raster = FALSE, cols = colors[,3]) 

# Call Average Peaks and Average Expression
DefaultAssay(sorghum_multiome_RNA_combined) <- "RNA"
Average_Expression<-as.data.frame(x = AverageExpression(object = sorghum_multiome_RNA_combined, verbose = FALSE)$RNA)

DefaultAssay(sorghum_multiome_RNA_combined) <- "ATAC"
Average_Peaks<-as.data.frame(AverageExpression(object = sorghum_multiome_RNA_combined)[[3]])

setwd("~/Desktop/Salk/Project - Light/R files")
save(sorghum_multiome_RNA_combined, file = "L1_sorghum_multiome_RNA_combined_with_ATAC.RData")
save(Average_Peaks, file = "L1_sorghum_multiome_Average_Peaks.RData")
save(Average_Expression, file = "L1_sorghum_multiome_Average_Expression.RData")
