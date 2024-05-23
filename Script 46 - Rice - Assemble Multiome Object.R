# This script will combine assembled RNA and ATAC multiome data into one object, and identify cell types. 

setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_rice_multiome_RNA_combined.RData")
# Pause
load("L1_rice_multiome_atac_combined.RData")

summary(colnames(rice_multiome_atac_combined) == colnames(rice_multiome_RNA_combined))
rice_multiome_RNA_combined[['ATAC']] <- rice_multiome_atac_combined[['ATAC']] # Transfer Information
rice_multiome_RNA_combined[['ATAC_umap']] <- rice_multiome_atac_combined[['umap']] # Transfer Information
rice_multiome_RNA_combined[['dataset']] <- rice_multiome_atac_combined[['dataset']] # Transfer Information
DefaultAssay(rice_multiome_RNA_combined) <- "ATAC"
remove(rice_multiome_atac_combined)

DimPlot(object = rice_multiome_RNA_combined, reduction = "ATAC_umap", label=TRUE, raster = FALSE) + NoLegend()
DimPlot(object = rice_multiome_RNA_combined, reduction = "umap", label=TRUE, raster = FALSE) + NoLegend()
DimPlot(object = rice_multiome_RNA_combined, reduction = "umap", split.by = 'dataset', label=TRUE, raster = FALSE) + NoLegend()

# Make your GTF file 
gene_coords <- rtracklayer::import('~/Desktop/Salk/Project - Light/R files/MSU_rice_onlygene_edit.gtf')
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
Annotation(rice_multiome_RNA_combined) <- gene_coords

# Designate Cell Types
rice_multiome_RNA_combined <-RenameIdents(rice_multiome_RNA_combined,                   `0`   = "mesophyll",
                                                                                        `1`   = "unknown",
                                                                                        `2`   = "unknown",
                                                                                        `3`   = "fiber",
                                                                                        `4`   = "unknown",
                                                                                        `5`   = "unknown",
                                                                                        `6`   = "mesophyll",
                                                                                        `7`   = "unknown",
                                                                                        `8`   = "unknown",
                                                                                        `9`   = "guard", 
                                                                                        `10`   = "unknown",
                                                                                        `11`   = "epidermis",
                                                                                        `12`   = "xylem",
                                                                                        `13`   = "bundle_sheath",
                                                                                        `14`   = "epidermis",
                                                                                        `15`   = "unknown",
                                                                                        `16`   = "phloem",
                                                                                        `17`   = "unknown")
# L1 UMAP 
colors<-as.data.frame(table(Idents(rice_multiome_RNA_combined)))
colors[,3]<-c("springgreen4","grey90", "grey90","tan3","navajowhite3","turquoise3","steelblue1","deepskyblue4")
DimPlot(object = rice_multiome_RNA_combined, reduction = "umap", label=FALSE, cols = colors[,3]) #+ theme(legend.position="none")

# Call Average Peaks and Average Expression
DefaultAssay(rice_multiome_RNA_combined) <- "RNA"
Average_Expression<-as.data.frame(x = AverageExpression(object = rice_multiome_RNA_combined, verbose = FALSE)$RNA)

DefaultAssay(rice_multiome_RNA_combined) <- "ATAC"
Average_Peaks<-as.data.frame(AverageExpression(object = rice_multiome_RNA_combined)[[3]])

setwd("~/Desktop/Salk/Project - Light/R files/")
save(rice_multiome_RNA_combined, file = "L1_rice_multiome_RNA_combined_with_ATAC.RData")
save(Average_Peaks, file = "L1_rice_multiome_Average_Peaks.RData")
save(Average_Expression, file = "L1_rice_multiome_Average_Expression.RData")
