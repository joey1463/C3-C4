# This script will compute over-represented cis-regulatory elements within each cell type of the rice multiome data. 
# It will also compute over-represented cis-regulatory elements that are responsive to light in each cell type 

setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_rice_multiome_RNA_combined_with_ATAC.RData")
DefaultAssay(rice_multiome_RNA_combined) <- "ATAC"

# Get a list of motif position frequency matrices from the JASPAR database, and add motif information
pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'plants', all_versions = FALSE))
rice_multiome_RNA_combined <- AddMotifs(object = rice_multiome_RNA_combined, genome = BSgenome.Osativa.MSU.1, pfm = pfm)

# Compute Cis 
register(SerialParam())
rice_multiome_RNA_combined <- RunChromVAR(object = rice_multiome_RNA_combined ,genome = BSgenome.Osativa.MSU.1) # takes a long time to run
rice_multiome_RNA_combined <- RegionStats(rice_multiome_RNA_combined, genome = BSgenome.Osativa.MSU.1) # first compute the GC content for each peak

# Save Object
save(rice_multiome_RNA_combined, file = "L1_rice_multiome_RNA_combined_with_ATAC_with_ChromVAR.RData")
setwd("~/Desktop/Salk/Project - Light/R files")
load("L1_rice_multiome_RNA_combined_with_ATAC_with_ChromVAR.RData")
JASPAR_Plants<-read.table("JASPAR_Plant_TFs.txt",header=TRUE, sep="\t")
rownames(JASPAR_Plants)<-JASPAR_Plants[,1]

# Visualize Activity in Feature Plot
DefaultAssay(rice_multiome_RNA_combined) <- 'chromvar'
FeaturePlot(object = rice_multiome_RNA_combined,features = "MA1183.1",min.cutoff = 'q10',max.cutoff = 'q90',pt.size = 0.1)
rice_multiome_RNA_combined$treatment<-substr(rice_multiome_RNA_combined$dataset,1,4)
FeaturePlot(object = rice_multiome_RNA_combined,features = "MA1190.1",min.cutoff = 'q50',max.cutoff = 'q95',pt.size = 0.5, split.by='treatment')

# call cell-type specific cis elements
DefaultAssay(rice_multiome_RNA_combined) <- 'chromvar'
motifs_bundle_sheath<- FindMarkers(object = rice_multiome_RNA_combined, ident.1 = 'bundle_sheath',only.pos = TRUE, mean.fxn = rowMeans,fc.name = "avg_diff")
motifs_mesophyll <- FindMarkers(object = rice_multiome_RNA_combined, ident.1 = 'mesophyll',only.pos = TRUE, mean.fxn = rowMeans,fc.name = "avg_diff")
motifs_phloem <- FindMarkers(object = rice_multiome_RNA_combined, ident.1 = 'phloem',only.pos = TRUE, mean.fxn = rowMeans,fc.name = "avg_diff")
motifs_epidermis <- FindMarkers(object = rice_multiome_RNA_combined, ident.1 = 'epidermis',only.pos = TRUE, mean.fxn = rowMeans,fc.name = "avg_diff")
motifs_xylem <- FindMarkers(object = rice_multiome_RNA_combined, ident.1 = 'xylem',only.pos = TRUE, mean.fxn = rowMeans,fc.name = "avg_diff")
motifs_guard <- FindMarkers(object = rice_multiome_RNA_combined, ident.1 = 'guard',only.pos = TRUE, mean.fxn = rowMeans,fc.name = "avg_diff")
rice_cis_elements<-list(motifs_guard, motifs_epidermis, motifs_xylem, motifs_phloem, motifs_bundle_sheath, motifs_mesophyll) # order matters for qplot
#save(rice_cis_elements, file="rice_cis_elements.RData")

# call light responsive cis elements
Idents(rice_multiome_RNA_combined)<-paste(Idents(rice_multiome_RNA_combined), str_remove(substr(rice_multiome_RNA_combined$dataset,1,5), "_"), sep="_") #  label identities 
motifs_bundle_sheath_light<- FindMarkers(object = rice_multiome_RNA_combined, ident.1 = 'bundle_sheath_Light', ident.2 = 'bundle_sheath_Dark',only.pos = TRUE, mean.fxn = rowMeans,fc.name = "avg_diff")
motifs_mesophyll_light<- FindMarkers(object = rice_multiome_RNA_combined, ident.1 = 'mesophyll_Light', ident.2 = 'mesophyll_Dark',only.pos = TRUE, mean.fxn = rowMeans,fc.name = "avg_diff")
motifs_phloem_light<- FindMarkers(object = rice_multiome_RNA_combined, ident.1 = 'phloem_Light', ident.2 = 'phloem_Dark',only.pos = TRUE, mean.fxn = rowMeans,fc.name = "avg_diff")
motifs_epidermis_light<- FindMarkers(object = rice_multiome_RNA_combined, ident.1 = 'epidermis_Light', ident.2 = 'epidermis_Dark',only.pos = TRUE, mean.fxn = rowMeans,fc.name = "avg_diff")
motifs_xylem_light<- FindMarkers(object = rice_multiome_RNA_combined, ident.1 = 'xylem_Light', ident.2 = 'xylem_Dark',only.pos = TRUE, mean.fxn = rowMeans,fc.name = "avg_diff")
motifs_guard_light<- FindMarkers(object = rice_multiome_RNA_combined, ident.1 = 'guard_Light', ident.2 = 'guard_Dark',only.pos = TRUE, mean.fxn = rowMeans,fc.name = "avg_diff")
rice_cis_elements_light<-list(motifs_guard_light, motifs_epidermis_light, motifs_xylem_light, motifs_phloem_light, motifs_bundle_sheath_light, motifs_mesophyll_light) # order matters for qplot
#save(rice_cis_elements_light, file="rice_cis_elements_light.RData")

# Supplementary Table Export
load(file="rice_cis_elements.RData")
motifs_bundle_sheath<-merge(JASPAR_Plants,rice_cis_elements[[5]],by='row.names',all.y = TRUE)
motifs_mesophyll <- merge(JASPAR_Plants,rice_cis_elements[[6]],by='row.names',all.y = TRUE)
motifs_phloem <- merge(JASPAR_Plants,rice_cis_elements[[4]],by='row.names',all.y = TRUE)
motifs_epidermis <- merge(JASPAR_Plants,rice_cis_elements[[2]],by='row.names',all.y = TRUE)
motifs_xylem <- merge(JASPAR_Plants,rice_cis_elements[[3]],by='row.names',all.y = TRUE)
motifs_guard <- merge(JASPAR_Plants,rice_cis_elements[[1]],by='row.names',all.y = TRUE)

load(file="sorghum_cis_elements_light.RData")
motifs_bundle_sheath_light<-merge(JASPAR_Plants,sorghum_cis_elements_light[[5]],by='row.names',all.y = TRUE)
motifs_mesophyll_light <- merge(JASPAR_Plants,sorghum_cis_elements_light[[6]],by='row.names',all.y = TRUE)
motifs_phloem_light <- merge(JASPAR_Plants,sorghum_cis_elements_light[[4]],by='row.names',all.y = TRUE)
motifs_epidermis_light <- merge(JASPAR_Plants,sorghum_cis_elements_light[[2]],by='row.names',all.y = TRUE)
motifs_xylem_light <- merge(JASPAR_Plants,sorghum_cis_elements_light[[3]],by='row.names',all.y = TRUE)
motifs_guard_light <- merge(JASPAR_Plants,sorghum_cis_elements_light[[1]],by='row.names',all.y = TRUE)
