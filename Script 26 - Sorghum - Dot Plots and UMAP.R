# This script will take the Sorghum Atlas and label cell types. It will also generate a dot plot using published markers. 

library(Seurat)
library(ggplot2)
load("L1_sorghum_combined.RData") 

# Dot Plot of Selected Markers
Sorghum_Marker_Genes<-c("Sobic.001G168400", # guard marker
                        "Sobic.006G186500", # guard marker
                        "Sobic.010G102800", # guard marker
                        "Sobic.004G342000", # epidermis marker
                        "Sobic.003G443200", # xylem marker
                        "Sobic.001G488700", # phloem marker
                        "Sobic.007G054500", # phloem marker
                        "Sobic.003G002600", # bundle sheath marker
                        "Sobic.008G039900", # bundle sheath marker
                        "Sobic.003G209800", # mesophyll marker
                        "Sobic.010G160700", # mesophyll marker
                        "Sobic.003G234200") # mesophyll marker 

sorghum_combined_labelled <-sorghum_combined
sorghum_combined_labelled <-RenameIdents(sorghum_combined_labelled, `0`   = "mesophyll",
                                                                    `1`   = "mesophyll",
                                                                    `2`   = "unknown",
                                                                    `3`   = "unknown",
                                                                    `4`   = "bundle_sheath",
                                                                    `5`   = "unknown",
                                                                    `6`   = "epidermis",
                                                                    `7`   = "unknown",
                                                                    `8`   = "epidermis",
                                                                    `9`   = "phloem",
                                                                    `10`   = "xylem",
                                                                    `11`   = "unknown",
                                                                    `12`   = "phloem",
                                                                    `13`   = "phloem",
                                                                    `14`   = "guard",
                                                                    `15`   = "bundle_sheath",
                                                                    `16`   = "unknown",
                                                                    `17`   = "epidermis",
                                                                    `18`   = "unknown")
sorghum_combined_labelled$cell_types<-Idents(sorghum_combined_labelled)

sorghum_combined_CI_labelled<-sorghum_combined_labelled[,sorghum_combined_labelled$assay_type %in% c("CI")]
sorghum_combined_CI_labelled<-sorghum_combined_CI_labelled[,substr(colnames(sorghum_combined_CI_labelled),6,7) == '48'] 

# Save It
save(sorghum_combined_labelled, file="L1_sorghum_combined_labelled.RData")
save(sorghum_combined_CI_labelled, file="L1_sorghum_combined_CI_labelled.RData")
#load("L1_sorghum_combined_labelled.RData")

# L1 UMAP 
colors<-as.data.frame(table(Idents(sorghum_combined_labelled)))
colors[,3]<-c("springgreen4","grey90", "steelblue1","navajowhite3","deepskyblue4","turquoise3","tan3")
pdf('sorghum_cluster.pdf')
DimPlot(object = sorghum_combined_labelled, reduction = "umap", label=FALSE, raster = FALSE, cols = colors[,3]) + theme(legend.position="none")
dev.off()

# proportions
proportions<-table(Idents(sorghum_combined_labelled))
proportions<-proportions/sum(proportions)

# Dot Plot
Idents(sorghum_combined_labelled) <- factor(sorghum_combined_labelled@active.ident, rev(c("mesophyll","bundle_sheath","phloem","xylem","epidermis","guard")))
sorghum_combined_dot_plot<-sorghum_combined_labelled[,WhichCells(sorghum_combined_labelled, idents = c("guard","epidermis","phloem","xylem","bundle_sheath","mesophyll"))]
pdf('sorghum_dot.pdf')
DotPlot(sorghum_combined_dot_plot, features = Sorghum_Marker_Genes, col.max = 2, scale.max=10) + RotatedAxis() + scale_colour_gradient2(low = "white", mid = "grey70", high = "hotpink3") +  theme_bw()  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()
