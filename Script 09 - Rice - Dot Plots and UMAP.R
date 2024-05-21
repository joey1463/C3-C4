# This script will take the Rice Atlas and label cell types. It will also generate a dot plot using published markers. 

library(Seurat)
library(ggplot2)
load("L1_rice_combined.RData") 

Rice_Marker_Genes<-c(
  "LOC-Os06g14030", # Guard
  "LOC-Os02g56920", # Epidermis
  "LOC-Os11g32650", # Epidermis
  "LOC-Os01g73980", # Xylem
  "LOC-Os06g41090", # Phloem
  "LOC-Os01g41710", # Mesophyll
  "LOC-Os12g19470", # Mesophyll
  "LOC-Os11g26160", # Mesophyll
  "LOC-Os08g06630") # Mesophyll

rice_combined_labelled <-rice_combined
rice_combined_labelled <-RenameIdents(rice_combined_labelled, `0`   = "mesophyll",
                                                              `1`   = "unknown",
                                                              `2`   = "unknown",
                                                              `3`   = "epidermis",
                                                              `4`   = "fiber",
                                                              `5`   = "mesophyll",
                                                              `6`   = "phloem",
                                                              `7`   = "epidermis",
                                                              `8`   = "unknown",
                                                              `9`   = "unknown",
                                                              `10`   = "bundle_sheath",
                                                              `11`   = "unknown",
                                                              `12`   = "xylem",
                                                              `13`   = "fiber",
                                                              `14`   = "epidermis",
                                                              `15`   = "guard",
                                                              `16`   = "mestome sheath",
                                                              `17`   = "unknown",
                                                              `18`   = "unknown")
rice_combined_labelled$cell_types<-Idents(rice_combined_labelled) 

# Save It
save(rice_combined_labelled, file="L1_rice_combined_labelled.RData")
#load("L1_rice_combined_labelled.RData")

# L1 UMAP Coloring 
colors<-as.data.frame(table(Idents(rice_combined_labelled)))
colors[,3]<-c("springgreen4","grey90", "navajowhite3","peachpuff4","deepskyblue4","steelblue1","turquoise3", "tan3","royalblue3")
pdf('rice_cluster.pdf')
DimPlot(object = rice_combined_labelled, reduction = "umap", label=FALSE, raster = FALSE, cols = colors[,3]) + theme(legend.position="none")
dev.off()

# Proportions
proportions<-table(Idents(rice_combined_labelled))
proportions<-proportions/sum(proportions)

# Dot Plot
Idents(rice_combined_labelled) <- factor(rice_combined_labelled@active.ident, rev(c("mesophyll","bundle_sheath","phloem","xylem","epidermis","guard")))
rice_combined_dot_plot<-rice_combined_labelled[,WhichCells(rice_combined_labelled, idents = c("guard","epidermis","phloem","xylem","bundle_sheath","mesophyll"))]
pdf('rice_dot.pdf')
DotPlot(rice_combined_dot_plot, features = Rice_Marker_Genes, col.max = 2, scale.max=10) + RotatedAxis() + scale_colour_gradient2(low = "white", mid = "grey70", high = "hotpink3") +  theme_bw()  + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
