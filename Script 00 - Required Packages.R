# Packages
library(Seurat)
library(Signac)
library(DESeq2)
require(ggplot2)
library(hdf5r)
library(EDASeq)
library(dplyr)
library(data.table)
library(plyranges)
library(MotifDb)
library(cowplot)
library(viridis)
library(VennDiagram)
library(TFBSTools)
library(memes)
library(JASPAR2020)
library(motifmatchr)
library(ggseqlogo)
library(rtracklayer)
library(plyr)
library(ggalluvial)
require(reshape2)
require(gridExtra)
library(gplots)
library(plotly)
library(splines)
library(igraph)
library(qlcMatrix)
library(ggforce)
library(Matrix)
library(devtools)
library(stringr)
library(DoubletFinder)
library(SingleCellExperiment)
library(clusterExperiment)

# BS Genomes
library(BSgenome.Osativa.MSU.1) 
library(BSgenome.SBicolor.454.2) 
library(Claxum.version.1.0)
library(Hvulgare.version.1.0)
library(BDistachyon.version.1.0)

# Read in Barcodes
setwd("~/Desktop/Salk/Protocols/Combo Indexing Oligos")

# RT Barcodes - all barcodes are of length 10
RT_Barcodes<-read.table("sc_RT_plate.txt",header=TRUE, sep="\t")
RT_Barcode<-NA
for (i in 1:384) {RT_Barcode[i]<-substring(as.character(RT_Barcodes[i,4]),22,31)}
RT_Barcodes[,5]<-RT_Barcode 

# LIG Barcodes - some LIG barcodes are 9 bp, others are 10bp
LIG_Barcodes<-read.table("sc_ligation_plate.txt",header=TRUE, sep="\t")
LIG_Barcode<-NA
for (i in 1:384) {if (nchar(as.character(LIG_Barcodes[i,4])) == 53) {LIG_Barcode[i]<-substring(as.character(LIG_Barcodes[i,4]),44,53)} 
                  if (nchar(as.character(LIG_Barcodes[i,4])) == 51) {LIG_Barcode[i]<-substring(as.character(LIG_Barcodes[i,4]),43,51)}}
LIG_Barcodes[,5]<-LIG_Barcode
