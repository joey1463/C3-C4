# Exaptation of ancestral cell identity networks enables C4 photosynthesis

This repository contains the code for the paper "Exaptation of ancestral cell identity networks enables C4 photosynthesis" published in XXXX.

## Environment Setup

Ensure that R version 4.0.5 is installed on your system. Download it from [CRAN](https://cran.r-project.org/mirrors.html).

### Installing Required Libraries
```R
# Read the requirements file
requirements <- read.table("R_requirements.txt", header = FALSE, skip = 1, sep = "\t", stringsAsFactors = FALSE, col.names = c("Package", "Version"))
# Check for missing packages and install them
missing_packages <- setdiff(requirements$Package, installed.packages()[, "Package"])
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
# Verify installation
# Re-check installed packages
installed <- installed.packages()[, "Package"]
# Identify any packages that failed to install
failed_to_install <- setdiff(missing_packages, installed)
if (length(failed_to_install) > 0) {
  warning("The following packages failed to install:", paste(failed_to_install, collapse = ", "))
} else {
  message("All packages installed successfully.")
}
```
### Clone the Repository

Clone this repository to access the scripts on your local machine:
 ```bash
git clone https://github.com/joey1463/C3-C4.git
cd C3-C4
```
# R Scripts for single nuclei expression and chromatin analyses

This repository contains a series of R scripts designed to handle various aspects of RNA-seq data analysis, particularly focusing on sci-RNA-seq3 demultiplexing, clustering, and visualization for Rice and Sorghum experiments. Each script is tailored to specific stages of the data analysis pipeline.

## Overview of Scripts

### Script 00 - Required Packages

- **Purpose**: Initializes the environment by loading necessary R packages.

### Script 01 - sci-RNA-seq3 Demultiplexing 0 - 12hr Time Points
- **Purpose**: Processes read counts from sci-RNA-seq3 experiments for 0-12 hour time points, separating data by species and compiling into counts matrices suitable for Seurat analysis.
- **Operations**:
  - **Barcode Reading**: Differentiates between RT and LIG barcodes, handling variations in barcode length (9 or 10 bp).
  - **Matrix Generation**: Constructs a matrix ensuring all rows are aligned and consistent.
  - **Data Cleaning**: Removes control wells and ensures no duplicate names in the datasets for Sorghum and Rice.
  - **Exporting Data**: Outputs a report detailing the barcode processing.

### Script 02 - sci-RNA-seq3 Demultiplexing 48hr Time Points
- **Purpose**: Similar to Script 01, but processes data from 48-hour time points.

### Script 03 - Rice - Clustering 10X and sci-RNA-seq3
- **Purpose**: Integrates sci-RNA-seq3 and 10X-RNA experiments for rice into a single Seurat object.

### Script 04 - Sorghum - Clustering 10X and sci-RNA-seq3
- **Purpose**: Similar to Script 03 but for sorghum.

### Script 05 - Rice + Sorghum Cross-Species Clustering
- **Purpose**: Performs cross-species clustering for rice and sorghum using ortholog matches.

### Script 06 - Rice - 10X & sci-RNA-seq3 Atlas Visualization
- **Purpose**: Visualizes the clustered rice atlas from RNA data, supporting subsampling for performance.
- **Data**: Loads `L1_rice_combined.RData` (full dataset).

### Script 07 - Rice - Barplot Composition
- **Purpose**: Generates barplots to compare cluster representation across assay types and time points.

### Script 08 - Rice - Find Cluster Markers
- **Purpose**: Identifies the top 200 specific markers for each cluster in the rice data.

### Script 09 - Rice - Dot Plots and UMAP
- **Purpose**: Labels cell types in the Rice Atlas and generates dot plots using published markers.
- **Data**: Loads `L1_rice_combined_labelled.RData`.

### Script 10 - Rice - mTurquoise Clustering
- **Purpose:** Processes and clusters 10X-RNA nuclei data from the rice mTurquoise line.
- **Features:**
  - Combines data from two replicates.
  - Identifies markers.
  - Generates plots.

### Script 11 - Rice -Identifying Bundle Sheath with mTurquoise Markers
- **Purpose:** Uses mTurquoise line specific markers to identify bundle sheath cells, comparing with the Rice atlas.
- **Features:**
  - Analysis of marker overlap.
  - Module score plot generation.

### Script 12 - Rice - Cell Class Clustering
- **Purpose:** Re-cluster various cell types including mesophyll, epidermal, and vasculature clusters.
- **Cell Classes:** Mesophyll, Vasculature, Epidermis

### Script 13 - Rice - Cell Class Investigating
- **Purpose:** Visualize clustering results for different rice cell classes.
- **Visualized Cell Classes:** Vasculature, Mesophyll, Epidermis

### Script 14 - Rice - Compute Differential Expression
- **Purpose:** Calculate differentially expressed genes across identified cell types.

### Script 15 - Rice - Pseudo-Bulk Expression Profiles and T0-T12 Comparison
- **Purpose:** Analyze transcriptional profiles and differential expression under etiolated conditions.
- **Features:**
  - Generate and analyze pseudo-bulked expression profiles.
  - Compare T0 and T12 responses.

### Script 16 - Rice - Recluster Mesophyll using DE Genes
- **Purpose:** Subcluster mesophyll cells using differentially expressed genes.
- **Features:**
  - Recluster using variable features.

### Script 17 - Rice - Load Data Sets
- **Purpose:** Import pseudo-bulked expression profiles and other relevant data sets.
- **Features:**
  - Handle time-series data.
  - Export gene lists for further analysis.

### Script 18 - Rice - Cell Type Specific Expression of Photosynthesis Genes
- **Purpose:** Analyze light-responsive genes involved in photosynthesis within specific cell types.
- **Features:**
  - Subset and analyze differentially expressed genes.
  - Generate heatmaps.

### Script 19 - Rice - Volcano Plot
- **Purpose:** Construct volcano plots to visualize differentially expressed genes between cell types over time.
- **Features:**
  - Generate volcano plots for 0h and 12h time points.
  - Additional scatterplot visualization.

### Script 20: Rice - Plotting Individual Genes
- **Description**: Generates cell-type specific gene expression profiles for selected candidate genes.
- **Key Operations**: Gene expression plotting.
- **Figures**: HY5s and PIFs, supplementary figures.

### Script 21: Rice - Z-score Plots
- **Description**: Clusters and visualizes expression patterns of genes differentially expressed in response to light.
- **Key Operations**: Differential expression clustering followed by Z-score plotting.

### Script 22: Rice - Pairwise Differential Expression of Cell Type Pairs
- **Description**: Performs ANCOVA analysis on pseduobulked transcriptional profiles to analyze differential expression between cell type pairs.

### Script 23: Sorghum - 10X & sci-RNA-seq3 Atlas Visualization
- **Description**: Visualizes clustered Sorghum Atlas using RNA data.
- **Key Operations**: Data subsampling for local processing, feature plotting.

### Script 24: Sorghum - Barplot Composition
- **Description**: Creates barplots to compare cluster representation across different assay types and time points.

### Script 25: Sorghum - Find Cluster Markers
- **Description**: Identifies top cluster-specific markers for each cluster.
- **Key Operations**: Marker identification.

### Script 26: Sorghum - Dot Plots and UMAP
- **Description**: Labels cell types in Sorghum Atlas and generates dot plots using published markers.
- **Key Operations**: UMAP visualization, dot plot creation.

### Script 27: Sorghum - Cell Class Clustering
- **Description**: Clusters specific cell classes including mesophyll, epidermal, and vasculature clusters.

### Script 28: Sorghum - Cell Class Investigating
- **Description**: Visualizes clustering of cell classes like vasculature, mesophyll, and epidermis.

### Script 29: Sorghum - Compute Differential Expression
- **Description**: Computes differentially expressed genes for each identified cell type.

### Script 30: Sorghum - Pseudo-Bulk Expression Profiles and T0-T12 Comparison
- **Purpose**: Computes pseudo-bulk transcriptional profiles for each cell type and identifies differentially expressed genes between mesophyll and bundle sheath cell types under etiolated conditions.
- **Key Functions**: Average and aggregate expression computation, T0 and T12 response analysis.

### Script 31: Sorghum - Recluster Bundle Sheath using DE Genes
- **Purpose**: Subclusters bundle sheath 10X nuclei from sorghum using differentially expressed genes as variable features.
- **Data**: Loads reclustered data from `L3_sorghum_bundle_sheath_reclustered.RData`.

### Script 32: Sorghum - Load Data Sets
- **Purpose**: Reads in pseudo-bulked expression profiles, names of photosynthesis-related genes, and transcription factors.
- **Operations**: Time averaging, gene list export for supplementary materials.

### Script 33: Sorghum - Cell Type Specific Expression of Photosynthesis Genes
- **Purpose**: Analyzes light-responsive differentially expressed genes in specific cell types, focusing on photosynthesis genes.
- **Output**: Generates heatmaps for visual analysis and saves them to `heatmap_sup_sorghum.txt`.

### Script 34: Sorghum - Volcano Plots
- **Purpose**: Constructs volcano plots to visualize differentially expressed genes between mesophyll and bundle sheath cell types at T0 and T12.
- **Features**: Includes extended data scatterplot.

### Script 35: Sorghum - Plotting Individual Genes
- **Purpose**: Creates cell-type specific expression profiles for candidate genes.

### Script 36: Sorghum - Z-score Plots
- **Purpose**: Clusters differentially expressed genes by light response to identify dominant expression trends, visualized through Z-score plots.
- **Figures**: Includes figures for both upregulated and downregulated genes.

### Script 37: Sorghum - Pairwise Differential Expression of Cell Type Pairs
- **Purpose**: Performs ANCOVA analysis on unnormalized pseudo-bulked transcriptional profiles to evaluate differential expression due to cell type and light response.
- **Details**: Outputs tissue or time factor coefficients and p-values for significant genes.

### Script 38: RNA Orthology - Calculate Cell Type Markers For Both Species
- **Purpose**: Identifies cell type-specific marker genes across species.

### Script 39: RNA Orthology - Overlapping Cell Type Specific Genes Across Species
- **Purpose**: Analyzes overlap of cell type-specific genes across species using orthology datasets.
- **Outputs**: Heatmaps of significance and gene overlap counts.

### Script 40 - RNA Orthology - Sankey Plots
- **Purpose**: Generates Sankey plots to display comparisons of cell-type marker genes across species.
- **Features**:
  - Plotting Sankey diagrams.
  - Exporting specific lists related to the plots.

### Script 41 - RNA Orthology - Cross Species UMAP
- **Purpose**: Processes cross-species clustered data to label nuclei based on their annotation from the rice or sorghum RNA atlas.
- **Features**:
  - Visualizing unlabelled and labelled UMAPs.
  - Conditional labeling based on cell type source.

### Script 42 - RNA Orthology - C4 Bundle Sheath's Gain and Loss of C3 Genes
- **Purpose**: Identifies differential expression of cell type markers in the bundle sheath cells of rice and sorghum.
- **Features**:
  - Calculation and visualization of gained and lost genes.
  - Saving and exporting data subsets for further analysis.

### Script 43 - RNA Orthology - Consistent and Differential Partitioning
- **Purpose**: Compares differentially partitioned genes across species.
- **Features**:
  - Handling data partitions and merging by orthogroups.
  - Visualizing results through heatmaps and bar plots.

### Script 44 - RNA Orthology - Heatmaps of M and BS Partition Patterns
- **Purpose**: Examines partitioning patterns specifically in mesophyll vs. bundle sheath across species.
- **Features**:
  - Generating heatmaps for differentially and consistently partitioned gene pairs.

### Script 45 - Rice + Sorghum Calling MACS2 Peaks
- **Purpose**: Calls peaks using MACS2 on different 10X-multiome libraries for rice and sorghum.
- **Features**:
  - Peak calling with replication management.

### Script 46 - Rice - 10X-Multiome RNA Modality Clustering
- **Purpose**: Clusters Rice 10X-Multiome RNA data.
- **Features**:
  - Analysis of multiple replicates.
  - Assessment of clustering metrics.

### Script 47 - Rice - 10X-Multiome ATAC Modality Clustering
- **Purpose**: Clusters Rice 10X-Multiome ATAC data.
- **Features**:
  - Data merging and object creation.
  - Visualization of clustering results.

### Script 48 - Rice - Assemble Multiome Object
- **Purpose**: Combines assembled RNA and ATAC multiome data into a single object to identify cell types.
- **Features**:
  - Comprehensive data assembly and cell type identification.

### Script 49 - Rice - Overlap Multiome Atlas with RNA Atlas
- **Purpose**: Compares cell type-specific gene expression markers from multiome data with RNA atlas data.
- **Features**:
  - Identification and comparison of markers.
  - Visualization via heatmaps.

### Script 50 - Rice - Compute M & BS Specific Genes, and DOF Expression Patterns
- **Features**:
  - Computes mesophyll-specific and bundle sheath-specific genes.
  - Visualizes DOF family gene expression patterns in these cell types.

### Script 51 - Rice - Compute Cis-Element Enrichment
- **Features**:
  - Identifies over-represented cis-regulatory elements responsive to light in each cell type.
  - Uses JASPAR database for motif information.


### Script 52 - Sorghum - 10X-Multiome RNA Modality Clustering
- **Features**:
  - Clusters Sorghum 10X-Multiome RNA data.

### Script 53 - Sorghum - 10X-Multiome ATAC Modality Clustering
- **Features**:
  - Clusters Sorghum 10X-Multiome ATAC data, including peak calling and data merging.

### Script 54 - Sorghum - Assemble Multiome Object
- **Features**:
  - Combines RNA and ATAC data into one object and identifies cell types.

### Script 55 - Sorghum - Overlap Multiome Atlas with RNA Atlas
- **Features**:
  - Compares cell type-specific gene expression markers with RNA Atlas data.

### Script 56 - Sorghum - Compute M & BS Specific Genes, and DOF Expression Patterns
- **Features**:
  - Similar to Script 50, adapted for Sorghum.

### Script 57 - Sorghum - Compute Cis-Element Enrichment
- **Features**:
  - Similar to Script 51, adapted for Sorghum.


### Script 58 - Chromatin Orthology - Partitioning RNA for Multiome
- **Features**:
  - Analyzes gene expression partitioning across species and identifies differentially and consistently partitioned genes.

### Script 59 - Chromatin Orthology - Cis Element Overlap
- **Features**:
  - Overlaps cis-regulatory elements across species, assessing statistical significance.
  
### Script 60 - Chromatin Orthology - Comparing Number of DOF Sites
- **Purpose**: Compares DOF site counts associated with differentially partitioned genes across species using a binomial test.
- **Visualizations**: Includes scatterplots for differential and consistent comparisons.
- **Output**: Exports analysis results.

### Script 61 - Rice - Chromatin Heatmaps and Plotting Tracks
- **Purpose**: Generates heatmaps of cell type-specific accessible chromatin from a multiome dataset.
- **Key Features**: Visualizes patterns of accessibility and outputs gene names for GO terms analysis.

### Script 62 - Rice - Light Responsive Changes in Chromatin Accessibility
- **Purpose**: Analyzes light-responsive changes in chromatin accessibility for selected genes.
- **Output**: Plots changes and reports statistical significance.

### Script 63 - Rice - Enriched Motifs in Differentially Partitioned Genes
- **Purpose**: Identifies and ranks enriched motifs within accessible chromatin of differentially partitioned genes.
- **Data Management**: Includes steps to save and load results for further iteration and analysis.

### Script 64 - Rice - Counting DOFs in Differentially Partitioned Genes
- **Purpose**: Counts Dof.2 motifs in differentially partitioned genes using motif scanning.
- **Verification**: Compares motif count to manual counts as a sanity check.

### Script 65 - Rice - Plotting Chromatin Accessibility Tracks of Individual Genes
- **Purpose**: Plots accessible chromatin tracks within specific cell types for candidate genes.
- **Key Feature**: Focuses on genes like GAPDH and NADP-ME.

### Script 66 - Sorghum - Chromatin Heatmaps and Plotting Tracks
- **Purpose**: Similar to the Rice script 61, adapted for Sorghum.
- **Data Handling**: Includes commands to set working directories and export results.

### Script 67 - Sorghum - Light Responsive Changes in Chromatin Accessibility
- **Purpose**: Assesses changes in chromatin accessibility in response to light conditions.
- **Statistical Analysis**: Reports p-values and conducts pairwise t-tests.

### Script 68 - Sorghum - Enriched Motifs in Differentially Partitioned Genes
- **Purpose**: Analyzes enriched motifs in Sorghum, similar to the corresponding Rice script.

### Script 69 - Sorghum - Counting DOFs in Differentially Partitioned Genes
- **Purpose**: Counts Dof.2 motifs in Sorghum, ensuring consistency with similar analyses in Rice.

### Script 70 - Sorghum - Plotting Chromatin Accessibility Tracks of Individual Genes

- **Purpose**: This script plots accessible chromatin tracks within mesophyll and bundle sheath cell types for candidate genes in Sorghum. It is designed to help visualize differences in chromatin accessibility between these cell types.
- **Main Functions**:
  - Read in data
  - Call peaks
  - Plot for genes like GAPDH and NADPME

### Script 71: C3 and C4 Grasses Dof2 Motif Enrichment

- **Purpose**: Computes and analyzes motif enrichment for genes that are consistently partitioned into the bundle sheath in both rice and sorghum. This script extends the analysis to include orthologs from other C3 grasses.
- **Main Functions**:
  - Read in required data for multiple species (Chasmanthium, Rice, Sorghum, Barley, Brachypodium)
  - Use AME from the MEME suite to assess cis-regulatory enrichment within their promoters
  - Plot and export outcomes
