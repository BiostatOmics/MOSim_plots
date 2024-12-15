# MOSim: bulk and single-cell multi-layer regulatory network simulator

Carolina Monzó, Carlos Martínez-Mira, Ángeles Arzalluz-Luque, Maider Aguerralde-Martin, Ana Conesa, Sonia Tarazona

The MOSim R package is available at  
<https://www.bioconductor.org/packages/release/bioc/html/MOSim.html>.

This repository contains the code used to generate the analyses and figures in the paper.  
The repository is organized in 3 directories, one for each of the figures depicting results of the manuscript.  
Below, we describe the files included in this repository.  

#### Figure_2

[2BC_plots](Figure_2/2BC_plots.R) performs Principal Component Analysis and creates a boxplot of absolute Pearson's correlation values from interactions of RNA-seq genes in Group 1 with ATAC-seq regulators. Regulation is TRUE when the regulator has been simulated to activate/repress gene expression, and FALSE if there is no simulated relationship between the gene and regulator.  
  
[2D_plotProfile](Figure_2/2D_plotProfile.R) uses MOSim's plotProfile function to depict the relationship between a gene and it's activator, and a gene and it's repressor. It does so for both Group 1 and 2, showing the relationship between RNA-seq and ATAC-seq simulated data.  

#### Figure_3

[3B_RNA_PCA_figure](Figure_3/3B_RNA_PCA_figure.R) performs Principal Component analysis and visualizations for single-cell RNA-seq data, showing clustering by cell type, experimental group and replicate.  

[3B_PCA_ATAC_figure](Figure_3/3B_PCA_ATAC_figure.R) performs Principal Component analysis and visualizations for single-cell ATAC-seq data, showing clustering by cell type, experimental group and replicate.  

[3C_SCcorrelations_figure](Figure_3/3C_SCcorrelations_figure.R) creates boxplots of absolute Kendall correlation values from interactions of scRNA-seq genes in Group 1 with scATAC-seq regulators. Regulation is TRUE when the regulator has been simulated to activate/repress gene expression or FALSE otherwise. 

[3D_ActRep](Figure_3/3D_ActRep.R) plots two examples of gene-regulator single-cell simulated profiles in each group. Left Y-axis shows gene expression values, while right Y-axis shows counts for scATAC-seq regions. Vertical bars at each time point show the standard error of the mean of the cells for the 3 simulated replicates. 

[make_many](Figure_3/make_many.R) simulates scRNA-seq and scATAC-seq data for six cell types in 2 groups and 3 replicates.  

[run](Figure_3/run.sh) calls script make_many.R to simulate the single-cell GRN multi-omics dataset in a high performance computing cluster.

#### Figure_4  

[4ABC_multilayer_regulatory_networks](Figure_4/4ABC_multilayer_regulatory_networks.R) creates a representation of multi-layer regulatory networks with bulk data simulated by MOSim. The file first simulates bulk data for RNA-seq, miRNA-seq and TF. Takes de first 100 DEGs and creates the Gene Regulatory Networks (GRN) for the two simulated groups. Finally it creates heatmaps for the expression profiles of the genes, miRNAs and TFs in the GRN of Groups 1 and 2.
