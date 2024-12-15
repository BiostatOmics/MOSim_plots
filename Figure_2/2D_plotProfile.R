# Load necessary libraries
library(MOSim)
library(igraph)
library(RCy3)
library(ComplexHeatmap)
library(ggplot2)

# Read simulated data

omics_list <- c("RNA-seq", "miRNA-seq", "DNase-seq")

rnaseq_simulation <- readRDS("~/workspace/mosim_paper/bulk_rnaseq_simulation.rds")

# Get the count tables
rnaseq_simulated <- omicResults(rnaseq_simulation, omics_list)
# get the settings used to generate each count table
all_settings <- omicSettings(rnaseq_simulation)
design_matrix <- experimentalDesign(rnaseq_simulation)

# Plot activator
profile_plot <- plotProfile(multi_simulation,
            omics = c("RNA-seq", "DNase-seq"),
            featureIDS = list(
                "RNA-seq" = "ENSMUSG00000052726",
                "DNase-seq" = "1_140257767_140257897"
            ))

# print
profile_plot +
    theme_bw() +
    theme(legend.position="top")

# Plot repressor

profile_plot <- plotProfile(multi_simulation,
            omics = c("RNA-seq", "DNase-seq"),
            featureIDS = list(
                "RNA-seq" = "ENSMUSG00000020434",
                "DNase-seq" = "11_3879784_3880304"
            ))

# print
profile_plot +
    theme_bw() +
    theme(legend.position="top")
