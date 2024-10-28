setwd("~/workspace/1_conesalab/MOSim")
devtools::load_all()

suppressPackageStartupMessages({
  library(tidyverse)
  library(scater)
  library(WGCNA)
  library(Seurat)
  library(RColorConesa)
  library(acorde)
})

set.seed(123)
# ## Plots for paper
# # Create a list of omics data types (e.g., scRNA-seq and scATAC-seq)
omicsList <- sc_omicData(list("scRNA-seq", "scATAC-seq"))

cell_types <- list('Treg' = c(1:10),'cDC' = c(11:20),'CD4_TEM' = c(21:30),
                   'Memory_B' = c(31:40))

sim <- scMOSim(omicsList, cell_types, numberReps = 3, numberGroups = 2,
               diffGenes = list(c(0.3, 0.2)), noiseRep = 0.1,
               noiseGroup = 0.3, clusters = 8, feature_no = 8000,
               regulatorEffect = list(c(0.2, 0.1), c(0.1, 0.2)), TF = TRUE)


saveRDS(object = sim, file = '~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_scMOSim_few_2602.rds')