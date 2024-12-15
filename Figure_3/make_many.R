library(argparser)

parser <- arg_parser("parse")
parser <- add_argument(parser, "--seed", help = "seed")

args <- parse_args(parser)

s <- args$seed
set.seed(s)

devtools::load_all("/home/cmonzo/workspace/MOSim")
library(Seurat)
library(tidyverse)
library(pbmcMultiome.SeuratData)

## MAKING BIG DATASET FOR EXAMPLE IN PAPER
omics_types <- c("scRNA-seq", "scATAC-seq")

omicsList <- list()

for (omics in omics_types){
      # Load data from seurat
      if(omics == "scRNA-seq"){
        dat <- pbmcMultiome.SeuratData::pbmc.rna
        dat <- subset(x = dat, subset = seurat_annotations %in% c("CD4 TEM", 
                      "cDC", "Memory B", "Treg", "NK", "CD16 Mono"))
        unique_cell_types <- unique(dat@meta.data$seurat_annotations)
        extracted_cells <- list()
        cellnames <- c()
        for (cell_type in unique_cell_types) {
          type_cells <- subset(dat, subset = seurat_annotations %in% cell_type)
          counts <- as.matrix(type_cells@assays[["RNA"]]@counts)
          extracted_cells[[cell_type]] <- counts
          cellnames <- append(cellnames, replicate(length(colnames(counts)), cell_type))
        }
        dato <- Reduce(cbind, extracted_cells)
        counts <- as.matrix(dato)
        
      } else if (omics == "scATAC-seq"){
        dat <- pbmcMultiome.SeuratData::pbmc.atac
        dat <- subset(x = dat, subset = seurat_annotations %in% c("CD4 TEM", 
                      "cDC", "Memory B", "Treg", "NK", "CD16 Mono"))
        unique_cell_types <- unique(dat@meta.data$seurat_annotations)
        extracted_cells <- list()
        for (cell_type in unique_cell_types) {
          type_cells <- subset(dat, subset = seurat_annotations %in% cell_type)
          counts <- as.matrix(type_cells@assays[["ATAC"]]@counts)
          extracted_cells[[cell_type]] <- counts
        }
        dato <- Reduce(cbind, extracted_cells)
        counts <- as.matrix(dato)
      }
      
      meta <- dat@meta.data$seurat_annotations
      
      metadf <- data.frame("meta" = meta, "cell" = colnames(counts))
      # Sort the metadata according to cell_type
      metadf <- metadf[order(metadf$meta),]
      # Sort by celltype
      counts <- counts[, metadf$cell]
      omicsList[[omics]] <- counts
    
} 

cell_types <- list()
for (char in unique(cellnames)){
  char_positions <- which(cellnames == char)
  cell_types[[char]] <- char_positions
}

sim <- scMOSim(omicsList, cell_types, numberReps = 3, numberGroups = 2,
               diffGenes = list(c(0.3, 0.2)), noiseRep = 0.1,
               noiseGroup = 0.3, clusters = 8, feature_no = 10000,
               regulatorEffect = list(c(0.2, 0.1), c(0.1, 0.2)), TF = TRUE)


saveRDS(object = sim, file = paste0('~/workspace/mosim_paper/sim_6cells8clus10000_scMOSim_2groups_', s, '.rds'))
