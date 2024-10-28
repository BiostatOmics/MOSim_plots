devtools::load_all("~/workspace/1_conesalab/MOSim")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(dplyr)
})

set.seed(000)

sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_6cells8clus8000_scMOSim_2groups_022.rds")
sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_6cells8clus1600_scMOSim_2groups_222.rds")

cell_types <- sim$cellTypes
associationMatrix <- sim$AssociationMatrices$AssociationMatrix_Group_2

res1 <- match_gene_regulator_cluster(sim$Group_2$Rep_1$`sim_scRNA-seq`$counts, 
                                     sim$Group_2$Rep_1$`sim_scATAC-seq`$counts, 
                                     cell_types, associationMatrix)

res2 <- match_gene_regulator_cluster(sim$Group_2$Rep_2$`sim_scRNA-seq`$counts, 
                                     sim$Group_2$Rep_2$`sim_scATAC-seq`$counts, 
                                     cell_types, associationMatrix)

res3 <- match_gene_regulator_cluster(sim$Group_2$Rep_3$`sim_scRNA-seq`$counts, 
                                     sim$Group_2$Rep_3$`sim_scATAC-seq`$counts, 
                                     cell_types, associationMatrix)




# Extract expression matrices from sim
mat_orange1 <- res1$rna
mat_orange2 <- res2$rna
mat_orange3 <- res3$rna

mat_green1 <- res1$atac
mat_green2 <- res2$atac
mat_green3 <- res3$atac

mat_orange1 <- sim$Group_2$Rep_1$`sim_scRNA-seq`$counts
mat_orange2 <- sim$Group_2$Rep_2$`sim_scRNA-seq`$counts
mat_orange3 <- sim$Group_2$Rep_3$`sim_scRNA-seq`$counts

mat_green1 <- sim$Group_2$Rep_1$`sim_scATAC-seq`$counts
mat_green2 <- sim$Group_2$Rep_2$`sim_scATAC-seq`$counts
mat_green3 <- sim$Group_2$Rep_3$`sim_scATAC-seq`$counts

# Vector indicating cell types
cell_type_vector <- rep(names(sim$cellTypes), times = lengths(sim$cellTypes))

# Function to aggregate data by cell type within each replicate
aggregate_by_celltype <- function(mat, cell_types) {
  sapply(unique(cell_types), function(ct) rowMeans(mat[, cell_types == ct], na.rm = TRUE))
}

# Aggregating data for each replicate
agg_orange1 <- aggregate_by_celltype(mat_orange1, cell_type_vector)
agg_orange2 <- aggregate_by_celltype(mat_orange2, cell_type_vector)
agg_orange3 <- aggregate_by_celltype(mat_orange3, cell_type_vector)

agg_green1 <- aggregate_by_celltype(mat_green1, cell_type_vector)
agg_green2 <- aggregate_by_celltype(mat_green2, cell_type_vector)
agg_green3 <- aggregate_by_celltype(mat_green3, cell_type_vector)

# Combining the aggregated data across replicates and calculating SEM
combine_replicates <- function(agg1, agg2, agg3) {
  unique_ct <- colnames(agg1)
  combined_mean <- sapply(unique_ct, function(ct) rowMeans(cbind(agg1[, ct], agg2[, ct], agg3[, ct]), na.rm = TRUE))
  combined_sem <- sapply(unique_ct, function(ct) apply(cbind(agg1[, ct], agg2[, ct], agg3[, ct]), 1, function(x) sd(x, na.rm = TRUE) / sqrt(3)))
  list(mean = combined_mean, sem = combined_sem)
}

agg_orange_combined <- combine_replicates(agg_orange1, agg_orange2, agg_orange3)
agg_green_combined <- combine_replicates(agg_green1, agg_green2, agg_green3)

# Select specific gene and peak
selected_gene <- "CYB5RL"
selected_peak <- "chr9-107963188-107963693"
agg_orange_selected <- data.frame(
  CellType = colnames(agg_orange_combined$mean), 
  Mean = agg_orange_combined$mean[selected_gene, ], 
  SEM = agg_orange_combined$sem[selected_gene, ],
  Omic = "Orange"
)
agg_green_selected <- data.frame(
  CellType = colnames(agg_green_combined$mean), 
  Mean = agg_green_combined$mean[selected_peak, ], 
  SEM = agg_green_combined$sem[selected_peak, ],
  Omic = "Green"
)

# Plotting
p1 <- ggplot(agg_orange_selected, aes(x = CellType, y = Mean, group = 1)) +
  ylim(0, 0.18) +
  geom_line(color = "#15918A") +
  geom_point(color = "#15918A") +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, color = "#15918A") +
  theme_classic() +
  labs(x = "Cell Type", y = "RBP7") +
  theme(axis.title.y.left = element_text(color = "black"), axis.text.y.left = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(agg_green_selected, aes(x = CellType, y = Mean, group = 1)) +
  geom_line(color = "#F58A53") +
  geom_point(color = "#F58A53") +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, color = "#F58A53") +
  theme_classic() +
  labs(y = "chr3-101753518-101753798") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "chr3-101753518-101753798")) +
  ylim(0, 0.08) +
  theme(axis.title.y.right = element_text(color = "black"), axis.text.y.right = element_text(color = "black"), 
        axis.ticks.y.right = element_line(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1))

# Combine the plots with patchwork
p1 + p2 + plot_layout(ncol = 1) + plot_annotation(title = "Gene and Peak Expression Across Cell Types")

ggsave(filename = "./Repressor_G2.pdf", height = 4, width = 3)



## Plot also an activator
# Select specific gene and peak
selected_gene <- "PTPN22"
selected_peak <- "chr12-31742761-31743451"
agg_orange_selected <- data.frame(
  CellType = colnames(agg_orange_combined$mean), 
  Mean = agg_orange_combined$mean[selected_gene, ], 
  SEM = agg_orange_combined$sem[selected_gene, ],
  Omic = "Orange"
)
agg_green_selected <- data.frame(
  CellType = colnames(agg_green_combined$mean), 
  Mean = agg_green_combined$mean[selected_peak, ], 
  SEM = agg_green_combined$sem[selected_peak, ],
  Omic = "Green"
)

# Plotting
p1 <- ggplot(agg_orange_selected, aes(x = CellType, y = Mean, group = 1)) +
  ylim(0, 1.5) +
  geom_line(color = "#15918A") +
  geom_point(color = "#15918A") +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, color = "#15918A") +
  theme_classic() +
  labs(x = "Cell Type", y = "PTPN22") +
  theme(axis.title.y.left = element_text(color = "black"), axis.text.y.left = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(agg_green_selected, aes(x = CellType, y = Mean, group = 1)) +
  geom_line(color = "#F58A53") +
  geom_point(color = "#F58A53") +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2, color = "#F58A53") +
  theme_classic() +
  labs(y = "chr12-31742761-31743451") +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "chr12-31742761-31743451")) +
  ylim(0, 1.5) +
  theme(axis.title.y.right = element_text(color = "black"), axis.text.y.right = element_text(color = "black"), 
        axis.ticks.y.right = element_line(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1))

library(patchwork)

# Combine the plots with patchwork
p1 + p2 + plot_layout(ncol = 1) + plot_annotation(title = "Gene and Peak Expression Across Cell Types")

ggsave(filename = "./Activator_G2.pdf", height = 4, width = 3)
