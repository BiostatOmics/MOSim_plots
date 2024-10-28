library(ggplot2)
library(tidyr)
library(dplyr)

sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_6cells8clus8000_scMOSim_2groups_022.rds")
associationMatrix <- sim$AssociationMatrices$AssociationMatrix_Group_2

cell_types <- sim$cellTypes

rna <- sim$Group_1$Rep_1$`sim_scRNA-seq`@counts
atac <- sim$Group_1$Rep_1$`sim_scATAC-seq`@counts

rna[is.na(rna)] <- 0
atac[is.na(atac)] <- 0

rna <- as.data.frame(rna)
atac <- as.data.frame(atac)

means_celltype_rna <- calculate_mean_per_list_df(rna, cell_types)
means_celltype_atac <- calculate_mean_per_list_df(atac, cell_types)

# Comparisons to do according to them being regulated
trues <- associationMatrix[associationMatrix$RegulatorEffect %in% "Activator",]
trues <- trues[, c("Gene_ID", "Peak_ID")]

comp_genes <- means_celltype_rna[match(trues$Gene_ID, rownames(means_celltype_rna)),]
comp_peaks <- means_celltype_atac[match(trues$Peak_ID, rownames(means_celltype_atac)),]

# Function to compute Kendall correlation between rows
compute_kendall <- function(gene_row, peak_row) {
  return(cor(gene_row, peak_row, method = "kendall"))
}

correl <- c()

for (i in 1:nrow(trues)) {
  correl[i] <- compute_kendall(as.numeric(comp_genes[i, ]), as.numeric(comp_peaks[i, ]))
}
trues$correlation <- correl

trues <- na.omit(trues)
trues$regulation <- TRUE


## Now get the falses

tmp <- associationMatrix[associationMatrix$Gene_cluster == 0,] %>% drop_na(Gene_ID)
tmp2 <- associationMatrix[associationMatrix$Peak_cluster == 0,] %>% drop_na(Peak_ID)

falses <- as.data.frame(sample(tmp$Gene_ID, 20000))
colnames(falses) <- c("Gene_ID")
falses$Peak_ID <- sample(tmp2$Peak_ID, 20000)

# extract our genes and peaks of interest
comp_genes <- means_celltype_rna[match(falses$Gene_ID, rownames(means_celltype_rna)),]
comp_peaks <- means_celltype_atac[match(falses$Peak_ID, rownames(means_celltype_atac)),]

correl <- c()

for (i in 1:nrow(falses)) {
  correl[i] <- compute_kendall(as.numeric(comp_genes[i, ]), as.numeric(comp_peaks[i, ]))
}
falses$correlation <- correl

falses <- falses %>% drop_na(correlation)
falses$regulation <- "FALSE"

df <- rbind(as.data.frame(trues), as.data.frame(falses))
df$correlation <- abs(df$correlation)

ggplot(df, aes(x = as.factor(regulation), y = correlation)) +
  geom_boxplot(aes(fill = as.factor(regulation)), notch = FALSE) +
  scale_x_discrete(labels = c("FALSE", "TRUE")) +
  labs(x = "Regulation", y = "Absolute correlation") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  scale_fill_manual(values=c("#FDC659", 
                             "#d0b0d4"))

ggsave(width = 3, height = 3.5, filename = "~/workspace/1_conesalab/test_scMOSim/paper_plots/3C_absCorrelation.pdf")

