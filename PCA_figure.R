# Concatting samples for PCA
devtools::load_all("/home/cmonzo/workspace/MOSim")
suppressPackageStartupMessages({
  library(tidyverse)
  library(scater)
  library(WGCNA)
  library(Seurat)
  library(RColorConesa)
  library(acorde)
})
library(pbmcMultiome.SeuratData)

set.seed(000)

sim <- readRDS("~/workspace/mosim_paper/sim_6cells8clus8000_scMOSim_2groups_022.rds")

numcells <- "RNA"

# Normalize using Seurat's method
datG1R1 <- Seurat::CreateSeuratObject(sim$Group_1$Rep_1$`sim_scRNA-seq`@counts)
datG1R1 <- Seurat::SCTransform(datG1R1, verbose = FALSE, return.only.var.genes = FALSE)
datG1R1 <- as.data.frame(datG1R1[["SCT"]]$data)
names(datG1R1) <- paste0(names(datG1R1), "_G1R1")

datG1R2 <- Seurat::CreateSeuratObject(sim$Group_1$Rep_2$`sim_scRNA-seq`@counts)
datG1R2 <- Seurat::SCTransform(datG1R2, verbose = FALSE, return.only.var.genes = FALSE)
datG1R2 <- as.data.frame(datG1R2[["SCT"]]$data)
names(datG1R2) <- paste0(names(datG1R2), "_G1R2")

datG1R3 <- Seurat::CreateSeuratObject(sim$Group_1$Rep_3$`sim_scRNA-seq`@counts)
datG1R3 <- Seurat::SCTransform(datG1R3, verbose = FALSE, return.only.var.genes = FALSE)
datG1R3 <- as.data.frame(datG1R3[["SCT"]]$data)
names(datG1R3) <- paste0(names(datG1R3), "_G1R3")

datG2R1 <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)
datG2R1 <- Seurat::SCTransform(datG2R1, verbose = FALSE, return.only.var.genes = FALSE)
datG2R1 <- as.data.frame(datG2R1[["SCT"]]$data)
names(datG2R1) <- paste0(names(datG2R1), "_G2R1")

datG2R2 <- Seurat::CreateSeuratObject(sim$Group_2$Rep_2$`sim_scRNA-seq`@counts)
datG2R2 <- Seurat::SCTransform(datG2R2, verbose = FALSE, return.only.var.genes = FALSE)
datG2R2 <- as.data.frame(datG2R2[["SCT"]]$data)
names(datG2R2) <- paste0(names(datG2R2), "_G2R2")

datG2R3 <- Seurat::CreateSeuratObject(sim$Group_2$Rep_3$`sim_scRNA-seq`@counts)
datG2R3 <- Seurat::SCTransform(datG2R3, verbose = FALSE, return.only.var.genes = FALSE)
datG2R3 <- as.data.frame(datG2R3[["SCT"]]$data)
names(datG2R3) <- paste0(names(datG2R3), "_G2R3")

dat <- merge(datG1R1, datG1R2, by = 0, all = TRUE)
row.names(dat) <- dat$Row.names
dat$Row.names <- NULL
dat <- merge(dat, datG1R3, by = 0, all = TRUE)
row.names(dat) <- dat$Row.names
dat$Row.names <- NULL
dat <- merge(dat, datG2R1, by = 0, all = TRUE)
row.names(dat) <- dat$Row.names
dat$Row.names <- NULL
dat <- merge(dat, datG2R2, by = 0, all = TRUE)
row.names(dat) <- dat$Row.names
dat$Row.names <- NULL
dat <- merge(dat, datG2R3, by = 0, all = TRUE)
row.names(dat) <- dat$Row.names
dat$Row.names <- NULL

dat[is.na(dat)] <- 0

# Create cell-tocelltype ID table
ct <- tibble::tibble(Cell = colnames(datG1R1),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

# Now make metadata table
repeated_ct <- do.call(rbind, rep(list(ct), times = 6))
repeated_ct$Group <- rep(1:2, each = nrow(ct)*3)
repeated_ct$Replicate <- c(rep(1:3, each = nrow(ct)), rep(1:3, each = nrow(ct)))
# Now update the Cell part so it corresponds to the new names
repeated_ct$Cell <- colnames(dat)

library(factoextra)
pca_sc <- prcomp(t(dat))
fviz_eig(pca_sc)
ggsave(width = 3, height = 2, filename = paste0("./paper_plots/", numcells, "PCA_scree.pdf"))

# Plot by celltype

pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type)) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC1 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC2 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
#theme_classic(base_size = 10) +
  scale_color_conesa(palette = "main")

ggsave(width = 3.05, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC1_PC2_celltype.pdf"))


pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC3, colour = cell_type)) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC1 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC3 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
#theme_classic(base_size = 10) +
  scale_color_conesa(palette = "main")
ggsave(width = 3.05, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC1_PC3_celltype.pdf"))


pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC4, colour = cell_type)) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC1 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC4 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
#theme_classic(base_size = 10) +
  scale_color_conesa(palette = "main")
ggsave(width = 3.05, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC1_PC4_celltype.pdf"))


# Since the third component is also important to visualize it
pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC2, y = PC3, colour = cell_type)) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC2 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC3 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
#theme_classic(base_size = 10) +
  scale_color_conesa(palette = "main")

ggsave(width = 3.05, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC2_PC3_celltype.pdf"))

pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC3, y = PC4, colour = cell_type)) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC3 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC4 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
#theme_classic(base_size = 10) +
  scale_color_conesa(palette = "main")

ggsave(width = 3.05, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC3_PC4_celltype.pdf"))

# Plot by group

pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC2, colour = as.factor(Group))) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC1 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC2 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
  #theme_classic(base_size = 10) +
  scale_color_manual(values = c("#d0b0d4", "#64c298"))

ggsave(width = 3.235, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC1_PC2_group.pdf"))


pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC3, colour = as.factor(Group))) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC1 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC3 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
  #theme_classic(base_size = 10) +
  scale_color_manual(values = c("#d0b0d4", "#64c298"))
ggsave(width = 3.235, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC1_PC3_group.pdf"))


pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC4, colour = as.factor(Group))) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC1 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC4 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
  #theme_classic(base_size = 10) +
  scale_color_manual(values = c("#d0b0d4", "#64c298"))
ggsave(width = 3.235, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC1_PC4_group.pdf"))


# Since the third component is also important to visualize it
pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC2, y = PC3, colour = as.factor(Group))) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC2 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC3 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
  #theme_classic(base_size = 10) +
  scale_color_manual(values = c("#d0b0d4", "#64c298"))

ggsave(width = 3.235, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC2_PC3_group.pdf"))

pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC3, y = PC4, colour = as.factor(Group))) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC3 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC4 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
  #theme_classic(base_size = 10) +
  scale_color_manual(values = c("#d0b0d4", "#64c298"))

ggsave(width = 3.235, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC3_PC4_group.pdf"))


## Plot by replicate

pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC2, colour = as.factor(Replicate))) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC1 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC2 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
  #theme_classic(base_size = 10) +
  scale_color_manual(values = c("#713f98", "#d11d57", "#efefa1"))

ggsave(width = 3.45, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC1_PC2_replicate.pdf"))


pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC3, colour = as.factor(Replicate))) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC1 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC3 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
  #theme_classic(base_size = 10) +
  scale_color_manual(values = c("#713f98", "#d11d57", "#efefa1"))
ggsave(width = 3.45, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC1_PC3_replicate.pdf"))


pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC4, colour = as.factor(Replicate))) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC1 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC4 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
  #theme_classic(base_size = 10) +
  scale_color_manual(values = c("#713f98", "#d11d57", "#efefa1"))
ggsave(width = 3.45, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC1_PC4_replicate.pdf"))


# Since the third component is also important to visualize it
pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC2, y = PC3, colour = as.factor(Replicate))) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC2 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC3 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
  #theme_classic(base_size = 10) +
  scale_color_manual(values = c("#713f98", "#d11d57", "#efefa1"))

ggsave(width = 3.45, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC2_PC3_replicate.pdf"))

pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(repeated_ct, by = "Cell") %>% 
  ggplot(aes(x = PC3, y = PC4, colour = as.factor(Replicate))) + 
  geom_point(alpha = 0.5, size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  xlab(paste0("PC3 (", round(summary(pca_sc)$importance[2,1]*100, 1),"%)")) +
  ylab(paste0("PC4 (", round(summary(pca_sc)$importance[2,2]*100, 1),"%)")) +
  #theme_classic(base_size = 10) +
  scale_color_manual(values = c("#713f98", "#d11d57", "#efefa1"))

ggsave(width = 3.45, height = 1.7, filename = paste0("./paper_plots/", numcells, "PCA_PC3_PC4_replicate.pdf"))



