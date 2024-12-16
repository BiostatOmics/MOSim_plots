##########   Briefings in Bioinformatics  ##############

## By Maider Aguerralde & Sonia Tarazona
## Created: Dec-2024


# install.packages("Seurat")
# install.packages("SeuratObject")
library(SeuratObject)
library(Seurat)


# Preparing data -----------------------------------------------------------

sim = readRDS("C:/Users/Sonia/Documents/sim_6cells8clus10000_scMOSim_2groups_byHand2.rds")

cell_types <- sim$cellTypes


### COUNTS group 1 replicate 1

rna = sim$Group_1$Rep_1$`sim_scRNA-seq`@counts
rna <- as.data.frame(rna)


### Cells & cell types
ct <- data.frame(Cell = colnames(sim$Group_1$Rep_1$`sim_scRNA-seq`@counts),
                 cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))


### Genes and profiles

perfils <- sim$AssociationMatrices$AssociationMatrix_Group_1[sim$AssociationMatrices$AssociationMatrix_Group_1$Gene_cluster %in% c(1:8), c("Gene_ID", "Gene_cluster")]
perfils = as.data.frame(perfils)

rna = rna[perfils$Gene_ID,]



### Averaging cells per cell type
pseudobulk = t(rna)
pseudobulk = pseudobulk[ct$Cell,]
pseudobulk = aggregate(pseudobulk, by = list(ct = ct$cell_type), mean)
rownames(pseudobulk) = pseudobulk$ct
pseudobulk = t(pseudobulk[,-1])


## Removing flat genes
pseudobulk = pseudobulk[perfils$Gene_cluster != 7,]
perfils = perfils[perfils$Gene_cluster != 7,]







# Optimal number of clusters ----------------------------------------------

library(cluster)
library(factoextra)

distS = stats::cor(t(pseudobulk), method = "spearman")

fviz_nbclust(x = pseudobulk, FUNcluster = pam, diss = 1-distS,
             method = "silhouette", k.max = 15, verbose = FALSE)
fviz_nbclust(x = pseudobulk, FUNcluster = pam, diss = 1-distS,
             method = "wss", k.max = 15, verbose = FALSE)


# k=10
pamSopt = pam(1-distS, diss = TRUE, k = 10, nstart = 10)




             
# Figure 3A ---------------------------------------------------------------

clusters = pamSopt$clustering
clusters = factor(clusters, levels = c(2,3,5,6,1,  7,4,8:10),
                  labels = c("1A", "1B", "2A", "2B", "3",
                             "4", "5", "6", "7A", "7B"))

perfilClust = aggregate(pseudobulk, by = list(clusters), mean)

library(RColorBrewer)
colores = brewer.pal(6, "Set3")[-2]

par(mar = c(5,2,1,1), mfrow = c(1,2))

matplot(1:6, t(perfilClust[1:5,-1]), type = "l", col = colores, lty = 1, lwd = 3,
        xlab = "", ylab = "", main = "", ylim = c(0,1.1), xaxt = "n")
axis(side = 1, at = 1:6, labels = colnames(perfilClust)[-1], 
     las = 2, cex.axis = 0.8)
legend("top", legend = levels(clusters)[1:5], col = colores, lty = 1, lwd = 3, 
       ncol = 3, bty = "n")

matplot(1:6, t(perfilClust[6:10,-1]), type = "l", col = colores, lty = 1, lwd = 3,
        xlab = "", ylab = "", main = "", ylim = c(0,1.8), xaxt = "n")
axis(side = 1, at = 1:6, labels = colnames(perfilClust)[-1], 
     las = 2, cex.axis = 0.8)
legend("topright", legend = levels(clusters)[6:10], col = colores, lty = 1, lwd = 3, 
       ncol = 2, bty = "n")



table(perfils$Gene_cluster, clusters)
100*(155+57+2+16)/7000  # 3.3%
