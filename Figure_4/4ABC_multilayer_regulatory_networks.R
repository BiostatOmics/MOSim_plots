# Load necessary libraries
library(MOSim)
library(igraph)
library(RCy3)
library(ComplexHeatmap)

# Do a small simulation of RNA-seq, miRNA-seq and TF data, ask for 3 replicates and 5 time points

set.seed(123)

nrep = 3
omics_list <- c("RNA-seq", "miRNA-seq")
omics_options <- c(omicSim('miRNA-seq', regulatorEffect = list('activator' = 0,'repressor' = 0.4, 'NE' = 0.6)))

ti = c(0:5) 
rnaseq_simulation <- mosim(omics = omics_list, omicsOptions = omics_options ,times = ti, depth = 30, numberReps = nrep, numberGroups = 2, diffGenes=0.05,TFtoGene = TRUE)
# Get the count tables
rnaseq_simulated <- omicResults(rnaseq_simulation, omics_list)
# get the settings used to generate each count table
all_settings <- omicSettings(rnaseq_simulation)
design_matrix <- experimentalDesign(rnaseq_simulation)


#Take only the first 100 DEGs
GE = rnaseq_simulated$`RNA-seq`
GE_alls = all_settings$`RNA-seq`[which(all_settings$`RNA-seq`$DE),]
genes = GE_alls$ID[1:100]

# Which are the regulators that present regulations
associations = all_settings[-1]
associations$`miRNA-seq` = associations$`miRNA-seq`[which(rowSums(is.na(associations$`miRNA-seq`[,c(3,4)]))!=2),]
associations$TF = associations$TF[which(rowSums(is.na(associations$TF[,c(3,4)]))!=2),]

#Take only the regulators that present regulations with the selected genes
sub_assoc = associations
sub_assoc$`miRNA-seq` = sub_assoc$`miRNA-seq`[which(genes %in% sub_assoc$`miRNA-seq`$Gene),]
sub_assoc$TF = sub_assoc$TF[which(genes %in% sub_assoc$TF$Gene),]

#Diferenciate per grupos

sub_assoc_G1 = sub_assoc
sub_assoc_G1$`miRNA-seq` = sub_assoc_G1$`miRNA-seq`[!is.na(sub_assoc_G1$`miRNA-seq`$Effect.Group1) & sub_assoc_G1$`miRNA-seq`$Effect.Group1 == 'repressor',]
sub_assoc_G1$TF = sub_assoc_G1$TF[!is.na(sub_assoc_G1$TF$Effect.Group1),]

sub_assoc_G2 = sub_assoc
sub_assoc_G2$`miRNA-seq` = sub_assoc_G2$`miRNA-seq`[!is.na(sub_assoc_G2$`miRNA-seq`$Effect.Group2) & sub_assoc_G2$`miRNA-seq`$Effect.Group2 == 'repressor',]
sub_assoc_G2$TF = sub_assoc_G2$TF[!is.na(sub_assoc_G2$TF$Effect.Group2),]

assoc_G1 = as.data.frame(rbind(sub_assoc_G1$`miRNA-seq`[,c(1,2,3)],sub_assoc_G1$TF[,c(1,2,3)]))
assoc_G2 = as.data.frame(rbind(sub_assoc_G2$`miRNA-seq`[,c(1,2,4)],sub_assoc_G2$TF[,c(1,2,4)]))

df_G1 = GE[unique(assoc_G1$Gene),,drop=FALSE]
df_G1 = rbind(df_G1, rnaseq_simulated$`miRNA-seq`[unique(sub_assoc_G1$`miRNA-seq`[,c(1,2,3)]$ID),,drop=FALSE])
df_G1 = rbind(df_G1, rnaseq_simulated$TF[unique(sub_assoc_G1$TF[,c(1,2,3)]$ID),,drop=FALSE])

# Create the mean across the replicates so it is not difficult to visualize the Figures
mean_columns <- list()

# Loop over columns in steps of 3 (since every 3 columns represent replicates)
for (i in seq(1, ncol(df_G1), by = 3)) {
  # Calculate the mean of every set of 3 columns (replicates)
  mean_col <- rowMeans(df_G1[, i:(i+2)], na.rm = TRUE)
  
  # Extract the column name for the first replicate to use it as the new column name
  new_col_name <- gsub("Rep\\d+", "Mean", colnames(df_G1)[i])
  
  # Store the result in the list, naming it appropriately
  mean_columns[[new_col_name]] <- mean_col
}

# Convert the list of mean columns into a dataframe
df_G1 <- as.data.frame(mean_columns)
df_G1 = df_G1[,grep('Group1', colnames(df_G1)),drop=FALSE]
df_G1 = scale(t(df_G1), center =TRUE, scale = TRUE)
data_matrix1 = t(as.matrix(df_G1))

df_G1 = as.data.frame(rbind(sub_assoc_G1$`miRNA-seq`[,c(1,2,3)],sub_assoc_G1$TF[,c(1,2,3)]))


df_G2 = GE[unique(assoc_G2$Gene),,drop=FALSE]
df_G2 = rbind(df_G2, rnaseq_simulated$`miRNA-seq`[unique(sub_assoc_G2$`miRNA-seq`[,c(1,2,4)]$ID),,drop=FALSE])
df_G2 = rbind(df_G2, rnaseq_simulated$TF[unique(sub_assoc_G2$TF[,c(1,2,4)]$ID),,drop=FALSE])

# Create the mean across the replicates so it is not difficult to visualize the Figures
mean_columns <- list()

# Loop over columns in steps of 3 (since every 3 columns represent replicates)
for (i in seq(1, ncol(df_G2), by = 3)) {
  # Calculate the mean of every set of 3 columns (replicates)
  mean_col <- rowMeans(df_G2[, i:(i+2)], na.rm = TRUE)
  
  # Extract the column name for the first replicate to use it as the new column name
  new_col_name <- gsub("Rep\\d+", "Mean", colnames(df_G2)[i])
  
  # Store the result in the list, naming it appropriately
  mean_columns[[new_col_name]] <- mean_col
}

# Convert the list of mean columns into a dataframe
df_G2 <- as.data.frame(mean_columns)
df_G2 = df_G2[,grep('Group2', colnames(df_G2)),drop=FALSE]
df_G2 = scale(t(df_G2), center =TRUE, scale = TRUE)
data_matrix2 = t(as.matrix(df_G2))

df_G2 = as.data.frame(rbind(sub_assoc_G2$`miRNA-seq`[,c(1,2,4)],sub_assoc_G2$TF[,c(1,2,4)]))


##### Figure 4 A ######
nodes_G1 = data.frame(id = c(unique(df_G1[,'Gene']),unique(df_G1[,'ID'])),
                      omic = c(rep('gene',length(unique(df_G1[,'Gene']))),rep('miRNA', length(unique(sub_assoc_G1$`miRNA-seq`$ID))),rep('TF',length(unique(sub_assoc_G1$TF$ID)))))
interactions_G1 = data.frame(from = df_G1[,'Gene'], to = df_G1[,'ID'], inter= df_G1[,3])
ig_G1 = igraph::graph_from_data_frame(interactions_G1, vertices = nodes_G1, directed = FALSE)
mem_G1 = igraph::components(ig_G1)$membership

RCy3::createNetworkFromIgraph(ig_G1)

##### Figure 4 B ######

nodes_G2 = data.frame(id = c(unique(df_G2[,'Gene']),unique(df_G2[,'ID'])),
                      omic = c(rep('gene',length(unique(df_G2[,'Gene']))),rep('miRNA', length(unique(sub_assoc_G2$`miRNA-seq`$ID))),rep('TF',length(unique(sub_assoc_G2$TF$ID)))))
interactions_G2 = data.frame(from = df_G2[,'Gene'], to = df_G2[,'ID'], inter= df_G2[,3])
ig_G2 = igraph::graph_from_data_frame(interactions_G2, vertices = nodes_G2, directed = FALSE)
mem_G2 = igraph::components(ig_G2)$membership

RCy3::createNetworkFromIgraph(ig_G2)

##### Figure 4 C ######

data_matrix = data_matrix1[names(mem_G1[order(mem_G1)]),,drop=FALSE]
colnames(data_matrix) = gsub('.Mean','',colnames(data_matrix))

ra = rowAnnotation(
  subnetwork = mem_G1[order(mem_G1)],
  omic = c(rep('gene', 41),rep('miRNA', 2), 'TF', rep("gene", 44),'miRNA',rep("TF", 2)),
  col = list(subnetwork = c("1" = "#74CDF0", "2" = "#FDA3D1"),
             omic = c("gene" = "#15918A", "miRNA" = "#FDC659","TF" = "#F58A53"))
)

Heatmap(data_matrix,  
        right_annotation = ra, 
        col = colorRampPalette(c("blue", "white", "red"))(50), 
        cluster_rows = FALSE, cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 5),
        column_title  = "Heatmap of Group1", name='Expression')

# Do KNN to order the genes and visualiza better the Heatmap
net2_gen = names(mem_G1[which(mem_G1==2)])
dt_matrix = data_matrix[net2_gen,,drop=FALSE]

net1_gen = names(mem_G1[which(mem_G1==1)])
dt_matrix1 = data_matrix[net1_gen,,drop=FALSE]

ra_net2 = rowAnnotation(
  omic = c( rep("gene", 44),'miRNA',rep("TF", 2)),
  col = list(omic = c("gene" = "#15918A", "miRNA" = "#FDC659","TF" = "#F58A53"))
)

ra_net1 = rowAnnotation(
  omic = c(rep('gene', 41),rep('miRNA', 2), 'TF'),
  col = list(
    omic = c("gene" = "#15918A", "miRNA" = "#FDC659","TF" = "#F58A53"))
)

ht<-Heatmap(dt_matrix,  
            right_annotation = ra_net2, 
            col = colorRampPalette(c("blue", "white", "red"))(50), 
            cluster_rows = TRUE, cluster_columns = FALSE, 
            row_names_gp = gpar(fontsize = 5),
            column_title  = "Heatmap of Group1", name='Expression')

ht1<-Heatmap(dt_matrix1,  
             right_annotation = ra_net1, 
             col = colorRampPalette(c("blue", "white", "red"))(50), 
             cluster_rows = TRUE, cluster_columns = FALSE, 
             row_names_gp = gpar(fontsize = 5),
             column_title  = "Heatmap of Group1", name='Expression')

#Extract the re-ordering of the data_matrix to create a heatmap that is easier to look at

#Extraido por row_order(ht) pero aquí escribo a mano el orden por el hecho de que el knn no pueda dar resultados distintos

my_order = c( 25, 35, 30, 33, 39, 36, 29, 27, 43, 42, 28, 9, 31, 44, 8, 32, 40, 
              11, 12, 37, 3, 17, 6, 41, 7, 23, 19, 2, 14, 26, 20, 34, 21, 24, 4, 13, 5, 1, 
              10, 22, 15, 16, 18, 38, 45, 46, 47)

dt_matrix = dt_matrix[my_order,,drop=FALSE]

my_order1 = c(14, 37, 15, 26, 17, 39, 13, 38, 16, 18, 23, 30, 12, 8, 32, 35, 10, 
              29, 21, 25, 5, 40, 34, 27, 33, 11, 19, 20, 24, 28, 22, 9, 1, 4, 6, 2, 7, 36, 31, 
              3, 41, 42, 43, 44)

dt_matrix1 = dt_matrix1[my_order1,,drop=FALSE]

data_matrix = rbind(dt_matrix1,dt_matrix)

ra2 = rowAnnotation(
  subnetwork = c(rep('1',44),rep('2',47)),
  omic = c(rep('gene', 41),rep('miRNA', 2), 'TF', rep("gene", 44),'miRNA',rep("TF", 2)),
  col = list(subnetwork = c("1" = "#74CDF0", "2" = "#FDA3D1"),
             omic = c("gene" = "#15918A", "miRNA" = "#FDC659","TF" = "#F58A53"))
)

Heatmap(data_matrix,  
        right_annotation = ra2, 
        col = colorRampPalette(c("blue", "white", "red"))(50), 
        cluster_rows = FALSE, cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 5),
        column_title  = "Heatmap of Group1", name='Expression')

Heatmap(data_matrix,  
        right_annotation = ra2, 
        col = colorRampPalette(c("blue", "white", "red"))(50), 
        cluster_rows = FALSE, cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),  # Adjust the font size of column names if necessary
        width = unit(5, "cm"),  # Adjust the width of the heatmap to make the cells narrower
        column_title = "Heatmap of Group1", name='Expression')

#Do the same for Group 2

data_matrix = data_matrix2[names(mem_G2[order(mem_G2)]),,drop=FALSE]
colnames(data_matrix) = gsub('.Mean','',colnames(data_matrix))

ra = rowAnnotation(
  subnetwork = mem_G2[order(mem_G2)],
  omic = c(rep('gene', 70),rep('miRNA', 2), rep('TF',3), rep("gene", 22),'miRNA'),
  col = list(network = c("1" = "#74CDF0", "2" = "#FDA3D1"),
             omic = c("gene" = "#15918A", "miRNA" = "#FDC659","TF" = "#F58A53"))
)

Heatmap(data_matrix,  
        right_annotation = ra, 
        col = colorRampPalette(c("blue", "white", "red"))(50), 
        cluster_rows = FALSE, cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 5),
        column_title  = "Heatmap of GroupB", name='Expression')

net1G2_gen = names(mem_G2[which(mem_G2==1)])
net2G2_gen = names(mem_G2[which(mem_G2==2)])
dt_matrix1 = data_matrix[net1G2_gen,,drop=FALSE]
dt_matrix2 = data_matrix[net2G2_gen,,drop=FALSE]

ra1 = rowAnnotation(
  omic = c(rep('gene', 70),rep('miRNA', 2), rep('TF',3)),
  col = list(
    omic = c("gene" = "#15918A", "miRNA" = "#FDC659","TF" = "#F58A53"))
)

ra2 = rowAnnotation(
  omic = c( rep("gene", 22),'miRNA'),
  col = list(
    omic = c("gene" = "#15918A", "miRNA" = "#FDC659","TF" = "#F58A53"))
)

ht1<-Heatmap(dt_matrix1,  
             right_annotation = ra1, 
             col = colorRampPalette(c("blue", "white", "red"))(50), 
             cluster_rows = TRUE, cluster_columns = FALSE, 
             row_names_gp = gpar(fontsize = 5),
             column_title  = "Heatmap of GroupB", name='Expression')

ht2<-Heatmap(dt_matrix2,  
             right_annotation = ra2, 
             col = colorRampPalette(c("blue", "white", "red"))(50), 
             cluster_rows = TRUE, cluster_columns = FALSE, 
             row_names_gp = gpar(fontsize = 5),
             column_title  = "Heatmap of GroupB", name='Expression')

#Extract the re-ordering of the data_matrix to create a heatmap that is easier to look at

#Extraido por row_order(ht) pero aquí escribo a mano el orden por el hecho de que el knn no pueda dar resultados distintos

my_order1 = c(15, 19, 61, 37, 43, 68, 40, 58, 30, 53, 33, 49, 51, 56, 38, 67, 57, 55, 
              25, 32, 28, 70, 27, 18, 48, 66, 42, 34, 26, 14, 69, 44, 45, 7, 24, 20, 
              62, 13, 46, 36, 12, 2, 11, 5, 41, 35, 8, 3, 31, 47, 63, 10, 23, 39, 22, 29, 
              65, 64, 50, 52, 6, 16, 1, 54, 60, 9, 4, 21, 17, 59, 71, 72, 73, 74, 75)

dt_matrix1 = dt_matrix1[my_order1,,drop=FALSE]

my_order2 = c(6, 17, 14, 15, 2, 16, 4, 1, 18, 20, 7, 21, 13, 8, 11, 3, 9, 5, 22, 10, 19, 12,23)

dt_matrix2 = dt_matrix2[my_order2,,drop=FALSE]

data_matrix= rbind(dt_matrix1,dt_matrix2)


ra = rowAnnotation(
  subnetwork = c(rep('1',75),rep('2',23)),
  omic =  c(rep('gene', 70),rep('miRNA', 2), rep('TF',3), rep("gene", 22),'miRNA'),
  col = list(subnetwork = c("1" = "#74CDF0", "2" = "#FDA3D1"),
             omic = c("gene" = "#15918A", "miRNA" = "#FDC659","TF" = "#F58A53"))
)

Heatmap(data_matrix,  
        right_annotation = ra, 
        col = colorRampPalette(c("blue", "white", "red"))(50), 
        cluster_rows = FALSE, cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),  # Adjust the font size of column names if necessary
        width = unit(5, "cm"),
        column_title  = "Heatmap of Group2", name='Expression')



