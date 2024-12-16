##########   Briefings in Bioinformatics  ##############

## By Maider Aguerralde & Sonia Tarazona
## Created: Dec-2024



# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("MOSim")
# BiocManager::install("NOISeq")


library(MOSim)
library(NOISeq)



# Preparing data -----------------------------------------------------------

load("latest_simu.RData")  


# RNA-seq
datos = simu.object@simulators$SimRNAseq@simData
disseny = do.call("rbind", strsplit(colnames(datos), split = ".", fixed = TRUE))
disseny = as.data.frame(disseny)
colnames(disseny) = c("Group", "Time", "Rep")
rownames(disseny) = colnames(datos)
disseny$GroupTime = apply(disseny[,1:2], 1, paste, collapse = "_")
datosFilt = filtered.data(datos, factor = disseny$GroupTime, norm = FALSE, 
                          method = 1, cv.cutoff = 500, cpm = 1)

DEGs = simu.object@simSettings$geneProfiles$SimRNAseq


# Group1: ALL + NOCENTER + SCALED

noised.nondes <- simu.object@simSettings$featureSamples[[paste0("Sim", gsub("-", "", "RNAseq"))]]$noiseNonDE
simu.raw <- datosFilt[! rownames(datosFilt) %in% noised.nondes, ]
simu.data.group1.raw <- simu.raw[,grep("Group1", colnames(simu.raw))]

temps = disseny$Time[disseny$Group == "Group1"]
simu.data.rep <- aggregate(t(simu.data.group1.raw), by = list(temps), mean)
rownames(simu.data.rep) = simu.data.rep$Group.1
simu.data.rep = t(simu.data.rep[,-1])

grupo1 <- t(scale(t(simu.data.rep), center = FALSE, scale = TRUE))
grupo1 = na.omit(grupo1)





# Optimal number of clusters ----------------------------------------------

library(factoextra)

fviz_nbclust(x = grupo1, FUNcluster = kmeans, 
             method = "silhouette", k.max = 20, verbose = FALSE,
             iter.max = 20) + labs(title = "Num. optimo clusters") 
fviz_nbclust(x = grupo1, FUNcluster = kmeans, 
                  method = "wss", k.max = 20, verbose = FALSE,
                  iter.max = 20) + labs(title = "Num. optimo clusters")

kmeans2 <- kmeans(x = grupo1, centers = 7, nstart = 20, iter.max = 20)




# Figure 2A ---------------------------------------------------------------

perfiles2 = aggregate(grupo1, by = list(kmeans2$cluster), mean)

library(RColorBrewer)
colores = brewer.pal(8, "Set3")[-2]
nombres = c("Cont.Repr.", "Trans.Induct.", "Trans.Induct.", 
            "Trans.Repr.", "Trans.Repr.", "Cont.Induct.", "Flat")
par(mar = c(4,2,1,1))
matplot(0:5, t(perfiles2[,-1]), type = "l", col = colores, lty = 1, lwd = 3,
        xlab = "Time points", ylab = "", main = "")
legend("topright", legend = nombres, col = colores[1:7], lty = 1, lwd = 3, ncol = 3, bty = "n")





# Classification error ----------------------------------------------------

DEGs2 = DEGs[DEGs$ID %in% rownames(grupo1), ]

table(DEGs2$Group1, kmeans2$cluster)

# errores: 0.35%
100*(11+4+13+16+3)/13337  



