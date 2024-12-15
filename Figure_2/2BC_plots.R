##########   Generating plots for MOSim paper  ##############

## By Sonia Tarazona
## Created: 4-oct-2022



# Set-up ------------------------------------------------------------------

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.15")
# 
# BiocManager::install("MOSim")
# BiocManager::install("NOISeq")

library(MOSim)
library(NOISeq)
library(ggplot2)
library(grid)
library(gridExtra)



# Explorando datos -------------------------------------------------

archivos = dir()
archivos = archivos[grep("RData", archivos)]

objetos = vector("list", length = length(archivos))
names(objetos) = archivos
for (a in archivos) {
  objetos[[a]] = load(a, verbose = TRUE)
}

simulaciones = archivos[grep("simu", archivos)]

i = 7
archivos[i]
load(archivos[i], verbose = TRUE)
names(simu.object@simulators)   # count matrices for RNAseq, miRNA, ATAC, chipseq
simu.object@totalGenes  # 38293
simu.object@diffGenes   # 5744
simu.object@numberReps  # 3
simu.object@numberGroups # 2
simu.object@times  # 0-5
# simu.object@geneNames
# simu.object@simSettings   ### NECESARIO!!! contiene muchas cosas
# simu.object@profiles
simu.object@profileProbs
simu.object@noiseParams
simu.object@depth
simu.object@minMaxQuantile
simu.object@minMaxFC
simu.object@TFtoGene  # No TFs




# Figure 2: PCA -----------------------------------------------------------

i = 1
archivos[i]
load(archivos[i], verbose = TRUE)

# RNA-seq
datos = simu.object@simulators$SimRNAseq@simData
disseny = do.call("rbind", strsplit(colnames(datos), split = ".", fixed = TRUE))
disseny = as.data.frame(disseny)
colnames(disseny) = c("Group", "Time", "Rep")
rownames(disseny) = colnames(datos)
disseny$GroupTime = apply(disseny[,1:2], 1, paste, collapse = "_")
datosFilt = filtered.data(datos, factor = disseny$GroupTime, norm = FALSE, 
                          method = 1, cv.cutoff = 500, cpm = 1)
miPCA = PCA.GENES(scale(t(log(datosFilt+1))))

# DEinfo = simu.object@simSettings$geneProfiles$SimRNAseq
# datosDE = datos[DEinfo$DE,]
# miPCA = PCA.GENES(scale(t(log(datosDE+1))))

dat4plot = data.frame(miPCA$scores[,1:2], disseny)
colnames(dat4plot)[1:2] = paste0("PC", 1:2)
p1 <- ggplot(data = dat4plot, aes(x = PC1, y = PC2, color = Time, shape = Group)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(alpha = 1.3) + ggtitle("RNA-seq") + 
  xlab(paste0("PC1 (", round(miPCA$var.exp[1,1]*100,1),"%)")) +
  ylab(paste0("PC2 (", round(miPCA$var.exp[2,1]*100,1),"%)"))


# ATAC-seq
datos = simu.object@simulators$SimDNaseseq@simData
disseny = do.call("rbind", strsplit(colnames(datos), split = ".", fixed = TRUE))
disseny = as.data.frame(disseny)
colnames(disseny) = c("Group", "Time", "Rep")
rownames(disseny) = colnames(datos)
disseny$GroupTime = apply(disseny[,1:2], 1, paste, collapse = "_")
datosFilt = filtered.data(datos, factor = disseny$GroupTime, norm = FALSE, 
                          method = 1, cv.cutoff = 500, cpm = 1)
miPCA = PCA.GENES(scale(t(log(datosFilt+1))))

dat4plot = data.frame(miPCA$scores[,1:2], disseny)
colnames(dat4plot)[1:2] = paste0("PC", 1:2)
p2 <- ggplot(data = dat4plot, aes(x = PC1, y = PC2, color = Time, shape = Group)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(alpha = 1.3) + ggtitle("ATAC-seq") + 
  xlab(paste0("PC1 (", round(miPCA$var.exp[1,1]*100,1),"%)")) +
  ylab(paste0("PC2 (", round(miPCA$var.exp[2,1]*100,1),"%)"))


## Final figure

pdf("Figure2B_PCA.pdf", height = 4, width = 10)
grid.arrange(p1,p2, nrow = 1)
dev.off()




# Figura 4: Box plot for correlation -----------------------------------

genes = simu.object@simSettings$geneProfiles$SimRNAseq
head(genes)
rownames(genes) = genes$ID
regulations = simu.object@simSettings$geneProfiles$SimDNaseseq
regulations = regulations[!is.na(regulations$Gene),]
head(regulations)
regulations$G1gene = genes[regulations$Gene, "Group1"]
table(regulations$G1gene, useNA = "i")
table(regulations$G1gene, regulations$Group1, useNA = "i")
flat2 = (regulations$G1gene == "flat")*(regulations$Group1 == "flat")
table(flat2, useNA = "i")
regulations = regulations[flat2 != 1,]
table(regulations$Effect.Group1, useNA = "i")
regulations$EffectG1 = !is.na(regulations$Effect.Group1)
table(regulations$EffectG1)

rnaseqG1 = simu.object@simulators$SimRNAseq@simData
rnaseqG1 = rnaseqG1[,grep("Group1", colnames(rnaseqG1))]
rnaseqG1 =t(aggregate(t(rnaseqG1), by = list(time = rep(0:5, each = 3)), mean))
rnaseqG1 = rnaseqG1[-1,]
atacseqG1 = simu.object@simulators$SimDNaseseq@simData
atacseqG1 = atacseqG1[,grep("Group1", colnames(atacseqG1))]
atacseqG1 =t(aggregate(t(atacseqG1), by = list(time = rep(0:5, each = 3)), mean))
atacseqG1 = atacseqG1[-1,]

corre = sapply(1:nrow(regulations), 
               function(i) {
                 cor(as.numeric(rnaseqG1[regulations$Gene[i],]),
                     as.numeric(atacseqG1[regulations$ID[i],]))
                 })
regulations$corre = corre
regulations$absCor = abs(corre)

quitar = (regulations$Group1 == "flat")*(regulations$EffectG1)
regulations = regulations[quitar == 0,]

ggplot(regulations, aes(x=EffectG1, y=absCor)) + 
  geom_violin(trim=FALSE)

p = ggplot(regulations, aes(x=EffectG1, y=absCor, fill = EffectG1)) + 
  geom_boxplot(outlier.shape=16, outlier.size=2, notch=TRUE) +
   scale_fill_brewer(palette="Dark2") + theme_classic() + theme(legend.position="none") +
  labs(title = "Gene regulation by ATAC-seq", x = "Regulation", y = "Absolute correlation")
p

pdf("Figure2C_Correlations.pdf", height = 5, width = 3.7)
p
dev.off()

png("Figure2C_Correlations.png", height = 480, width = 370, pointsize = 30)
p
dev.off()

correTRUE = regulations$absCor[regulations$EffectG1]
100*mean(correTRUE > 0.9)
100*mean(correTRUE <= 0.9)
