
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This is the second step in processing the scRNA-seq  
# primary breast datasets.
# You will need to have each dataset processed as Seurat objects
# and filtered from "Individual_Dataset_Loading.R"

# The output of this R script will be the input 
# of "Integration.R"

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PrelimObj_Dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/PrimaryBreast_Processing/Individual_Dataset_Objects/Unfiltered_and_PrelimQC"
DoubletResult_Dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/PrimaryBreast_Processing/Individual_Dataset_Objects/withDoubletMeta"

# Libraries -------------------------------------------------------------
library(BiocManager)
library(DoubletFinder)
library(GEOquery) 
library(plyr)
library(dplyr) 
library(Matrix)
library(devtools)
library(Seurat)
library(ggplot2)
library(cowplot) 
library(SAVER) 
library(metap)
library(multtest)
library(msigdbr)
library(fgsea)
library(monocle3)
library(velocyto.R)
library(loomR)
library(clustree)
library(tibble)
library(devtools)
library(harmony)
library(SeuratData)
library(UCell)
library(glmGamPoi)
library(sctransform)
library(matrixStats)
library(sparseMatrixStats)
library(DESeq2)
library(genefu)




#AziziPrim ====================

setwd(PrelimObj_Dir)
Azizi_prim <- readRDS(file = "AziziPRIMimmune_62622.rds")
Azizi_prim$orig.ident <- "Aziziimmune"
DefaultAssay(Azizi_prim) <- "RNA"

hsp.genes <- grep(pattern = "^HSP", x = rownames(Azizi_prim@assays[["RNA"]]), value = TRUE)
print(hsp.genes)
Azizi_prim <- PercentageFeatureSet(Azizi_prim, "^HSP", col.name = "percent.heatshock")


combo <- Azizi_prim
combo <- SCTransform(combo, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet", "percent.heatshock"), verbose = TRUE)
combo <- RunPCA(combo, npcs = 100, ndims.print = 1:10, nfeatures.print = 10)
DimHeatmap(combo, dims = 25:30, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 30:40, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 40:50, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 50:60, cells = 500, balanced = TRUE)

combo <- FindNeighbors(combo, reduction = "pca", dims = 1:60)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
combo <- FindClusters(combo, resolution = resolution.range)
clustree(combo, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")


combo <- FindNeighbors(combo, reduction = "pca", dims = 1:60, nn.method = "rann") 
combo <- FindClusters(combo, resolution = 0.7)
combo <- RunUMAP(object = combo, reduction = "pca", dims = 1:60, seed.use = 123)

combo.list <- SplitObject(combo, split.by = "samples")

pkchoose.list <- list()
for (i in 1:length(combo.list)) {
  seu_temp <- combo.list[[i]]
  
  sweep.res.list <- paramSweep_v3(seu_temp, PCs = 
                                    seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  #sweep.stats.list[[i]] <- sweep.stats
  sweep.stats.pk <- find.pK(sweep.stats)
  pK=as.numeric(as.character(sweep.stats.pk$pK))
  BCmetric=sweep.stats.pk$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  nExp_poi <- 0.0059*nrow(seu_temp@meta.data) #azizi 0.59% doublet rate
  seu_temp <- doubletFinder_v3(seu_temp, PCs = 
                                 seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE, pN = 0.25, pK = as.numeric(pK_choose), nExp = nExp_poi, reuse.pANN = FALSE)
  combo.list[[i]] <- seu_temp
}

doublet.cell.list <- list()
doublet.df <- list()
for (i in names(combo.list)) {
  seu_temp <- combo.list[[i]]
  doublet.sub <- seu_temp@meta.data[ ,grepl("DF", colnames(seu_temp@meta.data)), drop = F]
  head(doublet.sub)
  doublet.df[[i]] <- doublet.sub
  doublet.cell.list[[i]] <- rownames(doublet.sub) 
  Azizi_prim <- AddMetaData(Azizi_prim, doublet.sub, paste0("Doublet.", i))
}

colnames(Azizi_prim@meta.data)
Azizi_prim$Doublet.Call <- apply(Azizi_prim@meta.data[,c(35:62)], 1, function(x) x[!is.na(x)][1])
Azizi_prim@meta.data <- Azizi_prim@meta.data[,-c(35:62)]
table(Azizi_prim$Doublet.Call)

setwd(DoubletResult_Dir)
saveRDS(Azizi_prim, "AziziPrim_withDoublet_62622.rds")


#AziziT =====================================================

setwd(PrelimObj_Dir)
AziziT <- readRDS(file = "AziziT_filtered_62622.rds")

AziziT$orig.ident <- "AziziT"
AziziT$celltype_minor <- "T Cell"
DefaultAssay(AziziT) <- "RNA"

hsp.genes <- grep(pattern = "^HSP", x = rownames(AziziT@assays[["RNA"]]), value = TRUE)
print(hsp.genes)
AziziT <- PercentageFeatureSet(AziziT, "^HSP", col.name = "percent.heatshock")

combo <- AziziT
combo <- SCTransform(combo, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet","percent.heatshock"), verbose = TRUE)
combo <- RunPCA(combo, npcs = 100, ndims.print = 1:10, nfeatures.print = 10)
DimHeatmap(combo, dims = 25:30, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 30:40, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 40:50, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 50:60, cells = 500, balanced = TRUE)

combo <- FindNeighbors(combo, reduction = "pca", dims = 1:60)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
combo <- FindClusters(combo, resolution = resolution.range)
clustree(combo, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")


combo <- FindNeighbors(combo, reduction = "pca", dims = 1:60, nn.method = "rann")
combo <- FindClusters(combo, resolution = 0.7)
combo <- RunUMAP(object = combo, reduction = "pca", dims = 1:60, seed.use = 123)

combo.list <- SplitObject(combo, split.by = "samples")

pkchoose.list <- list()
for (i in 1:length(combo.list)) {
  seu_temp <- combo.list[[i]]
  
  sweep.res.list <- paramSweep_v3(seu_temp, PCs = 
                                    seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  #sweep.stats.list[[i]] <- sweep.stats
  sweep.stats.pk <- find.pK(sweep.stats)
  pK=as.numeric(as.character(sweep.stats.pk$pK))
  BCmetric=sweep.stats.pk$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  nExp_poi <- 0.054*nrow(seu_temp@meta.data) ## 5.4% doublet rate
  seu_temp <- doubletFinder_v3(seu_temp, PCs = 
                                 seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE, pN = 0.25, pK = as.numeric(pK_choose), nExp = nExp_poi, reuse.pANN = FALSE)
  combo.list[[i]] <- seu_temp
}

doublet.cell.list <- list()
doublet.df <- list()
for (i in names(combo.list)) {
  seu_temp <- combo.list[[i]]
  doublet.sub <- seu_temp@meta.data[ ,grepl("DF", colnames(seu_temp@meta.data)), drop = F]
  head(doublet.sub)
  doublet.df[[i]] <- doublet.sub
  doublet.cell.list[[i]] <- rownames(doublet.sub) 
  AziziT <- AddMetaData(AziziT, doublet.sub, paste0("Doublet.", i))
}

colnames(AziziT@meta.data)
AziziT$Doublet.Call <- apply(AziziT@meta.data[,c(32:36)], 1, function(x) x[!is.na(x)][1])
AziziT@meta.data <- AziziT@meta.data[,-c(32:36)]
table(AziziT$Doublet.Call)

setwd(DoubletResult_Dir)
saveRDS(AziziT, "AziziT_withDoublet_62622.rds")


#Karaayvaz ================================

setwd(PrelimObj_Dir)
Karaayvaz <- readRDS(file = "Karaayvaz_filtered_62622.rds")
Karaayvaz$orig.ident <- "Karaayvaz"
Karaayvaz <- subset(x = Karaayvaz, subset = Patient == "PT084", invert = TRUE)
DefaultAssay(Karaayvaz) <- "RNA"

hsp.genes <- grep(pattern = "^HSP", x = rownames(Karaayvaz@assays[["RNA"]]), value = TRUE)
print(hsp.genes)
Karaayvaz <- PercentageFeatureSet(Karaayvaz, "^HSP", col.name = "percent.heatshock")

combo <- Karaayvaz
combo <- SCTransform(combo, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet", "percent.heatshock"), verbose = TRUE)
combo <- RunPCA(combo, npcs = 100, ndims.print = 1:10, nfeatures.print = 10)
DimHeatmap(combo, dims = 10:25, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 25:30, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 30:40, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 40:50, cells = 500, balanced = TRUE)

combo <- FindNeighbors(combo, reduction = "pca", dims = 1:30)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
combo <- FindClusters(combo, resolution = resolution.range)
clustree(combo, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")


combo <- FindNeighbors(combo, reduction = "pca", dims = 1:30, nn.method = "rann") 
combo <- FindClusters(combo, resolution = 1.6)
combo <- RunUMAP(object = combo, reduction = "pca", dims = 1:30, seed.use = 123)

PT039_P1 <- subset(combo, subset = samples == "PT039_P1")
PT081_P3 <- subset(combo, subset = samples == "PT081_P3")
PT081_P5 <- subset(combo, subset = samples == "PT081_P5")
PT126_P7 <- subset(combo, subset = samples == "PT126_P7")
combo <- subset(combo, subset = samples == "PT039_P1", invert = T)
combo <- subset(combo, subset = samples == "PT081_P3", invert = T)
combo <- subset(combo, subset = samples == "PT081_P5", invert = T)
combo <- subset(combo, subset = samples == "PT126_P7", invert = T)

combo.list <- SplitObject(combo, split.by = "samples")

pkchoose.list <- list()
for (i in 1:length(combo.list)) {
  seu_temp <- combo.list[[i]]
  
  sweep.res.list <- paramSweep_v3(seu_temp, PCs = 
                                    seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  #sweep.stats.list[[i]] <- sweep.stats
  sweep.stats.pk <- find.pK(sweep.stats)
  pK=as.numeric(as.character(sweep.stats.pk$pK))
  BCmetric=sweep.stats.pk$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  nExp_poi <- 0.01*nrow(seu_temp@meta.data) #1% doublet rate
  seu_temp <- doubletFinder_v3(seu_temp, PCs = 
                                 seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE, pN = 0.25, pK = as.numeric(pK_choose), nExp = nExp_poi, reuse.pANN = FALSE)
  combo.list[[i]] <- seu_temp
}

doublet.cell.list <- list()
doublet.df <- list()
for (i in names(combo.list)) {
  seu_temp <- combo.list[[i]]
  doublet.sub <- seu_temp@meta.data[ ,grepl("DF", colnames(seu_temp@meta.data)), drop = F]
  head(doublet.sub)
  doublet.df[[i]] <- doublet.sub
  doublet.cell.list[[i]] <- rownames(doublet.sub) 
  Karaayvaz <- AddMetaData(Karaayvaz, doublet.sub, paste0("Doublet.", i))
}


colnames(Karaayvaz@meta.data)
Karaayvaz$Doublet.Call <- apply(Karaayvaz@meta.data[,c(41:50)], 1, function(x) x[!is.na(x)][1])
Karaayvaz@meta.data <- Karaayvaz@meta.data[,-c(41:50)]

setwd(DoubletResult_Dir)
saveRDS(Karaayvaz, "Karaayvaz_withDoublet_62622.rds")


#Pal ===============================

setwd(PrelimObj_Dir)
Pal <- readRDS(file  = "Pal_PRIM_filtered_62622.rds")
Pal$orig.ident <- "Pal_Prim"
DefaultAssay(Pal) <- "RNA"

hsp.genes <- grep(pattern = "^HSP", x = rownames(Pal@assays[["RNA"]]), value = TRUE)
print(hsp.genes)
Pal <- PercentageFeatureSet(Pal, "^HSP", col.name = "percent.heatshock")


combo <- Pal
combo <- SCTransform(combo, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet", "percent.heatshock"), verbose = TRUE)
combo <- RunPCA(combo, npcs = 100, ndims.print = 1:10, nfeatures.print = 10)
ElbowPlot(combo, ndims = 100)
DimHeatmap(combo, dims = 10:25, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 25:30, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 30:40, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 40:50, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 50:60, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 60:70, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 70:80, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 80:90, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 90:100, cells = 500, balanced = TRUE)



combo <- FindNeighbors(combo, reduction = "pca", dims = 1:60)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
combo <- FindClusters(combo, resolution = resolution.range)
clustree(combo, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")



combo <- FindNeighbors(combo, reduction = "pca", dims = 1:60, nn.method = "rann") 
combo <- FindClusters(combo, resolution = 0.7)
combo <- RunUMAP(object = combo, reduction = "pca", dims = 1:60, seed.use = 123)

combo.list <- SplitObject(combo, split.by = "samples")

pkchoose.list <- list()
for (i in 1:length(combo.list)) {
  seu_temp <- combo.list[[i]]
  
  sweep.res.list <- paramSweep_v3(seu_temp, PCs = 
                                    seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  #sweep.stats.list[[i]] <- sweep.stats
  sweep.stats.pk <- find.pK(sweep.stats)
  pK=as.numeric(as.character(sweep.stats.pk$pK))
  BCmetric=sweep.stats.pk$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  nExp_poi <- 0.031*nrow(seu_temp@meta.data) ##3.1 doublet rate
  seu_temp <- doubletFinder_v3(seu_temp, PCs = 
                                 seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE, pN = 0.25, pK = as.numeric(pK_choose), nExp = nExp_poi, reuse.pANN = FALSE)
  combo.list[[i]] <- seu_temp
}

doublet.cell.list <- list()
doublet.df <- list()
for (i in names(combo.list)) {
  seu_temp <- combo.list[[i]]
  doublet.sub <- seu_temp@meta.data[ ,grepl("DF", colnames(seu_temp@meta.data)), drop = F]
  head(doublet.sub)
  doublet.df[[i]] <- doublet.sub
  doublet.cell.list[[i]] <- rownames(doublet.sub) 
  Pal <- AddMetaData(Pal, doublet.sub, paste0("Doublet.", i))
}

colnames(Pal@meta.data)
Pal$Doublet.Call <- apply(Pal@meta.data[,38:65], 1, function(x) x[!is.na(x)][1])
Pal@meta.data <- Pal@meta.data[,-c(38:65)]

setwd(DoubletResult_Dir)
saveRDS(Pal, "Pal_withDoublet_62622.rds")



#Qian ===============================

setwd(PrelimObj_Dir)
Qian <- readRDS(file = "Qian_filtered_62622.rds")
Qian$orig.ident <- "Qian"
Qian <- subset(x = Qian, subset = Patient == "BC_3", invert = TRUE)
Qian$samples <- Qian$Patient
DefaultAssay(Qian) <- "RNA"

hsp.genes <- grep(pattern = "^HSP", x = rownames(Qian@assays[["RNA"]]), value = TRUE)
print(hsp.genes)
Qian <- PercentageFeatureSet(Qian, "^HSP", col.name = "percent.heatshock")

combo <- Qian
combo <- SCTransform(combo, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet", "percent.heatshock"), verbose = TRUE)
combo <- RunPCA(combo, npcs = 100, ndims.print = 1:10, nfeatures.print = 10)
ElbowPlot(combo, ndims = 100)
DimHeatmap(combo, dims = 10:25, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 25:30, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 30:40, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 40:50, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 50:60, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 60:70, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 70:80, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 80:90, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 90:100, cells = 500, balanced = TRUE)



combo <- FindNeighbors(combo, reduction = "pca", dims = 1:60)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
combo <- FindClusters(combo, resolution = resolution.range)
clustree(combo, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")



combo <- FindNeighbors(combo, reduction = "pca", dims = 1:60, nn.method = "rann")
combo <- FindClusters(combo, resolution = 1)
combo <- RunUMAP(object = combo, reduction = "pca", dims = 1:60, seed.use = 123)

combo.list <- SplitObject(combo, split.by = "samples")

pkchoose.list <- list()
for (i in 1:length(combo.list)) {
  seu_temp <- combo.list[[i]]
  
  sweep.res.list <- paramSweep_v3(seu_temp, PCs = 
                                    seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  #sweep.stats.list[[i]] <- sweep.stats
  sweep.stats.pk <- find.pK(sweep.stats)
  pK=as.numeric(as.character(sweep.stats.pk$pK))
  BCmetric=sweep.stats.pk$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  nExp_poi <- 0.039*nrow(seu_temp@meta.data) ##3.9 doublet rate
  seu_temp <- doubletFinder_v3(seu_temp, PCs = 
                                 seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE, pN = 0.25, pK = as.numeric(pK_choose), nExp = nExp_poi, reuse.pANN = FALSE)
  combo.list[[i]] <- seu_temp
}

doublet.cell.list <- list()
doublet.df <- list()
for (i in names(combo.list)) {
  seu_temp <- combo.list[[i]]
  doublet.sub <- seu_temp@meta.data[ ,grepl("DF", colnames(seu_temp@meta.data)), drop = F]
  head(doublet.sub)
  doublet.df[[i]] <- doublet.sub
  doublet.cell.list[[i]] <- rownames(doublet.sub) 
  Qian <- AddMetaData(Qian, doublet.sub, paste0("Doublet.", i))
}

colnames(Qian@meta.data)
Qian$Doublet.Call <- apply(Qian@meta.data[,c(43:52)], 1, function(x) x[!is.na(x)][1])
Qian@meta.data <- Qian@meta.data[,-c(43:52)]

setwd(DoubletResult_Dir)
saveRDS(Qian, "Qian_withDouble_62622t.rds")



#Savas ====================

setwd(PrelimObj_Dir)
Savas <- readRDS(file = "Savas_filtered_62622.rds")

Savas$orig.ident <- "Savas" 
Savas$samples <- Savas$Patient
DefaultAssay(Savas) <- "RNA"

hsp.genes <- grep(pattern = "^HSP", x = rownames(Savas@assays[["RNA"]]), value = TRUE)
print(hsp.genes)
Savas <- PercentageFeatureSet(Savas, "^HSP", col.name = "percent.heatshock")


combo <- Savas
combo <- SCTransform(combo, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet", "percent.heatshock"), verbose = TRUE)
combo <- RunPCA(combo, npcs = 100, ndims.print = 1:10, nfeatures.print = 10)
DimHeatmap(combo, dims = 25:30, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 30:40, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 40:50, cells = 500, balanced = TRUE)

combo <- FindNeighbors(combo, reduction = "pca", dims = 1:50)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
combo <- FindClusters(combo, resolution = resolution.range)
clustree(combo, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")


combo <- FindNeighbors(combo, reduction = "pca", dims = 1:50, nn.method = "rann") 
combo <- FindClusters(combo, resolution = 1.6)
combo <- RunUMAP(object = combo, reduction = "pca", dims = 1:50, seed.use = 123)

combo.list <- SplitObject(combo, split.by = "samples")

pkchoose.list <- list()
for (i in 1:length(combo.list)) {
  seu_temp <- combo.list[[i]]
  
  sweep.res.list <- paramSweep_v3(seu_temp, PCs = 
                                    seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  #sweep.stats.list[[i]] <- sweep.stats
  sweep.stats.pk <- find.pK(sweep.stats)
  pK=as.numeric(as.character(sweep.stats.pk$pK))
  BCmetric=sweep.stats.pk$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  nExp_poi <- 0.031*nrow(seu_temp@meta.data) ##3.1% doublet rate
  seu_temp <- doubletFinder_v3(seu_temp, PCs = 
                                 seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE, pN = 0.25, pK = as.numeric(pK_choose), nExp = nExp_poi, reuse.pANN = FALSE)
  combo.list[[i]] <- seu_temp
}

doublet.cell.list <- list()
doublet.df <- list()
for (i in names(combo.list)) {
  seu_temp <- combo.list[[i]]
  doublet.sub <- seu_temp@meta.data[ ,grepl("DF", colnames(seu_temp@meta.data)), drop = F]
  head(doublet.sub)
  doublet.df[[i]] <- doublet.sub
  doublet.cell.list[[i]] <- rownames(doublet.sub) 
  Savas <- AddMetaData(Savas, doublet.sub, paste0("Doublet.", i))
}


colnames(Savas@meta.data)
Savas$Doublet.Call <- apply(Savas@meta.data[,37:38], 1, function(x) x[!is.na(x)][1])
Savas@meta.data <- Savas@meta.data[,-c(37:38)]
table(Savas$Doublet.Call)

setwd(DoubletResult_Dir)
saveRDS(Savas, "Savas_withDoublet_62622.rds")

#wu (2020) ===============================

setwd(PrelimObj_Dir)
Wu <- readRDS(file = "OldWu_filtered_62622.rds")
Wu$orig.ident <- "Wu"
Wu <- subset(x = Wu, subset = Patient == "Patient 4", invert = TRUE)
DefaultAssay(Wu) <- "RNA"

hsp.genes <- grep(pattern = "^HSP", x = rownames(Wu@assays[["RNA"]]), value = TRUE)
print(hsp.genes)
Wu <- PercentageFeatureSet(Wu, "^HSP", col.name = "percent.heatshock")


combo <- Wu
combo <- SCTransform(combo, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet", "percent.heatshock"), verbose = TRUE)
combo <- RunPCA(combo, npcs = 100, ndims.print = 1:10, nfeatures.print = 10)
ElbowPlot(combo, ndims = 100)
DimHeatmap(combo, dims = 10:25, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 25:30, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 30:40, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 40:50, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 50:60, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 60:70, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 70:80, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 80:90, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 90:100, cells = 500, balanced = TRUE)

# #https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/77
# mat <- Seurat::GetAssayData(combo, assay = "RNA", slot = "scale.data")
# pca <- combo[["pca"]]
# # Get the total variance:
# total_variance <- combo@reductions$pca@misc$total.variance
# eigValues = (pca@stdev)^2  ## EigenValues
# varExplained = eigValues / total_variance
# total = cumsum(varExplained)


combo <- FindNeighbors(combo, reduction = "pca", dims = 1:75)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
combo <- FindClusters(combo, resolution = resolution.range)
clustree(combo, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")



combo <- FindNeighbors(combo, reduction = "pca", dims = 1:75, nn.method = "rann") #dims = 1:50, nn.method = "rann") 
combo <- FindClusters(combo, resolution = 0.7)
combo <- RunUMAP(object = combo, reduction = "pca", dims = 1:75, seed.use = 123)#dims = 1:50,seed.use = 123)

combo.list <- SplitObject(combo, split.by = "samples")

pkchoose.list <- list()
for (i in 1:length(combo.list)) {
  seu_temp <- combo.list[[i]]
  
  sweep.res.list <- paramSweep_v3(seu_temp, PCs = 
                                    seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  #sweep.stats.list[[i]] <- sweep.stats
  sweep.stats.pk <- find.pK(sweep.stats)
  pK=as.numeric(as.character(sweep.stats.pk$pK))
  BCmetric=sweep.stats.pk$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  nExp_poi <- 0.023*nrow(seu_temp@meta.data) ## 2.3% doublet rate
  seu_temp <- doubletFinder_v3(seu_temp, PCs = 
                                 seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE, pN = 0.25, pK = as.numeric(pK_choose), nExp = nExp_poi, reuse.pANN = FALSE)
  combo.list[[i]] <- seu_temp
}

doublet.cell.list <- list()
doublet.df <- list()
for (i in names(combo.list)) {
  seu_temp <- combo.list[[i]]
  doublet.sub <- seu_temp@meta.data[ ,grepl("DF", colnames(seu_temp@meta.data)), drop = F]
  head(doublet.sub)
  doublet.df[[i]] <- doublet.sub
  doublet.cell.list[[i]] <- rownames(doublet.sub) 
  Wu <- AddMetaData(Wu, doublet.sub, paste0("Doublet.", i))
}

colnames(Wu@meta.data)
Wu$Doublet.Call <- apply(Wu@meta.data[,39:42], 1, function(x) x[!is.na(x)][1])
Wu@meta.data <- Wu@meta.data[,-c(39:42)]

setwd(DoubletResult_Dir)
saveRDS(Wu, "OldWu_withDoublet_62622.rds")


#wu (2021) ===============================

setwd(PrelimObj_Dir)
Wu2021 <- readRDS(file = "NewWu_PRIM_filtered_62622.rds")
Wu2021$orig.ident <- "Wu2021prim"
Wu2021 <- subset(x = Wu2021, subset = Patient == "CID3963", invert = TRUE)
Wu2021 <- subset(x = Wu2021, subset = Patient == "CID4066", invert = TRUE)
Wu2021 <- subset(x = Wu2021, subset = Patient == "CID4398", invert = TRUE)
Wu2021 <- subset(x = Wu2021, subset = Patient == "CID44991", invert = TRUE)
Wu2021 <- subset(x = Wu2021, subset = Patient == "CID4513", invert = TRUE)
Wu2021 <- subset(x = Wu2021, subset = Patient == "CID4523", invert = TRUE)
DefaultAssay(Wu2021) <- "RNA"

hsp.genes <- grep(pattern = "^HSP", x = rownames(Wu2021@assays[["RNA"]]), value = TRUE)
print(hsp.genes)
Wu2021 <- PercentageFeatureSet(Wu2021, "^HSP", col.name = "percent.heatshock")

combo <- Wu2021
combo <- SCTransform(combo, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet", "percent.heatshock"), verbose = TRUE)
combo <- RunPCA(combo, npcs = 100, ndims.print = 1:10, nfeatures.print = 10)
ElbowPlot(combo, ndims = 100)
DimHeatmap(combo, dims = 10:25, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 25:30, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 30:40, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 40:50, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 50:60, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 60:70, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 70:80, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 80:90, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 90:100, cells = 500, balanced = TRUE)


combo <- FindNeighbors(combo, reduction = "pca", dims = 1:60)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
combo <- FindClusters(combo, resolution = resolution.range)
clustree(combo, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")



combo <- FindNeighbors(combo, reduction = "pca", dims = 1:60, nn.method = "rann") 
combo <- FindClusters(combo, resolution = 0.7)
combo <- RunUMAP(object = combo, reduction = "pca", dims = 1:60, seed.use = 123)

combo.list <- SplitObject(combo, split.by = "samples")

pkchoose.list <- list()
for (i in 1:length(combo.list)) {
  seu_temp <- combo.list[[i]]
  
  sweep.res.list <- paramSweep_v3(seu_temp, PCs = 
                                    seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  #sweep.stats.list[[i]] <- sweep.stats
  sweep.stats.pk <- find.pK(sweep.stats)
  pK=as.numeric(as.character(sweep.stats.pk$pK))
  BCmetric=sweep.stats.pk$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  nExp_poi <- 0.023*nrow(seu_temp@meta.data) ##2.3 doublet rate
  seu_temp <- doubletFinder_v3(seu_temp, PCs = 
                                 seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE, pN = 0.25, pK = as.numeric(pK_choose), nExp = nExp_poi, reuse.pANN = FALSE)
  combo.list[[i]] <- seu_temp
}

doublet.cell.list <- list()
doublet.df <- list()
for (i in names(combo.list)) {
  seu_temp <- combo.list[[i]]
  doublet.sub <- seu_temp@meta.data[ ,grepl("DF", colnames(seu_temp@meta.data)), drop = F]
  head(doublet.sub)
  doublet.df[[i]] <- doublet.sub
  doublet.cell.list[[i]] <- rownames(doublet.sub) 
  Wu2021 <- AddMetaData(Wu2021, doublet.sub, paste0("Doublet.", i))
}

colnames(Wu2021@meta.data)
Wu2021$Doublet.Call <- apply(Wu2021@meta.data[,43:62], 1, function(x) x[!is.na(x)][1])
Wu2021@meta.data <- Wu2021@meta.data[,-c(43:62)]

setwd(DoubletResult_Dir)
saveRDS(Wu2021, "Wu2021_withDoublet_62622.rds")



#Xu ================================

setwd(PrelimObj_Dir)
Xu <- readRDS(file = "Xu_PRIM_filtered_62622.rds")
Xu$orig.ident <- "Xu"
DefaultAssay(Xu) <- "RNA"

hsp.genes <- grep(pattern = "^HSP", x = rownames(Xu@assays[["RNA"]]), value = TRUE)
print(hsp.genes)
Xu <- PercentageFeatureSet(Xu, "^HSP", col.name = "percent.heatshock")

combo <- Xu
combo <- SCTransform(combo, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet", "percent.heatshock"), verbose = TRUE)
combo <- RunPCA(combo, npcs = 100, ndims.print = 1:10, nfeatures.print = 10)
DimHeatmap(combo, dims = 25:30, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 30:40, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 40:50, cells = 500, balanced = TRUE)

combo <- FindNeighbors(combo, reduction = "pca", dims = 1:40)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
combo <- FindClusters(combo, resolution = resolution.range)
clustree(combo, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")


combo <- FindNeighbors(combo, reduction = "pca", dims = 1:40, nn.method = "rann") 
combo <- FindClusters(combo, resolution = 1.6)
combo <- RunUMAP(object = combo, reduction = "pca", dims = 1:40, seed.use = 123)

combo.list <- SplitObject(combo, split.by = "samples")

pkchoose.list <- list()
for (i in 1:length(combo.list)) {
  seu_temp <- combo.list[[i]]
  
  sweep.res.list <- paramSweep_v3(seu_temp, PCs = 
                                    seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  #sweep.stats.list[[i]] <- sweep.stats
  sweep.stats.pk <- find.pK(sweep.stats)
  pK=as.numeric(as.character(sweep.stats.pk$pK))
  BCmetric=sweep.stats.pk$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  nExp_poi <- 0.0126*nrow(seu_temp@meta.data) #1.26% doublet rate
  seu_temp <- doubletFinder_v3(seu_temp, PCs = 
                                 seu_temp@commands$RunUMAP.SCT.pca$dims, sct = TRUE, pN = 0.25, pK = as.numeric(pK_choose), nExp = nExp_poi, reuse.pANN = FALSE)
  combo.list[[i]] <- seu_temp
}

doublet.cell.list <- list()
doublet.df <- list()
for (i in names(combo.list)) {
  seu_temp <- combo.list[[i]]
  doublet.sub <- seu_temp@meta.data[ ,grepl("DF", colnames(seu_temp@meta.data)), drop = F]
  head(doublet.sub)
  doublet.df[[i]] <- doublet.sub
  doublet.cell.list[[i]] <- rownames(doublet.sub) 
  Xu <- AddMetaData(Xu, doublet.sub, paste0("Doublet.", i))
}

colnames(Xu@meta.data)
Xu$Doublet.Call <- apply(Xu@meta.data[,24:28], 1, function(x) x[!is.na(x)][1])
Xu@meta.data <- Xu@meta.data[,-c(24:28)]

setwd(DoubletResult_Dir)
saveRDS(Xu, "Xu_withDoublet_62622.rds")

