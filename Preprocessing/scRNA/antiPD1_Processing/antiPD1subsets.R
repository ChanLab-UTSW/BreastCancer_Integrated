setwd("/work/InternalMedicine/s204665")
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/antiPD1")

library(TFEA.ChIP)
library(org.Hs.eg.db)
# Libraries -------------------------------------------------------------
#install.packages("future.batchtools")
library(BiocManager)
library(GEOquery) 
#install.packages("plyr")
library(plyr)#, lib.loc="/opt/R/4.0.2/lib64/R/library")
library(dplyr) 
library(Matrix)
library(devtools)
#install.packages("Seurat")
library(Seurat)#, lib.loc="/opt/R/4.0.2/lib64/R/library") 
#install.packages("ggplot2")
library(ggplot2)#, lib.loc="/opt/R/4.0.2/lib64/R/library") 
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

#If you need to find what package a function is in:
#install.packages("sos")
#library("sos")
#findFn("Heatmap")



library(devtools)
#install_github("immunogenomics/harmony")
library(harmony)
library(SeuratData)
library(UCell)

#install.packages("sctransform")
#remotes::install_github("ChristophH/sctransform@develop")
#devtools::install_github("const-ae/sparseMatrixStats")
#BiocManager::install("glmGamPoi")
library(glmGamPoi)
library(sctransform)#, lib.loc="/opt/R/4.0.2/lib64/R/library")
library(matrixStats)#, lib.loc="/opt/R/4.0.2/lib64/R/library")
library(sparseMatrixStats)
library(DESeq2)
library(genefu)

#Load in full objects =============================



sobj1 <- readRDS(file = "antiPD1nochemo_withpam50_nozallgenedem_71322.rds")
sobj2 <- readRDS(file = "antiPD1YESchemo_withpam50_nozallgenedem_71322")

#subset antiPD1 groups if need ==========================================

nochemo_pre <- subset(x = sobj1, subset = Treatment.Status == "Naive (Pre antiPD1)")
nochemo_post <- subset(x = sobj1, subset = Treatment.Status == "Naive (On antiPD1)")

chemo_pre <- subset(x = sobj2, subset = Treatment.Status == "Neoadjuvant_chemo (Pre antiPD1)")
chemo_post <- subset(x = sobj2, subset = Treatment.Status == "Neoadjuvant_chemo (On antiPD1)")


#subset (no chemo, prePD1) ==============================================

DefaultAssay(nochemo_pre) <- "RNA"
nochemo_pre <- NormalizeData(nochemo_pre, assay = "RNA")
table(Idents(nochemo_pre))
nochemo_pre$celltype_final <- Idents(nochemo_pre)


nochemo_pre <- SCTransform(nochemo_pre, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                            "percent.platelet", "percent.heatshock"), verbose = TRUE)
nochemo_pre <- RunPCA(nochemo_pre, npcs = 100, verbose = TRUE)
pdf("test.pdf")
DimHeatmap(nochemo_pre, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(nochemo_pre, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(nochemo_pre, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(nochemo_pre, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(nochemo_pre, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()

nochemo_pre <- FindNeighbors(nochemo_pre, reduction = "pca", dims = 1:30, nn.method = "rann")


resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6)
nochemo_pre <- FindClusters(nochemo_pre, resolution = resolution.range)
clustree(nochemo_pre, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")

DefaultAssay(nochemo_pre) <- "SCT"
nochemo_pre <- FindNeighbors(nochemo_pre, reduction = "pca", dims = 1:30, nn.method = "rann")
nochemo_pre <- FindClusters(nochemo_pre, resolution = 1)
nochemo_pre <- RunUMAP(nochemo_pre, reduction = "pca", dims = 1:30, verbose = TRUE, seed.use = 123)

Idents(nochemo_pre) <- nochemo_pre$celltype_final
pdf("test.pdf", width = 10, height = 10)
DimPlot(nochemo_pre, reduction = "umap", label = T, repel = TRUE, raster = FALSE) #+ ggtitle(label = "NK Cells from lumA samples (PC = 15)")DimPlot(nochemo_pre, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "BC.Subtype")
dev.off()

DefaultAssay(nochemo_pre) <- "RNA"

nochemo_pre <- NormalizeData(nochemo_pre, assay = "RNA")

saveRDS(nochemo_pre, file = "nochemo_PREtx_noZallgenedem_71322.rds")
nochemo_pre <- readRDS(file = "nochemo_PREtx_noZallgenedem_71322.rds")
# saveRDS(nochemo_pre, file = "nochemo_PREtx_FINAL_7522.rds")
# nochemo_pre <- readRDS(file = "nochemo_PREtx_FINAL_7522.rds")

#subset (no chemo, postPD1) ==============================================

DefaultAssay(nochemo_post) <- "RNA"
nochemo_post <- NormalizeData(nochemo_post, assay = "RNA")
table(Idents(nochemo_post))
nochemo_post$celltype_final <- Idents(nochemo_post)


nochemo_post <- SCTransform(nochemo_post, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                            "percent.platelet", "percent.heatshock"), verbose = TRUE)
nochemo_post <- RunPCA(nochemo_post, npcs = 100, verbose = TRUE)
pdf("test.pdf")
DimHeatmap(nochemo_post, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(nochemo_post, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(nochemo_post, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(nochemo_post, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(nochemo_post, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()

nochemo_post <- FindNeighbors(nochemo_post, reduction = "pca", dims = 1:30, nn.method = "rann")


resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6)
nochemo_post <- FindClusters(nochemo_post, resolution = resolution.range)
clustree(nochemo_post, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")

DefaultAssay(nochemo_post) <- "SCT"

nochemo_post <- FindNeighbors(nochemo_post, reduction = "pca", dims = 1:40, nn.method = "rann")
nochemo_post <- FindClusters(nochemo_post, resolution = 1)
nochemo_post <- RunUMAP(nochemo_post, reduction = "pca", dims = 1:40, verbose = TRUE, seed.use = 123)

Idents(nochemo_post) <- nochemo_post$celltype_final
pdf("test.pdf", width = 10, height = 10)
DimPlot(nochemo_post, reduction = "umap", label = T, repel = TRUE, raster = FALSE) #+ ggtitle(label = "NK Cells from lumA samples (PC = 15)")DimPlot(nochemo_pre, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "BC.Subtype")
dev.off()

DefaultAssay(nochemo_post) <- "RNA"

nochemo_post <- NormalizeData(nochemo_post, assay = "RNA")


saveRDS(nochemo_post, file = "nochemo_POSTtx_noZallgenedem_71322.rds")
nochemo_post <- readRDS(file = "nochemo_POSTtx_noZallgenedem_71322.rds")
# saveRDS(nochemo_post, file = "nochemo_POSTtx_FINAL_7522.rds")
# nochemo_post <- readRDS(file = "nochemo_POSTtx_SCTnormRNA_61022.rds")

#subset (YES chemo, prePD1) ==============================================

DefaultAssay(chemo_pre) <- "RNA"
table(Idents(chemo_pre))
chemo_pre$celltype_final <- Idents(chemo_pre)


chemo_pre <- SCTransform(chemo_pre, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                            "percent.platelet", "percent.heatshock"), verbose = TRUE)
chemo_pre <- RunPCA(chemo_pre, npcs = 100, verbose = TRUE)
pdf("test.pdf")
DimHeatmap(chemo_pre, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(chemo_pre, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(chemo_pre, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(chemo_pre, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(chemo_pre, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()

chemo_pre <- FindNeighbors(chemo_pre, reduction = "pca", dims = 1:25, nn.method = "rann")


resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6)
chemo_pre <- FindClusters(chemo_pre, resolution = resolution.range)
clustree(chemo_pre, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")

DefaultAssay(chemo_pre) <- "SCT"
chemo_pre <- FindNeighbors(chemo_pre, reduction = "pca", dims = 1:35, nn.method = "rann")
chemo_pre <- FindClusters(chemo_pre, resolution = 1)
chemo_pre <- RunUMAP(chemo_pre, reduction = "pca", dims = 1:35, verbose = TRUE, seed.use = 123)

Idents(chemo_pre) <- chemo_pre$celltype_final
pdf("test.pdf", width = 10, height = 10)
DimPlot(chemo_pre, reduction = "umap", label = T, repel = TRUE, raster = FALSE) #+ ggtitle(label = "NK Cells from lumA samples (PC = 15)")DimPlot(nochemo_pre, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "BC.Subtype")
dev.off()

DefaultAssay(chemo_pre) <- "RNA"

chemo_pre <- NormalizeData(chemo_pre, assay = "RNA")


saveRDS(chemo_pre, file = "YESchemo_PREtx_noZallgenedem_71322.rds")
chemo_pre <- readRDS(file = "YESchemo_PREtx_noZallgenedem_71322.rds")
# saveRDS(chemo_pre, file = "YESchemo_PREtx_FINAL_7522.rds")
# chemo_pre <- readRDS(file = "YESchemo_PREtx_FINAL_7522.rds")

#subset (YES chemo, postPD1) ==============================================

DefaultAssay(chemo_post) <- "RNA"
table(Idents(chemo_post))
chemo_post$celltype_final <- Idents(chemo_post)


chemo_post <- SCTransform(chemo_post, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                              "percent.platelet", "percent.heatshock"), verbose = TRUE)
chemo_post <- RunPCA(chemo_post, npcs = 100, verbose = TRUE)
pdf("test.pdf")
DimHeatmap(chemo_post, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(chemo_post, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(chemo_post, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(chemo_post, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(chemo_post, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()

chemo_post <- FindNeighbors(chemo_post, reduction = "pca", dims = 1:40, nn.method = "rann")


resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6)
chemo_post <- FindClusters(chemo_post, resolution = resolution.range)
clustree(chemo_post, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")

DefaultAssay(chemo_post) <- "SCT"
chemo_post <- FindNeighbors(chemo_post, reduction = "pca", dims = 1:40, nn.method = "rann")
chemo_post <- FindClusters(chemo_post, resolution = 1)
chemo_post <- RunUMAP(chemo_post, reduction = "pca", dims = 1:40, verbose = TRUE, seed.use = 123)

Idents(chemo_post) <- chemo_post$celltype_final
pdf("test.pdf", width = 10, height = 10)
DimPlot(chemo_post, reduction = "umap", label = T, repel = TRUE, raster = FALSE) #+ ggtitle(label = "NK Cells from lumA samples (PC = 15)")DimPlot(nochemo_pre, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "BC.Subtype")
dev.off()

DefaultAssay(chemo_post) <- "RNA"

chemo_post <- NormalizeData(chemo_post, assay = "RNA")


saveRDS(chemo_post, file = "YESchemo_POSTtx_noZallgenedem_71322.rds")
chemo_post <- readRDS(file = "YESchemo_POSTtx_noZallgenedem_71322.rds")
# saveRDS(chemo_post, file = "YESchemo_POSTtx_FINAL_7522.rds")
# chemo_post <- readRDS(file = "YESchemo_POSTtx_FINAL_7522.rds")

#NK (nochemo all) ==========================

NK_nochemo <- subset(sobj1, idents = "NK Cells")
DefaultAssay(NK_nochemo) <- "RNA"
table(Idents(NK_nochemo))
NK_nochemo$celltype_final <- Idents(NK_nochemo)


NK_nochemo <- SCTransform(NK_nochemo, vars.to.regress = c("percent.mt"), verbose = TRUE)
NK_nochemo <- RunPCA(NK_nochemo, npcs = 100, verbose = TRUE)
DimHeatmap(NK_nochemo, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(NK_nochemo, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(NK_nochemo, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(NK_nochemo, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(NK_nochemo, dims = 45:55, cells = 500, balanced = TRUE)

NK_nochemo <- FindNeighbors(NK_nochemo, reduction = "pca", dims = 1:25, nn.method = "rann")


resolution.range <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
NK_nochemo <- FindClusters(NK_nochemo, resolution = resolution.range)
clustree(NK_nochemo, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")

NK_nochemo <- FindNeighbors(NK_nochemo, reduction = "pca", dims = 1:25, nn.method = "rann")
NK_nochemo <- FindClusters(NK_nochemo, resolution = 2)
NK_nochemo <- RunUMAP(NK_nochemo, reduction = "pca", dims = 1:25, verbose = TRUE, seed.use = 123)

#Idents(NK_nochemo) <- NK_nochemo$celltype_final

DimPlot(NK_nochemo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE)


NK_nochemo <- RenameIdents(NK_nochemo, `1` = "1", `0` = "1",
                           `13` = "0", `2` = "0", `10` = "0", `12` = "0", 
                           `5` = "0", `6` = "0", `7` = "0", `8` = "0",
                           `3` = "2", `4` = "2", `9` = "2", `11` = "2")

saveRDS(NK_nochemo, file = "NK_nochemoboth_withNKsub.rds")
NK_nochemo <- readRDS(file = "NK_nochemoboth_withNKsub.rds")

DefaultAssay(NK_nochemo) <- "RNA"

NKdegs <- list(c("NR4A3", "CCL4", "RASGRP2",
                 "JUN", "DUSP1", "FOS", "NR4A1", "KLRG1", 
                 "NHSL2", "CX3CR1", "FOSB", "NR4A2", #upregulated reprog
                 "COX6A2-", "PYCR1-", "EXTL1-", "CHAC1-", "SLC6A9-",
                 "OSGIN1-", "OSBPL1A-", "PPP2R2C-", "CLBA1-",
                 "HMOX1-", "NQO1-", "CARS1-", "SSTR2-", #downregulated reprog
                 "GNLY", #NK specific
                 "PRF1",
                 "KLRD1",
                 "NKG7",
                 "HSPB1-",
                 "COX6C-",
                 "ADIRF-",
                 "CST3-",
                 "EPCAM-", #additions
                 "CD33-",
                 "CD3D-",
                 "MS4A1-",
                 #cancer up compared to healthy
                 "XCL1",
                 "MUCL1",
                 "CCL4L2",
                 "CCL3",
                 "CXCL13-",
                 "LAIR2-",
                 "CD8A-",
                 "MARCKSL1-"))
NK_nochemo <- AddModuleScore_UCell(NK_nochemo, features = NKdegs, name = "NKdegs", assay = "RNA")
FeaturePlot(object = NK_nochemo, features = "signature_1NKdegs", order = TRUE, label = FALSE, repel = TRUE, min.cutoff = 0, raster = FALSE) + 
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_colour_gradientn(colours = c("lightgrey",
                                     "#70bfed",
                                     "#414174")) + 
  ggtitle(label = " ")



NK_nochemo[["ClusterName"]] <- Idents(NK_nochemo)
NK_nochemo$cells <- colnames(NK_nochemo)
sobjlists <- FetchData(object = NK_nochemo, vars = c("signature_1NKdegs",
                                                     "cells", "BC.Subtype", "ClusterName",
                                                     "Expansion", "Treatment.Status"))
saveRDS(sobjlists, file = "anti_nochemoboth_tolaptop.rds")

#NK (nochemo, pre) ===========

NK_pre <- subset(x = NK_nochemo, subset = Treatment.Status == "Naive (Pre antiPD1)")
table(Idents(NK_pre))
NK_pre$celltype_NKsub <- Idents(NK_pre)

DefaultAssay(NK_pre) <- "RNA"
NK_pre <- SCTransform(NK_pre, vars.to.regress = c("percent.mt"), verbose = TRUE)
NK_pre <- RunPCA(NK_pre, npcs = 100, verbose = TRUE)
DimHeatmap(NK_pre, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(NK_pre, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(NK_pre, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(NK_pre, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(NK_pre, dims = 45:55, cells = 500, balanced = TRUE)

NK_pre <- FindNeighbors(NK_pre, reduction = "pca", dims = 1:25, nn.method = "rann")


resolution.range <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
NK_pre <- FindClusters(NK_pre, resolution = resolution.range)
clustree(NK_pre, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")

NK_pre <- FindNeighbors(NK_pre, reduction = "pca", dims = 1:25, nn.method = "rann")
NK_pre <- FindClusters(NK_pre, resolution = 2)
NK_pre <- RunUMAP(NK_pre, reduction = "pca", dims = 1:25, verbose = TRUE, seed.use = 123)

Idents(NK_pre) <- NK_pre$celltype_NKsub

DimPlot(NK_pre, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE)

saveRDS(NK_pre, file = "NKsubset_nochemoPRE_withNKsub.rds")

#NK (nochemo, post) ==========
NK_post <- subset(x = NK_nochemo, subset = Treatment.Status == "Naive (On antiPD1)")
DefaultAssay(NK_post) <- "RNA"
table(Idents(NK_post))
NK_post$celltype_NKsub <- Idents(NK_post)


DefaultAssay(NK_post) <- "RNA"
NK_post <- SCTransform(NK_post, vars.to.regress = c("percent.mt"), verbose = TRUE)
NK_post <- RunPCA(NK_post, npcs = 100, verbose = TRUE)
DimHeatmap(NK_post, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(NK_post, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(NK_post, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(NK_post, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(NK_post, dims = 45:55, cells = 500, balanced = TRUE)

NK_post <- FindNeighbors(NK_post, reduction = "pca", dims = 1:25, nn.method = "rann")


resolution.range <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
NK_post <- FindClusters(NK_post, resolution = resolution.range)
clustree(NK_post, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")

NK_post <- FindNeighbors(NK_post, reduction = "pca", dims = 1:25, nn.method = "rann")
NK_post <- FindClusters(NK_post, resolution = 2)
NK_post <- RunUMAP(NK_post, reduction = "pca", dims = 1:25, verbose = TRUE, seed.use = 123)

Idents(NK_post) <- NK_post$celltype_NKsub

DimPlot(NK_post, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE)


saveRDS(NK_post, file = "NKsubset_nochemoPOST_withNKsub.rds")



#Epi ALL GROUPS ===========================================

PD1_all <- merge(sobj1, sobj2)
PD1.all.Epi <- subset(PD1_all, idents = "Epithelial Cells")
DefaultAssay(PD1.all.Epi) <- "RNA"
table(Idents(PD1.all.Epi))
PD1.all.Epi$celltype_final <- Idents(PD1.all.Epi)


PD1.all.Epi <- SCTransform(PD1.all.Epi, vars.to.regress = c("percent.mt"), verbose = TRUE)
PD1.all.Epi <- RunPCA(PD1.all.Epi, npcs = 100, verbose = TRUE)


DimHeatmap(PD1.all.Epi, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(PD1.all.Epi, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(PD1.all.Epi, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(PD1.all.Epi, dims = 35:45, cells = 500, balanced = TRUE)
pdf("test.pdf")
DimHeatmap(PD1.all.Epi, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()

PD1.all.Epi <- FindNeighbors(PD1.all.Epi, reduction = "pca", dims = 1:50, nn.method = "rann")

resolution.range <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
PD1.all.Epi <- FindClusters(PD1.all.Epi, resolution = resolution.range)
p <- clustree(PD1.all.Epi, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")
ggsave("PD1all_Epiclustree.pdf", p, width = 14.22, height = 13.01)

PD1.all.Epi <- FindNeighbors(PD1.all.Epi, reduction = "pca", dims = 1:50, nn.method = "rann")
PD1.all.Epi <- FindClusters(PD1.all.Epi, resolution = 0.4)
PD1.all.Epi <- RunUMAP(PD1.all.Epi, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use = 123)

Idents(PD1.all.Epi) <- PD1.all.Epi$celltype_final
pdf("test.pdf")
DimPlot(PD1.all.Epi, reduction = "umap", label = TRUE, repel = TRUE,raster = FALSE)
dev.off()

#Subtyping (Epi everything) ===============

DefaultAssay(PD1.all.Epi) <- "RNA"

#Mydata <- subset(combo.reference, idents = c("Unspecified Epithelial Cells"))
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects")
sigdat <- read.csv("NatGen_Supplementary_table_S4.csv")
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/antiPD1")
temp_allgenes <- c(as.vector(sigdat[,"Basal_SC"]),
                   as.vector(sigdat[,"Her2E_SC"]),
                   as.vector(sigdat[,"LumA_SC"]),
                   as.vector(sigdat[,"LumB_SC"]))
temp_allgenes <- unique(temp_allgenes[!temp_allgenes == ""])

#Read in the single cell RDS object as 'Mydata'
# Mydata <- readRDS("path_to_seurat_object.Rdata")
Mydata <- ScaleData(PD1.all.Epi, features=temp_allgenes, assay = "RNA")

tocalc<-as.data.frame(Mydata@assays$RNA@scale.data)

#calculate mean scsubtype scores
outdat <- matrix(0,
                 nrow=ncol(sigdat),
                 ncol=ncol(tocalc),
                 dimnames=list(colnames(sigdat),
                               colnames(tocalc)))
for(i in 1:ncol(sigdat)){
  
  # sigdat[i,!is.na(sigdat[i,])]->module
  # row <- as.character(unlist(module))
  row <- as.character(sigdat[,i])
  row<-unique(row[row != ""])
  genes<-which(rownames(tocalc) %in% row)
  
  temp<-apply(tocalc[genes,],2,function(x){mean(as.numeric(x),na.rm=TRUE)})
  
  outdat[i,]<-as.numeric(temp)
  
}

final<-outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
final<-as.data.frame(final)
is.num <- sapply(final, is.numeric);final[is.num] <- lapply(final[is.num], round, 4)
finalm<-as.matrix(final)

##Scaling scores function before calling the highest Call
center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
  get_average <- function(v) sum(v * row.w)/sum(row.w)
  average <- apply(x, 2, get_average)
  sweep(x, 2, average)
}

##Obtaining the highest call
finalmt<-as.data.frame(t(finalm))
finalm.sweep.t<-center_sweep(finalmt)
Finalnames<-colnames(finalm.sweep.t)[max.col(finalm.sweep.t,ties.method="first")]
finalm.sweep.t$SCSubtypeCall <- Finalnames

##Writing out output files (rownames remain the same for both)

write.csv(finalm.sweep.t, "PD1Epi_everything_Scores.csv")

#write.table(finalm.sweep.t, "Mydata_Scores.txt", sep="\t")
#write.table(Finalnames, "Mydata_CALLS.txt", sep="\t")

Mydata <- AddMetaData(Mydata, finalm.sweep.t[,5], col.name = "sc50.Pred")
PD1.all.Epi <- AddMetaData(PD1.all.Epi, finalm.sweep.t[,5], col.name = "sc50.Pred")

cancer.epi <- Mydata
cancer.epi$samples <- paste(cancer.epi$samples, cancer.epi$BC.Subtype, sep = "_")
sobjlists <- FetchData(object = cancer.epi, vars = c("samples","Patient", "BC.Subtype", "sc50.Pred"))
library(reshape2)

sobjlists = melt(sobjlists, id.vars = c("samples","Patient", "BC.Subtype", "sc50.Pred"))
sobjlists$cells <- rownames(sobjlists)
head(sobjlists)

counts_call <- table(sobjlists$samples, sobjlists$sc50.Pred)
write.csv(counts_call, file = "PD1Epi_everything_call_count_persample.csv")


#labels =========

PD1.all.Epi$pam50.pred <- PD1.all.Epi$samples
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_1")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_10")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_11")] <- "Her2"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_12")] <- "LumB"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_13")] <- "Her2"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_14")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_15")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_16")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_17")] <- "LumB"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_18")] <- "LumA"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_19")] <- "Her2"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_2")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_20")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_21")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_22")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_23")] <- "LumA"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_24")] <- "LumA"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_25")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_26")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_27")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_28")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_29")] <- "LumA"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_3")] <- "Her2"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_30")] <- "LumB"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_31")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_32")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_33")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_34")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_35")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_36")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_37")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_38")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_39")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_4")] <- "LumB"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_40")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_41")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_42")] <- "LumB"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_5")] <- "LumB"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_6")] <- "LumA"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_7")] <- "Her2"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_8")] <- "Basal"
PD1.all.Epi$pam50.pred[which(PD1.all.Epi$pam50.pred == "BIOKEY_9")] <- "Basal"
table(PD1.all.Epi$pam50.pred)

saveRDS(PD1.all.Epi, file = "Epi_noyeschemo_prepost.rds")
prePD1 <- readRDS(file = "Epi_noyeschemo_prepost.rds")




#Epi no chemo pre post ==============
nochemoboth <- subset(PD1.all.Epi, (subset = Treatment.Status == "Naive (Pre antiPD1)") | (subset = Treatment.Status == "Naive (On antiPD1)"))

DefaultAssay(nochemoboth) <- "RNA"
table(Idents(nochemoboth))
nochemoboth$celltype_final <- Idents(nochemoboth)

nochemoboth <- SCTransform(nochemoboth, vars.to.regress = c("percent.mt"), verbose = TRUE)
nochemoboth <- RunPCA(nochemoboth, npcs = 100, verbose = TRUE)
pdf("test.pdf")
DimHeatmap(nochemoboth, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(nochemoboth, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(nochemoboth, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(nochemoboth, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(nochemoboth, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()

nochemoboth <- FindNeighbors(nochemoboth, reduction = "pca", dims = 1:50, nn.method = "rann")


resolution.range <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
nochemoboth <- FindClusters(nochemoboth, resolution = resolution.range)
pdf("test.pdf")
clustree(nochemoboth, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")
dev.off()

nochemoboth <- FindNeighbors(nochemoboth, reduction = "pca", dims = 1:50, nn.method = "rann")
nochemoboth <- FindClusters(nochemoboth, resolution = 0.7)
nochemoboth <- RunUMAP(nochemoboth, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use = 123)

Idents(nochemoboth) <- nochemoboth$celltype_final

pdf("test.pdf")
DimPlot(nochemoboth, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE)
dev.off()

saveRDS(nochemoboth, file = "Epi_nochemoboth_withsc50pam50.rds")

#Epi yes chemo pre post =============
yeschemoboth <- subset(PD1.all.Epi, (subset = Treatment.Status == "Neoadjuvant_chemo (Pre antiPD1)") | (subset = Treatment.Status == "Neoadjuvant_chemo (On antiPD1)"))

DefaultAssay(yeschemoboth) <- "RNA"
table(Idents(yeschemoboth))
yeschemoboth$celltype_final <- Idents(yeschemoboth)

yeschemoboth <- SCTransform(yeschemoboth, vars.to.regress = c("percent.mt"), verbose = TRUE)
yeschemoboth <- RunPCA(yeschemoboth, npcs = 100, verbose = TRUE)
pdf("test.pdf")
DimHeatmap(yeschemoboth, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(yeschemoboth, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(yeschemoboth, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(yeschemoboth, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(yeschemoboth, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()

yeschemoboth <- FindNeighbors(yeschemoboth, reduction = "pca", dims = 1:40, nn.method = "rann")


resolution.range <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
yeschemoboth <- FindClusters(yeschemoboth, resolution = resolution.range)
pdf("test.pdf")
clustree(yeschemoboth, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")
dev.off()

yeschemoboth <- FindNeighbors(yeschemoboth, reduction = "pca", dims = 1:40, nn.method = "rann")
yeschemoboth <- FindClusters(yeschemoboth, resolution = 0.7)
yeschemoboth <- RunUMAP(yeschemoboth, reduction = "pca", dims = 1:40, verbose = TRUE, seed.use = 123)

Idents(yeschemoboth) <- yeschemoboth$celltype_final

pdf("test.pdf")
DimPlot(yeschemoboth, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE)
dev.off()

saveRDS(yeschemoboth, file = "Epi_yeschemoboth_withsc50pam50.rds")

#Epi prePD1 both ============================

prePD1 <- subset(PD1.all.Epi, (subset = Treatment.Status == "Naive (Pre antiPD1)") | (subset = Treatment.Status == "Neoadjuvant_chemo (Pre antiPD1)"))

DefaultAssay(prePD1) <- "RNA"
table(Idents(prePD1))
prePD1$celltype_final <- Idents(prePD1)

prePD1 <- SCTransform(prePD1, vars.to.regress = c("percent.mt"), verbose = TRUE)
prePD1 <- RunPCA(prePD1, npcs = 100, verbose = TRUE)
DimHeatmap(prePD1, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(prePD1, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(prePD1, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(prePD1, dims = 35:45, cells = 500, balanced = TRUE)
pdf("test.pdf")
DimHeatmap(prePD1, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()

prePD1 <- FindNeighbors(prePD1, reduction = "pca", dims = 1:50, nn.method = "rann")


resolution.range <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
prePD1 <- FindClusters(prePD1, resolution = resolution.range)
pdf("test.pdf")
clustree(prePD1, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")
dev.off()

prePD1 <- FindNeighbors(prePD1, reduction = "pca", dims = 1:50, nn.method = "rann")
prePD1 <- FindClusters(prePD1, resolution = 0.7)
prePD1 <- RunUMAP(prePD1, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use = 123)

Idents(prePD1) <- prePD1$celltype_final

pdf("test.pdf")
DimPlot(prePD1, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE)
dev.off()

saveRDS(prePD1, file = "Epi_prePD1_both.rds")
#Epi postPD1 both ============================

postPD1 <- subset(PD1.all.Epi, (subset = Treatment.Status == "Naive (On antiPD1)") | (subset = Treatment.Status == "Neoadjuvant_chemo (On antiPD1)"))

DefaultAssay(postPD1) <- "RNA"
table(Idents(postPD1))
postPD1$celltype_final <- Idents(postPD1)

postPD1 <- SCTransform(postPD1, vars.to.regress = c("percent.mt"), verbose = TRUE)
postPD1 <- RunPCA(postPD1, npcs = 100, verbose = TRUE)
DimHeatmap(postPD1, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(postPD1, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(postPD1, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(postPD1, dims = 35:45, cells = 500, balanced = TRUE)
pdf("test.pdf")
DimHeatmap(postPD1, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()

postPD1 <- FindNeighbors(postPD1, reduction = "pca", dims = 1:50, nn.method = "rann")


resolution.range <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
postPD1 <- FindClusters(postPD1, resolution = resolution.range)
pdf("test.pdf")
clustree(postPD1, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")
dev.off()

postPD1 <- FindNeighbors(postPD1, reduction = "pca", dims = 1:50, nn.method = "rann")
postPD1 <- FindClusters(postPD1, resolution = 0.7)
postPD1 <- RunUMAP(postPD1, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use = 123)

Idents(postPD1) <- postPD1$celltype_final

pdf("test.pdf")
DimPlot(postPD1, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE)
dev.off()

saveRDS(postPD1, file = "Epi_postPD1_both.rds")

#Epi (nochemo, pre) ===========

Epi_pre <- subset(x = prePD1, subset = Treatment.Status == "Naive (Pre antiPD1)")
table(Idents(Epi_pre))
table(Epi_pre$Treatment.Status)

DefaultAssay(Epi_pre) <- "RNA"
Epi_pre <- SCTransform(Epi_pre, vars.to.regress = c("percent.mt"), verbose = TRUE)
Epi_pre <- RunPCA(Epi_pre, npcs = 100, verbose = TRUE)
DimHeatmap(Epi_pre, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(Epi_pre, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(Epi_pre, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(Epi_pre, dims = 35:45, cells = 500, balanced = TRUE)
pdf("test.pdf")
DimHeatmap(Epi_pre, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()
Epi_pre <- FindNeighbors(Epi_pre, reduction = "pca", dims = 1:50, nn.method = "rann")


resolution.range <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
Epi_pre <- FindClusters(Epi_pre, resolution = resolution.range)
pdf("test.pdf")
clustree(Epi_pre, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")
dev.off()

Epi_pre <- FindNeighbors(Epi_pre, reduction = "pca", dims = 1:50, nn.method = "rann")
Epi_pre <- FindClusters(Epi_pre, resolution = 0.7)
Epi_pre <- RunUMAP(Epi_pre, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use = 123)


pdf("test.pdf")
DimPlot(Epi_pre, reduction = "umap", label = F, repel = TRUE, raster = FALSE) + ggtitle ("Epithelial Clusters (no chemo, pre antiPD1)")
dev.off()
p1 <- DimPlot(Epi_pre, reduction = "umap", label = F, repel = TRUE, raster = FALSE, group.by = "BC.Subtype") 
ggsave("Epi_prenochemo.pdf", p1, width = 7.11, height = 6.5)

saveRDS(Epi_pre, file = "Epi_nochemoPRE_finalwithsubtypes.rds")
Epi_pre <- readRDS(file = "Epi_nochemoPRE_finalwithsubtypes.rds")


#Epi (nochemo, post) ==========
Epi_post <- subset(x = postPD1, subset = Treatment.Status == "Naive (On antiPD1)")
DefaultAssay(Epi_post) <- "RNA"
table(Idents(Epi_post))


DefaultAssay(Epi_post) <- "RNA"
Epi_post <- SCTransform(Epi_post, vars.to.regress = c("percent.mt"), verbose = TRUE)
Epi_post <- RunPCA(Epi_post, npcs = 100, verbose = TRUE)
pdf("test.pdf")
DimHeatmap(Epi_post, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(Epi_post, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(Epi_post, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(Epi_post, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(Epi_post, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()

Epi_post <- FindNeighbors(Epi_post, reduction = "pca", dims = 1:50, nn.method = "rann")


resolution.range <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
Epi_post <- FindClusters(Epi_post, resolution = resolution.range)
pdf("test.pdf")
clustree(Epi_post, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")
dev.off()

Epi_post <- FindNeighbors(Epi_post, reduction = "pca", dims = 1:50, nn.method = "rann")
Epi_post <- FindClusters(Epi_post, resolution = 1)
Epi_post <- RunUMAP(Epi_post, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use = 123)


pdf("test.pdf")
DimPlot(Epi_post, reduction = "umap", label = F, repel = TRUE, raster = FALSE) + ggtitle ("Epithelial Clusters (no chemo, post antiPD1)")
dev.off()
p1 <- DimPlot(Epi_post, reduction = "umap", label = F, repel = TRUE, raster = FALSE, group.by = "BC.Subtype") 

saveRDS(Epi_post, file = "Epi_nochemoPOST_finalwithsubtypes.rds")

Epi_post <- readRDS(file = "Epi_nochemoPOST_finalwithsubtypes.rds")

#Epi (yeschemo, pre) ===========

Epi_pre <- subset(x = prePD1, subset = Treatment.Status == "Neoadjuvant_chemo (Pre antiPD1)")
table(Idents(Epi_pre))
table(Epi_pre$Treatment.Status)

DefaultAssay(Epi_pre) <- "RNA"
Epi_pre <- SCTransform(Epi_pre, vars.to.regress = c("percent.mt"), verbose = TRUE)
Epi_pre <- RunPCA(Epi_pre, npcs = 100, verbose = TRUE)

pdf("test.pdf")
DimHeatmap(Epi_pre, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(Epi_pre, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(Epi_pre, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(Epi_pre, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(Epi_pre, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()

Epi_pre <- FindNeighbors(Epi_pre, reduction = "pca", dims = 1:40, nn.method = "rann")

resolution.range <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
Epi_pre <- FindClusters(Epi_pre, resolution = resolution.range)
p <- clustree(Epi_pre, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")
ggsave("Epi_prechemo_clustree.pdf", p, width = 14.22, height = 13.01)

Epi_pre <- FindNeighbors(Epi_pre, reduction = "pca", dims = 1:50, nn.method = "rann")
Epi_pre <- FindClusters(Epi_pre, resolution = 0.4)
Epi_pre <- RunUMAP(Epi_pre, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use = 123)

pdf("test.pdf")
DimPlot(Epi_pre, reduction = "umap", label = F, repel = TRUE, raster = FALSE) + ggtitle ("Epithelial Clusters (no chemo, pre antiPD1)")
dev.off()
p1 <- DimPlot(Epi_pre, reduction = "umap", label = F, repel = TRUE, raster = FALSE, group.by = "BC.Subtype") 

saveRDS(Epi_pre, file = "Epi_yeschemoPRE_finalwithsubtypes")

#Epi (yeschemo, post) ===========

Epi_post <- subset(x = postPD1, subset = Treatment.Status == "Neoadjuvant_chemo (On antiPD1)")
table(Idents(Epi_post))
table(Epi_post$Treatment.Status)

DefaultAssay(Epi_post) <- "RNA"
Epi_post <- SCTransform(Epi_post, vars.to.regress = c("percent.mt"), verbose = TRUE)
Epi_post <- RunPCA(Epi_post, npcs = 100, verbose = TRUE)
pdf("test.pdf")
DimHeatmap(Epi_post, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(Epi_post, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(Epi_post, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(Epi_post, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(Epi_post, dims = 45:55, cells = 500, balanced = TRUE)
dev.off()

Epi_post <- FindNeighbors(Epi_post, reduction = "pca", dims = 1:40, nn.method = "rann")

resolution.range <- c(0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
Epi_post <- FindClusters(Epi_post, resolution = resolution.range)
p <- clustree(Epi_post, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")
ggsave("Epi_postchemo_clustree.pdf", p, width = 14.22, height = 13.01)

Epi_post <- FindNeighbors(Epi_post, reduction = "pca", dims = 1:40, nn.method = "rann")
Epi_post <- FindClusters(Epi_post, resolution = 0.4)
Epi_post <- RunUMAP(Epi_post, reduction = "pca", dims = 1:40, verbose = TRUE, seed.use = 123)

pdf("test.pdf")
DimPlot(Epi_post, reduction = "umap", label = F, repel = TRUE, raster = FALSE) + ggtitle ("Epithelial Clusters (no chemo, pre antiPD1)")
dev.off()
p1 <- DimPlot(Epi_post, reduction = "umap", label = F, repel = TRUE, raster = FALSE, group.by = "BC.Subtype") 

saveRDS(Epi_post, "Epi_yeschemoPOST_finalwithsubtypes.rds")
