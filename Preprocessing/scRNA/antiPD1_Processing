

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This script loads and processes Bassez's antiPD1 dataset
# up until the point of labeling cell types,
# which is found in CellTypeLabeling.R

# Dataset Link: https://lambrechtslab.sites.vib.be/en/single-cell
# Paper Link: https://www.nature.com/articles/s41591-021-01323-8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


antiPD1_origDir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/PrimscRNA_Orig_Dataset_Files/antiPD1_dataset"
antiPD1_files <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/antiPD1_Processing/antiPD1_Intermobjects/"

# Libraries -------------------------------------------------------------

library(TFEA.ChIP)
library(org.Hs.eg.db)
library(DoubletFinder)
library(BiocManager)
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


#non-treated (cohort1) loading ========================================

setwd(antiPD1_origDir)
antiPD1 <- readRDS(file = "1863-counts_cells_cohort1.rds")

antiPD1.meta <- read.csv(file = '1872-BIOKEY_metaData_cohort1_web.csv')


head(antiPD1.meta)
colnames(antiPD1.meta) <- c("cells", "nCount_RNA", "nFeature_RNA", "Patient", 
                            "Treatment.Status", "Expansion", "BC.Subtype", "celltype_minor", "Treatment.Status")

antiPD1.meta <- antiPD1.meta[,-9]
#Please note that patient numbers have been randomized not to match those highlighted in Supplementary Table 1 from Bassez et al. (Nature Medicine 2021).
antiPD1.meta$Treatment.Status[which(antiPD1.meta$Treatment.Status == "On")] <- "Naive (On antiPD1)"
antiPD1.meta$Treatment.Status[which(antiPD1.meta$Treatment.Status == "Pre")] <- "Naive (Pre antiPD1)"
antiPD1.meta$BC.Subtype[which(antiPD1.meta$BC.Subtype == "ER+")] <- "HR+"

antiPD1.meta$BRCA.Status <- "NA"
antiPD1.meta$Ethnicity <- "NA"
antiPD1.meta$Grade <- "NA"
antiPD1.meta$Stage <- "NA"
antiPD1.meta$Status <- "Primary"
antiPD1.meta$Tissue.Source <- "Primary Tumor"
antiPD1.meta$BMI <- "NA"
antiPD1.meta$RNA.Type <-"mRNA"
antiPD1.meta$Gene.Coverage <- "5'"
antiPD1.meta$Library.Preparation <-"10X Genomics Chromium"
antiPD1.meta$Capture.Method <- "10x Genomics Chromium"
antiPD1.meta$Sequencer <- "Illumina HiSeq 4000"  
antiPD1.meta$Gender <- "Female"


#gene conversion ====================================

library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(antiPD1)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- antiPD1[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- antiPD1[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  #k <- data.frame(t(j))
  k <- Matrix(t(j), sparse = T)
  
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}

antiPD1 <- new_mat

test <- antiPD1[usegenes_final == "MCM2",]
test[1:7,1:5]
test_trans <- t(test)
test_trans <- as.data.frame(test_trans)
test_trans[1:5,1:6]
test_genes <- summary(test_trans$MCM2)



#create antiPD1 nochemo object ======================

fdat_naive <-as.data.frame(rownames(antiPD1))
common_colnames <- "gene_short_name"
colnames(fdat_naive) <- common_colnames

sobj_naive <- CreateSeuratObject(counts = antiPD1)
rownames(antiPD1.meta) <- antiPD1.meta$cells
sobj_naive <-AddMetaData(sobj_naive,metadata=antiPD1.meta)
sobj_naive[["RNA"]]@meta.features<-fdat_naive
head(sobj_naive[["RNA"]]@meta.features)
slotNames(sobj_naive[["RNA"]])

setwd(antiPD1_files)
saveRDS(sobj_naive, "unfiltered_Renamed_nochemo_7422.rds")

#treated(chemo) (cohort2) loading ========================================

setwd(antiPD1_origDir)
chemo <- readRDS(file = "1867-counts_cells_cohort2.rds")

chemo.meta <- read.csv(file = '1871-BIOKEY_metaData_cohort2_web (1).csv')
head(chemo.meta)

colnames(chemo.meta) <- c("cells", "nCount_RNA", "nFeature_RNA", "Patient", 
                          "Treatment.Status", "Expansion", "BC.Subtype", "celltype_minor", "Treatment.Status")

chemo.meta <- chemo.meta[,-9]
#Please note that patient numbers have been randomized not to match those highlighted in Supplementary Table 1 from Bassez et al. (Nature Medicine 2021).
chemo.meta$Treatment.Status[which(chemo.meta$Treatment.Status == "On")] <- "Neoadjuvant_chemo (On antiPD1)"
chemo.meta$Treatment.Status[which(chemo.meta$Treatment.Status == "Pre")] <- "Neoadjuvant_chemo (Pre antiPD1)"
chemo.meta$BC.Subtype[which(chemo.meta$BC.Subtype == "ER+")] <- "HR+"

chemo.meta$BRCA.Status <- "NA"
chemo.meta$Ethnicity <- "NA"
chemo.meta$Grade <- "NA"
chemo.meta$Stage <- "NA"
chemo.meta$Status <- "Primary"
chemo.meta$Tissue.Source <- "Primary Tumor"
chemo.meta$BMI <- "NA"
chemo.meta$RNA.Type <-"mRNA"
chemo.meta$Gene.Coverage <- "5'"
chemo.meta$Library.Preparation <-"10X Genomics Chromium"
chemo.meta$Capture.Method <- "10x Genomics Chromium"
chemo.meta$Sequencer <- "Illumina HiSeq 4000"  
chemo.meta$Gender <- "Female"




#gene conversion ====================================

library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(chemo)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- chemo[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- chemo[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  #k <- data.frame(t(j))
  k <- Matrix(t(j), sparse = T)
  
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}

chemo <- new_mat

test <- chemo[usegenes_final == "MCM2",]
test[1:7,1:5]
test_trans <- t(test)
test_trans <- as.data.frame(test_trans)
test_trans[1:5,1:6]
test_genes <- summary(test_trans$MCM2)




# chemo object creation ===========================

fdat_chemo <-as.data.frame(rownames(chemo))
common_colnames <- "gene_short_name"
colnames(fdat_chemo) <- common_colnames

sobj_chemo <- CreateSeuratObject(counts = chemo)
rownames(chemo.meta) <- chemo.meta$cells
sobj_chemo <-AddMetaData(sobj_chemo,metadata=chemo.meta)
sobj_chemo[["RNA"]]@meta.features<-fdat_chemo
head(sobj_chemo[["RNA"]]@meta.features)
slotNames(sobj_chemo[["RNA"]])

setwd(antiPD1_files)
saveRDS(sobj_chemo, "unfiltered_Renamed_YESchemo_7422.rds")


#cell cycle ========================================================================
#from: https://satijalab.org/seurat/articles/cell_cycle_vignette.html

sobj1 <- sobj_naive
sobj2 <- sobj_chemo

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sobj1 <- CellCycleScoring(sobj1, s.features = s.genes, g2m.features = g2m.genes)
sobj1$CC.Difference <- sobj1$S.Score - sobj1$G2M.Score
# view cell cycle scores and phase assignments
head(sobj1[[]])

sobj2 <- CellCycleScoring(sobj2, s.features = s.genes, g2m.features = g2m.genes)
sobj2$CC.Difference <- sobj2$S.Score - sobj2$G2M.Score
# view cell cycle scores and phase assignments
head(sobj2[[]])

# Preprocessing and Filtering (SCTransform) --------------------------------------------------------------------

#if print(mito.genes) returns "character(0)", try "^MT.\\"
mito.genes <- grep(pattern = "^MT-", x = rownames(sobj1@assays[["RNA"]]), value = TRUE)
print(mito.genes)
sobj1 <- PercentageFeatureSet(sobj1, pattern = "^MT-", col.name = "percent.mt")

#blood cell contamination
hemo.genes <- grep(pattern = "^HB[^(P)]", x = rownames(sobj1@assays[["RNA"]]), value = TRUE)
print(hemo.genes)
sobj1 <- PercentageFeatureSet(sobj1, "^HB[^(P)]", col.name = "percent.hb")
platelet.genes <- grep(pattern = "^PECAM1|^PF4\\b", x = rownames(sobj1@assays[["RNA"]]), value = TRUE)
print(platelet.genes)
sobj1 <- PercentageFeatureSet(sobj1, "^PECAM1|^PF4\\b", col.name = "percent.platelet")

#heat shock protein
hsp.genes <- grep(pattern = "^HSP", x = rownames(sobj1@assays[["RNA"]]), value = TRUE)
print(hsp.genes)
sobj1 <- PercentageFeatureSet(sobj1, "^HSP", col.name = "percent.heatshock")



nCounthi <- quantile(sobj1@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
nCountlo <- quantile(sobj1@meta.data$nCount_RNA, 0.05)#05)  
nFeathi <- quantile(sobj1@meta.data$nFeature_RNA, 0.95) #0.9) 
nFeatlo <- quantile(sobj1@meta.data$nFeature_RNA, 0.05) #0.1) 
ptMtlo <- quantile(sobj1@meta.data$percent.mt, 0.05)
ptMthi <- quantile(sobj1@meta.data$percent.mt, 0.99)
ptHbhi <- quantile(sobj1@meta.data$percent.hb, 0.95)

sobj1 <- subset(x = sobj1, subset = nFeature_RNA > nFeatlo & nFeature_RNA < nFeathi)
sobj1 <- subset(x = sobj1, subset = percent.mt > ptMtlo & percent.mt < ptMthi)
sobj1 <- subset(x = sobj1, subset = nCount_RNA > nCountlo & nCount_RNA < nCounthi)
sobj1 <- subset(x = sobj1, subset = percent.hb < ptHbhi)

summary(sobj1@meta.data$percent.mt) #max MSUT be between 10 and 15 at most
summary(sobj1@meta.data$percent.hb) #as low as possible
dim(sobj1)
dim(sobj_naive)



mito.genes <- grep(pattern = "^MT-", x = rownames(sobj2@assays[["RNA"]]), value = TRUE)
print(mito.genes)
sobj2 <- PercentageFeatureSet(sobj2, pattern = "^MT-", col.name = "percent.mt")

#blood cell contamination
hemo.genes <- grep(pattern = "^HB[^(P)]", x = rownames(sobj2@assays[["RNA"]]), value = TRUE)
print(hemo.genes)
sobj2 <- PercentageFeatureSet(sobj2, "^HB[^(P)]", col.name = "percent.hb")
platelet.genes <- grep(pattern = "^PECAM1|^PF4\\b", x = rownames(sobj2@assays[["RNA"]]), value = TRUE)
print(platelet.genes)
sobj2 <- PercentageFeatureSet(sobj2, "^PECAM1|^PF4\\b", col.name = "percent.platelet")

#heatshock
hsp.genes <- grep(pattern = "^HSP", x = rownames(sobj2@assays[["RNA"]]), value = TRUE)
print(hsp.genes)
sobj2 <- PercentageFeatureSet(sobj2, "^HSP", col.name = "percent.heatshock")



nCounthi <- quantile(sobj2@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
nCountlo <- quantile(sobj2@meta.data$nCount_RNA, 0.05)#05)  
nFeathi <- quantile(sobj2@meta.data$nFeature_RNA, 0.95) #0.9) 
nFeatlo <- quantile(sobj2@meta.data$nFeature_RNA, 0.05) #0.1) 
ptMtlo <- quantile(sobj2@meta.data$percent.mt, 0.05)
ptMthi <- quantile(sobj2@meta.data$percent.mt, 0.99)
ptHbhi <- quantile(sobj2@meta.data$percent.hb, 0.95)

sobj2 <- subset(x = sobj2, subset = nFeature_RNA > nFeatlo & nFeature_RNA < nFeathi)
sobj2 <- subset(x = sobj2, subset = percent.mt > ptMtlo & percent.mt < ptMthi)
sobj2 <- subset(x = sobj2, subset = nCount_RNA > nCountlo & nCount_RNA < nCounthi)
sobj2 <- subset(x = sobj2, subset = percent.hb < ptHbhi)

summary(sobj2@meta.data$percent.mt) #max MSUT be between 10 and 15 at most
summary(sobj2@meta.data$percent.hb) #as low as possible
dim(sobj2)
dim(sobj_chemo)



#doublet finder (nochemo) =======================================

sobj1$samples <- sobj1$Patient
combo <- sobj1
combo <- SCTransform(combo, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet", "percent.heatshock"), verbose = TRUE)
combo <- RunPCA(combo, npcs = 100, ndims.print = 1:10, nfeatures.print = 10)
DimHeatmap(combo, dims = 10:25, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 25:30, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 30:40, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 40:50, cells = 500, balanced = TRUE)

combo <- FindNeighbors(combo, reduction = "pca", dims = 1:45)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
combo <- FindClusters(combo, resolution = resolution.range)
clustree(combo, prefix = "SCT_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")


combo <- FindNeighbors(combo, reduction = "pca", dims = 1:45, nn.method = "rann") 
combo <- FindClusters(combo, resolution = 0.4)
combo <- RunUMAP(object = combo, reduction = "pca", dims = 1:45, seed.use = 123)


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
  nExp_poi <- 0.023*nrow(seu_temp@meta.data) #2.3% doublet rate
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
  sobj1 <- AddMetaData(sobj1, doublet.sub, paste0("Doublet.", i))
}



colnames(sobj1@meta.data)
sobj1$Doublet.Call <- apply(sobj1@meta.data[,32:62], 1, function(x) x[!is.na(x)][1])
sobj1@meta.data <- sobj1@meta.data[,-c(32:62)]

setwd(antiPD1_files)
saveRDS(sobj1, "nochemo_withDoublet_7422.rds")

#doublet finder (yeschemo) =======================================


sobj2$samples <- sobj2$Patient
combo <- sobj2
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
combo <- FindClusters(combo, resolution = 0.7)
combo <- RunUMAP(object = combo, reduction = "pca", dims = 1:30, seed.use = 123)


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
  nExp_poi <- 0.023*nrow(seu_temp@meta.data) #2.3% doublet rate
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
  sobj2 <- AddMetaData(sobj2, doublet.sub, paste0("Doublet.", i))
}



colnames(sobj2@meta.data)
sobj2$Doublet.Call <- apply(sobj2@meta.data[,c(32:42)], 1, function(x) x[!is.na(x)][1])
sobj2@meta.data <- sobj2@meta.data[,-c(32:42)]
table(sobj2$Doublet.Call)

setwd(antiPD1_files)
saveRDS(sobj2, "yeschemo_withDoublet_7422.rds")


#PCA and UMAP -----------------------------------------

sobj1 <- subset(sobj1, subset = Doublet.Call == "Doublet", invert = T)
sobj2 <- subset(sobj2, subset = Doublet.Call == "Doublet", invert = T)

sobj1 <- SCTransform(sobj1, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet", "percent.heatshock"), verbose = TRUE)

sobj2 <- SCTransform(sobj2, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet", "percent.heatshock"), verbose = TRUE)


sobj1 <- RunPCA(sobj1, npcs = 100, ndims.print = 1:10, nfeatures.print = 5)
DimHeatmap(sobj1, dims = 10:25, cells = 500, balanced = TRUE)
DimHeatmap(sobj1, dims = 25:30, cells = 500, balanced = TRUE)
DimHeatmap(sobj1, dims = 30:40, cells = 500, balanced = TRUE)
DimHeatmap(sobj1, dims = 40:50, cells = 500, balanced = TRUE)

DefaultAssay(sobj1) <- "SCT"

sobj1 <- FindNeighbors(sobj1, reduction = "pca", dims = 1:35, nn.method = "rann") 
sobj1 <- FindClusters(sobj1, resolution = 5)
sobj1 <- RunUMAP(object = sobj1, reduction = "pca", dims = 1:35, seed.use = 123)

DimPlot(sobj1, reduction = "umap", label = FALSE, raster = FALSE, 
        order = c("NK Cells", "Epithelial Cells", "T Cells", "Stroma", "B Cells",
                  "Endothelial Cells", "Plasma Cells", "Myeloid Cells")) 
p <- DimPlot(sobj1, reduction = "umap", group.by = "Treatment.Status", raster = FALSE) + ggtitle(label = " ")
p <- DimPlot(sobj1, reduction = "umap", group.by = "samples", raster = FALSE) +                                          # Apply guides function
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1)) + ggtitle(label = " ")
p <- DimPlot(sobj1, reduction = "umap", group.by = "Expansion", raster = FALSE) + ggtitle(label = " ")
sobj1$samples <- sobj1$Patient

DimPlot(sobj1, reduction = "umap", label = T, raster = FALSE) + ggtitle(label = " ")
DimPlot(sobj1, reduction = "umap", group.by = "celltype_minor", raster = FALSE) + ggtitle(label = " ")

#chemo ______________________________________

sobj2 <- RunPCA(sobj2, npcs = 100, ndims.print = 1:10, nfeatures.print = 5)
DimHeatmap(sobj2, dims = 10:25, cells = 500, balanced = TRUE)
DimHeatmap(sobj2, dims = 25:30, cells = 500, balanced = TRUE)
DimHeatmap(sobj2, dims = 30:40, cells = 500, balanced = TRUE)
DimHeatmap(sobj2, dims = 40:50, cells = 500, balanced = TRUE)

DefaultAssay(sobj2) <- "SCT"

sobj2 <- FindNeighbors(sobj2, reduction = "pca", dims = 1:40, nn.method = "rann") 
sobj2 <- FindClusters(sobj2, resolution = 5)
sobj2 <- RunUMAP(object = sobj2, reduction = "pca", dims = 1:40, seed.use = 123)

DimPlot(sobj2, reduction = "umap", label = T, raster = FALSE) + ggtitle(label = " ")
DimPlot(sobj2, reduction = "umap", group.by = "celltype_minor", raster = FALSE) + ggtitle(label = " ")

DimPlot(sobj2, reduction = "umap", pt.size = 0.1,label=TRUE, raster = FALSE) + ggtitle(label = "antiPD1") + theme(legend.position = "None")
DimPlot(sobj2, reduction = "umap", group.by = "Treatment.Status", raster = FALSE) + ggtitle(label = " ")
DimPlot(sobj2, reduction = "umap", group.by = "celltype_minor", raster = FALSE) + ggtitle(label = " ")

sobj2$samples <- sobj2$Patient


