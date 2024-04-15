# 
# DESCRIPTION 
#
# Loads and processes pembrolizumab-treated scRNA-seq dataset from Bassez et al.
# up to the point of cell type annotation (Cell_Type_Annotation.R)
#
# Bassez et al.: https://www.nature.com/articles/s41591-021-01323-8
# Data accessed: https://lambrechtslab.sites.vib.be/en/single-cell
#
# ------------------------------------------------
# ------------------------------------------------
# SETUP ------------------------------------------

antiPD1_origDir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/PrimscRNA_Orig_Dataset_Files/antiPD1_dataset"
antiPD1_files <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/antiPD1_Processing/antiPD1_Intermobjects/"

# Libraries --------------------------------------

library(plyr)
library(dplyr) 
library(Matrix)
library(Seurat)
library(ggplot2)
library(cowplot) 
library(SeuratData)
library(UCell)
library(sctransform)
library(matrixStats)
library(sparseMatrixStats)
library(limma)
library(org.Hs.eg.db)

# ------------------------------------------------
# ------------------------------------------------
# Load in raw data non-chemo-treated samples (cohort 1) ----

setwd(antiPD1_origDir)

# Load in counts matrix
antiPD1 <- readRDS(file = "1863-counts_cells_cohort1.rds")

# Load in metadata
antiPD1.meta <- read.csv(file = '1872-BIOKEY_metaData_cohort1_web.csv')
# Note: Patient IDs were randomized to unmatch those in Table S1 from Bassez et al. 

# Harmonize metadata
colnames(antiPD1.meta) <- c("cells", 
                            "nCount_RNA", 
                            "nFeature_RNA", 
                            "Patient", 
                            "Treatment.Status", 
                            "Expansion", 
                            "BC.Subtype", 
                            "celltype_minor", 
                            "Treatment.Status")

antiPD1.meta <- antiPD1.meta[,-9]
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

# Convert gene names -----------------------------

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
  k <- Matrix(t(j), sparse = T)
  
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}

antiPD1 <- new_mat

# Create Seurat object of non-chemo-treated samples ----

fdat_naive <-as.data.frame(rownames(antiPD1))
common_colnames <- "gene_short_name"
colnames(fdat_naive) <- common_colnames

sobj_naive <- CreateSeuratObject(counts = antiPD1)
rownames(antiPD1.meta) <- antiPD1.meta$cells
sobj_naive <-AddMetaData(sobj_naive,metadata=antiPD1.meta)
sobj_naive[["RNA"]]@meta.features<-fdat_naive
head(sobj_naive[["RNA"]]@meta.features)
slotNames(sobj_naive[["RNA"]])

# Cell cycling QC --------------------------------

# Adapted from: https://satijalab.org/seurat/articles/cell_cycle_vignette.html

sobj1 <- sobj_naive
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sobj1 <- CellCycleScoring(sobj1, s.features = s.genes, g2m.features = g2m.genes)
sobj1$CC.Difference <- sobj1$S.Score - sobj1$G2M.Score

# View cell cycle scores and phase assignments
head(sobj1[[]])

# Preprocessing and filtering --------------------

# Filter out mitochondrial genes
# If print(mito.genes) returns "character(0)", try "^MT.\\"
mito.genes <- grep(pattern = "^MT-", x = rownames(sobj1@assays[["RNA"]]), value = TRUE)
print(mito.genes)
sobj1 <- PercentageFeatureSet(sobj1, pattern = "^MT-", col.name = "percent.mt")

# Filter out blood cell contamination
hemo.genes <- grep(pattern = "^HB[^(P)]", x = rownames(sobj1@assays[["RNA"]]), value = TRUE)
print(hemo.genes)
sobj1 <- PercentageFeatureSet(sobj1, "^HB[^(P)]", col.name = "percent.hb")
platelet.genes <- grep(pattern = "^PECAM1|^PF4\\b", x = rownames(sobj1@assays[["RNA"]]), value = TRUE)
print(platelet.genes)
sobj1 <- PercentageFeatureSet(sobj1, "^PECAM1|^PF4\\b", col.name = "percent.platelet")

# Filter out heat shock protein
hsp.genes <- grep(pattern = "^HSP", x = rownames(sobj1@assays[["RNA"]]), value = TRUE)
print(hsp.genes)
sobj1 <- PercentageFeatureSet(sobj1, "^HSP", col.name = "percent.heatshock")

# Set thresholds for filtering 
nCounthi <- quantile(sobj1@meta.data$nCount_RNA, 0.95)
nCountlo <- quantile(sobj1@meta.data$nCount_RNA, 0.05)  
nFeathi <- quantile(sobj1@meta.data$nFeature_RNA, 0.95)  
nFeatlo <- quantile(sobj1@meta.data$nFeature_RNA, 0.05) 
ptMtlo <- quantile(sobj1@meta.data$percent.mt, 0.05)
ptMthi <- quantile(sobj1@meta.data$percent.mt, 0.99)
ptHbhi <- quantile(sobj1@meta.data$percent.hb, 0.95)

# Filter out cells that fall outside of threshold
sobj1 <- subset(x = sobj1, subset = nFeature_RNA > nFeatlo & nFeature_RNA < nFeathi)
sobj1 <- subset(x = sobj1, subset = percent.mt > ptMtlo & percent.mt < ptMthi)
sobj1 <- subset(x = sobj1, subset = nCount_RNA > nCountlo & nCount_RNA < nCounthi)
sobj1 <- subset(x = sobj1, subset = percent.hb < ptHbhi)

# Check filtering
summary(sobj1@meta.data$percent.mt) # Max should be between 10-15 
summary(sobj1@meta.data$percent.hb) # Value should be as low as possible
dim(sobj1)
dim(sobj_naive)

# Remove doublets using DoubletFinder ------------

sobj1$samples <- sobj1$Patient
combo <- sobj1
combo <- SCTransform(combo, 
                     vars.to.regress = c("percent.mt", 
                                         "nCount_RNA",
                                         "nFeature_RNA", 
                                         "percent.hb", 
                                         "percent.platelet", 
                                         "percent.heatshock"), 
                     verbose = TRUE)

combo <- RunPCA(combo, 
                npcs = 100, 
                ndims.print = 1:10, 
                nfeatures.print = 10)

DimHeatmap(combo, dims = 1:50, cells = 500, balanced = TRUE)

combo <- FindNeighbors(combo, 
                       reduction = "pca", 
                       dims = 1:45)

resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)

combo <- FindClusters(combo, 
                      resolution = resolution.range)

combo <- FindNeighbors(combo, 
                       reduction = "pca", 
                       dims = 1:45, 
                       nn.method = "rann") 

combo <- FindClusters(combo, resolution = 0.4)

combo <- RunUMAP(object = combo, 
                 reduction = "pca", 
                 dims = 1:45, 
                 seed.use = 123)

combo.list <- SplitObject(combo, split.by = "samples")

pkchoose.list <- list()
for (i in 1:length(combo.list)) {
  seu_temp <- combo.list[[i]]
  
  sweep.res.list <- paramSweep_v3(seu_temp, 
                                  PCs = seu_temp@commands$RunUMAP.SCT.pca$dims, 
                                  sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats.pk <- find.pK(sweep.stats)
  
  pK = as.numeric(as.character(sweep.stats.pk$pK))
  BCmetric = sweep.stats.pk$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  nExp_poi = 0.023 * nrow(seu_temp@meta.data) # 2.3% doublet rate
  
  seu_temp <- doubletFinder_v3(seu_temp, 
                               PCs = seu_temp@commands$RunUMAP.SCT.pca$dims, 
                               sct = TRUE, 
                               pN = 0.25, 
                               pK = as.numeric(pK_choose), 
                               nExp = nExp_poi, 
                               reuse.pANN = FALSE)
  combo.list[[i]] <- seu_temp
}

doublet.cell.list <- list()
doublet.df <- list()
for (i in names(combo.list)) {
  seu_temp <- combo.list[[i]]
  doublet.sub <- seu_temp@meta.data[ , grepl("DF", 
                                             colnames(seu_temp@meta.data)), 
                                     drop = F]
  head(doublet.sub)
  doublet.df[[i]] <- doublet.sub
  doublet.cell.list[[i]] <- rownames(doublet.sub) 
  sobj1 <- AddMetaData(sobj1, doublet.sub, paste0("Doublet.", i))
}

colnames(sobj1@meta.data)
sobj1$Doublet.Call <- apply(sobj1@meta.data[,32:62], 1, function(x) x[!is.na(x)][1])
sobj1@meta.data <- sobj1@meta.data[,-c(32:62)]

# Perform PCA, UMAP, unsupervised clustering -----

sobj1 <- subset(sobj1, subset = Doublet.Call == "Doublet", invert = T)

sobj1 <- SCTransform(sobj1, 
                     vars.to.regress = c("percent.mt", 
                                         "nCount_RNA",
                                         "nFeature_RNA", 
                                         "percent.hb", 
                                         "percent.platelet", 
                                         "percent.heatshock"), 
                     verbose = TRUE)

sobj1 <- RunPCA(sobj1, 
                npcs = 100, 
                ndims.print = 1:10, 
                nfeatures.print = 5)
DimHeatmap(sobj1, dims = 1:50, cells = 500, balanced = TRUE)

DefaultAssay(sobj1) <- "SCT"
sobj1 <- FindNeighbors(sobj1, 
                       reduction = "pca", 
                       dims = 1:35, 
                       nn.method = "rann") 

sobj1 <- FindClusters(sobj1, resolution = 5)

sobj1 <- RunUMAP(object = sobj1, 
                 reduction = "pca", 
                 dims = 1:35, 
                 seed.use = 123)

DimPlot(sobj1, 
        reduction = "umap", 
        label = FALSE, 
        raster = FALSE, 
        order = c("NK Cells", "Epithelial Cells", "T Cells", 
                  "Stroma", "B Cells", "Endothelial Cells", 
                  "Plasma Cells", "Myeloid Cells")) 

DimPlot(sobj1, 
        reduction = "umap", 
        group.by = "Treatment.Status", 
        raster = FALSE) + 
  ggtitle(label = " ")

DimPlot(sobj1, 
        reduction = "umap", 
        group.by = "samples", 
        raster = FALSE) + 
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1)) + 
  ggtitle(label = " ")

DimPlot(sobj1, 
        reduction = "umap", 
        group.by = "Expansion", 
        raster = FALSE) + 
  ggtitle(label = " ")
sobj1$samples <- sobj1$Patient
