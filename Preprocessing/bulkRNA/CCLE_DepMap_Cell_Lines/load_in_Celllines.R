# SETUP ------------------------------------------

# Libraries --------------------------------------

setwd("~/Preprocessing/bulkRNA/CCLE_DepMap_Cell_Lines")

library(BiocManager)
library(GEOquery) 
library(plyr)
library(dplyr) 
library(Matrix)
library(devtools)
library(ggplot2)
library(multtest)
library(tibble)
library(UCell)
library(glmGamPoi)
library(sctransform)
library(matrixStats)
library(sparseMatrixStats)
library(SeuratObject, lib.loc = "/cm/shared/apps/R/gcc/4.1.1/lib64/R/library")
library(Seurat, lib.loc = "/cm/shared/apps/R/gcc/4.1.1/lib64/R/library")

# ------------------------------------------------
# ------------------------------------------------
# Load in counts matrix --------------------------

# from DepMap public 22Q2 files: RNAseq TPM gene expression data for all genes using RSEM; log2 transformed, using a pseudo-count of 1 
fullmat <- read.csv("CCLE_expression_full.csv")

rownames(fullmat) <- fullmat$X
fullmat <- fullmat[,-1]
colnames(fullmat) <- sub('\\.\\..*', '', colnames(fullmat))
fullmat <- t(fullmat)
fullmat[1:5,1:5]

# Load in metadata -------------------------------

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Cell_Lines/")

# from DepMap public 22Q2 files: metadata for all cell lines
sample.info <- read.csv("sample_info.csv")
sample.info <- sample.info[,-1]
rownames(sample.info) <- sample.info$DepMap_ID
head(sample.info)
cell.lines.wanted <- colnames(fullmat)
sample.info <- sample.info[rownames(sample.info) %in% cell.lines.wanted, ]
dim(sample.info)

# Convert gene names -----------------------------

library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(fullmatTPM)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- fullmatTPM[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- fullmatTPM[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}

dim(new_mat)
dim(fullmatTPM)

fullmatTPM <- new_mat

# Create Seurat object of DepMap cell lines ---------------------------

cell_lines <- CreateSeuratObject(fullmat, assay = "RNA", meta.data = sample.info)
DefaultAssay(cell_lines) <- "RNA"
cell_lines <- NormalizeData(cell_lines)
