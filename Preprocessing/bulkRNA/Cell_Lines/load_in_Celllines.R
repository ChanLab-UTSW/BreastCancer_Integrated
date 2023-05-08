# SETUP ------------------------------------------

# Libraries --------------------------------------

library(BiocManager)
library(GEOquery) 
library(plyr)
library(dplyr) 
library(Matrix)
library(devtools)
library(Seurat)
library(ggplot2)
library(multtest)
library(tibble)
library(SeuratData)
library(UCell)
library(glmGamPoi)
library(sctransform)
library(matrixStats)
library(sparseMatrixStats)

# ------------------------------------------------
# ------------------------------------------------
# Load in counts matrix --------------------------

fullmat <- read.csv("CCLE_expression_full.csv")
rownames(fullmat) <- fullmat$X
fullmat <- fullmat[,-1]
colnames(fullmat) <- sub('\\.\\..*', '', colnames(fullmat))
fullmat <- t(fullmat)
fullmat[1:5,1:5]

# Load in metadata -------------------------------

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Cell_Lines/")

sample.info <- read.csv("sample_info.csv")
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

# Create Seurat object ---------------------------

cell_lines <- CreateSeuratObject(fullmat, assay = "RNA", meta.data = sample.info)
DefaultAssay(cell_lines) <- "RNA"
cell_lines <- NormalizeData(cell_lines)
