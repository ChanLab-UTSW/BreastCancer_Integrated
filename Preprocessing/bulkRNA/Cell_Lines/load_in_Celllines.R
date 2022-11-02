


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



#Load in matrix ================

#TPM good but log2 not ideal: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3247-x#:~:text=Based%20on%20these%20results%2C%20TPM,association%20between%20sample%20and%20gene

fullmatTPM <- read.csv("CCLE_expression_full.csv")
rownames(fullmatTPM) <- fullmatTPM$X
fullmatTPM <- fullmatTPM[,-1]
colnames(fullmatTPM) <- sub('\\..*', '', colnames(fullmatTPM))
fullmatTPM <- t(fullmatTPM)

#convert genes ==============

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


#metadata ===========


sample.info <- read.csv("sample_info.csv")
rownames(sample.info) <- sample.info$DepMap_ID
head(sample.info)
dim(sample.info)
cell.lines.wanted <- colnames(fullmatTPM)
sample.info.subset <- sample.info[rownames(sample.info) %in% cell.lines.wanted, ]
dim(sample.info.subset)

cell.lines.missing <- cell.lines.wanted[!(cell.lines.wanted %in% rownames(sample.info.subset)) ]
dim(sample.info.subset)
colnames(sample.info.subset)

cell.lines.missing.df <- data.frame(row.names = cell.lines.missing)
cell.lines.missing.df$DepMap_ID <- cell.lines.missing
cell.lines.missing.df$cell_line_name <- NA
cell.lines.missing.df$stripped_cell_line_name <- NA
cell.lines.missing.df$CCLE_Name <- NA
cell.lines.missing.df$alias <- NA
cell.lines.missing.df$COSMICID <- NA
cell.lines.missing.df$sex <- NA
cell.lines.missing.df$source <- NA
cell.lines.missing.df$RRID <- NA
cell.lines.missing.df$WTSI_Master_Cell_ID <- NA
cell.lines.missing.df$sample_collection_site <- NA
cell.lines.missing.df$primary_or_metastasis <- NA
cell.lines.missing.df$primary_disease <- NA
cell.lines.missing.df$Subtype <- NA
cell.lines.missing.df$age <- NA
cell.lines.missing.df$Sanger_Model_ID <- NA
cell.lines.missing.df$depmap_public_comments <- NA
cell.lines.missing.df$lineage <- NA
cell.lines.missing.df$lineage_subtype <- NA
cell.lines.missing.df$lineage_sub_subtype <- NA
cell.lines.missing.df$lineage_molecular_subtype <- NA
cell.lines.missing.df$culture_type <- NA

head(cell.lines.missing.df)

sample.full <- rbind(sample.info.subset,cell.lines.missing.df)
head(sample.full)


sample.full <- sample.full[match(colnames(fullmatTPM), rownames(sample.full)), ] 
head(sample.full)


#DGElist object for easier transport =========================

#https://github.com/yunshun/HumanBreast10X/blob/main/RCode/BulkRNAseq.R

library(edgeR)
group <- factor(colnames(fullmatTPM))
Bulk.Object <- DGEList(counts = fullmatTPM,
                       samples = sample.full,
                       group = group)

Bulk.Object$counts[1:5,1:5] #raw expression matrix
head(Bulk.Object$samples) #metadata
Bulk.Object$samples$group <- group #sample names

saveRDS(Bulk.Object, "cell_lines_TPM_genesconvert_withmeta_72822.rds")


