


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This will generate all of the sub-figures for Figure 2,
# and is the exploratory analysis for the NK subsets
# as well as the analysis around the rNK population.


#_____________________________________________________________________

##note: for extra code check files 

##ReprogSCT.R, NKAnalysis.R, similboxplot.R in
## /project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only 

##and simil_reprogNK2_s204665.R, scRNA_seq_Functions.R in
## /project/InternalMedicine/Chan_lab/shared

##and Prim_Subset_Objects.R, Reprogrammed_ID.R in 
## /project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA

##and SurvivalFinal.R in /project/InternalMedicine/Chan_lab/shared/FinalObjects/TCGA_Analysis

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


sharedDir <- "/project/InternalMedicine/Chan_lab/shared/"

PrimDir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/PrimaryBreast_Processing/Main_Prim_Object"
NKsubDir <- "/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only"
DEGdir <- "/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only"

silhouette_dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/ReviewerAdditions/silhouette_test"
subcluster_dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/subclustering"

rNKdir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2"
rNKcell_dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/rNK_labels/"

nichenet_dir <- "/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/NichenetNKsub/"

TCGA_dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/bulkRNA/TCGA/TCGA_Objects"
TCGAresults_Dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2"
Metabric_Dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/bulkRNA/Metabric/Metabric_Objects"


ReprogSig.final = c("CALD1", "CLU", "ALOX12", "LTBP1", "CAVIN2", "PARVB", 
                    "GP6", "SCD", "ITGAX", "NR4A3", "CCL4", "CR2", "HEATR9", 
                    "XDH", "RASGRP2", "MID1", "JUN", "CMKLR1", "DUSP1", "FOS", 
                    "ABCA1", "TNFAIP3", "NR4A1", "KLRG1", "DTX1", "NHSL2", 
                    "GFRA2", "FAM81A", "CX3CR1", "RHPN1", "HES1", "F5", "GAS2L1", 
                    "THBS1", "MYLK", "TMTC1", "FOSB", "NR4A2", "MPIG6B", "SLC6A4", 
                    "PLXNA4", "VWF", "TUBB1", #up
                    "BCAT1-", "ALDH1L2-", "COX6A2-", "PYCR1-", "LHFPL2-", "AHRR-", 
                    "EXTL1-", "ASNS-", "CHAC1-", "MTHFD2-", "NEK6-", "SLC6A9-", 
                    "FMNL2-", "ASB2-", "SLC7A3-", "AVIL-", "CDH1-", "CISH-", "LGALS3-", 
                    "GPT2-", "CXCR6-", "TRIB3-", "CDKN1A-", "ATF5-", "SLC7A5-",
                    "SLC1A4-", "PMEPA1-", "CEMIP2-", "OSGIN1-", "ZNF503-", "ITGA1-", "ISG20-", 
                    "PACSIN1-", "TBC1D16-", "RN7SL1-", "SH3PXD2B-", "SCN3B-", "OSBPL1A-", 
                    "ME1-", "HPGDS-", "PPP2R2C-", "CLBA1-", "HMOX1-", "NQO1-", 
                    "CARS1-", "SSTR2-", "SNORA23-") #down



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

library(TFEA.ChIP)
library(org.Hs.eg.db)
library(TCGAbiolinks)
library(survminer)
library(survival)
library(openxlsx)
library(edgeR)
library(RTCGA)
library(RTCGA.clinical)
library(pheatmap)
library(cetcolor)
library(readxl)
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






# ================================================================== ======

# FUNCTIONS ----- ------------------------------

SCT_Ref_Integrate_Seurat <- function(seurat_obj, 
                                     batch = "Capture.Method", 
                                     regress.vars = c("percent.mt", "nCount_RNA",
                                                      "nFeature_RNA", "percent.hb","percent.platelet", "percent.heatshock"),
                                     num_integration_features = 3000,
                                     reference_datasets,
                                     project_name, 
                                     output_dir,
                                     save_integration_anchors = TRUE,
                                     save_integrated_object = TRUE)
{
  
  
  combo.list <- SplitObject(seurat_obj, split.by = batch)
  
  #https://github.com/satijalab/sctransform/issues/94 <- sctransform for read counts (smart-seq2)
  for (i in 1:length(combo.list)) {
    combo.list[[i]] <- SCTransform(combo.list[[i]], verbose = T, vars.to.regress = regress.vars)
  }
  
  
  features <- SelectIntegrationFeatures(object.list = combo.list, nfeatures = num_integration_features)
  combo.list <- PrepSCTIntegration(object.list = combo.list, anchor.features = features)
  
  reference.all <- c()
  for (ref in reference_datasets)
  {
    reference.1 <-  which(names(combo.list) == c(ref))
    reference.all <- c(reference.all, reference.1)
  }
  
  dimension.list <- c()
  for (i in 1:length(combo.list))
  {
    dimension.list <- c(dimension.list, dim(combo.list[[i]])[2])
  }
  
  if (any(dimension.list < 30))
  {
    
    combo.anchors <- FindIntegrationAnchors(object.list = combo.list, normalization.method = "SCT",
                                            anchor.features = features, reference = reference.all,
                                            k.score = (min(dimension.list)-1), dims = 1:(min(dimension.list)-1))
    setwd(output_dir)
    if (save_integration_anchors == T){
      saveRDS(combo.anchors, file = paste0(project_name,"_anchors_",Sys.Date(),".rds"))
    }
    
    combo.reference <- IntegrateData(anchorset = combo.anchors, normalization.method = "SCT", k.weight = (min(dimension.list)-1))
    
    if (save_integrated_object == T){
      saveRDS(combo.reference, file = paste0(project_name,"_integrated_",Sys.Date(),".rds"))
    }
    
  }
  
  else
  {
    combo.anchors <- FindIntegrationAnchors(object.list = combo.list, normalization.method = "SCT",
                                            anchor.features = features, reference = reference.all)
    
    
    setwd(output_dir)
    if (save_integration_anchors == T){
      saveRDS(combo.anchors, file = paste0(project_name,"_anchors_",Sys.Date(),".rds"))
    }
    
    combo.reference <- IntegrateData(anchorset = combo.anchors, normalization.method = "SCT")
    
    if (save_integrated_object == T){
      saveRDS(combo.reference, file = paste0(project_name,"_integrated_",Sys.Date(),".rds"))
    }
    
  }
  
  DefaultAssay(combo.reference) <- "integrated"
  
  return(combo.reference)
}

DEG_Remove_mito <- function (df){
  df_rm_mito <- df[!grepl("^MT-|^MT.",rownames(df)),]
  return(df_rm_mito)
}


human_mouse_conversion_genesymbol <- function(x)
{
  library(biomaRt)
  #searchAttributes(mart = ensembl, pattern = "hgnc")
  human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  
  ## use the two lines below if the function crashes and it says to try the ensembl server:
  #human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  #mouse <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  mousegenes <- x
  #mousegenes <- x$gene
  hggenes = getLDS(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id","description"),
                   filters = "external_gene_name", values = mousegenes,
                   mart = mouse, attributesL = c("hgnc_symbol","external_gene_name","entrezgene_id","ensembl_gene_id"),
                   martL = human, uniqueRows=T)
  library(data.table)
  hggenes <- as.data.table(hggenes)
  result <- setkey(hggenes, HGNC.symbol)
  return(result)
}

human_mouse_conversion_ensembl <- function(x)
{
  library(biomaRt)
  #searchAttributes(mart = ensembl, pattern = "hgnc")
  human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  
  ## use the two lines below if the function crashes and it says to try the ensembl server:
  #human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  #mouse <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  
  mousegenes <- x
  #mousegenes <- x$gene
  hggenes = getLDS(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id","description"),
                   filters = "ensembl_gene_id_version", values = mousegenes,
                   mart = mouse, attributesL = c("hgnc_symbol","external_gene_name","entrezgene_id","ensembl_gene_id"),
                   martL = human, uniqueRows=T)
  library(data.table)
  hggenes <- as.data.table(hggenes)
  result <- setkey(hggenes, HGNC.symbol)
  return(result)
}

simil <- function(df, drop, file, method) {
  # drop the input gene list from analysis
  if (length(drop) > 0) {
    jc <- df[-drop, , drop = TRUE]
  }
  else {
    jc <- df
  }
  #jc[jc > 0] <- 1
  #jc[jc <= 0] <- -1
  jc <- as.matrix(jc)
  # calculate jaccard similarity matrix if method == "jaccard"
  if (method == "jaccard") {
    jc <- prabclus::jaccard(jc)
    jc <- 1 - jc
  }
  # calculate correlation matrix if method == "corr"
  else if (method == "corr") {
    jc <- Rfast::cora(jc)
  }
  # save file
  saveRDS(jc, file = file)
  # return quantiles and mean similarity value
  v <- c(quantile(as.vector(jc), na.rm = TRUE), mean(as.vector(jc), na.rm = TRUE))
  names(v) <- c('0%','25%', '50%', '75%', '100%', 'mean')
  print(v)
  return(v)
}



# ================================================================== ======
# ================================================================== ======


#2A NK Clustering ========================
# ================================================================== ======
# load in Prim object ------

# setwd(PrimDir)
# combo.reference <- readRDS("PrimObject_withreprog_noZallgenedem_71322.rds")
# DefaultAssay(combo.reference) <- "RNA"
# combo.reference <- NormalizeData(combo.reference, assay = "RNA")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2")
combo.reference <- readRDS("PrimObject_withCORRECTreprog_nootherchange_82422.rds") ##newnewnew
setwd("/project/InternalMedicine/Chan_lab/shared/")
combo.inferCNV <- readRDS("nichenetobj_111122.rds")

combo.reference@meta.data <- combo.reference@meta.data[,-c(82,85)] 
head(combo.inferCNV@meta.data)

unique(colnames(combo.inferCNV@meta.data))
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"totalCNV", drop = F], "totalCNV")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"corCNV", drop = F], "corCNV")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"cancer", drop = F], "cancer")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"ESR1", drop = F], "ESR1")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"ERBB2", drop = F], "ERBB2")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"PIK3CA", drop = F], "PIK3CA")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"NTRK", drop = F], "NTRK")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"CD274", drop = F], "CD274")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"ERBB3", drop = F], "ERBB3")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"EGFR", drop = F], "EGFR")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"FGFR", drop = F], "FGFR")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"TACSTD2", drop = F], "TACSTD2")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"CDK", drop = F], "CDK")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"AR", drop = F], "AR")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"NECTIN2", drop = F], "NECTIN2")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"LAG3", drop = F], "LAG3")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"raw_GM1", drop = F], "raw_GE1")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"raw_GM2", drop = F], "raw_GE2")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"raw_GM3", drop = F], "raw_GE3")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"raw_GM4", drop = F], "raw_GE4")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"raw_GM5", drop = F], "raw_GE5")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"raw_GM6", drop = F], "raw_GE6")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"raw_GM7", drop = F], "raw_GE7")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"raw_GM8", drop = F], "raw_GE8")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"raw_GM9", drop = F], "raw_GE9")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"raw_GM10", drop = F], "raw_GE10")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"maxGM", drop = F], "maxGE")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"maxZscore", drop = F], "maxZscore")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM1", drop = F], "GE1")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM2", drop = F], "GE2")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM3", drop = F], "GE3")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM4", drop = F], "GE4")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM5", drop = F], "GE5")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM6", drop = F], "GE6")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM7", drop = F], "GE7")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM8", drop = F], "GE8")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM9", drop = F], "GE9")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM10", drop = F], "GE10")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM1_idents", drop = F], "GE1_idents")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM2_idents", drop = F], "GE2_idents")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM3_idents", drop = F], "GE3_idents")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM4_idents", drop = F], "GE4_idents")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM5_idents", drop = F], "GE5_idents")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM6_idents", drop = F], "GE6_idents")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM7_idents", drop = F], "GE7_idents")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM8_idents", drop = F], "GE8_idents")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM9_idents", drop = F], "GE9_idents")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"GM10_idents", drop = F], "GE10_idents")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"celltype_final", drop = F], "celltype_final")
combo.reference <- AddMetaData(combo.reference, combo.inferCNV@meta.data[,"celltype_withreprog", drop = F], "celltype_withreprog")

head(combo.reference@meta.data)
table(combo.reference$celltype_final)
table(combo.reference$celltype_withreprog)

setwd("/project/InternalMedicine/Chan_lab/shared")
#saveRDS(combo.reference, "PrimObj_withGEmeta_with_correctassays_withfinalreprog_111622.rds")


combo.reference <- readRDS("PrimObj_withGEmeta_with_correctassays_withfinalreprog_111622.rds")
DefaultAssay(combo.reference) <- "RNA"
combo.reference <- NormalizeData(combo.reference, assay = "RNA")

# NK clustering ------------------

table(Idents(combo.reference))
Idents(combo.reference) <- combo.reference$celltype_final

NKsub <- subset(combo.reference, idents = "NK Cells")
DefaultAssay(NKsub) <- "RNA"
NKsub <- NormalizeData(NKsub, assay = "RNA")


table(NKsub$Capture.Method)
NKsub.list <- SplitObject(NKsub, split.by = "Capture.Method")

for (i in 1:length(NKsub.list)) {
  NKsub.list[[i]] <- SCTransform(NKsub.list[[i]], verbose = T, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                                                   "percent.platelet", "percent.heatshock"))
}

NK.features <- SelectIntegrationFeatures(object.list = NKsub.list, nfeatures = 3000)
NKsub.list <- PrepSCTIntegration(object.list = NKsub.list, anchor.features = NK.features)


reference.1 <-  which(names(NKsub.list) == c("10X Genomics Chromium"))
reference.2 <-  which(names(NKsub.list) == c("10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"))
reference.3 <-  which(names(NKsub.list) == c("10X Genomics Chromium v2 5'"))
reference.4 <-  which(names(NKsub.list) == c("10X Genomics Single Cell 3' v2"))

reference.list <- c(reference.1, reference.2, reference.3, reference.4)

NKsub.anchors <- FindIntegrationAnchors(object.list = NKsub.list, normalization.method = "SCT",
                                        anchor.features = NK.features, reference = reference.list, 
                                        k.score = 27, dims = 1:27)

NK.all.combo <- IntegrateData(anchorset = NKsub.anchors, normalization.method = "SCT", k.weight = 27)
setwd(subcluster_dir)
saveRDS(NK.all.combo, "NK_subintegrated_noclustering_fromCNVobj_111622.rds")

DefaultAssay(NK.all.combo) <- "integrated"
#https://github.com/satijalab/seurat/issues/1963
NK.all.combo <- RunPCA(NK.all.combo, npcs = 300, verbose = FALSE, approx=FALSE)

DimPlot(NK.all.combo, reduction = "pca", raster = F, group.by = "Capture.Method",
        order = c("Smart-Seq2", "inDrop v2", "Singleron GEXSCOPE Single Cell RNAseq Library Kit",
                  "10X Genomics Single Cell 3' v2", "10X Genomics Chromium v2 5'",
                  "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library",
                  "10X Genomics Chromium"))
DimPlot(NK.all.combo, reduction = "pca", raster = F, group.by = "orig.ident",
        order =c("Karaayvaz", "Savas", "Wu", "Aziziimmune", "Xu",
                 "AziziT", "Qian", "Wu2021prim", "Pal_Prim"))

DimHeatmap(NK.all.combo, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(NK.all.combo, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(NK.all.combo, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(NK.all.combo, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(NK.all.combo, dims = 45:55, cells = 500, balanced = TRUE)
DimHeatmap(NK.all.combo, dims = 55:65, cells = 500, balanced = TRUE)


NK.all.combo <- FindNeighbors(NK.all.combo, reduction = "pca", dims = 1:18)
resolution.range <- c(0.01, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6, 1.8, 2.0)
NK.all.combo <- FindClusters(NK.all.combo, resolution = resolution.range)
clustree(NK.all.combo, prefix = "integrated_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")

DefaultAssay(NK.all.combo) <- "integrated"
NK.all.combo <- FindNeighbors(NK.all.combo, reduction = "pca", dims = 1:18, annoy.metric = "manhattan")
NK.all.combo <- FindClusters(NK.all.combo, resolution = 0.3) 
NK.all.combo <- RunUMAP(NK.all.combo, reduction = "pca", dims = 1:18, verbose = TRUE, seed.use=123)

NK.all.combo$celltype_NKsub <- Idents(NK.all.combo)


library(viridis)
p <- DimPlot(NK.all.combo, reduction = "umap", label = F, 
             repel = TRUE, raster = FALSE, pt.size = 1.5,
             cols = c(turbo(6))) + ggtitle(label = " ")


pdf("test.pdf")
p + theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 25)) +
  theme(axis.text = element_text(size = 20))
dev.off()

DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")

setwd(NKsubDir)
saveRDS(NK.all.combo, "NKall_PC17manhattan0.4res_71322.rds")


# silhouette (REVIEWER ADDITION) +++++++++++++++++++++++++++++++++++++++++++++ ---------

library(Seurat)
library(cluster)
library(tidyverse)
library(viridis)
library(crayon)

#https://www.biorxiv.org/content/10.1101/2022.05.31.494081v1.full

setwd(subcluster_dir)
#start with just integrated object, no clustering/PC calc on it
#setwd(silhouette_dir)
#NK.all.combo <- readRDS("NK_subintegrated_noclustering_11822.rds")
NK.all.combo <- readRDS("NK_subintegrated_noclustering_fromCNVobj_111622.rds")

DefaultAssay(NK.all.combo) <- "integrated"
#https://github.com/satijalab/seurat/issues/1963
NK.all.combo <- RunPCA(NK.all.combo, npcs = 300, verbose = FALSE, approx=FALSE)

DimHeatmap(NK.all.combo, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(NK.all.combo, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(NK.all.combo, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(NK.all.combo, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(NK.all.combo, dims = 45:55, cells = 500, balanced = TRUE)
DimHeatmap(NK.all.combo, dims = 55:65, cells = 500, balanced = TRUE)

##this works as long as you don't rename the columns!!!!

# testPCs <- c(18, seq(from = 5, to = 100, by = 5))
# resolutions <- c(seq(from = 0.5, to = 6.5, by = 0.5), seq(from = 7.5, to = 10, by = 0.5)) #for some reason, 7 doesn't work??

resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6, 1.8, 2.0)
testPCs <- 18
mean_silhouette_score_list <- list()
#loopy
for (pc in testPCs)
{
  for (res in resolutions)
  {
    
    cat(blue$bold("Clustering for: \n"))
    print(paste0("PC number: ", pc))
    print(paste0("Resolution number: ", res))
    
    DefaultAssay(NK.all.combo) <- "integrated"
    NK.all.combo <- FindNeighbors(NK.all.combo, reduction = "pca", dims = 1:pc, annoy.metric = "manhattan")
    NK.all.combo <- FindClusters(NK.all.combo, resolution = res) 
    NK.all.combo <- RunUMAP(NK.all.combo, reduction = "pca", dims = 1:pc, verbose = TRUE, seed.use=123)
    
    ## THESE LINES ARE IFFY TBH __________
    # print(ncol(NK.all.combo@meta.data))
    # print(colnames(NK.all.combo@meta.data)[c((ncol(NK.all.combo@meta.data) - 1):(ncol(NK.all.combo@meta.data)))])
    # colnames(NK.all.combo@meta.data)[ncol(NK.all.combo@meta.data)] <- paste0(colnames(NK.all.combo@meta.data)[ncol(NK.all.combo@meta.data)],"_Dim.", pc)
    # print(colnames(NK.all.combo@meta.data)[c((ncol(NK.all.combo@meta.data) - 1):(ncol(NK.all.combo@meta.data)))])
    ##after this is fine, though if you don't run above section, input_res object will fail _________________
    
    
    #silhouette _________
    cat(blue$bold("Calculating silhouette score for: \n"))
    print(paste0("PC number: ", pc))
    print(paste0("Resolution number: ", res))
    
    distance_matrix <- dist(Embeddings(NK.all.combo[['pca']])[, 1:pc])
    ##if don't run iffy block, this next immediate line will fail because it assumes the metadata columns were renamed
    input_res <- NK.all.combo@meta.data[[paste0("integrated_snn_res.",res)]]#, "_Dim.", pc)]]
    silhouette <- silhouette(as.numeric(input_res), dist = distance_matrix)
    NK.all.combo@meta.data$silhouette_score <- silhouette[,3]
    add_meta <- as.data.frame(NK.all.combo@meta.data$silhouette_score)
    rownames(add_meta) <- row.names(NK.all.combo@meta.data)
    
    meta_colname <- paste0("silhouette_score_res.", res, "_Dim.", pc)
    cat(blue$bold("meta_colname: \n"))
    print(meta_colname)
    
    NK.all.combo <- AddMetaData(NK.all.combo, metadata = add_meta, col.name = meta_colname)
    
    mean_silhouette_score_list[[meta_colname]] <- mean(NK.all.combo@meta.data[[meta_colname]])
  }
  
  # cat(blue$bold("Saving object for: \n"))
  # print(paste0("PC number: ", pc))
  # 
  # setwd(paste0(silhouette_dir,"/NK_objects/"))
  # saveRDS(NK.all.combo, paste0("NK_clustered_silhouette_PC", pc, "_res", res, Sys.Date(), ".rds"))
}

setwd(subcluster_dir)
saveRDS(mean_silhouette_score_list, "mean_NKsihouette_PC18_res0.1to2_111622.rds")

# saveRDS(NK.all.combo, "NK_with_silhouette_PCs5to100_res0.5to10_11822.rds")
# saveRDS(mean_silhouette_score_list2, "mean_NKsihouette_PCs5to100_res0.5to10_11822.rds")

colnames(NK.all.combo@meta.data)

# NK.all.combo <- readRDS("NK_with_silhouette_PCs5to100_res0.5to10_11822.rds")
# mean_silhouette_score_list <- readRDS("mean_NKsihouette_PCs5to100_res0.5to10_11822.rds")


colnames(NK.all.combo@meta.data)[seq(88, 124, by = 2)]
setwd(paste0(subcluster_dir, "/Silhouette_pdfs"))

for (indsilette in colnames(NK.all.combo@meta.data)[c(seq(88, 124, by = 2), 886)])
{
  p <- FeaturePlot(object = NK.all.combo, features = indsilette, 
                   order = TRUE, label = FALSE, repel = TRUE, 
                   min.cutoff = 0, raster = FALSE, pt.size = 1.5) +
    #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    scale_color_viridis(option = "D")+
    ggtitle(label = " ")
  
  pdf(paste0(indsilette,"_featureplot_", Sys.Date(),".pdf"), width = 7.11, height = 6.5)
  print(
    p + theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
      theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
      theme(text = element_text(size = 25)) +
      theme(axis.text = element_text(size = 20))
  )
  dev.off()
}


setwd(subcluster_dir)
# for (col in colnames(NK.all.combo@meta.data)[seq(88, 124, by = 2)])#143)])
# {
#   for (indmean in mean_silhouette_score_list)
#   {
#     col_num <- which(colnames(NK.all.combo@meta.data)==col )
#     pdf(paste0(col, "_Silouhette_", Sys.Date(), ".pdf"), height = 6,width=7)
#     print(
#       test <- NK.all.combo@meta.data %>%
#         mutate(barcode = rownames(.)) %>%
#         arrange(celltype_NKsub,-col_num) %>%
#         mutate(barcode = factor(barcode, levels = barcode)) %>%
#         ggplot() + 
#         geom_col(aes(barcode, col, fill = celltype_NKsub), show.legend = TRUE) +
#         geom_hline(yintercept = indmean, color = 'red', linetype = 'dashed') +
#         scale_x_discrete(name = 'Cells') +
#         scale_y_continuous(name = 'Silhouette score') +
#         scale_fill_manual(values = c(turbo(6))) +
#         theme_bw() + 
#         theme(
#           axis.title.x = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank()
#         )
#     )
#     dev.off()
#   }
#   
# }


# #no loopy - settings for original clustering (PC = 18, cell_label = celltype_NKsub)
# resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6, 1.8, 2.0)
# 
# for(res in resolutions)
# {
#   DefaultAssay(NK.all.combo) <- "integrated"
#   NK.all.combo <- FindNeighbors(NK.all.combo, reduction = "pca", dims = 1:18, annoy.metric = "manhattan")
#   NK.all.combo <- FindClusters(NK.all.combo, resolution = res) 
#   NK.all.combo <- RunUMAP(NK.all.combo, reduction = "pca", dims = 1:18, verbose = TRUE, seed.use=123)
# 
# }

colnames(NK.all.combo@meta.data)[c(133, seq(from = 136, to = 156, by = 2))] <- c("integrated_snn_res.0.1_Dim.18",
                                                                                 "integrated_snn_res.0.2_Dim.18",
                                                                                 "integrated_snn_res.0.3_Dim.18", "integrated_snn_res.0.4_Dim.18",
                                                                                 "integrated_snn_res.0.5_Dim.18", "integrated_snn_res.0.6_Dim.18",
                                                                                 "integrated_snn_res.0.7_Dim.18", "integrated_snn_res.1_Dim.18",
                                                                                 "integrated_snn_res.1.3_Dim.18", "integrated_snn_res.1.6_Dim.18", 
                                                                                 "integrated_snn_res.1.8_Dim.18", "integrated_snn_res.2_Dim.18")



# distance_matrix <- dist(Embeddings(NK.all.combo[['pca']])[, 1:18])
# clusters <- NK.all.combo@meta.data$integrated_snn_res.2_Dim.18
# silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
# NK.all.combo@meta.data$silhouette_score <- silhouette[,3]
# add_meta <- as.data.frame(NK.all.combo@meta.data$silhouette_score)
# rownames(add_meta) <- row.names(NK.all.combo@meta.data)
# NK.all.combo <- AddMetaData(NK.all.combo, metadata = add_meta, col.name = "silhouette_score_res.2_Dim.18")

c("integrated_snn_res.0.1_Dim.18",
  "integrated_snn_res.0.2_Dim.18",
  "integrated_snn_res.0.3_Dim.18", "integrated_snn_res.0.4_Dim.18",
  "integrated_snn_res.0.5_Dim.18", "integrated_snn_res.0.6_Dim.18",
  "integrated_snn_res.0.7_Dim.18", "integrated_snn_res.1_Dim.18",
  "integrated_snn_res.1.3_Dim.18", "integrated_snn_res.1.6_Dim.18", 
  "integrated_snn_res.1.8_Dim.18", "integrated_snn_res.2_Dim.18")

p <- NK.all.combo@meta.data %>%
  mutate(barcode = rownames(.)) %>%
  arrange(integrated_snn_res.0.1_Dim.18,-silhouette_score_res.0.1_Dim.18) %>%
  mutate(barcode = factor(barcode, levels = barcode)) %>%
  ggplot() +
  geom_col(aes(barcode, silhouette_score_res.0.1_Dim.18, fill = integrated_snn_res.0.1_Dim.18), show.legend = TRUE) +
  geom_hline(yintercept = mean_silhouette_score_list$silhouette_score_res.0.1_Dim.18, color = 'red', linetype = 'dashed') +
  scale_x_discrete(name = 'Cells') +
  scale_y_continuous(name = 'Silhouette score') +
  #scale_fill_manual(values = c(turbo(6))) +
  #scale_fill_manual(values = c(turbo(12))) +
  scale_fill_manual(values = c(turbo(6))) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

setwd(paste0(subcluster_dir, "/Silhouette_pdfs"))
pdf("NKsub_silouhette_silhouette_score_res.0.1_Dim.18_111122.pdf", width = 10, height =5)
#pdf("NKsub_silouhette_silhouette_score_res.0.2_Dim.18_111122.pdf", width = 10, height =5)
#pdf("NKsub_silouhette_silhouette_score_res.0.3_Dim.18_111122.pdf", width = 10, height =5)
#pdf("NKsub_silouhette_silhouette_score_res.0.4_Dim.18_111122.pdf", width = 10, height =5)
#pdf("NKsub_silouhette_silhouette_score_res.0.5_Dim.18_111122.pdf", width = 10, height =5)
#pdf("NKsub_silouhette_silhouette_score_res.0.6_Dim.18_111122.pdf", width = 10, height =5)
#pdf("NKsub_silouhette_silhouette_score_res.0.7_Dim.18_111122.pdf", width = 10, height =5)
#pdf("NKsub_silouhette_silhouette_score_res.1_Dim.18_111122.pdf", width = 10, height =5)
#pdf("NKsub_silouhette_silhouette_score_res.1.3_Dim.18_111122.pdf", width = 10, height =5)
#pdf("NKsub_silouhette_silhouette_score_res.1.6_Dim.18_111122.pdf", width = 10, height =5)
#pdf("NKsub_silouhette_silhouette_score_res.1.8_Dim.18_111122.pdf", width = 10, height =5)
#pdf("NKsub_silouhette_silhouette_score_res.2_Dim.18_111122.pdf", width = 10, height =5)
p
dev.off()

# mean_silhouette_score<- list()
# mean_silhouette_score$silhouette_score_res.0.1_Dim.18 <- mean(NK.all.combo@meta.data$silhouette_score_res.0.1_Dim.18)
# mean_silhouette_score$silhouette_score_res.0.2_Dim.18 <- mean(NK.all.combo@meta.data$silhouette_score_res.0.2_Dim.18)
# mean_silhouette_score$silhouette_score_res.0.3_Dim.18 <- mean(NK.all.combo@meta.data$silhouette_score_res.0.3_Dim.18)
# mean_silhouette_score$silhouette_score_res.0.4_Dim.18 <- mean(NK.all.combo@meta.data$silhouette_score_res.0.4_Dim.18)
# mean_silhouette_score$silhouette_score_res.0.5_Dim.18 <- mean(NK.all.combo@meta.data$silhouette_score_res.0.5_Dim.18)
# mean_silhouette_score$silhouette_score_res.0.6_Dim.18 <- mean(NK.all.combo@meta.data$silhouette_score_res.0.6_Dim.18)
# mean_silhouette_score$silhouette_score_res.0.7_Dim.18 <- mean(NK.all.combo@meta.data$silhouette_score_res.0.7_Dim.18)
# mean_silhouette_score$silhouette_score_res.1_Dim.18 <- mean(NK.all.combo@meta.data$silhouette_score_res.1_Dim.18)
# mean_silhouette_score$silhouette_score_res.1.3_Dim.18 <- mean(NK.all.combo@meta.data$silhouette_score_res.1.3_Dim.18)
# mean_silhouette_score$silhouette_score_res.1.6_Dim.18 <- mean(NK.all.combo@meta.data$silhouette_score_res.1.6_Dim.18)
# mean_silhouette_score$silhouette_score_res.1.8_Dim.18 <- mean(NK.all.combo@meta.data$silhouette_score_res.1.8_Dim.18)
# mean_silhouette_score$silhouette_score_res.2_Dim.18 <- mean(NK.all.combo@meta.data$silhouette_score_res.2_Dim.18)
# 
# setwd(silhouette_dir)
# 
# saveRDS(mean_silhouette_score, "mean_silhouette_pc18_res0.1to0.1to2_dim18_111122.rds")
# saveRDS(NK.all.combo, "NKobj_with_dim18_res0.1to2_111122.rds")
# 

mean_silhouette_score_df <- as.data.frame(unlist(mean_silhouette_score_list))
colnames(mean_silhouette_score_df) <- "mean_silhouette_score"
rownames(mean_silhouette_score_df) <- gsub("silhouette_score_","", rownames(mean_silhouette_score_df))

write.csv(mean_silhouette_score_df, "mean_silhouette_score_NK_fromCNVobj_dim18_111622.csv")

# ================================================================== ======
#2B NKsub Markers ========================
# ================================================================== ======

# load NKsubset =================================

setwd(NKsubDir)
NK.all.combo <- readRDS("NKall_PC17manhattan0.4res_71322.rds")
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")


# FindMarkers (NK) ======================================================

Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub

NK.mark <- FindAllMarkers(NK.all.combo, test.use = "MAST")
NK.mark.nomit <- DEG_Remove_mito(NK.mark)

head(NK.mark.nomit)
NK.mark.nomit <- NK.mark.nomit[NK.mark.nomit$p_vabl_adj <= 0.05,]
NK.mark.nomit <- NK.mark.nomit[abs(NK.mark.nomit$avg_log2FC) >= 0.56,]
NK.mark.nomit$diff.pct = abs(NK.mark.nomit$pct.1 - NK.mark.nomit$pct.2)


setwd(DEGdir)
write.csv(NK.mark.nomit, "NKclustDEGs_pvaladj0.05log056_71322.csv")
write.csv(NK.mark, "NKclustDEGs_NOTHRESH_71322.csv")


# NK subtype markers (c0) =========================================

setwd(DEGdir)
NK.mark.nomit <- read.csv("NKclustDEGs_pvaladj0.05log056_71322.csv")


c0.keep <- NK.mark.nomit
head(c0.keep)
c0.keep <- c0.keep[c0.keep$cluster == "0",]


c0.up <- c0.keep[c0.keep$avg_log2FC >0,]
# c0.up <- c0.up[c0.up$pct.1 >= 0.2,] #newcutoff
c0.up.final <- c0.up$gene

c0.down <- c0.keep[c0.keep$avg_log2FC < 0,]
# c0.down <- c0.down[c0.down$diff.pct >= 0.1,]
# c0.down <- c0.down[abs(c0.down$avg_log2FC) >= 0.7,]
c0.down.final <- c0.down$gene
c0.down.final <- paste(c0.down.final, "-", sep="")

c0.mark <- list(c(c0.up.final, c0.down.final))

saveRDS(c0.mark, "NK.c0markers_71322.rds")




# NK subtype markers (c1) =========================================

c1.keep <- NK.mark.nomit
head(c1.keep)
c1.keep <- c1.keep[c1.keep$cluster == "1",]
#c1.keep <- c1.keep[abs(c1.keep$avg_log2FC) >= 0.56,]

c1.up <- c1.keep[c1.keep$avg_log2FC >0,]
c1.up <- c1.up[c1.up$diff.pct >= 0.218,]
c1.up.final <- c1.up$gene

c1.down <- c1.keep[c1.keep$avg_log2FC < 0,]
c1.down <- c1.down[abs(c1.down$avg_log2FC) >= 0.65,]
c1.down <- c1.down[c1.down$diff.pct >= 0.1,]
c1.down.final <- c1.down$gene
c1.down.final <- paste(c1.down.final, "-", sep="")

c1.mark <- list(c(c1.up.final, c1.down.final))
saveRDS(c1.mark, "NK.c1markers_71322.rds")


# NK subtype markers (c2) =========================================

c2.keep <- NK.mark.nomit
c2.keep <- c2.keep[c2.keep$cluster == "2",]
head(c2.keep)

c2.up <- c2.keep[c2.keep$avg_log2FC >0,]
c2.up <- c2.up[c2.up$diff.pct >0.21,]
c2.up <- c2.up[abs(c2.up$avg_log2FC) >0.75,]
c2.up.final <- c2.up$gene

c2.down <- c2.keep[c2.keep$avg_log2FC < 0,]
c2.down <- c2.down[c2.down$diff.pct >0.1,]
c2.down <- c2.down[abs(c2.down$avg_log2FC) >0.75,]
c2.down.final <- c2.down$gene
c2.down.final <- paste(c2.down.final, "-", sep="")

c2.mark <- list(c(c2.up.final, c2.down.final))
saveRDS(c2.mark, "NK.c2markers_71322.rds")


# NK subtype markers (c3) =========================================

c3.keep <- NK.mark.nomit
c3.keep <- c3.keep[c3.keep$cluster == "3",]
head(c3.keep)

c3.up <- c3.keep[c3.keep$avg_log2FC >0,]
c3.up <- c3.up[c3.up$avg_log2FC >0.8,]
c3.up <- c3.up[c3.up$diff.pct >0.1,]
c3.up.final <- c3.up$gene

c3.down <- c3.keep[c3.keep$avg_log2FC < 0,]
c3.down <- c3.down[abs(c3.down$avg_log2FC) >0.8,]
c3.down <- c3.down[c3.down$diff.pct >0.15,]
c3.down.final <- c3.down$gene
c3.down.final <- paste(c3.down.final, "-", sep="")

c3.mark <- list(c(c3.up.final, c3.down.final))
saveRDS(c3.mark, "NK.c3markers_71322.rds")

# NK subtype markers (c4) =========================================

c4.keep <- NK.mark.nomit
c4.keep <- c4.keep[c4.keep$cluster == "4",]
head(c4.keep)

c4.up <- c4.keep[c4.keep$avg_log2FC >0,]
c4.up <- c4.up[c4.up$diff.pct >0.28,]
c4.up.final <- c4.up$gene

c4.down <- c4.keep[c4.keep$avg_log2FC < 0,] ##none
# c4.down.final <- c4.down$gene
# c4.down.final <- paste(c4.down.final, "-", sep="")

c4.mark <- list(c(c4.up.final))
saveRDS(c4.mark, "NK.c4markers_71322.rds")


# NK subtype markers (c5) =========================================

c5.keep <- NK.mark.nomit
c5.keep <- c5.keep[c5.keep$cluster == "5",]
head(c5.keep)

c5.up <- c5.keep[c5.keep$avg_log2FC > 0,]
c5.up <- c5.up[c5.up$diff.pct > 0.25,]
c5.up.final <- c5.up$gene

c5.down <- c5.keep[c5.keep$avg_log2FC < 0,]
#c5.down <- c5.down[c5.down$diff.pct > 0.1,]
c5.down.final <- c5.down$gene
c5.down.final <- paste(c5.down.final, "-", sep="")

c5.mark <- list(c(c5.up.final, c5.down.final))
saveRDS(c5.mark, "NK.c5markers_71322.rds")


# NK subset Bubble heatmap -----------------


keep <- NK.mark.nomit[NK.mark.nomit$cluster == "5",]

keep[keep$gene == "HLA-DRB1", ]




setwd(DEGdir)
c0.mark <- unlist(readRDS("NK.c0markers_71322.rds"))
c0.up <- c("FCGR3A", "PRF1", "KLRC2", c0.mark[c(1,2,4)]) #not in DEGs: FCGR3A,KLRC2,PRF1
c1.mark <- unlist(readRDS("NK.c1markers_71322.rds"))
c1.up <- c1.mark[c(11,22,4,13,2,8)]
#c1.up <- c1.mark[1:23] ##preprint
c2.mark <- unlist(readRDS("NK.c2markers_71322.rds"))
c2.up <- c("FGFBP2", "GZMA", "GZMB", "CXCR1",c2.mark[c(1,4,2,3,5)]) #not in DEGs: GZMB, CXCR1, CXC3R1, S1PR5
c3.mark <- unlist(readRDS("NK.c3markers_71322.rds"))
c3.up <- c("GZMK", c3.mark[c(5,7,10)]) #not in DEGs: GZMK
c4.mark <- unlist(readRDS("NK.c4markers_71322.rds"))
c4.up <- c4.mark[c(2,1,4,6)]
#c4.up <- c4.mark ##preprint
c5.mark <- unlist(readRDS("NK.c5markers_71322.rds"))
c5.up <- c("CCL5", "HLA-DRB1", c5.mark[c(12,1,10)])


features <- list("NK0" = c0.up,
                 "NK1" = c1.up,
                 "NK2" = c2.up,
                 "NK3" = c3.up,
                 "NK4" = c4.up,
                 "NK5" = c5.up)

myLevels <- c("0", "1", "2", "3",
              "4", "5", "6")
Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub
factor(Idents(NK.all.combo), levels= myLevels)
Idents(NK.all.combo) <- factor(Idents(NK.all.combo), levels= rev(myLevels))

table(Idents(NK.all.combo))

features2 <- unique(c(features$NK0, features$NK1, features$NK2, 
               features$NK3, features$NK4, features$NK5))
a <- DotPlot(object = NK.all.combo, features=features2,
             dot.scale = 10) + theme(axis.text.x = element_text(angle = 90)) +
  scale_colour_gradient2(low="steelblue", mid="lightgrey", high="red")

#pdf("NKsubMarkers_71322.pdf", width = 26.3, height = 4.9)
setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2")
pdf("NKsubMarkers2B_112522.pdf", width = 26.3, height = 4.9)
a+ theme(axis.line = element_line(colour = 'black', size = 1.5)) +
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(size = 18)) +
  theme(strip.text = element_text(size=20))
dev.off()

# ================================================================== ======
#2D rNK sig per clust  ========================
# ================================================================== ======
# converting mouse rNK to human rNK gene names ==================================

library(biomaRt)

reproglist_mus_all <- as.character(c("Cald1", "Clu", "Alox12", "Ltbp1", "Fhl1", "Cavin2", "Parvb", "Gp6",
                                     "Scd2","Itgax", "Nr4a3", "Ccl4", "Cr2", "Heatr9", "Xdh", "Rasgrp2", "Mid1",  
                                     "Jun", "Cmklr1", "Dusp1", "Fos", "Abca1", "Tnfaip3", "Nr4a1","Klrg1","Dtx1", "Nhsl2","Gfra2", "Fam81a",
                                     "Cx3cr1","Rhpn1", "Hes1", "F5", "Gas2l1", "Thbs1", "Mylk", "Tmtc1", "Fosb", "Nr4a2",
                                     "Mpig6b", "Slc6a4", "Plxna4", "Vwf", "Tubb1",
                                     "ENSMUSG00000045065", "ENSMUSG00000085180", "ENSMUSG00000098292", "ENSMUSG00000108090",
                                     "Bcat1", "Aldh1l2", "Cox6a2", "Pycr1", "Lhfpl2", "Ahrr", "Extl1", "Asns", "Chac1",
                                     "Mthfd2", "Nek6", "Slc6a9", "Fmnl2", "Asb2", "Slc7a3", "Avil", "Cdh1", "Cish", "Trp53cor1",
                                     "Lgals3", "1700017B05Rik", "Gpt2", "Cxcr6", "Trib3", "Cdkn1a","Atf5", "Slc7a5", "Slc1a4", "Pmepa1",
                                     "Cemip2", "Osgin1","Zfp503", "Itga1", "Isg20", "Pacsin1", "Tbc1d16", "Rn7s1", "Sh3pxd2b",
                                     "Scn3b", "Osbpl1a", "Me1", "Hpgds", "Ppp2r2c", "Clba1", "Hmox1", "Nqo1", "Cars", "Sstr2", "Snora23",
                                     "ENSMUSG00000108053", "ENSMUSG00000087113"))
searchFilters(mart = mouse, pattern = "gene") #check what format your genes are in

reproglist_human_all <- human_mouse_conversion_genesymbol(reproglist_mus_all)

write.csv(reproglist_human_all, file = "allreprogupanddownconverted.csv")


# load NKsubset =================================

setwd(NKsubDir)
NK.all.combo <- readRDS("NKall_PC17manhattan0.4res_71322.rds")
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")

# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/")
# NK.all.combo <- readRDS("NK_withNKpathways_noZallgenedem_71422.rds") #newnewnew
# DefaultAssay(NK.all.combo) <- "RNA"
# NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")

# separate into subtypes =======================

#NK.all.combo$celltype_NKsub <- Idents(NK.all.combo)
Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub

LumA_HR <- subset(NK.all.combo, (subset = pam50.pred == "HR+") | (subset = pam50.pred == "LumA"))
LumB_HR <- subset(NK.all.combo, (subset = pam50.pred == "HR+") | (subset = pam50.pred == "LumB"))
Basal_TNBC <- subset(NK.all.combo, (subset = pam50.pred == "TNBC") | (subset = pam50.pred == "Basal"))

LumA_HR@meta.data$celltype_NKsub <- Idents(LumA_HR)
LumB_HR@meta.data$celltype_NKsub <- Idents(LumB_HR)
Basal_TNBC@meta.data$celltype_NKsub <- Idents(Basal_TNBC)


# NK.LumAHR =========================================

#https://github.com/satijalab/seurat/issues/1883
table(LumA_HR$orig.ident)
DefaultAssay(LumA_HR) <- "RNA"
LumA_HR <- NormalizeData(LumA_HR, assay = "RNA")

#just to get Xu and Pal clustered
LumA_HR$Capture.Method[which(LumA_HR$Capture.Method == "Singleron GEXSCOPE Single Cell RNAseq Library Kit")] <- "inDrop v2"
LumA_HR$Capture.Method[which(LumA_HR$Capture.Method == "10X Genomics Chromium")] <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"

LumA_HR.list <- SplitObject(LumA_HR, split.by = "Capture.Method")
for (i in 1:length(LumA_HR.list)) {
  LumA_HR.list[[i]] <- SCTransform(LumA_HR.list[[i]], verbose = T, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                                                       "percent.platelet", "percent.heatshock"))
}

LumA_HR.features <- SelectIntegrationFeatures(object.list = LumA_HR.list, nfeatures = 3000)
LumA_HR.list <- PrepSCTIntegration(object.list = LumA_HR.list, anchor.features = LumA_HR.features)

reference.2 <-  which(names(LumA_HR.list) == c("10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"))

LumA_HR.anchors <- FindIntegrationAnchors(object.list = LumA_HR.list, normalization.method = "SCT",
                                          anchor.features = LumA_HR.features, reference = reference.2)

NK.lumAHR.combo <- IntegrateData(anchorset = LumA_HR.anchors, normalization.method = "SCT", k.weight = 69)

DefaultAssay(NK.lumAHR.combo) <- "integrated"
#https://github.com/satijalab/seurat/issues/1963
NK.lumAHR.combo <- RunPCA(NK.lumAHR.combo, npcs = 300, verbose = FALSE, approx=FALSE)

DimPlot(NK.lumAHR.combo, reduction = "pca", raster = F, group.by = "Capture.Method")
DimPlot(NK.lumAHR.combo, reduction = "pca", raster = F, group.by = "orig.ident")

DimHeatmap(NK.lumAHR.combo, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(NK.lumAHR.combo, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(NK.lumAHR.combo, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(NK.lumAHR.combo, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(NK.lumAHR.combo, dims = 45:55, cells = 500, balanced = TRUE)
DimHeatmap(NK.lumAHR.combo, dims = 55:65, cells = 500, balanced = TRUE)


NK.lumAHR.combo <- FindNeighbors(NK.lumAHR.combo, reduction = "pca", dims = 1:15)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
NK.lumAHR.combo <- FindClusters(NK.lumAHR.combo, resolution = resolution.range)
clustree(NK.lumAHR.combo, prefix = "integrated_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")


DefaultAssay(NK.lumAHR.combo) <- "integrated"
NK.lumAHR.combo <- FindNeighbors(NK.lumAHR.combo, reduction = "pca", dims = 1:15, annoy.metric = "manhattan")
NK.lumAHR.combo <- FindClusters(NK.lumAHR.combo, resolution = 0.4)
NK.lumAHR.combo <- RunUMAP(NK.lumAHR.combo, reduction = "pca", dims = 1:15, verbose = TRUE, seed.use=123)

DimPlot(NK.lumAHR.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE) + ggtitle(label = " ")
DimPlot(NK.lumAHR.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "celltype_NKsub") 
DimPlot(NK.lumAHR.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "Capture.Method") 
DimPlot(NK.lumAHR.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "orig.ident")


NK.lumAHR.combo$Capture.Method[which(NK.lumAHR.combo$orig.ident == "Xu")] <- "Singleron GEXSCOPE Single Cell RNAseq Library Kit"
NK.lumAHR.combo$Capture.Method[which(NK.lumAHR.combo$orig.ident == "Pal_Prim")] <- "10X Genomics Chromium"

DefaultAssay(NK.lumAHR.combo) <- "RNA"
NK.lumAHR.combo <- NormalizeData(NK.lumAHR.combo, assay = "RNA")


#All NK ____________________________________________________________________

ReprogSig <- list(ReprogSig.final)
NK.lumAHR.combo <- AddModuleScore_UCell(NK.lumAHR.combo, features = ReprogSig, name = "ReprogSig", assay = "RNA")
FeaturePlot(object = NK.lumAHR.combo, features = "signature_1ReprogSig", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, raster = FALSE) + ggtitle(label = "Reprogrammed NK Signature (Expression > 20%, PC = 40)")

nCount75 <- quantile(NK.lumAHR.combo$signature_1ReprogSig, 0.75)#95) # calculate value in the 95th percentile)

reprog <- WhichCells(NK.lumAHR.combo, expression = signature_1ReprogSig > nCount75)
setwd(rNKcell_dir)
#saveRDS(reprog, file = "reprog_lumAHR_SCTRNAnoZ_71322.rds")
saveRDS(reprog, file = "reprog_lumAHR_FIXEDrsig_82422.rds")

Idents(NK.lumAHR.combo, cells = reprog) <- "Reprogrammed NK Cells"

# NK.LumBHR =========================================

#https://github.com/satijalab/seurat/issues/1883
table(LumB_HR$orig.ident)
DefaultAssay(LumB_HR) <- "RNA"
LumB_HR <- NormalizeData(LumB_HR, assay = "RNA")

#just to get Wu and Xu and Qian cell clustered
LumB_HR$Capture.Method[which(LumB_HR$Capture.Method == "Singleron GEXSCOPE Single Cell RNAseq Library Kit")] <- "inDrop v2"
LumB_HR$Capture.Method[which(LumB_HR$Capture.Method == "10X Genomics Single Cell 3' v2")] <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
LumB_HR$Capture.Method[which(LumB_HR$Capture.Method == "10X Genomics Chromium v2 5'")] <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"

LumB_HR.list <- SplitObject(LumB_HR, split.by = "Capture.Method")
for (i in 1:length(LumB_HR.list)) {
  LumB_HR.list[[i]] <- SCTransform(LumB_HR.list[[i]], verbose = T, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                                                       "percent.platelet", "percent.heatshock"))
}

LumB_HR.features <- SelectIntegrationFeatures(object.list = LumB_HR.list, nfeatures = 3000)
LumB_HR.list <- PrepSCTIntegration(object.list = LumB_HR.list, anchor.features = LumB_HR.features)


reference.1 <-  which(names(LumB_HR.list) == c("10X Genomics Chromium"))
reference.2 <-  which(names(LumB_HR.list) == c("10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"))

reference.list <- c(reference.1, reference.2)

LumB_HR.anchors <- FindIntegrationAnchors(object.list = LumB_HR.list, normalization.method = "SCT",
                                          anchor.features = LumB_HR.features, reference = reference.list)

NK.lumBHR.combo <- IntegrateData(anchorset = LumB_HR.anchors, normalization.method = "SCT")#, k.weight = 48)

DefaultAssay(NK.lumBHR.combo) <- "integrated"
#https://github.com/satijalab/seurat/issues/1963
NK.lumBHR.combo <- RunPCA(NK.lumBHR.combo, npcs = 300, verbose = FALSE, approx=FALSE)

DimPlot(NK.lumBHR.combo, reduction = "pca", raster = F, group.by = "Capture.Method")
DimPlot(NK.lumBHR.combo, reduction = "pca", raster = F, group.by = "orig.ident")

DimHeatmap(NK.lumBHR.combo, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(NK.lumBHR.combo, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(NK.lumBHR.combo, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(NK.lumBHR.combo, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(NK.lumBHR.combo, dims = 45:55, cells = 500, balanced = TRUE)
DimHeatmap(NK.lumBHR.combo, dims = 55:65, cells = 500, balanced = TRUE)

NK.lumBHR.combo <- FindNeighbors(NK.lumBHR.combo, reduction = "pca", dims = 1:30)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
NK.lumBHR.combo <- FindClusters(NK.lumBHR.combo, resolution = resolution.range)
clustree(NK.lumBHR.combo, prefix = "integrated_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")


DefaultAssay(NK.lumBHR.combo) <- "integrated"
NK.lumBHR.combo <- FindNeighbors(NK.lumBHR.combo, reduction = "pca", dims = 1:30, annoy.metric = "manhattan")
NK.lumBHR.combo <- FindClusters(NK.lumBHR.combo, resolution = 0.4)
NK.lumBHR.combo <- RunUMAP(NK.lumBHR.combo, reduction = "pca", dims = 1:30, verbose = TRUE, seed.use=123)

DimPlot(NK.lumBHR.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE) + ggtitle(label = " ")
DimPlot(NK.lumBHR.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "celltype_NKsub") 
DimPlot(NK.lumBHR.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "Capture.Method") 
DimPlot(NK.lumBHR.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "orig.ident")


NK.lumBHR.combo$Capture.Method[which(NK.lumBHR.combo$orig.ident == "Xu")] <- "Singleron GEXSCOPE Single Cell RNAseq Library Kit"
NK.lumBHR.combo$Capture.Method[which(NK.lumBHR.combo$orig.ident == "Wu")] <- "10X Genomics Single Cell 3' v2"
NK.lumBHR.combo$Capture.Method[which(NK.lumBHR.combo$orig.ident == "Qian")] <- "10X Genomics Chromium v2 5'"


DefaultAssay(NK.lumBHR.combo) <- "RNA"
NK.lumBHR.combo <- NormalizeData(NK.lumBHR.combo, assay = "RNA")


#applying signature ______________________

ReprogSig <- list(ReprogSig.final)
NK.lumBHR.combo <- AddModuleScore_UCell(NK.lumBHR.combo, features = ReprogSig, name = "ReprogSig", assay = "RNA")
FeaturePlot(object = NK.lumBHR.combo, features = "signature_1ReprogSig", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, raster = FALSE) + ggtitle(label = "Reprogrammed NK Signature (Expression > 20%, PC = 40)")

nCount75 <- quantile(NK.lumBHR.combo$signature_1ReprogSig, 0.75)#95) # calculate value in the 95th percentile)

reprog <- WhichCells(NK.lumBHR.combo, expression = signature_1ReprogSig > nCount75)
setwd(rNKcell_dir)
#saveRDS(reprog, file = "reprog_lumBHR_SCTRNAnoZ_71322.rds")
saveRDS(reprog, file = "reprog_lumBHR_FIXEDrsig_82422.rds")

Idents(NK.lumBHR.combo, cells = reprog) <- "Reprogrammed NK Cells"


# NK.TNBCbasal =========================================

#https://github.com/satijalab/seurat/issues/1883
table(Basal_TNBC$orig.ident)
DefaultAssay(Basal_TNBC) <- "RNA"
Basal_TNBC <- NormalizeData(Basal_TNBC, assay = "RNA")

table(Basal_TNBC$orig.ident)
table(Basal_TNBC$Capture.Method)


Basal_TNBC.list <- SplitObject(Basal_TNBC, split.by = "Capture.Method")
for (i in 1:length(Basal_TNBC.list)) {
  Basal_TNBC.list[[i]] <- SCTransform(Basal_TNBC.list[[i]], verbose = T, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                                                             "percent.platelet", "percent.heatshock"))
}

Basal_TNBC.features <- SelectIntegrationFeatures(object.list = Basal_TNBC.list, nfeatures = 3000)
Basal_TNBC.list <- PrepSCTIntegration(object.list = Basal_TNBC.list, anchor.features = Basal_TNBC.features)


reference.1 <-  which(names(Basal_TNBC.list) == c("10X Genomics Chromium"))
reference.2 <-  which(names(Basal_TNBC.list) == c("10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"))
reference.3 <-  which(names(Basal_TNBC.list) == c("10X Genomics Chromium v2 5'"))
reference.4 <-  which(names(Basal_TNBC.list) == c("10X Genomics Single Cell 3' v2"))

reference.list <- c(reference.1, reference.2, reference.3, reference.4)

Basal_TNBC.anchors <- FindIntegrationAnchors(object.list = Basal_TNBC.list, normalization.method = "SCT",
                                             anchor.features = Basal_TNBC.features, reference = reference.list,
                                             k.score = 19, dims = 1:19)

NK.BasalTNBC.combo <- IntegrateData(anchorset = Basal_TNBC.anchors, normalization.method = "SCT", k.weight = 19)

DefaultAssay(NK.BasalTNBC.combo) <- "integrated"
#https://github.com/satijalab/seurat/issues/1963
NK.BasalTNBC.combo <- RunPCA(NK.BasalTNBC.combo, npcs = 300, verbose = FALSE, approx=FALSE)

DimPlot(NK.BasalTNBC.combo, reduction = "pca", raster = F, group.by = "Capture.Method")
DimPlot(NK.BasalTNBC.combo, reduction = "pca", raster = F, group.by = "orig.ident")

DimHeatmap(NK.BasalTNBC.combo, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(NK.BasalTNBC.combo, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(NK.BasalTNBC.combo, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(NK.BasalTNBC.combo, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(NK.BasalTNBC.combo, dims = 45:55, cells = 500, balanced = TRUE)
DimHeatmap(NK.BasalTNBC.combo, dims = 55:65, cells = 500, balanced = TRUE)



NK.BasalTNBC.combo <- FindNeighbors(NK.BasalTNBC.combo, reduction = "pca", dims = 1:35)

resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
NK.BasalTNBC.combo <- FindClusters(NK.BasalTNBC.combo, resolution = resolution.range)
clustree(NK.BasalTNBC.combo, prefix = "integrated_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")

DefaultAssay(NK.BasalTNBC.combo) <- "integrated"

NK.BasalTNBC.combo <- FindNeighbors(NK.BasalTNBC.combo, reduction = "pca", dims = 1:35, annoy.metric = "manhattan")
NK.BasalTNBC.combo <- FindClusters(NK.BasalTNBC.combo, resolution = 0.4)
NK.BasalTNBC.combo <- RunUMAP(NK.BasalTNBC.combo, reduction = "pca", dims = 1:35, verbose = TRUE, seed.use=123)

DimPlot(NK.BasalTNBC.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE) + ggtitle(label = " ")
DimPlot(NK.BasalTNBC.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "celltype_NKsub") 
DimPlot(NK.BasalTNBC.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "Capture.Method") 
DimPlot(NK.BasalTNBC.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "orig.ident")


DefaultAssay(NK.BasalTNBC.combo) <- "RNA"
NK.BasalTNBC.combo <- NormalizeData(NK.BasalTNBC.combo, assay = "RNA")


#applying signature ______________________

ReprogSig <- list(ReprogSig.final)
NK.BasalTNBC.combo <- AddModuleScore_UCell(NK.BasalTNBC.combo, features = ReprogSig, name = "ReprogSig", assay = "RNA")
FeaturePlot(object = NK.BasalTNBC.combo, features = "signature_1ReprogSig", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, raster = FALSE) + ggtitle(label = "Reprogrammed NK Signature (Expression > 20%, PC = 40)")

nCount75 <- quantile(NK.lumBHR.combo$signature_1ReprogSig, 0.75)#95) # calculate value in the 95th percentile)

reprog <- WhichCells(NK.BasalTNBC.combo, expression = signature_1ReprogSig > nCount75)
setwd(rNKcell_dir)
#saveRDS(reprog, file = "reprog_basalTNBC_SCTRNAnoZ_71322.rds")
saveRDS(reprog, file = "reprog_basalTNBC_FIXEDrsig_82422.rds")

Idents(NK.BasalTNBC.combo, cells = reprog) <- "Reprogrammed NK Cells"


# final reprogrammed labeling =====================================================

setwd(rNKcell_dir)
setwd(rNKcell_dir)
lumAHR_reprog <- readRDS(file = "reprog_lumAHR_FIXEDrsig_82422.rds")
lumBHR_reprog <- readRDS(file = "reprog_lumBHR_FIXEDrsig_82422.rds")
BasalTNBC_reprog <- readRDS(file = "reprog_basalTNBC_FIXEDrsig_82422.rds")

#ones used in paper:
# setwd(NKsubDir)
# lumAHR_reprog <- readRDS(file = "reprog_lumAHR_SCTRNAnoZ_71322.rds")
# lumBHR_reprog <- readRDS(file = "reprog_lumBHR_SCTRNAnoZ_71322.rds")
# BasalTNBC_reprog <- readRDS(file = "reprog_basalTNBC_SCTRNAnoZ_71322.rds")

allreprog <- c(lumAHR_reprog,lumBHR_reprog,BasalTNBC_reprog)
allreprog_unique <- unique(allreprog)
length(allreprog_unique)

Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub
Idents(NK.all.combo, cells = allreprog_unique) <- "Reprogrammed NK Cells"

Idents(combo.reference) <- combo.reference$celltype_final
Idents(combo.reference, cells = allreprog_unique) <- "Reprogrammed NK Cells"
table(Idents(combo.reference))
combo.reference$celltype_withreprog <- Idents(combo.reference)
setwd(rNKdir)
saveRDS(combo.reference, file = "PrimObject_withCORRECTreprog_nootherchange_82422.rds")
#saveRDS(combo.reference, file = "PrimObject_withreprog_noZallgenedem_71322.rds")



NK.all.combo$celltype_withreprog <- Idents(NK.all.combo)
NK.all.combo <- RenameIdents(NK.all.combo, `0` = "Non-Reprogrammed NK Cells",
                             `1` = "Non-Reprogrammed NK Cells",
                             `2` = "Non-Reprogrammed NK Cells",
                             `3` = "Non-Reprogrammed NK Cells",
                             `4` = "Non-Reprogrammed NK Cells",
                             `5` = "Non-Reprogrammed NK Cells",
                             `6` = "Non-Reprogrammed NK Cells")
table(Idents(NK.all.combo))
NK.all.combo$celltype_withreprog_simply <- Idents(NK.all.combo)
setwd(rNKdir)
saveRDS(NK.all.combo, "NKsubset_withCORRECTreprog_nootherchange_82422.rds")
#saveRDS(NK.all.combo, "NKsubset_withreprog_nozallgenedem_71322.rds")




# rNK per clust boxplot ------

sobjlists <- FetchData(object = NK.all.combo,
                       vars = c("samples",
                                #"Patient",
                                "BC.Subtype",
                                "celltype_NKsub",
                                "signature_1ReprogSig"))



sobjlists <- sobjlists %>% dplyr::group_by(samples,
                                           #Patient,
                                           BC.Subtype,
                                           celltype_NKsub,
                                           signature_1ReprogSig) %>%
  dplyr::summarise(Nb = n()) %>%
  dplyr::mutate(C = sum(Nb)) %>%
  dplyr::mutate(percent = Nb/C*100)


head(sobjlists)

# HR_sobj <- sobjlists[sobjlists$BC.Subtype == "HR+",]
# HER2_sobj <- sobjlists[sobjlists$BC.Subtype == "HER2+",]
# TNBC_sobj <- sobjlists[sobjlists$BC.Subtype == "TNBC",]
# summary(HR_sobj$signature_1ReprogSig)
# summary(HER2_sobj$signature_1ReprogSig)
# summary(TNBC_sobj$signature_1ReprogSig)

library(ggpubr)

colnames(sobjlists)

stat.test <- compare_means(signature_1ReprogSig~celltype_NKsub, sobjlists, 
                           method = "kruskal.test", 
                           p.adjust.method = "bonferroni")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2")
write.csv(stat.test, "correctreprog_Fig2D_kruskal_reprogsig_vNKsub_bonferronicorrec_111822.csv")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/correctrNK_preprint")
write.csv(stat.test, "preprint_reprog_Fig2D_kruskal_reprogsig_vNKsub_bonferronicorrec_112322.csv")


# setwd(rNKdir)
# write.csv(stat.test, "Fig2D_kruskal_reprogsig_vNKsub_bonferronicorrec_72822.csv")

library(rstatix)

colnames(sobjlists)

sobjlists_grouped <- sobjlists %>%
  group_by(celltype_NKsub) 

sobjlists_grouped <- as.data.frame(sobjlists_grouped)


#perform Dunn's Test with Bonferroni correction for p-values
posthoctest <- dunn_test(formula = signature_1ReprogSig~celltype_NKsub, 
                         data = sobjlists_grouped,
                         p.adjust.method="bonferroni")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2")
write.csv(posthoctest, "correct_rNK_perNKsub_posthoc_bonferronipost_111822.csv")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/correctrNK_preprint")
write.csv(posthoctest, "preprint_correct_rNK_perNKsub_posthoc_bonferronipost_112322.csv")

# setwd(rNKdir)
# write.csv(posthoctest, "rNK_perNKsub_posthoc_bonferronipost_102422.csv") ##preprint
#posthoctest <- read.csv("rNK_perNKsub_posthoc_nocorreconpost_72822.csv")

my_comparisons <- list( c("1", "0"),
                        c("1", "2"),
                        c("1", "3"),
                        c("1", "4"),
                        c("1", "5"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("****", "***", "**", "*", "ns"))

p <- ggplot(sobjlists, aes(x = celltype_NKsub, y = signature_1ReprogSig, fill = celltype_NKsub)) + 
  #geom_violin(trim = FALSE) + 
  geom_boxplot(outlier.shape = NA) + 
  ylab("ReNK Signature (UScore)") + 
  xlab("NK Cell Subset") +
  ylim(0,0.22) +
  theme_bw() + 
  scale_fill_manual(values=c(turbo(6))) +
  # scale_fill_manual(values=c("#B53E8E", "#56B1B7", "#E6B650", "#53A84C", "#6E8BE9", 
  #                            "#D97D72", "#65E65A")) +
  #theme(axis.text.x = element_text(angle = 270, hjust = 1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  NoLegend() + 
  
  stat_compare_means(comparisons = my_comparisons,
                     method="wilcox.test", label="..p.adj..", color="black",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")))

library(viridis)
library(viridisLite)
# pdf("reprogsig_perclust_violin_7122.pdf", width = 12, height = 4.05)
# pdf("reprogsig_perclust_violin_7522.pdf", width = 12, height = 4.05)
pdf("reprogsig_perclust_violin_71322.pdf", width = 14, height = 8.05) ##preprint
pdf("preprint_reprogsig_perclust_violin_112322.pdf", width = 14, height = 8.05) ##correction
pdf("correct_reprogsig_perclust_violin_111822.pdf", width = 14, height = 8.05)
p + theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 25)) +
  theme(axis.text = element_text(size = 20))
dev.off()


# ================================================================== ======
#2E rNK MA plot  ========================
# ================================================================== ======

# load NKsubset =================================

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2")
NK.all.combo <- readRDS("NKsubset_withCORRECTreprog_nootherchange_82422.rds") #newnewnew
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")

setwd("/project/InternalMedicine/Chan_lab/shared/")
combo.inferCNV <- readRDS("nichenetobj_111122.rds")


head(combo.inferCNV@meta.data)

unique(colnames(combo.inferCNV@meta.data))
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"totalCNV", drop = F], "totalCNV")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"corCNV", drop = F], "corCNV")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"cancer", drop = F], "cancer")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"ESR1", drop = F], "ESR1")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"ERBB2", drop = F], "ERBB2")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"PIK3CA", drop = F], "PIK3CA")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"NTRK", drop = F], "NTRK")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"CD274", drop = F], "CD274")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"ERBB3", drop = F], "ERBB3")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"EGFR", drop = F], "EGFR")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"FGFR", drop = F], "FGFR")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"TACSTD2", drop = F], "TACSTD2")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"CDK", drop = F], "CDK")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"AR", drop = F], "AR")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"NECTIN2", drop = F], "NECTIN2")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"LAG3", drop = F], "LAG3")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"raw_GM1", drop = F], "raw_GE1")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"raw_GM2", drop = F], "raw_GE2")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"raw_GM3", drop = F], "raw_GE3")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"raw_GM4", drop = F], "raw_GE4")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"raw_GM5", drop = F], "raw_GE5")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"raw_GM6", drop = F], "raw_GE6")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"raw_GM7", drop = F], "raw_GE7")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"raw_GM8", drop = F], "raw_GE8")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"raw_GM9", drop = F], "raw_GE9")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"raw_GM10", drop = F], "raw_GE10")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"maxGM", drop = F], "maxGE")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"maxZscore", drop = F], "maxZscore")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM1", drop = F], "GE1")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM2", drop = F], "GE2")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM3", drop = F], "GE3")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM4", drop = F], "GE4")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM5", drop = F], "GE5")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM6", drop = F], "GE6")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM7", drop = F], "GE7")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM8", drop = F], "GE8")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM9", drop = F], "GE9")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM10", drop = F], "GE10")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM1_idents", drop = F], "GE1_idents")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM2_idents", drop = F], "GE2_idents")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM3_idents", drop = F], "GE3_idents")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM4_idents", drop = F], "GE4_idents")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM5_idents", drop = F], "GE5_idents")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM6_idents", drop = F], "GE6_idents")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM7_idents", drop = F], "GE7_idents")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM8_idents", drop = F], "GE8_idents")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM9_idents", drop = F], "GE9_idents")
NK.all.combo <- AddMetaData(NK.all.combo, combo.inferCNV@meta.data[,"GM10_idents", drop = F], "GE10_idents")

head(NK.all.combo@meta.data)
table(NK.all.combo$celltype_final)
table(NK.all.combo$celltype_withreprog)
table(NK.all.combo$celltype_NKsub)

setwd("/project/InternalMedicine/Chan_lab/shared")
saveRDS(NK.all.combo, "NKsub_withGEmeta_withCORRECTreprog_111722.rds")


NK.all.combo <- readRDS("NKsub_withGEmeta_withCORRECTreprog_111722.rds")

# rNK DEGs -------

table(NK.all.combo$celltype_withreprog_simply)
Idents(NK.all.combo) <- NK.all.combo$celltype_withreprog_simply


DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")


Reprog.mark <- FindAllMarkers(NK.all.combo,
                              test.use = "MAST")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2")
write.csv(Reprog.mark, "rNKvnon_DEGs_NOTHRESH_111822.csv") ##reviewer

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/correctrNK_preprint")
write.csv(Reprog.mark, "preprint_rNKvnon_DEGs_NOTHRESH_112222.csv") ##correction

setwd(DEGdir)
write.csv(Reprog.mark, "rNKvnon_DEGs_NOTHRESH_71422.csv") ##preprint


# j.markers_DGE_filtered <- read.csv("rNKvnonDEG_ABSOLUTELYNOTHRESH.csv")
# j.markers_DGE_filtered <- DEG_Remove_mito(j.markers_DGE_filtered)
# write.csv(j.markers_DGE_filtered, "rNK_NOTHRESHnomito.csv")


# rNK ma plot -------
dfsample <- read.csv("rNKvnon_DEGs_NOTHRESH_111822.csv", header = T, row.names = 1)
head(dfsample)
dfsample <- read.csv("preprint_rNKvnon_DEGs_NOTHRESH_112222.csv", header = T, row.names = 1)
#dfsample <- read.csv("rNKvnon_DEGs_NOTHRESH_71422.csv", header = T, row.names = 1) ##preprint
dfsample <- DEG_Remove_mito(dfsample)

avgexp <- AverageExpression(object = NK.all.combo, features = unique(dfsample$gene), assays = c("RNA", "SCT"))
avgexp <- as.data.frame(avgexp[['RNA']])
avgexp <- avgexp[which(rownames(avgexp) %in% unique(dfsample$gene)),]

head(avgexp)

rNKDEGs <- dfsample[dfsample$cluster=="Reprogrammed NK Cells",]
avgexpcRNK <- avgexp[rownames(avgexp) %in% rNKDEGs$gene,"Reprogrammed NK Cells", drop = F]


data <- cbind(rNKDEGs, avgexpcRNK)
colnames(data)
data_input <- data[,c(8,2, 5, 7)]
colnames(data_input) <- c("baseMean", "log2FoldChange", "padj", "gene")


which(data_input$gene == "KLRG1")
which(data_input$gene == "TIGIT")
which(data_input$gene == "NR4A3")
which(data_input$gene == "NR4A2")


library(ggpubr)
options(ggrepel.max.overlaps = 15)#15)
#options(ggrepel.max.overlaps = 300)

p <- ggmaplot(data_input, 
              fdr = 0.05, fc = 1.5, 
              #fdr = 0.05, fc = 1.7, 
              size = 1,
              palette = c("#B31B21", "#1465AC", "darkgray"),
              top = 0,              
              #top = 100,
              label.select = c("NR4A1", "NR4A2", "NR4A3", "FOSB", #up
                               "JUN", "DNAJB1", "HSPA1B", "FOS", "DUSP1", "HSPA1A",
                               "HSPA6", "TNFAIP3",
                               "MUCL1", "CYBA"), #down ##correction
              # label.select = c("NR4A1", "NR4A2", "NR4A3", "RHOB", #up
              #                  "HSPA1A", "FOS", "DUSP1", "JUN", "DNAJB1", "HSPA1B",
              #                  "FOSB", "TNFAIP3",
              #                  "MUCL1", "CYBA"), #down ##preprint
              #label.select = c("NR4A1", "NR4A2", "NR4A3"),
              genenames = as.vector(data_input$gene),
              
              font.label = c("plain", 11),
              select.top.method = "fc",
              ggtheme = ggplot2::theme_classic()) +
  xlim(0,7) + ylim(-2.5,4) + NoLegend()

ggsave("test.pdf", plot = p, width = 4.3, height = 4)


ggsave("correctrNK_MAplot_112322.pdf", plot = p, width = 4.3, height = 4) ##correction
ggsave("rNK_MAplot_71422.pdf", plot = p, width = 4.3, height = 4) ##preprint
ggsave("rNKMAplot_NR4only.pdf", plot = p, width = 4.3, height = 4)

# ================================================================== ======
#2H rNK v nonrNK similarity  ========================
# ================================================================== ======
# load NKsubset =================================

setwd(NKsubDir)
NK.all.combo <- readRDS("NK_withNKpathways_noZallgenedem_71422.rds") #preprint
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")

setwd("/project/InternalMedicine/Chan_lab/shared")
NK.all.combo <- readRDS("NKsub_withGEmeta_withCORRECTreprog_111722.rds")
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")


# lily similarity boxplot =================

setwd(NKsubDir) ##preprint
setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/similarity")

NK.all.mat <- GetAssayData(NK.all.combo[rownames(NK.all.combo@assays$RNA@data) %in% ReprogSig.final, ], slot = "data", assay = "RNA")
NK.all.mat <- as.data.frame(NK.all.mat)
saveRDS(NK.all.mat, "NK.all.mat_justreproggenes_111822.rds")
saveRDS(NK.all.mat, "preprint_NK.all.mat_justreproggenes_112222.rds") ##correction
#saveRDS(NK.all.mat, "NK.all.mat_justreproggenes_71522.rds") ##preprint

# simil(GetAssayData(NK.all.combo, slot = "data", assay = "RNA"),
#       drop = NULL,
#       #housekeeping_genes,
#       "allNK_all_simil.rds",
#       "corr")

NK.reprog.mat <- GetAssayData(subset(NK.all.combo[rownames(NK.all.combo@assays$RNA@data) %in% ReprogSig.final, ], subset = celltype_withreprog_simply == "Reprogrammed NK Cells"), slot = "data", assay = "RNA")
NK.reprog.mat <- as.data.frame(NK.reprog.mat)
saveRDS(NK.reprog.mat, "NK.reprog.mat_justreproggenes_111822.rds")
saveRDS(NK.reprog.mat, "preprint_NK.reprog.mat_justreproggenes_112222.rds") ##correction
#saveRDS(NK.reprog.mat, "NK.reprog.mat_justreproggenes_71522.rds") ##preprint

# simil(GetAssayData(subset(NK.all.combo[rownames(NK.all.combo@assays$SCT@scale.data), ], subset = celltype_withreprog_simply == "Reprogrammed NK Cells"), slot = "data", assay = "RNA"),
#       #housekeeping_genes,
#       drop = NULL,
#       "reprogNK_all_simil.rds",
#       "corr")

NK.NONreprog.mat <- GetAssayData(subset(NK.all.combo[rownames(NK.all.combo@assays$RNA@data) %in% ReprogSig.final, ], subset = celltype_withreprog_simply == "Non-Reprogrammed NK Cells"), slot = "data", assay = "RNA")
NK.NONreprog.mat <- as.data.frame(NK.NONreprog.mat)
saveRDS(NK.NONreprog.mat, "NK.NONreprog.mat_justreproggenes_111822.rds")
saveRDS(NK.NONreprog.mat, "preprint_NK.NONreprog.mat_justreproggenes_111822.rds") ##correction
#saveRDS(NK.NONreprog.mat, "NK.NONreprog.mat_justreproggenes_71522.rds") ##preprint


simil(GetAssayData(subset(NK.all.combo, subset = celltype_withreprog_simply == "Non-Reprogrammed NK Cells"), slot = "data", assay = "RNA"),
      #housekeeping_genes,
      drop = NULL,
      "nonreprogNK_all_simil.rds",
      "corr")

# get list of files for similarity matrices
setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/similarity")
files <- list.files(pattern = "_all_simil_corr_111822.rds")
setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/correctrNK_preprint")
files <- list.files(pattern = "_all_simil_corr_preprint_112222.rds") ##correction
#files <- list.files(pattern = "_all_simil_corr_71522.rds")
files <- sort(files)

# read in all Jaccard similarity scores into dataframe
list <- data.frame()
# for (i in files) {
#   list <- rbind(list,
#                 data.frame(name = strsplit(i, ".rds")[[1]],
#                            value = as.vector(readRDS(i)[upper.tri(readRDS(i), diag = F)])))
# }

##newnew
list <- rbind(list,
              data.frame(name = "Reprogrammed NK vs. \nReprogrammed NK",
                         value = as.vector(readRDS("reprogNK_all_simil_corr_111822.rds")[upper.tri(readRDS("reprogNK_all_simil_corr_111822.rds"), diag = F)])))

##correction
list <- rbind(list,
              data.frame(name = "Reprogrammed NK vs. \nReprogrammed NK",
                         value = as.vector(readRDS("reprogNK_all_simil_corr_preprint_112222.rds")[upper.tri(readRDS("reprogNK_all_simil_corr_preprint_112222.rds"), diag = F)])))

##preprint (next three lines)
# list <- rbind(list,
#               data.frame(name = "Reprogrammed NK vs. \nReprogrammed NK",
#                          value = as.vector(readRDS("reprogNK_all_simil_corr_71522.rds")[upper.tri(readRDS("reprogNK_all_simil_corr_71522.rds"), diag = F)])))


# list <- rbind(list,
#               data.frame(name = "Non-reprogrammed NK vs. \nNon-reprogrammed NK",
#                          value = as.vector(readRDS("nonreprogNK_all_simil_corr.rds")[upper.tri(readRDS("nonreprogNK_all_simil_corr.rds"), diag = F)])))
allNK <- readRDS("allNK_all_simil_corr_111822.rds")
allNK <- readRDS("allNK_all_simil_corr_preprint_112222.rds") ##correction
#allNK <- readRDS("allNK_all_simil_corr_71522.rds") ##preprint
NK_idents <- subset(NK.all.combo,
                    subset = celltype_withreprog_simply == "Reprogrammed NK Cells")
NK_idents <- colnames(GetAssayData(NK_idents))
allNK <- allNK[which(rownames(allNK) %in% NK_idents),
               which(!(colnames(allNK) %in% NK_idents))]
list <- rbind(list,
              data.frame(name = "Reprogrammed NK vs. \nNon-Reprogrammed NK",
                         value = as.vector(allNK)))

list.reprog <- list[list$name == "Reprogrammed NK vs. \nReprogrammed NK",]
summary(list.reprog$value)
list.Nonreprog <- list[list$name == "Reprogrammed NK vs. \nNon-Reprogrammed NK",]
summary(list.Nonreprog$value)

# plot violin plot of all Jaccard scores
my_comparisons <- list( c("Reprogrammed NK vs. \nReprogrammed NK", "Reprogrammed NK vs. \nNon-Reprogrammed NK"))#,
#c("Reprogrammed NK vs. \nReprogrammed NK", "Non-reprogrammed NK vs. \nNon-reprogrammed NK"),
#c("Non-reprogrammed NK vs. \nNon-reprogrammed NK", "Reprogrammed NK vs. \nNon-reprogrammed NK"))

test <- compare_means(value~name, list, method = "wilcox.test", 
                      p.adjust.method = "BH")

test

#load library


library(ggpubr)
p <- ggplot(list, aes(x = name, y = value, fill = name)) +
  #geom_violin(trim = FALSE) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Pearson Correlation") +
  xlab(" ") +
  ylim(0,1.2) +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 270, hjust = 1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  NoLegend() +
  stat_compare_means(comparisons = my_comparisons,
                     method="wilcox.test", label="..p.adj..", color="black")
#stat_pvalue_manual(test, label = "p.adj")

# pdf("ReprogvallNK_similbox_7122.pdf", width = 5.5, height = 5)
# pdf("ReprogvallNK_similbox_7522.pdf", width = 5.5, height = 5)
pdf("ReprogvallNK_similbox_71522.pdf", width = 5.5, height = 5) ##preprint
pdf("correct_ReprogvallNK_similbox_112222.pdf", width = 5.5, height = 5)
p + theme(axis.line = element_line(colour = 'black', size = 1)) +
  theme(axis.ticks = element_line(colour = "black", size = 1)) +
  theme(text = element_text(size = 13)) +
  theme(axis.text = element_text(size = 13))
dev.off()
ggsave("reprogNK_all_simil_violin.pdf", plot = p, width = 4, height = 4)




# ================================================================== ======
#2C NKsub Proportion barchart ========================
# ================================================================== ======


# reduce samples to just those >10 cells (Addition should've done +++++++++++) -------

Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub

table(Idents(NK.all.combo))

sort(table(NK.all.combo$samples))

NK.reduced <- NK.all.combo

##removing samples with < 10 cells
NK.reduced <- subset(NK.reduced, subset = samples == "0029_9C_HR", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0040_HR", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "P2prim", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT039_P11", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT039_P121", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT081_P5", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC_14", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC_6", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "CID3946", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT089_P1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "P1prim", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT058_P1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0031_Her2", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC1_TUMOR_3", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT081_P3", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0135_TNBC", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT126_P3", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0069_Her2", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC_9", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC3_TUMOR_1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "P3prim", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT039_P10", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0043_HR", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0163_HR", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC2_TUMOR_4", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "Patient 1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC1_TUMOR_1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0151_HR", invert = T)

sort(table(NK.reduced$samples))


# Barplot of ReprogNK  by sample ------------

table(Idents(NK.all.combo))

#what did
sobjlists <- FetchData(object = NK.all.combo,
                       vars = c("samples",
                                "Patient",
                                "BC.Subtype",
                                "celltype_NKsub",
                                "celltype_withreprog_simply"))

#what should have done :(
sobjlists <- FetchData(object = NK.reduced,
                       vars = c("samples",
                                "Patient",
                                "BC.Subtype",
                                "celltype_NKsub",
                                "celltype_withreprog_simply"))


sobjlists <- sobjlists %>% dplyr::group_by(samples,
                                           Patient,
                                           BC.Subtype,
                                           celltype_NKsub) %>%
  dplyr::summarise(Nb = n()) %>%
  dplyr::mutate(C = sum(Nb)) %>%
  dplyr::mutate(percent = Nb/C*100)



sobjlists <- sobjlists[order(match(sobjlists$celltype_NKsub, c("0", "1", "2", "3",
                                                               "4", "5"))),]
#sobjlists <- sobjlists[order(sobjlists$percent_reprog, decreasing = T),] # order by HER2

sobjlists <- sobjlists[order(sobjlists$BC.Subtype),] # order by BC Subtype
#sobjlists$samples <- factor(sobjlists$samples, levels = unique(sobjlists$samples))
sobjlists$celltype_NKsub <- factor(sobjlists$celltype_NKsub, levels = c("0", "1", "2", "3",
                                                                        "4", "5"))

library(viridis)
library(ggh4x)
bp <- ggplot(sobjlists,
             aes(x = samples,
                 y = percent,
                 group = as.factor(BC.Subtype),
                 fill = celltype_NKsub)) +
  scale_fill_manual(name = "NK Subtype",
                    values = turbo(6),
                    unique(sobjlists$celltype_NKsub)) +
  
  # scale_fill_manual(name = "Reprogrammed Status",
  #                   values = rev(c("#440154FF", "light grey")),
  #                   rev(unique(sobjlists$celltype_withreprog_simply))) +
  geom_bar(stat = "identity", width = 0.93) +
  # geom_text(aes(label = paste0(round(sobjlists$percent_reprog), "%"),
  #               y = 120,
  #               fill = NULL),
  #           angle = 90,
  #           colour = "#666666",
  #           size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 105),#135),
                     breaks = seq(0,100,25)) +
  facet_nested( ~ BC.Subtype,
                scales = "free",
                space = "free",
                switch = "x"
  ) +
  theme(legend.text=element_text(size=20)) +
  theme(legend.title=element_text(size=20)) +
  ylab("% Cells") +
  xlab("Samples (IHC Subtype)") +
  theme(strip.text.x = element_text(size = 20, face = "italic"),
        strip.background = element_rect(colour = "#FFFFFF",
                                        size = 1.5,
                                        fill = "#EEEEEE"),
        panel.spacing.x = unit(-0.1, "lines"))


#pdf("percentNKsubpersamplebar_71322.pdf", width = 18.4, height = 3.4) ##preprint
setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2")
pdf("percentNKsubpersamplebar_GREATER10SAMPLES_111822.pdf", width = 18.4, height = 3.4)
pdf("preprint_percentNKsubpersamplebar_GREATER10SAMPLES_112222.pdf", width = 18.4, height = 3.4)
bp
dev.off()




# ================================================================== ======
#2F rNK UScore per BCsub ========================
# ================================================================== ======
# load NKsubset =================================

setwd(NKsubDir)
NK.all.combo <- readRDS("NK_withNKpathways_noZallgenedem_71422.rds") #preprint
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")

setwd("/project/InternalMedicine/Chan_lab/shared")
NK.all.combo <- readRDS("NKsub_withGEmeta_withCORRECTreprog_111722.rds")
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")




# Prepare data for correlations (>10 samples only) =======================

Idents(NK.all.combo) <- NK.all.combo$celltype_withreprog_simply 

sort(table(NK.all.combo$samples))

NK.reduced <- NK.all.combo

##removing samples with < 10 cells
NK.reduced <- subset(NK.reduced, subset = samples == "0029_9C_HR", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0040_HR", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "P2prim", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT039_P11", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT039_P121", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT081_P5", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC_14", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC_6", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "CID3946", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "P1prim", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT089_P1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT058_P1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0031_Her2", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC1_TUMOR_3", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT081_P3", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0135_TNBC", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT126_P3", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0069_Her2", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC_9", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC3_TUMOR_1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "P3prim", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT039_P10", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0043_HR", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0163_HR", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC2_TUMOR_4", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "Patient 1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC1_TUMOR_1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0151_HR", invert = T)


sort(table(NK.reduced$samples))



colnames(NK.reduced@meta.data)

Sig.DF <- NK.reduced@meta.data[,c(4,101, 87, 102:289)] #reprog


#means of UScores
data <- Sig.DF %>% dplyr::group_by(samples) %>%
  dplyr::summarise(across(everything(), mean))



data <- as.data.frame(data)
head(data)

sort(table(NK.all.combo$samples))

sobjlists <- FetchData(object = NK.reduced,
                       vars = c("samples",
                                "Patient",
                                "BC.Subtype",
                                "celltype_withreprog_simply",
                                "Age",
                                "Stage",
                                "Grade",
                                "Tumor.Size",
                                "nodal_involvement",
                                "TNM.Classification",
                                "Ki67"))



sobjlists <- sobjlists %>% dplyr::group_by(samples,
                                           Patient,
                                           Age,
                                           Stage,
                                           Grade,
                                           Tumor.Size,
                                           nodal_involvement,
                                           Ki67,
                                           BC.Subtype,
                                           TNM.Classification,
                                           celltype_withreprog_simply) %>%
  dplyr::summarise(Nb = n()) %>%
  dplyr::mutate(C = sum(Nb)) %>%
  dplyr::mutate(percent = Nb/C*100)



sobjlists$percent_reprog <- NA
for (i in unique(sobjlists$samples)) {
  if(identical(which((sobjlists$samples == i)
                     & (sobjlists$celltype_withreprog_simply == "Reprogrammed NK Cells")), integer(0)) == T) {
    sobjlists$percent_reprog[which(sobjlists$samples == i)] <- 0
  }
  else {
    sobjlists$percent_reprog[which(sobjlists$samples == i)] <- sobjlists$percent[which((sobjlists$samples == i)
                                                                                       & (sobjlists$celltype_withreprog_simply == "Reprogrammed NK Cells"))]
  }
}

which(sobjlists$samples == "0308_Her2")
sobjlists[9,8]

sobjlists <- sobjlists[order(match(sobjlists$celltype_withreprog_simply, c("Reprogrammed NK Cells", "Non-Reprogrammed NK Cells"))),]
sobjlists <- sobjlists[order(sobjlists$percent_reprog, decreasing = T),] # order by HER2

sobjlists <- sobjlists[order(sobjlists$BC.Subtype),] # order by BC Subtype
sobjlists$samples <- factor(sobjlists$samples, levels = unique(sobjlists$samples))


sobjlists$TsizeStatus <- sobjlists$TNM.Classification
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "PpT1c, pN1, Mi, Stage IIA")] <- "T1"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT(m)2, N2a")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT1c, N1a, Mx")] <- "T1"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT1c, pN0")] <- "T1"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT1cN0M0")] <- "T1"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT1N0M0")] <- "T1"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N0, Mx")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N1(sn), Mx")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N1a")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N1a, Mx")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N2a")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N2a (Stage IIIA)")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, pN0 (i+),Stage IIB")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2N0M0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2N1aM0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2N1M0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2NxM0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT3, N1, Mx")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT3, pN0 (i+)")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT3, pN3, pMx, Stage IIIA")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT3N0M0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT4b, Nx")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT3, N1a, Mx")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "T2N1M0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "T3N2M0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "T2N2M0")] <- "T2-4"
table(sobjlists$TsizeStatus)

sobjlists$N <- NA
sobjlists$N[which(grepl("N0", sobjlists$TNM.Classification))] <- "N0"
sobjlists$N[which(grepl("N1", sobjlists$TNM.Classification))] <- "N1"
sobjlists$N[which(grepl("N2", sobjlists$TNM.Classification))] <- "N2"
sobjlists$N[which(grepl("N3", sobjlists$TNM.Classification))] <- "N3"
sobjlists$N[which(grepl("no", sobjlists$nodal_involvement))] <- "N0"
sobjlists$N[which(grepl("1", sobjlists$nodal_involvement))] <- "N1"
sobjlists$N[which(sobjlists$N %in% c("N1", "N2", "N3"))] <- "N1-3"


head(sobjlists)

age.split.df <- sobjlists
age.split.df$Age <- as.numeric(age.split.df$Age)
age.split.df <- age.split.df[which(!is.na(age.split.df$Age)),]
age.split.df$Age_split <- ifelse(age.split.df$Age >45, ">45", "<45")
colnames(age.split.df)

sobjlists <- merge(sobjlists, age.split.df[,c(1,18)], by = "samples")
sobjlists <- merge(sobjlists, data, by = "samples")

head(sobjlists)
unique(sobjlists$samples)

# Boxplot ================================

colnames(sobjlists)


colnames(sobjlists) <- sub("signature_1", "", colnames(sobjlists))

library(ggpubr)
library(ggplot2)
library(rstatix)

my_comparisons <- list(c("HR+", "HER2+"), c("HR+", "TNBC"),
                       c("HER2+", "TNBC"))

head(sobjlists)

stat.test<- compare_means(
  ReprogSig ~ BC.Subtype, data = sobjlists, method = "kruskal.test",
  p.adjust.method = "bonferroni"
)
write.csv(stat.test, "reprogsigvBCsub_kruskal_bonferroni_111822.csv") ##newnew
write.csv(stat.test, "preprint_reprogsigvBCsub_kruskal_bonferroni_112222.csv") ##correction

stat.test

sobjlists_grouped <- sobjlists %>%
  group_by(BC.Subtype) 

sobjlists_grouped <- as.data.frame(sobjlists_grouped)


#perform Dunn's Test with Bonferroni correction for p-values
posthoctest <- dunn_test(formula = ReprogSig~BC.Subtype, 
                         data = sobjlists_grouped,
                         p.adjust.method="none")

setwd(NKsubDir)
write.csv(posthoctest, "preprint_reprogsigvBCsub_nocorrectposthoc_112222.csv") ##correction
write.csv(posthoctest, "reprogsigvBCsub_nocorrectposthoc_72822.csv") ##preprint

my_comparisons <- list(c("HR+", "HER2+"), c("HR+", "TNBC"),
                       c("HER2+", "TNBC"))
test <- compare_means(value~name, list, method = "wilcox.test", p.adjust.method = "BH")

test <- stat_compare_means(sobjlists, comparisons = my_comparisons,
                           method="wilcox.test", label="p.format", color="black")

library(ggpubr)
#p <- ggplot(unique(sobjlists[which(!is.na(sobjlists.subset$BC.Subtype)),c(1,9,15)]), 
p <- ggplot(unique(sobjlists[which(!is.na(sobjlists$BC.Subtype)),c(1,9,20)]), 
            aes(x = BC.Subtype, y = ReprogSig, fill = BC.Subtype)) +
  geom_boxplot(width = 0.5) +
  #ylab("% Reprogrammed NK Cells") +
  ylab("Reprogrammed NK Signature (UScore)") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  NoLegend() +
  stat_compare_means(comparisons = my_comparisons,
                     method="wilcox.test", label="..p.adj..", color="black") +
  #stat_compare_means(label.y = 0.18) +
  scale_y_continuous(breaks = c(0, 0.11, 0.18),
                     limits = c(0, 0.18)) #+ 
#stat_pvalue_manual(test, label = "p.adj")
# scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
#                    limits = c(0, 110))

#p$layers[[2]]$aes_params$textsize <- 5
# pdf("percentreprog_BCsubtype_62322.pdf", width = 3, height = 2.5)
# pdf("reprogNKsig_BCsubtype_62322.pdf", width = 3, height = 2.5)
# pdf("reprogNKsig_bcsub_7522_bigger.pdf", width = 5.5, height = 5)
pdf("reprogNKsig_bcsub_71522_bigger.pdf", width = 5.5, height = 5) ##preprint

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2")
pdf("correctreprogNKsig_bcsub_111822.pdf", width = 5.5, height = 5)
pdf("preprint_correctreprogNKsig_bcsub_112222.pdf", width = 5.5, height = 5)

p +
  theme(axis.line = element_line(colour = 'black', size = 1)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1)) +
  theme(text = element_text(size = 16)) +
  theme(axis.text = element_text(size = 16))
dev.off()

# ================================================================== ======
#2I %rNK v age ========================
# ================================================================== ======
# load NKsubset =================================

setwd(NKsubDir)
NK.all.combo <- readRDS("NK_withNKpathways_noZallgenedem_71422.rds") #preprint
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")

setwd("/project/InternalMedicine/Chan_lab/shared")
NK.all.combo <- readRDS("NKsub_withGEmeta_withCORRECTreprog_111722.rds")
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")


# Prepare data for correlations (>10 samples only) =======================

Idents(NK.all.combo) <- NK.all.combo$celltype_withreprog_simply 

sort(table(NK.all.combo$samples))

NK.reduced <- NK.all.combo

##removing samples with < 10 cells
NK.reduced <- subset(NK.reduced, subset = samples == "0029_9C_HR", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0040_HR", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "P2prim", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT039_P11", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT039_P121", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT081_P5", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC_14", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC_6", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "CID3946", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "P1prim", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT089_P1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT058_P1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0031_Her2", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC1_TUMOR_3", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT081_P3", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0135_TNBC", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT126_P3", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0069_Her2", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC_9", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC3_TUMOR_1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "P3prim", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "PT039_P10", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0043_HR", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0163_HR", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC2_TUMOR_4", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "Patient 1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "BC1_TUMOR_1", invert = T)
NK.reduced <- subset(NK.reduced, subset = samples == "0151_HR", invert = T)


sort(table(NK.reduced$samples))



colnames(NK.reduced@meta.data)

Sig.DF <- NK.reduced@meta.data[,c(4,101, 87, 102:289)] #reprog


#means of UScores
data <- Sig.DF %>% dplyr::group_by(samples) %>%
  dplyr::summarise(across(everything(), mean))



data <- as.data.frame(data)
head(data)

sort(table(NK.all.combo$samples))

sobjlists <- FetchData(object = NK.reduced,
                       vars = c("samples",
                                "Patient",
                                "BC.Subtype",
                                "celltype_withreprog_simply",
                                "Age",
                                "Stage",
                                "Grade",
                                "Tumor.Size",
                                "nodal_involvement",
                                "TNM.Classification",
                                "Ki67"))



sobjlists <- sobjlists %>% dplyr::group_by(samples,
                                           Patient,
                                           Age,
                                           Stage,
                                           Grade,
                                           Tumor.Size,
                                           nodal_involvement,
                                           Ki67,
                                           BC.Subtype,
                                           TNM.Classification,
                                           celltype_withreprog_simply) %>%
  dplyr::summarise(Nb = n()) %>%
  dplyr::mutate(C = sum(Nb)) %>%
  dplyr::mutate(percent = Nb/C*100)



sobjlists$percent_reprog <- NA
for (i in unique(sobjlists$samples)) {
  if(identical(which((sobjlists$samples == i)
                     & (sobjlists$celltype_withreprog_simply == "Reprogrammed NK Cells")), integer(0)) == T) {
    sobjlists$percent_reprog[which(sobjlists$samples == i)] <- 0
  }
  else {
    sobjlists$percent_reprog[which(sobjlists$samples == i)] <- sobjlists$percent[which((sobjlists$samples == i)
                                                                                       & (sobjlists$celltype_withreprog_simply == "Reprogrammed NK Cells"))]
  }
}

which(sobjlists$samples == "0308_Her2")
sobjlists[9,8]

sobjlists <- sobjlists[order(match(sobjlists$celltype_withreprog_simply, c("Reprogrammed NK Cells", "Non-Reprogrammed NK Cells"))),]
sobjlists <- sobjlists[order(sobjlists$percent_reprog, decreasing = T),] # order by HER2

sobjlists <- sobjlists[order(sobjlists$BC.Subtype),] # order by BC Subtype
sobjlists$samples <- factor(sobjlists$samples, levels = unique(sobjlists$samples))


sobjlists$TsizeStatus <- sobjlists$TNM.Classification
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "PpT1c, pN1, Mi, Stage IIA")] <- "T1"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT(m)2, N2a")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT1c, N1a, Mx")] <- "T1"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT1c, pN0")] <- "T1"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT1cN0M0")] <- "T1"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT1N0M0")] <- "T1"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N0, Mx")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N1(sn), Mx")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N1a")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N1a, Mx")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N2a")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, N2a (Stage IIIA)")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2, pN0 (i+),Stage IIB")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2N0M0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2N1aM0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2N1M0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT2NxM0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT3, N1, Mx")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT3, pN0 (i+)")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT3, pN3, pMx, Stage IIIA")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT3N0M0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT4b, Nx")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "pT3, N1a, Mx")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "T2N1M0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "T3N2M0")] <- "T2-4"
sobjlists$TsizeStatus[which(sobjlists$TsizeStatus == "T2N2M0")] <- "T2-4"
table(sobjlists$TsizeStatus)

sobjlists$N <- NA
sobjlists$N[which(grepl("N0", sobjlists$TNM.Classification))] <- "N0"
sobjlists$N[which(grepl("N1", sobjlists$TNM.Classification))] <- "N1"
sobjlists$N[which(grepl("N2", sobjlists$TNM.Classification))] <- "N2"
sobjlists$N[which(grepl("N3", sobjlists$TNM.Classification))] <- "N3"
sobjlists$N[which(grepl("no", sobjlists$nodal_involvement))] <- "N0"
sobjlists$N[which(grepl("1", sobjlists$nodal_involvement))] <- "N1"
sobjlists$N[which(sobjlists$N %in% c("N1", "N2", "N3"))] <- "N1-3"


head(sobjlists)

age.split.df <- sobjlists
age.split.df$Age <- as.numeric(age.split.df$Age)
age.split.df <- age.split.df[which(!is.na(age.split.df$Age)),]
age.split.df$Age_split <- ifelse(age.split.df$Age >45, ">45", "<45")
colnames(age.split.df)

sobjlists <- merge(sobjlists, age.split.df[,c(1,18)], by = "samples")
sobjlists <- merge(sobjlists, data, by = "samples")

head(sobjlists)
unique(sobjlists$samples)

# Lily correlation line ================================


##sobjlists$percent_neg <- 100 - sobjlists$percent_neg
# sobjlists.subset <- sobjlists[which(sobjlists$BC.Subtype == "TNBC"),]
# unique(sobjlists.subset$samples)
# 
# sobjlists.subset <- sobjlists[which(sobjlists$BC.Subtype == "HER2+"),]
# unique(sobjlists.subset$samples)
# 
# sobjlists.subset <- sobjlists[which(sobjlists$BC.Subtype == "HR+"),]
# unique(sobjlists.subset$samples)

sobjlists.subset <- sobjlists

colnames(sobjlists.subset)
age_linreg <- unique(sobjlists.subset[,c(1,3,15)])
age_linreg$Age <- as.numeric(age_linreg$Age)
age_linreg <- age_linreg[which(!is.na(age_linreg$Age)),]

library(ggpubr)
p <- ggscatter(age_linreg, x = "Age", 
               y = "percent_reprog",
               add = "reg.line",
               conf.int = TRUE,
               add.params = list(color = "blue", fill = "lightgray"),
               xlab = "Age",
               ylab = "% Reprogrammed NK Cells") + theme_classic() +
  stat_cor(method = "pearson", size = 6, mapping = aes(label = "..p.adj..")) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     limits = c(0, 120))

# pdf("ReprogvAge_TNBC_scatter.pdf", width = 3, height = 2.5)
# pdf("ReprogvAge_HER2+_scatter.pdf", width = 3, height = 2.5)
# pdf("ReprogvAge_HR+_scatter.pdf", width = 3, height = 2.5)
# pdf("ReprogvAge_ALLBCsub_scatter.pdf", width = 3, height = 2.5)
# 
# pdf("ReprogvAge_ALLBCsub_scatter_62922.pdf", width = 4, height = 3.5)
# pdf("ReprogvAge_ALLBCsub_scatter_7622.pdf", width = 4, height = 3.5)
pdf("ReprogvAge_ALLBCsub_scatter_71522.pdf", width = 4, height = 3.5) ##preprint
setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2")
pdf("correctReprogvAge_ALLBCsub_scatter_111822.pdf", width = 4, height = 3.5)
pdf("preprint_correctReprogvAge_ALLBCsub_scatter_112222.pdf", width = 4, height = 3.5)

p+ theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 16)) +
  theme(axis.text = element_text(size = 16))
dev.off()

# ================================================================== ======
#2G rNK NicheNet ========================
# ================================================================== ======
# load in Prim object ------

setwd(PrimDir)
combo.reference <- readRDS("PrimObject_withreprog_noZallgenedem_71322.rds")
DefaultAssay(combo.reference) <- "RNA"
combo.reference <- NormalizeData(combo.reference, assay = "RNA")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2")
combo.reference <- readRDS("PrimObject_withCORRECTreprog_nootherchange_82422.rds") ##correction
DefaultAssay(combo.reference) <- "RNA"
combo.reference <- NormalizeData(combo.reference, assay = "RNA")



# epiNK clustering =========================================

table(Idents(combo.reference))
epiNKsub <- subset(combo.reference, idents = c("Epithelial Cells",
                                               "Reprogrammed NK Cells",
                                               "NK Cells"))
DefaultAssay(epiNKsub) <- "RNA"
epiNKsub <- NormalizeData(epiNKsub, assay = "RNA")
table(epiNKsub$Capture.Method)
epiNKsub.list <- SplitObject(epiNKsub, split.by = "Capture.Method")


for (i in 1:length(epiNKsub.list)) {
  epiNKsub.list[[i]] <- SCTransform(epiNKsub.list[[i]], verbose = T, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                                                         "percent.platelet", "percent.heatshock"))
}

epiNK.features <- SelectIntegrationFeatures(object.list = epiNKsub.list, nfeatures = 3000)
epiNKsub.list <- PrepSCTIntegration(object.list = epiNKsub.list, anchor.features = epiNK.features)


reference.1 <-  which(names(epiNKsub.list) == c("10X Genomics Chromium"))
reference.2 <-  which(names(epiNKsub.list) == c("10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"))
reference.3 <-  which(names(epiNKsub.list) == c("10X Genomics Chromium v2 5'"))
reference.4 <-  which(names(epiNKsub.list) == c("10X Genomics Single Cell 3' v2"))

reference.list <- c(reference.1, reference.2, reference.3, reference.4)

epiNK.anchors <- FindIntegrationAnchors(object.list = epiNKsub.list, normalization.method = "SCT",
                                        anchor.features = epiNK.features, reference = reference.list)

epiNK.combo <- IntegrateData(anchorset = epiNK.anchors, normalization.method = "SCT")

DefaultAssay(epiNK.combo) <- "integrated"
epiNK.combo <- RunPCA(epiNK.combo, npcs = 100, verbose = FALSE)


DimPlot(epiNK.combo, reduction = "pca", raster = F, group.by = "Capture.Method",
        order = c("Smart-Seq2", "inDrop v2", "Singleron GEXSCOPE Single Cell RNAseq Library Kit",
                  "10X Genomics Single Cell 3' v2", "10X Genomics Chromium v2 5'",
                  "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library",
                  "10X Genomics Chromium"))
DimPlot(epiNK.combo, reduction = "pca", raster = F, group.by = "orig.ident",
        order =c("Karaayvaz", "Savas", "Wu", "Aziziimmune", "Xu",
                 "AziziT", "Qian", "Wu2021prim", "Pal_Prim"))

pdf("test.pdf")
DimHeatmap(epiNK.combo, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(epiNK.combo, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(epiNK.combo, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(epiNK.combo, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(epiNK.combo, dims = 45:55, cells = 500, balanced = TRUE)
DimHeatmap(epiNK.combo, dims = 55:65, cells = 500, balanced = TRUE)
DimHeatmap(epiNK.combo, dims = 65:75, cells = 500, balanced = TRUE)
DimHeatmap(epiNK.combo, dims = 75:85, cells = 500, balanced = TRUE)
dev.off()

epiNK.combo <- FindNeighbors(epiNK.combo, reduction = "pca", dims = 1:60)
resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
epiNK.combo <- FindClusters(epiNK.combo, resolution = resolution.range)
clustree(epiNK.combo, prefix = "integrated_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")

DefaultAssay(epiNK.combo) <- "integrated"

epiNK.combo <- FindNeighbors(epiNK.combo, reduction = "pca", dims = 1:30)
epiNK.combo <- FindClusters(epiNK.combo, resolution = 0.7)
epiNK.combo <- RunUMAP(epiNK.combo, reduction = "pca", dims = 1:30, verbose = TRUE, seed.use=123)

DimPlot(epiNK.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE) 
DimPlot(epiNK.combo, reduction = "umap", label = TRUE, repel = TRUE, group.by = "celltype_withreprog", raster = FALSE) +theme(legend.position = "None")
DimPlot(epiNK.combo, reduction = "umap", label = TRUE, repel = TRUE, group.by = "celltype_final", raster = FALSE) +theme(legend.position = "None")
DimPlot(epiNK.combo, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident", raster = FALSE) +theme(legend.position = "None")


DefaultAssay(epiNK.combo) <- "RNA"
epiNK.combo <- NormalizeData(epiNK.combo, assay = "RNA")

NK.all.combo <- readRDS("NKsubset_withCORRECTreprog_nootherchange_82422.rds") ##correction
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")

Idents(epiNK.combo) <- epiNK.combo$celltype_final
table(Idents(epiNK.combo))
Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub
Idents(NK.all.combo) <- NK.all.combo$celltype_withreprog

c0 <- WhichCells(NK.all.combo, idents = "0")
c1 <- WhichCells(NK.all.combo, idents = "1")
c2 <- WhichCells(NK.all.combo, idents = "2")
c3 <- WhichCells(NK.all.combo, idents = "3")
c4 <- WhichCells(NK.all.combo, idents = "4")
c5 <- WhichCells(NK.all.combo, idents = "5")
reprog <- WhichCells(NK.all.combo, idents = "Reprogrammed NK Cells")

Idents(epiNK.combo, cells = c0) <- "0"
Idents(epiNK.combo, cells = c1) <- "1"
Idents(epiNK.combo, cells = c2) <- "2"
Idents(epiNK.combo, cells = c3) <- "3"
Idents(epiNK.combo, cells = c4) <- "4"
Idents(epiNK.combo, cells = c5) <- "5"
Idents(epiNK.combo, cells = reprog) <- "Reprogrammed NK Cells"

table(Idents(epiNK.combo))
epiNK.combo$celltype_withreprog <- Idents(epiNK.combo)
epiNK.combo$celltype_NKsub <- Idents(epiNK.combo)

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/correctrNK_preprint")
saveRDS(epiNK.combo, "epiNK_withCORRECTreprog_nootherchange_112322.rds") ##correction
saveRDS(epiNK.combo, "epiNK_withreprog_noZallgenedem_71322.rds") ##preprint
# nichenet models =========================

library(nichenetr)
library(tidyverse)
library(circlize)

# Sys.setenv(http_proxy="http://proxy.swmed.edu:3128")
# Sys.setenv(https_proxy="http://proxy.swmed.edu:3128")
# Sys.getenv("http_proxy")
# Sys.getenv("https_proxy")

setwd(sharedDir)
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

lr_network = readRDS("lr_network.rds")
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)


weighted_networks = readRDS("weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network

head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

# nichenet (reprog) ===========================================

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/subclustering")
epiNK <- readRDS("CancerEpicnvNK_withreprog_111722.rds")
setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/correctrNK_preprint")
epiNK <- readRDS("epiNK_withCORRECTreprog_nootherchange_112322.rds") ##correction
# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only")
# epiNK <- readRDS("epiNK_withreprog_noZallgenedem_71322.rds") ##preprint
DefaultAssay(epiNK) <- "RNA"
epiNK <- NormalizeData(epiNK, assay = "RNA")

Idents(epiNK) <- epiNK$celltype_withreprog
table(Idents(epiNK))
epiNK <- RenameIdents(epiNK, `0` = "Non-Reprogrammed NK Cells", `1` = "Non-Reprogrammed NK Cells",
                      `2` = "Non-Reprogrammed NK Cells",
                      `3` = "Non-Reprogrammed NK Cells", `4` = "Non-Reprogrammed NK Cells",
                      `5` = "Non-Reprogrammed NK Cells")
table(Idents(epiNK))

temp_nichenet1 = nichenet_seuratobj_cluster_de(
  seurat_obj = epiNK,
  assay_oi = "RNA",
  receiver_affected = "Reprogrammed NK Cells",
  receiver_reference = "Non-Reprogrammed NK Cells",
  sender = "Cancer Epithelial Cells",#"Epithelial Cells",
  geneset = "up", top_n_ligands = 50, top_n_targets = 200,
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/Nichenet")
saveRDS(temp_nichenet1, "reprog_non_CANCERepi_nichenetresults_111722.rds")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/correctrNK_preprint")
saveRDS(temp_nichenet1, "correct_reprog_non_cancerepi_nichenetresults_NOOTHERCHANGE_112322.rds") ##correction


setwd(nichenet_dir)
saveRDS(temp_nichenet1, "reprog_non_epi_nichenetresults_SCTrna_71922.rds") ##preprint

#epiNK$celltype_mainniche <- Idents(epiNK)
# epiNK <- RenameIdents(epiNK, `0` = "Non-Reprogrammed NK Cells", `1` = "Non-Reprogrammed NK Cells",
#                       `2` = "Non-Reprogrammed NK Cells",
#                       `3` = "Non-Reprogrammed NK Cells", `4` = "Non-Reprogrammed NK Cells",
#                       `5` = "Non-Reprogrammed NK Cells")
table(Idents(epiNK))

temp_nichenet1 = nichenet_seuratobj_cluster_de(
  seurat_obj = epiNK,
  assay_oi = "RNA",
  receiver_affected = "Non-Reprogrammed NK Cells",
  receiver_reference = "Reprogrammed NK Cells",
  sender = "Cancer Epithelial Cells", #Epithelial Cells",
  geneset = "up", top_n_ligands = 50, top_n_targets = 200,
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/Nichenet")
saveRDS(temp_nichenet1, "NONreprog_reprogref_CANCERepi_nichenetresults_111722.rds")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/correctrNK_preprint")
saveRDS(temp_nichenet1, "correct_NONreprog_reprogref_cancerepi_nichenetresults_NOOTHERCHANGE_112322.rds") ##correction


# setwd(nichenet_dir)
# saveRDS(temp_nichenet1, "NON_reprogref_epi_nichenetresults_SCTrna_71922.rds") ##preprint
#saveRDS(temp_nichenet1, "reprog_non_epi_nichenetresults_SCTrna_7522.rds")


# circos plot plain ==============================

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/Nichenet")
temp_nichenet1 <- readRDS(file = "reprog_non_CANCERepi_nichenetresults_111722.rds")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/Nichenet")
write.csv(temp_nichenet1$ligand_receptor_df, "CANCEREPI_Reprog_nonReF_Ligandreceptor_FULLLIST_111722.csv")
write.csv(temp_nichenet1$ligand_target_df, "CANCEREPI_Reprog_nonReF_Ligand_target_FULLLIST_111722.csv")


setwd(nichenet_dir)
temp_nichenet1 <- readRDS(file = "reprog_non_epi_nichenetresults_SCTrna_71922.rds")
temp_nichenet1 <- readRDS(file = "NON_reprogref_epi_nichenetresults_SCTrna_71922.rds")

# write.csv(temp_nichenet1$ligand_receptor_df, "Reprog_nonReF_Ligandreceptor_FULLLIST_71922.csv")
# write.csv(temp_nichenet1$ligand_target_df, "Reprog_nonReF_Ligand_target_FULLLIST_71922.csv")
write.csv(temp_nichenet1$ligand_receptor_df, "NonReprog_reprogReF_Ligandreceptor_FULLLIST_71922.csv")
write.csv(temp_nichenet1$ligand_target_df, "NonReprog_reprogReF_Ligand_target_FULLLIST_71922.csv")


setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/NichenetEMTFinal/")

lr.df.top <- temp_nichenet1$ligand_receptor_df
lr.df.top <- temp_nichenet1$ligand_target_df
lr.df.top <- lr.df.top %>% top_n(5, weight)


#generate circos plot with top 20 ligands by weight
pdf("top20LR_reprog_non_epi_8522.pdf")
pdf("top20LR_NON_reprogref_epi_71922.pdf")


pdf("test.pdf")
circos.par(canvas.ylim=c(-1.5,1.5), # edit canvas size
           track.margin = c(0.01, 0)) # adjust bottom and top margin



chordDiagram(lr.df.top,
             directional = 1,
             link.sort = TRUE,
             link.decreasing = FALSE,
             link.visible = lr.df.top$weight, # >= weight.cutoff,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.05))



circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA)
circos.clear()



dev.off()
#
# Read in NicheNet results -------

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/Nichenet")
files <- list.files(pattern = "CANCERepi_nichenetresults_111722.rds")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/correctrNK_preprint")
files <- list.files(pattern = "_NOOTHERCHANGE_112322.rds") ##correct


##preprint next 3 lines
# setwd(nichenet_dir)
# files <- list.files(pattern = "nichenetresults_SCTrna_71922.rds")
# files <- files[c(7,12)] #for rNK and non-rNK
# #files <- files[c(1:6)] #for NK subsets

#files <- "reprog_non_epi_nichenetresults_SCTrna_8522.rds"

# Nichenet object subsets (NKsub) ----------

# setwd(PrimDir)
# combo.reference <- readRDS("PrimObject_withreprog_noZallgenedem_71322.rds") ##preprint

setwd("/project/InternalMedicine/Chan_lab/shared")
combo.reference <- readRDS("PrimObj_withGEmeta_with_correctassays_withfinalreprog_111622.rds")
DefaultAssay(combo.reference) <- "RNA"
combo.reference <- NormalizeData(combo.reference, assay = "RNA")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2")
combo.reference <- readRDS("PrimObject_withCORRECTreprog_nootherchange_82422.rds") ##correction
DefaultAssay(combo.reference) <- "RNA"
combo.reference <- NormalizeData(combo.reference, assay = "RNA")

table(Idents(combo.reference))


cell <- unique(combo.reference$celltype_withreprog)
cell <- cell[which(!(cell %in%c("NK Cells", "Reprogrammed NK Cells")))]
NKsubset <- subset(combo.reference, subset = celltype_withreprog %in% c("NK Cells", "Reprogrammed NK Cells"))
DefaultAssay(NKsubset) <- "RNA"
NKsubset <- NormalizeData(NKsubset, assay = "RNA")

Idents(NKsubset) <- NKsubset$celltype_withreprog

primsubset <- subset(combo.reference, subset = celltype_withreprog %in% cell)
primsubset$newidents <- "other"
primsubset$newidents[which(primsubset$celltype_withreprog == "Cancer Epithelial Cells")] <- "epi" #Epithelial Cells")] <- "epi"
Idents(primsubset) <- primsubset$newidents

episubset <- subset(primsubset, subset = celltype_withreprog == "Cancer Epithelial Cells") #Epithelial Cells")
DefaultAssay(episubset) <- "RNA"
episubset <- NormalizeData(episubset, assay = "RNA")

episubset.save <- episubset # #episubset.save is a temp object to store all BC subtypes
episubset <- subset(episubset, subset = BC.Subtype == "TNBC")
episubset <- subset(episubset, subset = BC.Subtype == "HR+")
episubset <- subset(episubset, subset = BC.Subtype == "HER2+")
episubset <- episubset.save ##run this after finish one subtype


# NicheNet Circos plot by grouping --------------

## for NK subset nichenet vvvvvvvvvvvvvvv
# setwd(NKsubDir)
# NK.all.combo <- readRDS("NK_withNKpathways_noZallgenedem_71422.rds") #preprint

# setwd("/project/InternalMedicine/Chan_lab/shared")
# NK.all.combo <- readRDS("NKsub_withGEmeta_withCORRECTreprog_111722.rds") #newnewnew
# 
# Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub
# table(Idents(NK.all.combo))
# 
# NK0 <- WhichCells(NK.all.combo, idents = "0")
# NK1 <- WhichCells(NK.all.combo, idents = "1")
# NK2 <- WhichCells(NK.all.combo, idents = "2")
# NK3 <- WhichCells(NK.all.combo, idents = "3")
# NK4 <- WhichCells(NK.all.combo, idents = "4")
# NK5 <- WhichCells(NK.all.combo, idents = "5")
# 
# Idents(NKsubset, cells = NK0) <- "NK0"
# Idents(NKsubset, cells = NK1) <- "NK1"
# Idents(NKsubset, cells = NK2) <- "NK2"
# Idents(NKsubset, cells = NK3) <- "NK3"
# Idents(NKsubset, cells = NK4) <- "NK4"
# Idents(NKsubset, cells = NK5) <- "NK5"

table(Idents(NKsubset))

all_receptors <- data.frame()
all_targets <- data.frame()

targetlist <- list()
receptorlist <- list()

for (j in files) {
  #setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/071322/Nichenet")
  

  #setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/Nichenet") ##newnew
  setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/correctrNK_preprint") ##correction
  #setwd(nichenet_dir)  ##preprint
  temp_nichenet <- readRDS(j)
  if (j == files[1]) { cell = "NK Cells" ; print(j); print(cell)} else { cell = "Reprogrammed NK Cells" ; print(j); print(cell)} #for rNK and non-rNK
  # if (j == files[1]) { cell = "NK0" ; print(j); print(cell)} 
  # else if (j == files[2]) { cell = "NK1" ; print(j); print(cell)}
  # else if (j == files[3]) { cell = "NK2" ; print(j); print(cell)}
  # else if (j == files[4]) { cell = "NK3" ; print(j); print(cell)}
  # else if (j == files[5]) { cell = "NK4" ; print(j); print(cell)}
  # else { cell = "NK5" ; print(j); print(cell)}
  
  LR_heatmap <- t(temp_nichenet$ligand_receptor_matrix)
  LR_matrix <- data.frame()
  for (i in 1:dim(LR_heatmap)[1]) {
    for (j in 1:dim(LR_heatmap)[2]) {
      temp <- c(rownames(LR_heatmap)[i], colnames(LR_heatmap)[j], LR_heatmap[i,j])
      LR_matrix <- rbind(LR_matrix, temp)
    }
  }
  colnames(LR_matrix) <- c("ligand", "receptor", "weight")
  LR_matrix$weight <- as.numeric(LR_matrix$weight)
  LR_matrix <- LR_matrix[which(LR_matrix$weight > 0),]
  
  #get top ligands from NicheNet
  keep_rec <- AverageExpression(NKsubset, features = gsub("\\.", "-", unique(LR_matrix$receptor)), slot = "data", assays = "RNA")$RNA
  keep_rec <- as.data.frame(keep_rec)
  keep_rec <- keep_rec[which(keep_rec[,cell] > keep_rec[,setdiff(colnames(keep_rec),cell)]),]
  keep_rec$receptor <- rownames(keep_rec)
  keep_rec$recexp <- keep_rec[,cell] / keep_rec[,setdiff(colnames(keep_rec[,1:2]),cell)] #rNK v non rNK
  #keep_rec$recexp <- keep_rec[,cell] / keep_rec[,setdiff(colnames(keep_rec[,1:6]),cell)] #nk clust
  keep_rec <- keep_rec[,-c(1:2)] #rnK v non rNK
  #keep_rec <- keep_rec[,-c(1:6)] #nk clust
  LR_matrix <- LR_matrix[which(LR_matrix$receptor %in% keep_rec$receptor),]
  
  cellsubset <- GetAssayData(episubset, assay = "RNA", slot = "data") 
  cellsubset <- as.data.frame(cellsubset[which(rownames(cellsubset) %in% gsub("\\.", "-", unique(LR_matrix$ligand))),])
  cellsubset[cellsubset > 0] <- 1
  cellsubset$percent <- apply(cellsubset, 1, function(x) sum(x))/dim(cellsubset)[2]
  cellsubset$ligand <- rownames(cellsubset)
  keep_lig <- as.data.frame(cbind(cellsubset$ligand, cellsubset$percent))
  colnames(keep_lig) <- c("ligand", "ligpercent")
  LR_matrix <- LR_matrix[which(LR_matrix$ligand %in% keep_lig$ligand),]
  
  LR_matrix <- left_join(LR_matrix, keep_lig)
  LR_matrix <- left_join(LR_matrix, keep_rec)
  LR_matrix$LRexp <- as.numeric(LR_matrix$ligpercent) * as.numeric(LR_matrix$recexp) #rNK v nonrNK
  # LR_matrix$LRexp <- as.numeric(LR_matrix$ligpercent) * as.numeric(LR_matrix$recexp[,1] +LR_matrix$recexp[,2] +
  #                                                                    LR_matrix$recexp[,3] + LR_matrix$recexp[,4] +
  #                                                                    LR_matrix$recexp[,5]) #nk clusters
  dim(LR_matrix)
  
  lr.df.top <- LR_matrix %>% top_n(20, LRexp) #paper cutoff
  #lr.df.top <- LR_matrix %>% top_n(50, LRexp)
  lr.df.top <- lr.df.top[,c(1:3)]
  #lr.df.top <- unique(rbind(lr.df.top, PD1, PD2))
  
  ##setwd(DEGdir)
  
  setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/Nichenet")
  #generate circos plot with top 20 LR pairs by recexp
  # pdf(paste(cell, "TNBCcancerEpi_rNK_top20_circos_111722.pdf", sep = "_"), width = 8, height = 8)
  # pdf(paste(cell, "HRcancerEpi_rNK_top20_circos_111722.pdf", sep = "_"), width = 8, height = 8)
  # pdf(paste(cell, "HER2cancerEpi_rNK_top20_circos_111722.pdf", sep = "_"), width = 8, height = 8)
  
  setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/correctrNK_preprint")
  #pdf(paste(cell, "correct_TNBCEpi_rNK_top20_circos_112322.pdf", sep = "_"), width = 8, height = 8)
  #pdf(paste(cell, "correct_HREpi_rNK_top20_circos_112322.pdf", sep = "_"), width = 8, height = 8)
  pdf(paste(cell, "correct_HER2Epi_rNK_top20_circos_112322.pdf", sep = "_"), width = 8, height = 8)
  
  
  # pdf(paste(cell, "allcancerepi_8522_top50_TNBConly_circos.pdf", sep = "_"), width = 8, height = 8)
  # pdf(paste(cell, "allcancerepi_8522_top50_HRonly_circos.pdf", sep = "_"), width = 8, height = 8)
  # pdf(paste(cell, "allcancerepi_8522_top50_HER2only_circos.pdf", sep = "_"), width = 8, height = 8)
  
  # pdf(paste(cell, "allcancerepi_72522_top20_ALLEpi_circos.pdf", sep = "_"), width = 8, height = 8)
  # pdf(paste(cell, "allcancerepi_72522_top20_TNBConly_circos.pdf", sep = "_"), width = 8, height = 8)
  # pdf(paste(cell, "allcancerepi_72522_top20_HER2only_circos.pdf", sep = "_"), width = 8, height = 8)
  #pdf(paste(cell, "allcancerepi_72522_top20_HRonly_circos.pdf", sep = "_"), width = 8, height = 8)
  
  circos.par(canvas.ylim=c(-1.5,1.5), # edit  canvas size
             track.margin = c(0.01, 0)) # adjust bottom and top margin

  chordDiagram(lr.df.top,
               directional = 1,
               link.sort = TRUE,
               link.decreasing = T,
               link.visible = lr.df.top$weight,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.05))

  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
  }, bg.border = NA)
  circos.clear()

  dev.off()

  
  # # create NicheNet output plots
  # p <- temp_nichenet$ligand_activity_target_heatmap
  # ggsave(paste0("GM", GM, "_", cell, "_RL_activity_heatmap.pdf"), plot = p, width = 15, height = 10)

  targets <- temp_nichenet$ligand_target_df[which(temp_nichenet$ligand_target_df$ligand %in% LR_matrix$ligand),]
  targets <- targets %>% group_by(target) %>% mutate(avg_weight = mean(weight)) %>%
  add_count(target)
  targets <- unique(targets[,-c(1,3)])
  # # p <- DotPlot(primsubset, assay = "RNA", features = unique(targets$target), cols = "RdBu") + RotatedAxis()
  # # ggsave(paste0("GM", GM, "_RL_toptarget_dotplot.pdf"), plot = p, width = 20, height = 5)
  targets <- as.data.frame(targets)
  targets$NK <- cell
  targets$cell <- "HER2+ Epithelial Cells" #"Epithelial Cells"
  targetlist[[cell]] <- targets
  all_targets <- rbind(all_targets, targets)
  
  receptors <- LR_matrix
  # # p <- DotPlot(primsubset, assay = "RNA", features = unique(receptors$receptor), cols = "RdBu") + RotatedAxis()
  # # ggsave(paste0("GM", GM, "_RL_topreceptor_dotplot.pdf"), plot = p, width = 20, height = 5)
  receptors <- as.data.frame(receptors)
  receptors$NK <- cell
  receptors$cell <- "HER2+ Epithelial Cells"#"Epithelial Cells"
  receptorlist[[cell]] <- receptors
  all_receptors <- rbind(all_receptors, receptors)
  
}

##revision _________________
write.csv(all_receptors, "rNK_TNBCcancerEpiONLY_receptors_111822.csv")
write.csv(all_targets, "rNK_TNBCcancerEpiONLY_targets_111822.csv")
write.csv(all_receptors, "rNK_HRcancerEpiONLY_receptors_111822.csv")
write.csv(all_targets, "rNK_HRcancerEpiONLY_targets_111822.csv")
write.csv(all_receptors, "rNK_HER2cancerEpiONLY_receptors_111822.csv")
write.csv(all_targets, "rNK_HER2cancerEpiONLY_targets_111822.csv")


##correction _________________
setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/correctrNK_preprint")
# write.csv(all_receptors, "preprint_rNK_TNBCepi_receptors_112322.csv")
# write.csv(all_targets, "preprint_rNK_TNBCepi_targets_112322.csv")
# write.csv(all_receptors, "preprint_rNK_HRepi_receptors_112322.csv")
# write.csv(all_targets, "preprint_rNK_HRepi_targets_112322.csv")
write.csv(all_receptors, "preprint_rNK_HER2epi_receptors_112322.csv")
write.csv(all_targets, "preprint_rNK_HER2epi_targets_112322.csv")


##preprint __________
saveRDS(receptorlist, "NK_TNBCepi_receptors_72622.rds")
saveRDS(targetlist, "NK_TNBCepi_targets_72622.rds")

saveRDS(receptorlist, "NK_HRepi_receptors_72622.rds")
saveRDS(targetlist, "NK_HRepi_targets_72622.rds")

saveRDS(receptorlist, "NK_HER2epi_receptors_72622.rds")
saveRDS(targetlist, "NK_HER2epi_targets_72622.rds")


write.csv(all_receptors, "NK_allcancerepi_receptorsTNBC.csv")
write.csv(all_targets, "NK_allcancerepi_targetsTNBC.csv")

write.csv(all_receptors, "NK_allcancerepi_receptorsHR.csv")
write.csv(all_targets, "NK_allcancerepi_targetsHR.csv")

write.csv(all_receptors, "NK_allcancerepi_receptorsHER2.csv")
write.csv(all_targets, "NK_allcancerepi_targetsHER2.csv")



# ================================================================== ======


#2J Survival (see TCGA_Preprocessing.R in /project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/bulkRNA/TCGA) ========================
# ================================================================== ======


# Filtering TCGA =======================================================================

setwd(TCGA_dir)
BRCA <- readRDS(file = "BRCA_SummExpObj.rds")

#https://rdrr.io/bioc/TCGAbiolinks/man/TCGAanalyze_SurvivalKM.html and
#https://costalab.ukaachen.de/open_data/Bioinformatics_Analysis_in_R_2019/BIAR_D3/handout.html
# clinical_patient_Cancer <- GDCquery_clinic("TCGA-BRCA","clinical")
# #also pick with high immune infiltrate
# clinical_BRCAsub <- clinical_patient_Cancer[clinical_patient_Cancer$gender=="female", #& clinical_patient_Cancer$prior_treatment=="No", 
#                                             colnames(clinical_patient_Cancer)]
# head(clinical_BRCAsub)
# table(clinical_BRCAsub$gender)
# table(clinical_BRCAsub$prior_treatment)

colnames(colData(BRCA)) #metadata

#BRCA.subset <- BRCA[,((BRCA$gender == "female") & is.na(BRCA$gender == "female") == F)]
BRCA.subset <- BRCA
dim(BRCA.subset)
dim(BRCA)
table(BRCA$gender)
table(BRCA.subset$gender)

head(colData(BRCA)$treatments)
table(BRCA.subset$gender)


#https://costalab.ukaachen.de/open_data/Bioinformatics_Analysis_in_R_2019/BIAR_D3/handout.html#6_gene_expression_and_survival
#dataBRCAcomplete <- log2(assay(BRCA))
dataBRCAcomplete <- assay(BRCA.subset)


geneinfo_df <- rowData(BRCA.subset)

dataBRCAcomplete <- as.data.frame(dataBRCAcomplete)
dataBRCAcomplete$ensembl_gene_id <- rownames(dataBRCAcomplete)
dataBRCAcomplete <- dataBRCAcomplete[,c(1223, 1:1222)]
dataBRCAcomplete[1:5,1:5]

head(geneinfo_df) 
geneinfo_df$ConvertedGenes <- NA
geneinfo_df<- geneinfo_df[,c(2, 4, 1, 3)]
head(geneinfo_df)


dataBRCAcomplete <- merge(geneinfo_df, dataBRCAcomplete, by = "ensembl_gene_id")
dataBRCAcomplete <- as.data.frame(dataBRCAcomplete)
dataBRCAcomplete[1:5,1:5]


# converting genes =====================================

library(limma)
library(org.Hs.eg.db)
oldgenes <- dataBRCAcomplete[,"external_gene_name"]
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

dataBRCAcomplete$ConvertedGenes <- usegenes
new_mat <- dataBRCAcomplete
new_mat[1:5,1:5]

#new_mat$ConvertedGenes <- usegenes
new_mat <- new_mat[which(!(new_mat$ConvertedGenes %in%
                             usegenes[duplicated(usegenes)])),]

# dup_mat <- mRNAmicro[which((mRNAmicro$ConvertedGenes %in%
#                              usegenes[duplicated(usegenes)])),]
# 
# dup_mat[1:5,1:5]
# dim(mRNAmicro)
# dim(new_mat)

which(grepl("SMAD5", new_mat$ConvertedGenes))
new_mat

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- dataBRCAcomplete[which(dataBRCAcomplete$ConvertedGenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  k$ConvertedGenes <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}

dim(new_mat)
dim(dataBRCAcomplete)

rownames(new_mat) <- new_mat[,"ConvertedGenes"]

dataBRCAcomplete <- new_mat
dataBRCAcomplete[1:5,1:5]

dataBRCAcomplete_genedf <- dataBRCAcomplete[,1:4]
head(dataBRCAcomplete_genedf)

rownames(dataBRCAcomplete) <- dataBRCAcomplete$ConvertedGenes
dataBRCAcomplete <- dataBRCAcomplete[,-c(1:4)]
dataBRCAcomplete[1:5,1:5]

dataBRCAcomplete <- as.data.frame(sapply(dataBRCAcomplete, as.numeric))
rownames(dataBRCAcomplete) <- dataBRCAcomplete_genedf$external_gene_name
colnames(dataBRCAcomplete) <- gsub(".", "-", colnames(dataBRCAcomplete), fixed=TRUE)
dataBRCAcomplete[1:5,1:5]

setwd(TCGA_dir)
saveRDS(dataBRCAcomplete, "TCGA_genesconverted_72722.rds")
saveRDS(dataBRCAcomplete_genedf, "TCGA_genesinfo_withconverted_72722.rds")


# normalization ========
library(EDASeq)

setwd(TCGA_dir)
dataBRCAcomplete <- readRDS("TCGA_genesconverted_72722.rds")
dataBRCAcomplete_genedf <- readRDS("TCGA_genesinfo_withconverted_72722.rds")

dataBRCAcomplete[1:5,1:5]
head(dataBRCAcomplete_genedf)

rownames(dataBRCAcomplete) <- dataBRCAcomplete_genedf$ensembl_gene_id
#names(dataBRCAcomplete) <- colnames(dataBRCAcomplete)
dataNorm <- TCGAanalyze_Normalization(tabDF = dataBRCAcomplete, geneInfo =  geneInfoHT) #use this one for ensembl IDs

#just primary
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataNorm),
                                   typesample = c("TP"))
#including other sample types
# samplesmet <- TCGAquery_SampleTypes(barcode = colnames(dataNorm),
#                                    typesample = c("TM"))

dataNorm<-dataNorm[,colnames(dataNorm) %in% samplesTP]
#met
#dataNorm<-dataNorm[,colnames(dataNorm) %in% samplesmet]
tcga.vst<-DESeq2::vst(as.matrix(round(dataNorm)))



# Find genes in object =================================================

# geneID.list <- character()
# genelist <- c("CALD1", "CLU", "ALOX12", "LTBP1", "CAVIN2", "PARVB", 
#               "GP6", "SCD", "ITGAX", "NR4A3", "CCL4", "CR2", "HEATR9", 
#               "XDH", "RASGRP2", "MID1", "JUN", "CMKLR1", "DUSP1", "FOS", 
#               "ABCA1", "TNFAIP3", "NR4A1", "KLRG1", "DTX1", "NHSL2", 
#               "GFRA2", "FAM81A", "CX3CR1", "RHPN1", "HES1", "F5", "GAS2L1", 
#               "THBS1", "MYLK", "TMTC1", "FOSB", "NR4A2", "MPIG6B", "SLC6A4", 
#               "PLXNA4", "VWF", "TUBB1", "SLC7A5", #up
#               "BCAT1", "ALDH1L2", "COX6A2", "PYCR1", "LHFPL2", "AHRR", 
#               "EXTL1", "ASNS", "CHAC1", "MTHFD2", "NEK6", "SLC6A9", 
#               "FMNL2", "ASB2", "SLC7A3", "AVIL", "CDH1", "CISH", "LGALS3", 
#               "GPT2", "CXCR6", "TRIB3", "CDKN1A", "ATF5", 
#               "SLC1A4", "PMEPA1", "CEMIP2", "OSGIN1", "ZNF503", "ITGA1", "ISG20", 
#               "PACSIN1", "TBC1D16", "RN7SL1", "SH3PXD2B", "SCN3B", "OSBPL1A", 
#               "OSBPL1A", "ME1", "HPGDS", "PPP2R2C", "CLBA1", "HMOX1", "NQO1", 
#               "CARS1", #not in TCGA object
#               "SSTR2", "SNORA23") ##down
# 
# ##find IDs for gene list
# for (gene in genelist) {
#   gene.query <- which(grepl(gene, dataBRCAcomplete_genedf$external_gene_name)) # 5311 *stagnates >4000 days but worse otherwise
#   gene_name = dataBRCAcomplete_genedf[gene.query, "external_gene_name"] #SNAIL = ; TP53 = 
#   gene.name.match <- grep(paste0("^",gene,"$"), gene_name)
#   gene.final.result = gene.query[gene.name.match]
#   geneID.list <- c(geneID.list, gene.final.result)
# }
# 
# geneID.list <- as.numeric(geneID.list)
# 
# 

tcga.vst[1:5,1:5]

geneID.list <- character()


genelist <- c("CALD1", "CLU", "ALOX12", "LTBP1", "CAVIN2", "PARVB", 
              "GP6", "SCD", "ITGAX", "NR4A3", "CCL4", "CR2", "HEATR9", 
              "XDH", "RASGRP2", "MID1", "JUN", "CMKLR1", "DUSP1", "FOS", 
              "ABCA1", "TNFAIP3", "NR4A1", "KLRG1", "DTX1", "NHSL2", 
              "GFRA2", "FAM81A", "CX3CR1", "RHPN1", "HES1", "F5", "GAS2L1", 
              "THBS1", "MYLK", "TMTC1", "FOSB", "NR4A2", "MPIG6B", "SLC6A4", 
              "PLXNA4", "VWF", "TUBB1")#, "SLC7A5") #up

##find IDs for gene list
for (gene in genelist) {
  gene.query <- which(grepl(gene, dataBRCAcomplete_genedf$ConvertedGenes)) # 5311 *stagnates >4000 days but worse otherwise
  gene_name = dataBRCAcomplete_genedf[gene.query, "ConvertedGenes"] #SNAIL = ; TP53 = 
  gene.name.match <- grep(paste0("^",gene,"$"), gene_name)
  gene.final.result = gene.query[gene.name.match]
  geneID.list <- c(geneID.list, gene.final.result)
}

geneID.list <- as.numeric(geneID.list)

# Plot survival for gene signature =================================================

#triplecheck correct gene names 
gene_name = dataBRCAcomplete_genedf[geneID.list, "ConvertedGenes"] 
gene_name

tcga.vst <- as.data.frame(tcga.vst)

#match geneID to ensembl row name
gene_id = dataBRCAcomplete_genedf[geneID.list, "ensembl_gene_id"] 
gene_id

#dataSignature<-dataNorm[rownames(dataNorm) %in% genesHs$HGNC.symbol,]
dataSignature<-tcga.vst[rownames(tcga.vst) %in% gene_id,]

dim(dataSignature)
rownames(dataSignature)


#order the samples
idx<-order(-colSums(log2(dataSignature[rownames(dataSignature),]+1)) + #%in% subset(NKsignatureReprogrammed,group=="down")$gene,]+1))+
             colSums(log2(dataSignature[rownames(dataSignature),]+1)), decreasing = TRUE)# %in% subset(NKsignatureReprogrammed,group=="up")$gene,]+1)), decreasing = TRUE)
idx<-order(colSums(log2(dataSignature[rownames(dataSignature),]+1)), decreasing = TRUE) #%in% subset(NKsignatureReprogrammed,group=="up")$gene,]+1)), decreasing = TRUE)
SigScore<- -colSums(log2(dataSignature[rownames(dataSignature),]+1))+ #%in% subset(NKsignatureReprogrammed,group=="down")$gene,]+1))+
  colSums(log2(dataSignature[rownames(dataSignature),]+1)) #%in% subset(NKsignatureReprogrammed,group=="up")$gene,]+1))
hist(SigScore)

# #original
# PatientHigh<-substr(colnames(dataSignature)[idx],1,12)[1:300]
# PatientLow<-substr(colnames(dataSignature)[idx],1,12)[301:dim(dataSignature)[2]]

#mine
dim(dataSignature)
PatientHigh<-substr(colnames(dataSignature)[idx],1,12)[1:300]
PatientLow<-substr(colnames(dataSignature)[idx],1,12)[301:dim(dataSignature)[2]]

#original
PatientQ4<-substr(colnames(dataSignature)[idx],1,12)[1:200]
PatientQ1<-substr(colnames(dataSignature)[idx],1,12)[900:dim(dataSignature)[2]]


BRCA.surv.select<-survivalTCGA(BRCA.clinical, extract.cols = c("patient.histological_type",
                                                               "patient.clinical_cqcf.tumor_type" ,
                                                               "patient.age_at_initial_pathologic_diagnosis",
                                                               "patient.days_to_death",
                                                               "patient.new_tumor_events.new_tumor_event_after_initial_treatment"))
##try filtering here for age stuff if bad ;'(
brca.subtype <- TCGAquery_subtype(tumor = "brca")
timelim = 3650#2000

#BRCA.surv.select<-BRCA.surv.select %>% #dplyr::filter(bcr_patient_barcode %in% subset(brca.subtype,BRCA_Subtype_PAM50 != "Basal")$patient) %>% 
#  dplyr::filter(times<timelim)
BRCA.surv.select$scoreCut<- ifelse(BRCA.surv.select$bcr_patient_barcode %in% PatientHigh, "High","Low")
BRCA.surv.select$scoreQ14<- ifelse(BRCA.surv.select$bcr_patient_barcode %in% PatientQ1, "Q1",
                                   ifelse(BRCA.surv.select$bcr_patient_barcode %in% PatientQ4,"Q4","Other"))


# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/TCGA_Analysis/ReprogGenes/noinfiltrate/")

# pdf(paste("ReprogSig_survplot_NOinfiltrate_61322.pdf", sep = ""), height = 6,width=7)
# pdf("test.pdf")
# RTCGA::kmTCGA(
#   BRCA.surv.select,
#   times = "times",
#   xlim=c(0,timelim),
#   status = "patient.vital_status",
#   explanatory.names = c("scoreCut"),
#   pval = TRUE,
#   conf.int = FALSE
# )
# dev.off()

library(readxl)
setwd(TCGAresults_Dir)
cibersort<-read_xlsx("Table_1_Landscape of Immune Microenvironment Under Immune Cell Infiltration Pattern in Breast Cancer.xlsx",
                     sheet="S I", skip = 1)
useSample<-substr(subset(cibersort,`NK cells activated`>0.015 | `NK cells resting` >0.015)$id,6,17)

BRCA.surv.select.sub<-BRCA.surv.select %>% dplyr::filter(bcr_patient_barcode %in% useSample)
BRCA.surv.select.sub$scoreCut<- ifelse(BRCA.surv.select.sub$bcr_patient_barcode %in% PatientHigh, "High","Low")



setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/TCGA_Analysis/ReprogGenes/yesinfiltrate/") ##preprint
setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/survival")

#pdf("ReprogSig_survplot_immuneinfiltrate_NKrestact0.015_72722.pdf", width = 5, height = 4.4)
#pdf("ReprogSig_survplot_immuneinfiltrate_NKrestact0.015_10years_72722.pdf", width = 5, height = 4.4) ##preprint
pdf("reprogSig_survplot_immuneinfiltrate_NKrestact0.015_10years_noSLC7A5_102022.pdf", width = 5, height = 4.4)
#pdf("reprogSig_survplot_immuneinfiltrate_NKrestact0.015_2000days_noSLC7A5_111722.pdf", width = 5, height = 4.4)
RTCGA::kmTCGA(
  BRCA.surv.select.sub,
  times = "times",
  xlim=c(0,timelim),
  break.time.by = 365,
  status = "patient.vital_status",
  explanatory.names = c("scoreCut"),
  pval = TRUE,
  conf.int = FALSE,
  risk.table = F
)
dev.off()

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/TCGA_Analysis/ReprogGenes/")


# survival for individual rNK genes (REVIEWER ADDITION) +++++++++++++++++++++ ========================

#triplecheck correct gene names 
gene_name = dataBRCAcomplete_genedf[geneID.list, "external_gene_name"] 
gene_name

tcga.vst <- as.data.frame(tcga.vst)

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/rNK_individual")

for (Indgene in geneID.list)
{
  
  
  #match geneID to ensembl row name
  gene_id = dataBRCAcomplete_genedf[c(Indgene), "ensembl_gene_id"] 
  gene_id
  
  #get gene symbol to label pdf files
  gene_name = dataBRCAcomplete_genedf[Indgene, "external_gene_name"] 
  gene_name
  
  dataSignature<-tcga.vst[rownames(tcga.vst) %in% gene_id,]
  
  
  #order the samples
  idx<-order(-colSums(log2(dataSignature[rownames(dataSignature),]+1)) + #%in% subset(NKsignatureReprogrammed,group=="down")$gene,]+1))+
               colSums(log2(dataSignature[rownames(dataSignature),]+1)), decreasing = TRUE)# %in% subset(NKsignatureReprogrammed,group=="up")$gene,]+1)), decreasing = TRUE)
  idx<-order(colSums(log2(dataSignature[rownames(dataSignature),]+1)), decreasing = TRUE) #%in% subset(NKsignatureReprogrammed,group=="up")$gene,]+1)), decreasing = TRUE)
  SigScore<- -colSums(log2(dataSignature[rownames(dataSignature),]+1))+ #%in% subset(NKsignatureReprogrammed,group=="down")$gene,]+1))+
    colSums(log2(dataSignature[rownames(dataSignature),]+1)) #%in% subset(NKsignatureReprogrammed,group=="up")$gene,]+1))
  hist(SigScore)
  
  
  dim(dataSignature)
  PatientHigh<-substr(colnames(dataSignature)[idx],1,12)[1:300]
  PatientLow<-substr(colnames(dataSignature)[idx],1,12)[301:dim(dataSignature)[2]]
  
  PatientQ4<-substr(colnames(dataSignature)[idx],1,12)[1:200]
  PatientQ1<-substr(colnames(dataSignature)[idx],1,12)[900:dim(dataSignature)[2]]
  
  
  BRCA.surv.select<-survivalTCGA(BRCA.clinical, extract.cols = c("patient.histological_type",
                                                                 "patient.clinical_cqcf.tumor_type" ,
                                                                 "patient.age_at_initial_pathologic_diagnosis",                                                          "patient.new_tumor_events.new_tumor_event_after_initial_treatment"))
  
  brca.subtype <- TCGAquery_subtype(tumor = "brca")
  timelim = 3650
  
  BRCA.surv.select$scoreCut<- ifelse(BRCA.surv.select$bcr_patient_barcode %in% PatientHigh, "High","Low")
  BRCA.surv.select$scoreQ14<- ifelse(BRCA.surv.select$bcr_patient_barcode %in% PatientQ1, "Q1",
                                     ifelse(BRCA.surv.select$bcr_patient_barcode %in% PatientQ4,"Q4","Other"))
  
  library(readxl)
  setwd(TCGAresults_Dir)
  cibersort<-read_xlsx("Table_1_Landscape of Immune Microenvironment Under Immune Cell Infiltration Pattern in Breast Cancer.xlsx",
                       sheet="S I", skip = 1)
  useSample<-substr(subset(cibersort,`NK cells activated`>0.015 | `NK cells resting` >0.015)$id,6,17)
  
  BRCA.surv.select.sub<-BRCA.surv.select %>% dplyr::filter(bcr_patient_barcode %in% useSample)
  BRCA.surv.select.sub$scoreCut<- ifelse(BRCA.surv.select.sub$bcr_patient_barcode %in% PatientHigh, "High","Low")
  
  setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/survival/rNK_individual")  
  #setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/survival/rNKind_2000days")  
  
  pdf(paste(gene_name, "_survplot_Immuneinfiltrate_NKrestact0.015_10years_101722.pdf", sep = ""), height = 6,width=7)
  #pdf(paste(gene_name, "_survplot_Immuneinfiltrate_NKrestact0.015_2000days_111722.pdf", sep = ""), height = 6,width=7)
  print(RTCGA::kmTCGA(
    BRCA.surv.select.sub,
    times = "times",
    xlim=c(0,timelim),
    break.time.by = 365,
    status = "patient.vital_status",
    explanatory.names = c("scoreCut"),
    pval = TRUE,
    conf.int = FALSE,
    risk.table = F
  ))
  dev.off()
  
}


## correlation age (REVIEWER ADDITION) ++++++++++++++++++++++++++++++++++++++ -----

head(colData(BRCA.surv.select))
colnames(BRCA.clinical)[grep("days_to_death", colnames(BRCA.clinical))]
colnames(BRCA.clinical)[grep("age", colnames(BRCA.clinical))]

BRCA.surv.select<-survivalTCGA(BRCA.clinical, extract.cols = c("patient.histological_type",
                                                               "patient.clinical_cqcf.tumor_type" ,
                                                               "patient.age_at_initial_pathologic_diagnosis",
                                                               "patient.days_to_death",
                                                               "patient.new_tumor_events.new_tumor_event_after_initial_treatment"))

head(BRCA.surv.select)

head(BRCA.surv.select$patient.age_at_initial_pathologic_diagnosis)
table(BRCA.surv.select$patient.days_to_death)

class(BRCA.surv.select$patient.age_at_initial_pathologic_diagnosis)
BRCA.surv.select$patient.age_at_initial_pathologic_diagnosis <- as.numeric(BRCA.surv.select$patient.age_at_initial_pathologic_diagnosis)
BRCA.surv.select <- BRCA.surv.select[which(!is.na(BRCA.surv.select$patient.age_at_initial_pathologic_diagnosis)),]

class(BRCA.surv.select$patient.days_to_death)
BRCA.surv.select$patient.days_to_death <- as.numeric(BRCA.surv.select$patient.days_to_death)
BRCA.surv.select <- BRCA.surv.select[which(!is.na(BRCA.surv.select$patient.days_to_death)),]

BRCA.surv.select <- BRCA.surv.select %>% 
  mutate(patient.years_to_death = patient.days_to_death / 365)
summary(BRCA.surv.select$patient.years_to_death)
summary(BRCA.surv.select$patient.age_at_initial_pathologic_diagnosis)

library(ggpubr)
p <- ggscatter(BRCA.surv.select, x = "patient.age_at_initial_pathologic_diagnosis", 
               y = "patient.years_to_death",
               add = "reg.line",
               conf.int = TRUE,
               add.params = list(color = "blue", fill = "lightgray"),
               xlab = "Age",
               ylab = "Years to Death") + theme_classic() +
  stat_cor(method = "pearson", size = 6, mapping = aes(label = "..p.adj..")) +
  scale_y_continuous(breaks = c(0, 5, 10, 15),
                     limits = c(0, 15))

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2/survival")
pdf("agevdeath_TCGAsurvival_112122.pdf", width = 4, height = 3.5)

p+ theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 16)) +
  theme(axis.text = element_text(size = 16))
dev.off()


##x = age at diagnosis/365 for years
## y = days to death/365 for years

##scatter plot with linear regression line and see if significant or not
## if not significant, then yay, if is significant, then do survival individuall on each age group


