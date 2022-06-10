


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





#kara =======

#test cells ========================================

test.df <- data.frame(
  NK.Total = 0,
  NK.genes = "",
  T.Total = 0,
  T.genes = "",
  B.Total = 0,
  B.genes = "",
  Plasma.Total = 0,
  Plasma.genes = "",
  Mye.Total = 0,
  Mye.genes = "",
  MyoEpi.Total = 5,
  MyoEpi.genes = c("KRT5","KRT14","ACTA2",
                   "TAGLN"),
  PVL.Total = 0,
  PVL.genes = "",
  Fibro.Total = 0,
  Fibro.genes = "",
  Endo.Total = 0,
  Endo.genes = "",
  Immune.Total = 0,
  Immune.genes = "",
  Epi.Total = 2,
  Epi.genes = c("KRT5","KRT14") ,
  StrongEpi.Total = 1,
  StrongEpi.genes = "KRT18")


# ____________________________________________________________________________
# ____________________________________________________________________________
## Reduce expression matrix to just celltype markers _________________________
# ____________________________________________________________________________
# ____________________________________________________________________________


CellType.upMark <- c("PTPRC","FGFBP2", "SPON2", "KLRF1", "FCGR3A", "KLRD1", "TRDC", "CD3E",
                     "NCAM1", #NK
                     "CD2","CD3D", "CD3E","CD3G","CD8A","CD8B",
                     "CD4",#T
                     "CD79A","CD79B","BLNK","CD19",
                     "MS4A1",#B
                     "CD27","IGHA1","SDC1","TNFRSF17","JCHAIN",
                     "MZB1","DERL3","CD38","IGHG1","IGHG3",
                     "IGHG4", #plasma
                     "ITGAM","ITGAX","CD14",
                     "CD1C","CD1A","CD68",
                     "CD33", #myeloid
                     "EPCAM", "EGFR","FZR1", "KRT14", "ITGA6",
                     "KRT5", "TP63", "KRT17", "MME", "KRT8",
                     "KRT18", "KRT19", "FOXA1", "GATA3", "MUC1",
                     "CD24", "GABRP", #Epi
                     "KRT5","KRT14","ACTA2",
                     "TAGLN", #myoEpi
                     "MCAM", "PDGFRB", #PVL
                     "FAP","THY1","DCN","COL1A1","COL1A2","COL6A1",
                     "COL6A2","COL6A3", #Fibro
                     "PECAM1","VWF","CDH5","SELE", #Endo
                    "CD8A","CD8B", #CD8
                       "CD4", "FOXP3",#CD4, Treg
                    "CD33", #MDSc
                    "ITGAM", 
                    "FUT4", 
                    "CEACAM8",
                    "IL4R",
                    "HLA-DRA",
                    "INHBA", "IL1RN", "CCL4", "NLRP3", #Macro
                    "EREG", "IL1B", "LYVE1", "PLTP",
                    "SELENOP", "C1QC", "C1QA", "APOE",
                    "FCN1", "S100A9", "S100A8", "FCGR3A", "LST1", "LILRB2", #Mono
                    "CXCR1",#Neut
                    "FUT4", 
                    "FCGR3A", 
                    "CSF3R", 
                    "S100A9","TNF","CD274",
                    "LILRA4","GZMB","IL3RA", #DC
                    "CLEC9A", "FLT3", "IDO1", 
                    "CD1C", "FCER1A","HLA-DQA1", 
                    "LAMP3", "CCR7", "FSCN1",
                    "KIT","TPSAB1","CPA4")   #Mast

#table(Idents(combo.reference))
##subset cell population from main primary object
#combo.NK <- subset(combo.reference, idents = "NK Cells")

#get expression matrix for marker genes


table(Idents(combo.reference))
DefaultAssay(combo.reference) <- "RNA"

combo.reference <- NormalizeData(combo.reference, assay = "RNA")


#NK cells ______________________________________________________________________

NKsub <- subset(combo.reference, idents = "NK Cells")
NKsub <- NormalizeData(NKsub, assay = "RNA")
NKsub.mat <- GetAssayData(NKsub[unique(CellType.upMark),], slot = "data") #because it's corrected for sequencing depth and such
NKsub.mat <- as.matrix(NKsub.mat)

#Epithelial cells ______________________________________________________________________

Episub <- subset(combo.reference, idents = "Epithelial Cells")
Episub <- NormalizeData(Episub, assay = "RNA")
Episub.mat <- GetAssayData(Episub[unique(CellType.upMark),], slot = "data") #because it's corrected for sequencing depth and such
Episub.mat <- as.matrix(Episub.mat)

#Stromal cells ______________________________________________________________________

Stromasub <- subset(combo.reference, idents = c("Endothelial Cells", "Fibroblasts",
                                                "Perivascular-like (PVL) Cells", "Myoepithelial Cells"))
Stromasub <- NormalizeData(Stromasub, assay = "RNA")
Stromasub.mat <- GetAssayData(Stromasub[unique(CellType.upMark),], slot = "data") #because it's corrected for sequencing depth and such
Stromasub.mat <- as.matrix(Stromasub.mat)

#B cells ______________________________________________________________________

Bsub <- subset(combo.reference, idents = c("B Cells", "Plasma Cells"))
Bsub <- NormalizeData(Bsub, assay = "RNA")
Bsub.mat <- GetAssayData(Bsub[unique(CellType.upMark),], slot = "data") #because it's corrected for sequencing depth and such
Bsub.mat <- as.matrix(Bsub.mat)

#Myeloid cells ______________________________________________________________________

Myesub <- subset(combo.reference, idents = c("Mast Cells",
                                             "Dendritic Cells", "Neutrophils",
                                             "Monocytes", "MDSCs"))
Myesub <- NormalizeData(Myesub, assay = "RNA")
Myesub.mat <- GetAssayData(Myesub[unique(CellType.upMark),], slot = "data") #because it's corrected for sequencing depth and such
Myesub.mat <- as.matrix(Myesub.mat)

#T cells ______________________________________________________________________

Tsub <- subset(combo.reference, idents = c("Regulatory T Cells",
                                           "CD4+ T Cells", "CD8+ T Cells"))
Tsub <- NormalizeData(Tsub, assay = "RNA")
Tsub.mat <- GetAssayData(Tsub[unique(CellType.upMark),], slot = "data") #because it's corrected for sequencing depth and such
Tsub.mat <- as.matrix(Tsub.mat)


Subset.Mat <- NKsub.mat
Subset.Mat <- Episub.mat
Subset.Mat <- Stromasub.mat
Subset.Mat <- Bsub.mat
Subset.Mat <- Myesub.mat

Subset.Mat <- Tsub.mat

# ____________________________________________________________________________
# ____________________________________________________________________________
## Make list, where each entry is individual cell __________________________________________________
## (breaking up expression matrix to where each cell has its own expression matrix)   ______________
## Remove genes that aren't expressed for each cell's expression matrix  ___________________________
# ____________________________________________________________________________
# ____________________________________________________________________________


cell.list <- list()
zero.for.everything <- list()
for (i in 1:length(colnames(Subset.Mat))){
  if (all(Subset.Mat[,i] == 0) == T){
    get.cell <- Subset.Mat[,i, drop = F]
    get.genes.expressed <- rownames(get.cell)
    cell.name <- colnames(Subset.Mat)[i]
    zero.for.everything[[cell.name]] <- "Zero Expression for All Markers"
    
  }
  else{
    get.cell <- Subset.Mat[,i, drop = F]
    get.genes.expressed <- rownames(get.cell)[which(get.cell != 0)]
    cell.name <- colnames(Subset.Mat)[i]
    cell.list[[cell.name]] <- get.genes.expressed
  }
  
}

length(cell.list)
head(zero.for.everything)
length(zero.for.everything)
# ____________________________________________________________________________
# ____________________________________________________________________________
## Find number of cell type markers expressed by each cell  ________________________________________
# ____________________________________________________________________________
# ____________________________________________________________________________


setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/GeneralAnalyses/")
Markers <- read.csv("primmarkers_onlyup6322.csv")
head(Markers)

Marker.list <- as.list(as.data.frame(Markers))
Marker.list <- lapply(Marker.list, function(z){ z[!is.na(z) & z != ""]})

Marker.list <- lapply(Marker.list, factor)
cell.list <- lapply(cell.list, factor)

#test cell: Wu2021prim_CID3586_AGGGATGCATGGAATA <- should be a T cell?
markersum.list <- list()
for (c in names(cell.list)){
  cell.test <- as.character(cell.list[[c]])
  #test NK marker content
  test.NK <- as.data.frame(table(cell.test[cell.test %in% Marker.list$NK_Cells]))
  test.NK$Freq <- as.numeric(test.NK$Freq)
  NK.sum <- sum(test.NK$Freq)
  NK.genes <- as.character(test.NK$Var1)
  #test T marker content
  test.T <- as.data.frame(table(cell.test[cell.test %in% Marker.list$T_Cells]))  
  test.T$Freq <- as.numeric(test.T$Freq)
  T.sum <- sum(test.T$Freq)
  T.genes <- as.character(test.T$Var1)
  #test B marker content
  test.B <- as.data.frame(table(cell.test[cell.test %in% Marker.list$B_Cells]))  
  test.B$Freq <- as.numeric(test.B$Freq)
  B.sum <- sum(test.B$Freq)
  B.genes <- as.character(test.B$Var1)
  #test Plasma marker content
  test.Plasma <- as.data.frame(table(cell.test[cell.test %in% Marker.list$Plasma_Cells]))  
  test.Plasma$Freq <- as.numeric(test.Plasma$Freq)
  Plasma.sum <- sum(test.Plasma$Freq)
  Plasma.genes <- as.character(test.Plasma$Var1)
  #test Myeloid marker content
  test.Mye <- as.data.frame(table(cell.test[cell.test %in% Marker.list$Myeloid_Cells]))  
  test.Mye$Freq <- as.numeric(test.Mye$Freq)
  Mye.sum <- sum(test.Mye$Freq)
  Mye.genes <- as.character(test.Mye$Var1)
  #test Myoepi marker content
  test.Myoepi <- as.data.frame(table(cell.test[cell.test %in% Marker.list$Myoepithelial_Cells]))  
  test.Myoepi$Freq <- as.numeric(test.Myoepi$Freq)
  Myoepi.sum <- sum(test.Myoepi$Freq)
  Myoepi.genes <- as.character(test.Myoepi$Var1)
  #test PVL marker content
  test.PVL <- as.data.frame(table(cell.test[cell.test %in% Marker.list$PVL]))  
  test.PVL$Freq <- as.numeric(test.PVL$Freq)
  PVL.sum <- sum(test.PVL$Freq)
  PVL.genes <- as.character(test.PVL$Var1)
  #test Fibro marker content
  test.Fibro <- as.data.frame(table(cell.test[cell.test %in% Marker.list$Fibroblasts]))  
  test.Fibro$Freq <- as.numeric(test.Fibro$Freq)
  Fibro.sum <- sum(test.Fibro$Freq)
  Fibro.genes <- as.character(test.Fibro$Var1)
  #test Endo marker content
  test.Endo <- as.data.frame(table(cell.test[cell.test %in% Marker.list$Endothelial_Cells]))  
  test.Endo$Freq <- as.numeric(test.Endo$Freq)
  Endo.sum <- sum(test.Endo$Freq)
  Endo.genes <- as.character(test.Endo$Var1)
  #test CD45 content
  test.Immune <- as.data.frame(table(cell.test[cell.test %in% Marker.list$Immune]))  
  test.Immune$Freq <- as.numeric(test.Immune$Freq)
  Immune.sum <- sum(test.Immune$Freq)
  Immune.genes <- as.character(test.Immune$Var1)
  #test Epithelial content
  test.Epi <- as.data.frame(table(cell.test[cell.test %in% Marker.list$Epithelial]))  
  test.Epi$Freq <- as.numeric(test.Epi$Freq)
  Epi.sum <- sum(test.Epi$Freq)
  Epi.genes <- as.character(test.Epi$Var1)
  #test EPCAM content
  test.StrongEpi <- as.data.frame(table(cell.test[cell.test %in% Marker.list$StrongEpi]))  
  test.StrongEpi$Freq <- as.numeric(test.StrongEpi$Freq)
  StrongEpi.sum <- sum(test.StrongEpi$Freq)
  StrongEpi.genes <- as.character(test.StrongEpi$Var1)
  
  
  
  SumTotal <- data.frame(
    NK.Total = NK.sum,
    NK.genes = paste(NK.genes, collapse = ", "),
    T.Total = T.sum,
    T.genes = paste(T.genes, collapse = ", "),
    B.Total = B.sum,
    B.genes = paste(B.genes, collapse = ", "),
    Plasma.Total = Plasma.sum,
    Plasma.genes = paste(Plasma.genes, collapse = ", "),
    Mye.Total = Mye.sum,
    Mye.genes = paste(Mye.genes, collapse = ", "),
    MyoEpi.Total = Myoepi.sum,
    MyoEpi.genes = paste(Myoepi.genes, collapse = ", "),
    PVL.Total = PVL.sum,
    PVL.genes = paste(PVL.genes, collapse = ", "),
    Fibro.Total = Fibro.sum,
    Fibro.genes = paste(Fibro.genes, collapse = ", "),
    Endo.Total = Endo.sum,
    Endo.genes = paste(Endo.genes, collapse = ", "),
    Immune.Total = Immune.sum,
    Immune.genes = paste(Immune.genes, collapse = ", "),
    Epi.Total = Epi.sum,
    Epi.genes = paste(Epi.genes, collapse = ", "),
    StrongEpi.Total = StrongEpi.sum,
    StrongEpi.genes = paste(StrongEpi.genes, collapse = ", ")
  )
  
  Totals.Only <- SumTotal[,c("NK.Total", "T.Total", "B.Total",
                             "Plasma.Total", "Mye.Total", "MyoEpi.Total", 
                             "PVL.Total", "Fibro.Total", "Endo.Total",
                             "Immune.Total", "Epi.Total", "StrongEpi.Total")]
  
  Totals.Only$tie <- apply(Totals.Only, 1, function(x) sum(x == max(x)) > 1)
  SumTotal$tie <- Totals.Only$tie
  markersum.list[[c]] <- SumTotal
  
  
}

head(markersum.list)

# ____________________________________________________________________________
# ____________________________________________________________________________
## Make decision about each cell type _________________________
## based on criteria from Karaayvaz (supplemental pg 36: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-06052-0/MediaObjects/41467_2018_6052_MOESM1_ESM.pdf)
## Karaayvaz's github: https://github.com/Michorlab/tnbc_scrnaseq/blob/master/code/funcs_markers.R
# ____________________________________________________________________________
# ____________________________________________________________________________


#these cells are for sure their cell type because 
#only have markers for exactly one cell type AND
#immune cells are PTPRC+/EPCAM-
#epi cells are EPCAM+/PTPRC-
#stromal cells are EPCAM-/PTPRC-

Extreme.Cell.Type.Testing <- function(df){
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Literally only has markers for one type  +++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  #when only has markers for that one cell type
  
  if (df$NK.Total >= 1 & df$Immune.Total == 0 &
      df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("NK Cells")
  
  #when only has markers for that one cell type
  if (df$T.Total >= 1 & df$Immune.Total == 0 &
      df$NK.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("T Cells")
  
  #when only has markers for that one cell type
  if (df$B.Total >= 1 & df$Immune.Total == 0 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("B Cells")
  
  #when only has markers for that one cell type
  if (df$Plasma.Total >= 1 & df$Immune.Total == 0 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Plasma Cells")
  
  #when only has markers for that one cell type
  if (df$Mye.Total >= 1 & df$Immune.Total == 0 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Myeloid Cells")
  
  #when only has markers for that one cell type
  if (df$MyoEpi.Total >= 1 & df$Immune.Total == 0 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 & 
      df$Mye.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Myoepithelial Cells")
  
  #when only has markers for that one cell type
  if (df$PVL.Total >= 1 & df$Immune.Total == 0 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 & 
      df$Mye.Total == 0 & df$MyoEpi.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Perivascular-like (PVL) Cells")
  
  #when only has markers for that one cell type
  if (df$Fibro.Total >= 1 & df$Immune.Total == 0 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 & 
      df$Mye.Total == 0 & df$MyoEpi.Total == 0 & 
      df$PVL.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Fibroblasts")
  
  #when only has markers for that one cell type
  if (df$Endo.Total >= 1 & df$Immune.Total == 0 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 & 
      df$Mye.Total == 0 & df$MyoEpi.Total == 0 & 
      df$PVL.Total == 0 & df$Fibro.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Endothelial")
  
  if (df$Epi.Total >= 0 & df$StrongEpi.Total >= 1 &
      df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$NK.Total == 0 & df$Immune.Total == 0)
    return("Epithelial Cells")
  
  if ((df$Epi.Total >= 1 & df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$NK.Total == 0 & df$Immune.Total == 0) |
      (df$StrongEpi.Total >= 1 &
      df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$NK.Total == 0 & df$Immune.Total == 0))
    return("Epithelial Cells")
  
  if (df$Epi.Total == 0 & df$T.Total == 0 & df$B.Total == 0 &
       df$Plasma.Total == 0 & df$Mye.Total == 0 & 
       df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
       df$Fibro.Total == 0 & df$Endo.Total == 0 & 
       df$NK.Total == 0 & df$Immune.Total == 0)
    return("Not extreme Case")
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # is immune cell (PTPRC and no EPCAM or stromal) +++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (df$Immune.Total == 1 && 
        #not epi
        df$Epi.Total == 0 & df$StrongEpi.Total == 0 &
        #not stromal
        df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
        df$Fibro.Total == 0 & df$Endo.Total == 0 &
        
        df$NK.Total >= 1 & df$B.Total == 0 &
        df$Plasma.Total == 0 & df$Mye.Total == 0 & 
        df$T.Total == 0) 
      return("NK Cells")
  
  if (df$Immune.Total == 1 && 
      #not epi
      df$Epi.Total == 0 & df$StrongEpi.Total == 0 &
      #not stromal
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 &
      df$T.Total >= 1 & df$B.Total == 0 &
         df$Plasma.Total == 0 & df$Mye.Total == 0 & 
         df$NK.Total == 0)
        return("T Cells")
  
  if (df$Immune.Total == 1 && 
      #not epi
      df$Epi.Total == 0 & df$StrongEpi.Total == 0 &
      #not stromal
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & df$B.Total >= 1 & df$T.Total == 0 &
         df$Plasma.Total == 0 & df$Mye.Total == 0 &
         df$NK.Total == 0)
        return("B Cells")
  
  if (df$Immune.Total == 1 && 
      #not epi
      df$Epi.Total == 0 & df$StrongEpi.Total == 0 &
      #not stromal
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & df$Plasma.Total >= 1 & df$T.Total == 0 &
         df$B.Total == 0 & df$Mye.Total == 0 & 
         df$NK.Total == 0)
        return("Plasma Cells")
  
  if (df$Immune.Total == 1 && 
      #not epi
      df$Epi.Total == 0 & df$StrongEpi.Total == 0 &
      #not stromal
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & df$Mye.Total >= 1 & df$T.Total == 0 &
         df$B.Total == 0 & df$Plasma.Total == 0 &
         df$NK.Total == 0)
        return("Myeloid Cells")
      
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # is Epithelial cell (no PTPRC or stromal) +++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ##one of "Strong Epi markers" : EPCAM, KRT8, KRT18, KRT19
  if ((df$Epi.Total == 0 & df$StrongEpi.Total >= 1 &
      #not immune
      df$Immune.Total == 0 & df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$NK.Total == 0 & 
      #not stromal
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0) | 
      
  ## at least 2 Epi markers
      (df$StrongEpi.Total == 0 & df$Epi.Total > 2 &
      #not immune
      df$Immune.Total == 0 &
      df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$NK.Total == 0 & 
      #not stromal
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0) | 
      
      (df$StrongEpi.Total >= 1 & df$Epi.Total >= 1 &
       #not immune
       df$Immune.Total == 0 &
       df$T.Total == 0 & df$B.Total == 0 &
       df$Plasma.Total == 0 & df$Mye.Total == 0 & 
       df$NK.Total == 0 & 
       #not stromal
       df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
       df$Fibro.Total == 0 & df$Endo.Total == 0)) 
    {
    return("Epithelial Cells")
    }
  
  if (df$StrongEpi.Total == 0 && df$MyoEpi.Total >= 1 &
      df$Immune.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & df$Epi.Total == 0) 
    {
    return("Myoepithelial Cells")
    }
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # is Stromal cell (no EPCAM and no PTPRC) +++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  if (df$Epi.Total == 0 && df$StrongEpi.Total == 0 &#not epi
      #not immune
      df$Immune.Total == 0 & df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$NK.Total == 0 &
      #not myoepithelial
      df$MyoEpi.Total == 0 & df$PVL.Total >= 1 & 
        df$Fibro.Total == 0 & df$Endo.Total == 0)
      return("Perivascular-like (PVL) Cells")
  
  if (df$Epi.Total == 0 && df$StrongEpi.Total == 0 &#not epi
      #not immune
      df$Immune.Total == 0 & df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$NK.Total == 0 &
      #not myoepithelial
      df$MyoEpi.Total == 0 &df$Fibro.Total >= 1 & 
        df$PVL.Total == 0 & df$Endo.Total == 0)
      return("Fibroblasts")
  
  if (df$Epi.Total == 0 && df$StrongEpi.Total == 0 &#not epi
      #not immune
      df$Immune.Total == 0 & df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$NK.Total == 0 &
      #not myoepithelial
      df$MyoEpi.Total == 0 & df$Endo.Total >= 1 & 
        df$PVL.Total == 0 & df$Fibro.Total == 0)
      return("Endothelial Cells")
  
  else
    return("Not extreme Case")
}


Extreme.Results <- sapply(markersum.list, Extreme.Cell.Type.Testing)
#see which ones are null #test: Xu_P4_GCGTGATGCCCC_P4_Prim
Extreme.Results[unlist(lapply(Extreme.Results,is.null))]

Extreme.dataframe <- as.data.frame(Extreme.Results)
head(Extreme.dataframe)
table(Extreme.dataframe$Extreme.Results)

# Unknown.list <- subset(Extreme.Results, subset = Extreme.Results == "Non-Specific Immune Cells")
# head(Unknown.list)
# Unknown.Cells <- rownames(Unknown.list)
# Unknown.Cells.list <- markersum.list[(names(markersum.list) %in% Unknown.Cells)]
# head(Unknown.Cells.list)

Not.Extreme <- subset(Extreme.dataframe, subset = Extreme.Results == "Not extreme Case")
table(Not.Extreme$Extreme.Results)

Not.Extreme.Cells <- rownames(Not.Extreme)

Notextreme.list <- markersum.list[(names(markersum.list) %in% Not.Extreme.Cells)]
head(Notextreme.list)

#these cells are weirder and have multiple markers for multiple cell types
## implementing Karaayvaz's cutoffs

NotExtreme.Cell.Type.Testing <- function(df){
  
  
  if (df$Immune.Total == 1 & df$NK.Total == 0 &
      df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Non-Specific Immune Cells")
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # NK refinement ++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # only one type with at least 2 markers (without PTPRC) --> that type
  if (df$NK.Total > 1 & df$Immune.Total == 0 &
      df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("NK Cells")
  # only one type with at least 1 marker, and PTPRC --> that type
  if (df$NK.Total >= 1 && df$Immune.Total == 1 &
      df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("NK Cells")
  
  # at least 3 markers of that type and at most 1 of another immune --> that type
  # ______________________________________________________________________________
  
  if (df$NK.Total > 2 && df$Immune.Total >= 0 & 
      df$B.Total == 1 & df$T.Total == 0 &
       df$Plasma.Total == 0 & df$Mye.Total == 0 & 
       df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
       df$Fibro.Total == 0 & df$Endo.Total == 0 & 
       df$Epi.Total == 0 & df$StrongEpi.Total == 0)
      return("NK Cells")

  if (df$NK.Total > 2 && df$Immune.Total >= 0 & 
      df$B.Total == 0 & df$T.Total == 1 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("NK Cells")
  
  if (df$NK.Total > 2 && df$Immune.Total >= 0 & 
      df$B.Total == 0 & df$T.Total == 0 &
      df$Plasma.Total == 1 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("NK Cells")
  
  if (df$NK.Total > 2 && df$Immune.Total >= 0 & 
      df$B.Total == 0 & df$T.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 1 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("NK Cells")
  
  # EPCAM+ NK Cells
  if (df$NK.Total >= 1 && df$Immune.Total >= 0 &
      df$T.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 &
      df$Fibro.Total == 0 & df$Endo.Total == 0 &
      df$Epi.Total == 0 & df$StrongEpi.Total == 1 & df$StrongEpi.genes == "EPCAM") 
      return("EPCAM+ NK Cells")

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # T refinement ++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # only one type with at least 2 markers (without PTPRC) --> that type
  if (df$T.Total > 1 & df$Immune.Total == 0 &
      df$NK.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("T Cells")
  # only one type with at least 1 marker, and PTPRC --> that type
  if (df$T.Total >= 1 && df$Immune.Total == 1 &
      df$NK.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("T Cells")
  
  # at least 3 markers of that type and at most 1 of another immune --> that type
  # ______________________________________________________________________________
  
  if (df$T.Total > 2 && df$Immune.Total >= 0 & 
      df$NK.Total == 1 & df$B.Total == 0 &
       df$Plasma.Total == 0 & df$Mye.Total == 0 & 
       df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
       df$Fibro.Total == 0 & df$Endo.Total == 0 & 
       df$Epi.Total == 0 & df$StrongEpi.Total == 0)
      return("T Cells")

  if (df$T.Total > 2 && df$Immune.Total >= 0 & 
      df$NK.Total == 0 & df$B.Total == 1 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("T Cells")
  
  if (df$T.Total > 2 && df$Immune.Total >= 0 & 
      df$NK.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 1 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("T Cells")
  
  if (df$T.Total > 2 && df$Immune.Total >= 0 & 
      df$NK.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 1 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("T Cells")
  
  # EPCAM+ T Cells
  if (df$T.Total >= 1 && df$Immune.Total >= 0 &
      df$NK.Total == 0 & df$B.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 &
      df$Fibro.Total == 0 & df$Endo.Total == 0 &
      df$Epi.Total == 0 & df$StrongEpi.Total == 1 & df$StrongEpi.genes == "EPCAM") 
    
      return("EPCAM+ T Cells")
      
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # B refinement ++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # only one type with at least 2 markers (without PTPRC) --> that type
  if (df$B.Total > 1 & df$Immune.Total == 0 &
      df$T.Total == 0 & df$NK.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("B Cells")
  # only one type with at least 1 marker, and PTPRC --> that type
  if (df$B.Total >= 1 && df$Immune.Total == 1 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("B Cells")
  
  # at least 3 markers of that type and at most 1 of another immune --> that type
  # _____________________________________________________________________________
  
  if (df$B.Total > 2 && df$Immune.Total >= 0 &
      df$NK.Total == 1 & df$T.Total == 0 &
       df$Plasma.Total == 0 & df$Mye.Total == 0 & 
       df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
       df$Fibro.Total == 0 & df$Endo.Total == 0 & 
       df$Epi.Total == 0 & df$StrongEpi.Total == 0)
      return("B Cells")

  if (df$B.Total > 2 && df$Immune.Total >= 0 &
      df$NK.Total == 0 & df$T.Total == 1 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("B Cells")
  
  if (df$B.Total > 2 && df$Immune.Total >= 0 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$Plasma.Total == 1 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("B Cells")
  
  if (df$B.Total > 2 && df$Immune.Total >= 0 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 1 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("B Cells")
  
  
  # EPCAM+ B Cells
  if (df$B.Total >= 1 && df$Immune.Total >= 0 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$Plasma.Total == 0 & df$Mye.Total == 0 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 &
      df$Fibro.Total == 0 & df$Endo.Total == 0 &
      df$Epi.Total == 0 & df$StrongEpi.Total == 1 & df$StrongEpi.genes == "EPCAM") 
      return("EPCAM+ B Cells")

  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Plasma refinement ++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # only one type with at least 2 markers (without PTPRC) --> that type
  if (df$Plasma.Total > 1 & df$Immune.Total == 0 &
      df$T.Total == 0 & df$NK.Total == 0 &
      df$B.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Plasma Cells")
  
  # only one type with at least 1 marker, and PTPRC --> that type
  if (df$Plasma.Total >= 1 && df$Immune.Total == 1 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Mye.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Plasma Cells")
  
  # at least 3 markers of that type and at most 1 of another immune --> that type
  #_______________________________________________________________________________
  
  if (df$Plasma.Total > 2 && df$Immune.Total >= 0 & 
      df$NK.Total == 1 & df$T.Total == 0 &
       df$B.Total == 0 & df$Mye.Total == 0 &
       df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
       df$Fibro.Total == 0 & df$Endo.Total == 0 & 
       df$Epi.Total == 0 & df$StrongEpi.Total == 0)
      return("Plasma Cells")

  if (df$Plasma.Total > 2 && df$Immune.Total >= 0 & 
      df$NK.Total == 0 & df$T.Total == 1 &
      df$B.Total == 0 & df$Mye.Total == 0 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Plasma Cells")
  
  if (df$Plasma.Total > 2 && df$Immune.Total >= 0 & 
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 1 & df$Mye.Total == 0 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Plasma Cells")
  
  if (df$Plasma.Total > 2 && df$Immune.Total >= 0 & 
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Mye.Total == 1 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Plasma Cells")
  

  # EPCAM+ Plasma Cells
  if (df$Plasma.Total >= 1 && df$Immune.Total >= 0 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Mye.Total == 0 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 &
      df$Fibro.Total == 0 & df$Endo.Total == 0 &
      df$Epi.Total == 0 & df$StrongEpi.Total == 1 & df$StrongEpi.genes == "EPCAM") 
      return("EPCAM+ Plasma Cells")

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Myeloid refinement ++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # only one type with at least 2 markers (without PTPRC) --> that type
  if (df$Mye.Total > 1 & df$Immune.Total == 0 &
      df$T.Total == 0 & df$NK.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Myeloid Cells")
  
  # only one type with at least 1 marker, and PTPRC --> that type
  if (df$Mye.Total >= 1 && df$Immune.Total == 1 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 & 
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Myeloid Cells")
  
  # at least 3 markers of that type and at most 1 of another immune --> that type
  # __________________________________________________________________________-
  
  if (df$Mye.Total > 2 && df$Immune.Total >= 0 & 
      df$NK.Total == 1 & df$T.Total == 0 &
       df$B.Total == 0 & df$Plasma.Total == 0 &
       df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
       df$Fibro.Total == 0 & df$Endo.Total == 0 & 
       df$Epi.Total == 0 & df$StrongEpi.Total == 0)
      return("Myeloid Cells")

  if (df$Mye.Total > 2 && df$Immune.Total >= 0 & 
      df$NK.Total == 0 & df$T.Total == 1 &
      df$B.Total == 0 & df$Plasma.Total == 0 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Myeloid Cells")
  
  if (df$Mye.Total > 2 && df$Immune.Total >= 0 & 
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 1 & df$Plasma.Total == 0 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Myeloid Cells")
  
  if (df$Mye.Total > 2 && df$Immune.Total >= 0 & 
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 1 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 & 
      df$Fibro.Total == 0 & df$Endo.Total == 0 & 
      df$Epi.Total == 0 & df$StrongEpi.Total == 0)
    return("Myeloid Cells")
  
  
  # EPCAM+ Myeloid Cells
  if (df$Mye.Total >= 1 && df$Immune.Total >= 0 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 &
      df$Fibro.Total == 0 & df$Endo.Total == 0 &
      df$Epi.Total == 0 & df$StrongEpi.Total == 1 & df$StrongEpi.genes == "EPCAM") 
    return("EPCAM+ Myeloid Cells")

  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Myoepi refinement ++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # PTPRC+ MyoEpi Cells
  if ((df$MyoEpi.Total >= 1 && df$Immune.Total == 0 && df$Epi.Total == 0 & 
       df$StrongEpi.Total == 1 & df$StrongEpi.genes == "EPCAM" &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 &
      df$Mye.Total == 0 & df$PVL.Total == 0 &
      df$Fibro.Total == 0 & df$Endo.Total == 0 &
      df$Immune.Total == 1) | 
      
      (df$MyoEpi.Total >= 1 && df$Immune.Total == 0 && df$Epi.Total == 0 & df$StrongEpi.Total == 0 &
       df$NK.Total == 0 & df$T.Total == 0 &
       df$B.Total == 0 & df$Plasma.Total == 0 &
       df$Mye.Total == 0 & df$PVL.Total == 0 &
       df$Fibro.Total == 0 & df$Endo.Total == 0 &
       df$Immune.Total == 1)) 
      return("PTPRC+ Myoepithelial Cells")

  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # EMT cells  +++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # at least 2 epi markers
  if ((df$Epi.Total > 1 & df$Immune.Total == 0 & df$StrongEpi.Total >= 0) | 
      (df$StrongEpi.Total > 1 & df$Immune.Total == 0 & df$Epi.Total >= 0))
    return("Epithelial Cells")
  
  #epi with strom markers
  if (((df$Epi.Total > 1 | df$StrongEpi.Total > 1) & df$MyoEpi.Total >= 1 & 
       #no immune markers
       df$Immune.Total == 0 & df$NK.Total == 0 & 
       df$T.Total == 0 & df$B.Total == 0 & 
       df$Plasma.Total == 0 & df$Mye.Total == 0) |
      
      ((df$Epi.Total >= 1 | df$StrongEpi.Total >= 1) & df$PVL.Total >= 1 & 
       #no immune markers
       df$Immune.Total == 0 & df$NK.Total == 0 & 
       df$T.Total == 0 & df$B.Total == 0 & 
       df$Plasma.Total == 0 & df$Mye.Total == 0) |
      
      ((df$Epi.Total >= 1 | df$StrongEpi.Total >= 1) & df$Fibro.Total >= 1 &
       #no immune markers
       df$Immune.Total == 0 & df$NK.Total == 0 & 
       df$T.Total == 0 & df$B.Total == 0 & 
       df$Plasma.Total == 0 & df$Mye.Total == 0) |
      
      ((df$Epi.Total >= 1 | df$StrongEpi.Total >= 1) & df$Endo.Total >= 1 &
       #no immune markers
       df$Immune.Total == 0 & df$NK.Total == 0 & 
       df$T.Total == 0 & df$B.Total == 0 & 
       df$Plasma.Total == 0 & df$Mye.Total == 0)) 
    return("Epithelial with Stromal Markers (EMT)")
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # PVL refinement ++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # at least some stroma and at most 1 endothelial
  if (df$PVL.Total >= 1 && df$Endo.Total <= 1 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 &
      df$Mye.Total == 0 & df$MyoEpi.Total == 0 &
      df$Fibro.Total == 0 & df$Epi.Total == 0 & df$StrongEpi.Total == 0 & df$Immune.Total == 0)
       return("Perivascular-like (PVL) Cells")
  
  
  # PTPRC+ PVL Cells
  if (df$PVL.Total >= 1 && df$Immune.Total == 1 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 &
      df$MyoEpi.Total == 0 & df$Mye.Total == 0 &
      df$Fibro.Total == 0 & df$Endo.Total == 0 &
      df$Epi.Total == 0 & df$StrongEpi.Total == 0) 
      return("PTPRC+ PVL Cells")

  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Fibro refinement ++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  # at least some stroma and at most 1 endothelial
  if (df$Fibro.Total >= 1 && df$Endo.Total <= 1 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 &
      df$Mye.Total == 0 & df$MyoEpi.Total == 0 &
      df$PVL.Total == 0 & df$Epi.Total == 0 & df$StrongEpi.Total == 0 & df$Immune.Total == 0)
    return("Fibroblasts")
  
  
  # PTPRC+ Fibroblasts
  if (df$Fibro.Total >= 1 && df$Immune.Total == 1 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 &
      df$Mye.Total == 0 & df$Endo.Total == 0 &
      df$Epi.Total == 0 & df$StrongEpi.Total == 0) 
    return("PTPRC+ Fibroblasts")
  
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Endothelial refinement ++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # at least some endo and at most 1 stroma
  if (df$Endo.Total >= 1 && 
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 &
      df$Mye.Total == 0 & df$Epi.Total == 0 & df$StrongEpi.Total == 0 & df$Immune.Total == 0 & 
      
      ((df$MyoEpi.Total <= 1 & !(df$PVL.Total <= 1 | df$Fibro.Total <= 1)) | 
         (df$PVL.Total <= 1 & !(df$MyoEpi.Total <= 1 | df$Fibro.Total <= 1)) | 
         (df$Fibro.Total <= 1 & !(df$MyoEpi.Total <= 1 | df$PVL.Total <= 1))))
        return("Endothelial Cells")
    
    
  

  # PTPRC+ Endothelial
  if (df$Endo.Total >= 1 && df$Immune.Total == 1 &
      df$NK.Total == 0 & df$T.Total == 0 &
      df$B.Total == 0 & df$Plasma.Total == 0 &
      df$MyoEpi.Total == 0 & df$PVL.Total == 0 &
      df$Mye.Total == 0 & df$Fibro.Total == 0 &
      df$Epi.Total == 0 & df$StrongEpi.Total == 0) 
    return("PTPRC+ Endothelial Cells")
  
  else
    return("Unknown")
}


Notextreme.Results <- sapply(Notextreme.list, NotExtreme.Cell.Type.Testing)
Notextreme.dataframe <- as.data.frame(Notextreme.Results)
head(Notextreme.dataframe)
table(Notextreme.dataframe)  

  
Unknown.list <- subset(Notextreme.dataframe, subset = Notextreme.Results == "Unknown")
head(Unknown.list)

Unknown.Cells <- rownames(Unknown.list)
Unknown.Cells.list <- markersum.list[(names(markersum.list) %in% Unknown.Cells)]
head(Unknown.Cells.list)


for (cell in names(Unknown.Cells.list)) {
  
  extract <- markersum.list[[cell]]
  
  #extract the genes that were expressed
  test.NK <- unlist(strsplit(extract$NK.genes,split=', ',fixed=TRUE))
  test.T <- unlist(strsplit(extract$T.genes,split=', ',fixed=TRUE))
  test.B <- unlist(strsplit(extract$B.genes,split=', ',fixed=TRUE))
  test.Plasma <- unlist(strsplit(extract$Plasma.genes,split=', ',fixed=TRUE))
  test.Mye <- unlist(strsplit(extract$Mye.genes,split=', ',fixed=TRUE))
  test.MyoEpi <- unlist(strsplit(extract$MyoEpi.genes,split=', ',fixed=TRUE))
  test.PVL <- unlist(strsplit(extract$PVL.genes,split=', ',fixed=TRUE))
  test.Fibro <- unlist(strsplit(extract$Fibro.genes,split=', ',fixed=TRUE))
  test.Endo <- unlist(strsplit(extract$Endo.genes,split=', ',fixed=TRUE))
  test.Epi <- unlist(strsplit(extract$Epi.genes,split=', ',fixed=TRUE))
  test.StrongEpi <- unlist(strsplit(extract$StrongEpi.genes,split=', ',fixed=TRUE))
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #get the expression matrix for this cell and those genes ++++++++++++++
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  #NK Expression Sum ++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  NK.test <- as.data.frame(Subset.Mat[c(test.NK, "PTPRC"), cell])
  colnames(NK.test) <- cell
  rownames(NK.test) <- c(test.NK, "PTPRC")
  NK.exp.sum <- colSums(NK.test)
  NK.exp.avg <- NK.exp.sum/69059  #average expression (denominator is #rows of dataset)
  #NK.exp.avg <- colMeans(NK.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #T Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Tcell.test <- as.data.frame(Subset.Mat[c(test.T, "PTPRC"), cell])
  colnames(Tcell.test) <- cell
  rownames(Tcell.test) <- c(test.T, "PTPRC")
  T.exp.sum <- colSums(Tcell.test)
  T.exp.avg <- T.exp.sum/69059  #average expression (denominator is #rows of dataset)
  #T.exp.avg <- colMeans(Tcell.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #B Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Bcell.test <- as.data.frame(Subset.Mat[c(test.B, "PTPRC"), cell])
  colnames(Bcell.test) <- cell
  rownames(Bcell.test) <- c(test.B, "PTPRC")
  B.exp.sum <- colSums(Bcell.test)
  B.exp.avg <- B.exp.sum/69059  #average expression (denominator is #rows of dataset)
  #B.exp.avg <- colMeans(Bcell.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #Plasma Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  Plasmacell.test <- as.data.frame(Subset.Mat[c(test.Plasma, "PTPRC"), cell])
  colnames(Plasmacell.test) <- cell
  rownames(Plasmacell.test) <- c(test.Plasma, "PTPRC")
  Plasma.exp.sum <- colSums(Plasmacell.test)
  Plasma.exp.avg <- Plasma.exp.sum/69059  #average expression (denominator is #rows of dataset)
  #Plasma.exp.avg <- colMeans(Plasmacell.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #Myeloid Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Myeloid.test <- as.data.frame(Subset.Mat[c(test.Mye, "PTPRC"), cell])
  colnames(Myeloid.test) <- cell
  rownames(Myeloid.test) <- c(test.Mye, "PTPRC")
  Myeloid.exp.sum <- colSums(Myeloid.test)
  Myeloid.exp.avg <- Myeloid.exp.sum/69059 #average expression (denominator is #rows of dataset)
  #Myeloid.exp.avg <- colMeans(Myeloid.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #Myoepi Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  MyoEpi.test <- as.data.frame(Subset.Mat[c(test.MyoEpi), cell])
  colnames(MyoEpi.test) <- cell
  rownames(MyoEpi.test) <- c(test.MyoEpi)
  MyoEpi.exp.sum <- colSums(MyoEpi.test)
  MyoEpi.exp.avg <- MyoEpi.exp.sum/69059 #average expression (denominator is #rows of dataset)
  #MyoEpi.exp.avg <- colMeans(MyoEpi.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #PVL Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  PVL.test <- as.data.frame(Subset.Mat[c(test.PVL), cell])
  colnames(PVL.test) <- cell
  rownames(PVL.test) <- c(test.PVL)
  PVL.exp.sum <- colSums(PVL.test)
  PVL.exp.avg <- PVL.exp.sum/69059 #average expression (denominator is #rows of dataset)
  #PVL.exp.avg <- colMeans(PVL.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  #Fibro Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Fibro.test <- as.data.frame(Subset.Mat[c(test.Fibro), cell])
  colnames(Fibro.test) <- cell
  rownames(Fibro.test) <- c(test.Fibro)
  Fibro.exp.sum <- colSums(Fibro.test)
  Fibro.exp.avg <- Fibro.exp.sum/69059 #average expression (denominator is #rows of dataset)
  #Fibro.exp.avg <- colMeans(Fibro.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  #Endothelial Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Endo.test <- as.data.frame(Subset.Mat[c(test.Endo), cell])
  colnames(Endo.test) <- cell
  rownames(Endo.test) <- c(test.Endo)
  Endo.exp.sum <- colSums(Endo.test)
  Endo.exp.avg <- Endo.exp.sum/69059 #average expression (denominator is #rows of dataset)
  #Endo.exp.avg <- colMeans(Endo.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  #Epithelial Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Epi.test <- as.data.frame(Subset.Mat[c(test.Epi, test.StrongEpi), cell])
  colnames(Epi.test) <- cell
  rownames(Epi.test) <- c(test.Epi, test.StrongEpi)
  Epi.exp.sum <- colSums(Epi.test)
  Epi.exp.avg <- Epi.exp.sum/69059 #average expression (denominator is #rows of dataset)
  #Epi.exp.avg <- colMeans(Epi.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Final Call ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  if ((NK.exp.avg > T.exp.avg) & (NK.exp.avg > B.exp.avg) &
      (NK.exp.avg > Plasma.exp.avg) & (NK.exp.avg > Myeloid.exp.avg) &
      (NK.exp.avg > MyoEpi.exp.avg) & (NK.exp.avg > PVL.exp.avg) & 
      (NK.exp.avg > Fibro.exp.avg) & (NK.exp.avg > Endo.exp.avg) &
      (NK.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "NK Cells"
    
    
  }
  
  else if ((T.exp.avg > NK.exp.avg) & (T.exp.avg > B.exp.avg) &
           (T.exp.avg > Plasma.exp.avg) & (T.exp.avg > Myeloid.exp.avg) &
           (T.exp.avg > MyoEpi.exp.avg) & (T.exp.avg > PVL.exp.avg) & 
           (T.exp.avg > Fibro.exp.avg) & (T.exp.avg > Endo.exp.avg) &
           (T.exp.avg > Epi.exp.avg)) {
    
    
    if(("CD8A" %in% rownames(Subset.Mat)) | ("CD8B" %in% rownames(Subset.Mat)))
    {
      CD8T.test <- as.data.frame(Subset.Mat[c("CD8A","CD8B"), cell])
      colnames(CD8T.test) <- cell
      rownames(CD8T.test) <-c("CD8A","CD8B")
      CD8.exp.sum <- colSums(CD8T.test)
      CD8.exp.avg <- CD8.exp.sum/69059 
      #CD8.exp.avg <- colMeans(CD8T.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      CD8.exp.avg = 0
    }
    
    if("CD4" %in% rownames(Subset.Mat))
    {
      CD4.test <- as.data.frame(Subset.Mat[c("CD4"), cell])
      colnames(CD4.test) <- cell
      rownames(CD4.test) <-c("CD4")
      CD4.exp.sum <- colSums(CD4.test)
      CD4.exp.avg <- CD4.exp.sum/69059 
      #CD4.exp.avg <- colMeans(CD4.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }  
    else {
      CD4.exp.avg = 0
    }
    if("FOXP3" %in% rownames(Subset.Mat))
    {
      Treg.test <- as.data.frame(Subset.Mat[c("FOXP3"), cell])
      colnames(Treg.test) <- cell
      rownames(Treg.test) <-c("FOXP3")
      Treg.exp.sum <- colSums(Treg.test)
      Treg.exp.avg <- Treg.exp.sum/69059 
      #Treg.exp.avg <- colMeans(Treg.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }  
    else {
      Treg.exp.avg = 0
    }
    
    if ((CD8.exp.avg > CD4.exp.avg) & (CD8.exp.avg > Treg.exp.avg))
    {
      cell.result[[cell]] <- "CD8+ T Cells"
    }
    
    else if ((CD4.exp.avg > CD8.exp.avg) & (CD4.exp.avg > Treg.exp.avg))
    {
      cell.result[[cell]] <- "CD4+ T Cells"
    }
    
    else if ((Treg.exp.avg > CD8.exp.avg) & (Treg.exp.avg > CD4.exp.avg))
    {
      cell.result[[cell]] <- "Regulatory T Cells"
    }
    
    else 
    {
      cell.result[[cell]] <- "Unspecific T Cells"
    }
    
  }
  
  else if ((B.exp.avg > NK.exp.avg) & (B.exp.avg > T.exp.avg) &
           (B.exp.avg > Plasma.exp.avg) & (B.exp.avg > Myeloid.exp.avg) &
           (B.exp.avg > MyoEpi.exp.avg) & (B.exp.avg > PVL.exp.avg) & 
           (B.exp.avg > Fibro.exp.avg) & (B.exp.avg > Endo.exp.avg) &
           (B.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "B Cells"
    
    
  }
  
  else if ((Plasma.exp.avg > NK.exp.avg) & (Plasma.exp.avg > T.exp.avg) &
           (Plasma.exp.avg > B.exp.avg) & (Plasma.exp.avg > Myeloid.exp.avg) &
           (Plasma.exp.avg > MyoEpi.exp.avg) & (Plasma.exp.avg > PVL.exp.avg) & 
           (Plasma.exp.avg > Fibro.exp.avg) & (Plasma.exp.avg > Endo.exp.avg) &
           (Plasma.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "Plasma Cells"
    
    
  }
  
  else if ((Myeloid.exp.avg > NK.exp.avg) & (Myeloid.exp.avg > T.exp.avg) &
           (Myeloid.exp.avg > B.exp.avg) & (Myeloid.exp.avg > Plasma.exp.avg) &
           (Myeloid.exp.avg > MyoEpi.exp.avg) & (Myeloid.exp.avg > PVL.exp.avg) & 
           (Myeloid.exp.avg > Fibro.exp.avg) & (Myeloid.exp.avg > Endo.exp.avg) &
           (Myeloid.exp.avg > Epi.exp.avg)) {
    
    
    if(("CD33" %in% rownames(Subset.Mat)) | ("ITGAM" %in% rownames(Subset.Mat)) | 
       ("FUT4" %in% rownames(Subset.Mat)) | ("CEACAM8" %in% rownames(Subset.Mat)) |
       ("IL4R" %in% rownames(Subset.Mat)) | ("HLA-DR" %in% rownames(Subset.Mat)))
    {
      MDSc.test <- as.data.frame(Subset.Mat[c("CD33", #https://www.nature.com/articles/ncomms12150; table 1
                                              "ITGAM", 
                                              "FUT4", #CD15
                                              "CEACAM8",
                                              "IL4R",
                                              "HLA-DRA"), cell])
      colnames(MDSc.test) <- cell
      rownames(MDSc.test) <-c("CD33", #https://www.nature.com/articles/ncomms12150; table 1
                              "ITGAM", 
                              "FUT4", #CD15
                              "CEACAM8",
                              "IL4R",
                              "HLA-DRA")
      MDSc.exp.sum <- colSums(MDSc.test)
      MDSc.exp.avg <- MDSc.exp.sum/69059 
      #MDSc.exp.avg <- colMeans(MDSc.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      MDSc.exp.avg = 0
    }
    if(("INHBA" %in% rownames(Subset.Mat)) | ("IL1RN" %in% rownames(Subset.Mat)) | 
       ("CCL4" %in% rownames(Subset.Mat)) | ("NLRP3" %in% rownames(Subset.Mat)) |
       ("EREG" %in% rownames(Subset.Mat)) | ("IL1B" %in% rownames(Subset.Mat)) |
       ("LYVE1" %in% rownames(Subset.Mat)) |  ("PLTP" %in% rownames(Subset.Mat)) |
       ("SELENOP" %in% rownames(Subset.Mat)) |  ("C1QC" %in% rownames(Subset.Mat)) |
       ("C1QA" %in% rownames(Subset.Mat)) |  ("APOE" %in% rownames(Subset.Mat)))
    {
      Macro.test <- as.data.frame(Subset.Mat[c("INHBA", "IL1RN", "CCL4", "NLRP3",
                                               "EREG", "IL1B", "LYVE1", "PLTP",
                                               "SELENOP", "C1QC", "C1QA", "APOE"), cell])
      colnames(Macro.test) <- cell
      rownames(Macro.test) <-c("INHBA", "IL1RN", "CCL4", "NLRP3",
                               "EREG", "IL1B", "LYVE1", "PLTP",
                               "SELENOP", "C1QC", "C1QA", "APOE")
      Macro.exp.sum <- colSums(Macro.test)
      Macro.exp.avg <- Macro.exp.sum/69059 
      #Macro.exp.avg <- colMeans(Macro.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      Macro.exp.avg = 0
    }
    
    if(("FCN1" %in% rownames(Subset.Mat)) | ("S100A9" %in% rownames(Subset.Mat)) | 
       ("S100A8" %in% rownames(Subset.Mat)) | ("FCGR3A" %in% rownames(Subset.Mat)) |
       ("LST1" %in% rownames(Subset.Mat)) | ("LILRB2" %in% rownames(Subset.Mat)))
    {
      Mono.test <- as.data.frame(Subset.Mat[c("FCN1", "S100A9", "S100A8", "FCGR3A", "LST1", "LILRB2"), cell])
      colnames(Mono.test) <- cell
      rownames(Mono.test) <-c("FCN1", "S100A9", "S100A8", "FCGR3A", "LST1", "LILRB2")
      Mono.exp.sum <- colSums(Mono.test)
      Mono.exp.avg <- Mono.exp.sum/69059 
      #Mono.exp.avg <- colMeans(Mono.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      Mono.exp.avg = 0
    }
    if(("CXCR1" %in% rownames(Subset.Mat)) | ("FUT4" %in% rownames(Subset.Mat)) | 
       ("FCGR3A" %in% rownames(Subset.Mat)) | ("CSF3R" %in% rownames(Subset.Mat)) |
       ("S100A9" %in% rownames(Subset.Mat)) | ("TNF" %in% rownames(Subset.Mat)) |
       ("CD274" %in% rownames(Subset.Mat)))
    {
      Neut.test <- as.data.frame(Subset.Mat[c("CXCR1", #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
                                              "FUT4", #https://www.rndsystems.com/resources/cell-markers/immune-cells/granulocytes/neutrophil-cell-markers
                                              "FCGR3A", #CD16
                                              "CSF3R", #TCGA, Chung
                                              "S100A9","TNF","CD274"), cell])
      colnames(Neut.test) <- cell
      rownames(Neut.test) <-c("CXCR1", #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
                              "FUT4", #https://www.rndsystems.com/resources/cell-markers/immune-cells/granulocytes/neutrophil-cell-markers
                              "FCGR3A", #CD16
                              "CSF3R", #TCGA, Chung
                              "S100A9","TNF","CD274")
      Neut.exp.sum <- colSums(Neut.test)
      Neut.exp.avg <- Neut.exp.sum/69059 
      #Neut.exp.avg <- colMeans(Neut.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      Neut.exp.avg = 0
    }
    if(("LILRA4" %in% rownames(Subset.Mat)) | ("GZMB" %in% rownames(Subset.Mat)) | 
       ("IL3RA" %in% rownames(Subset.Mat)) | ("CLEC9A" %in% rownames(Subset.Mat)) |
       ("FLT3" %in% rownames(Subset.Mat)) | ("IDO1" %in% rownames(Subset.Mat)) |
       ("CD1C" %in% rownames(Subset.Mat)) | ("FCER1A" %in% rownames(Subset.Mat)) |
       ("HLA-DQA1" %in% rownames(Subset.Mat)) | ("LAMP3" %in% rownames(Subset.Mat)) |
       ("CCR7" %in% rownames(Subset.Mat)) | ("FSCN1" %in% rownames(Subset.Mat)))
    {
      DC.test <- as.data.frame(Subset.Mat[c("LILRA4","GZMB","IL3RA",
                                            "CLEC9A", "FLT3", "IDO1", 
                                            "CD1C", "FCER1A","HLA-DQA1", 
                                            "LAMP3", "CCR7", "FSCN1"), cell])
      colnames(DC.test) <- cell
      rownames(DC.test) <-c("LILRA4","GZMB","IL3RA",
                            "CLEC9A", "FLT3", "IDO1", 
                            "CD1C", "FCER1A","HLA-DQA1", 
                            "LAMP3", "CCR7", "FSCN1")
      DC.exp.sum <- colSums(DC.test)
      DC.exp.avg <- DC.exp.sum/69059 
      #DC.exp.avg <- colMeans(DC.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      DC.exp.avg = 0
    }
    if(("KIT" %in% rownames(Subset.Mat)) | ("TPSAB1" %in% rownames(Subset.Mat)) | 
       ("CPA4" %in% rownames(Subset.Mat)))
    {
      Mast.test <- as.data.frame(Subset.Mat[c("KIT","TPSAB1","CPA4"), cell])
      colnames(Mast.test) <- cell
      rownames(Mast.test) <-c("KIT","TPSAB1","CPA4")
      Mast.exp.sum <- colSums(Mast.test)
      Mast.exp.avg <- Mast.exp.sum/69059 
      #Mast.exp.avg <- colMeans(Mast.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      Mast.exp.avg = 0
    }
    
    
    if ((MDSc.exp.avg > Macro.exp.avg) & (MDSc.exp.avg > Mono.exp.avg) & 
        (MDSc.exp.avg > Neut.exp.avg) & (MDSc.exp.avg > DC.exp.avg) &
        (MDSc.exp.avg > Mast.exp.avg))
    {
      cell.result[[cell]] <- "MDSCs"
    }
    
    else if ((Macro.exp.avg > MDSc.exp.avg) & (Macro.exp.avg > Mono.exp.avg) & 
             (Macro.exp.avg > Neut.exp.avg) & (Macro.exp.avg > DC.exp.avg) &
             (Macro.exp.avg > Mast.exp.avg))
    {
      cell.result[[cell]] <- "Macrophages"
    }
    
    else if ((Mono.exp.avg > MDSc.exp.avg) & (Mono.exp.avg > Macro.exp.avg) & 
             (Mono.exp.avg > Neut.exp.avg) & (Mono.exp.avg > DC.exp.avg) &
             (Mono.exp.avg > Mast.exp.avg))
    {
      cell.result[[cell]] <- "Monocytes"
    }
    
    else if ((Neut.exp.avg > MDSc.exp.avg) & (Neut.exp.avg > Macro.exp.avg) & 
             (Neut.exp.avg > Mono.exp.avg) & (Neut.exp.avg > DC.exp.avg) &
             (Neut.exp.avg > Mast.exp.avg))
    {
      cell.result[[cell]] <- "Neutrophils"
    }
    
    else if ((DC.exp.avg > MDSc.exp.avg) & (DC.exp.avg > Macro.exp.avg) & 
             (DC.exp.avg > Mono.exp.avg) & (DC.exp.avg > Neut.exp.avg) &
             (DC.exp.avg > Mast.exp.avg))
    {
      cell.result[[cell]] <- "Dendritic Cells"
    }
    
    else if ((Mast.exp.avg > MDSc.exp.avg) & (Mast.exp.avg > Macro.exp.avg) & 
             (Mast.exp.avg > Mono.exp.avg) & (Mast.exp.avg > Neut.exp.avg) &
             (Mast.exp.avg > DC.exp.avg))
    {
      cell.result[[cell]] <- "Mast Cells"
    }
    
    else 
    {
      cell.result[[cell]] <- "Unspecific Myeloid Cells"
    }
    
    
  }
  
  else if ((MyoEpi.exp.avg > NK.exp.avg) & (MyoEpi.exp.avg > T.exp.avg) &
           (MyoEpi.exp.avg > B.exp.avg) & (MyoEpi.exp.avg > Plasma.exp.avg) &
           (MyoEpi.exp.avg > Myeloid.exp.avg) & (MyoEpi.exp.avg > PVL.exp.avg) & 
           (MyoEpi.exp.avg > Fibro.exp.avg) & (MyoEpi.exp.avg > Endo.exp.avg) &
           (MyoEpi.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "Myoepithelial Cells"
    
    
  }
  
  else if ((PVL.exp.avg > NK.exp.avg) & (PVL.exp.avg > T.exp.avg) &
           (PVL.exp.avg > B.exp.avg) & (PVL.exp.avg > Plasma.exp.avg) &
           (PVL.exp.avg > Myeloid.exp.avg) & (PVL.exp.avg > MyoEpi.exp.avg) & 
           (PVL.exp.avg > Fibro.exp.avg) & (PVL.exp.avg > Endo.exp.avg) &
           (PVL.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "Perivascular-like (PVL) Cells"
    
    
  }
  
  else if ((Fibro.exp.avg > NK.exp.avg) & (Fibro.exp.avg > T.exp.avg) &
           (Fibro.exp.avg > B.exp.avg) & (Fibro.exp.avg > Plasma.exp.avg) &
           (Fibro.exp.avg > Myeloid.exp.avg) & (Fibro.exp.avg > MyoEpi.exp.avg) & 
           (Fibro.exp.avg > PVL.exp.avg) & (Fibro.exp.avg > Endo.exp.avg) &
           (Fibro.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "Fibroblasts"
    
    
  }
  
  else if ((Endo.exp.avg > NK.exp.avg) & (Endo.exp.avg > T.exp.avg) &
           (Endo.exp.avg > B.exp.avg) & (Endo.exp.avg > Plasma.exp.avg) &
           (Endo.exp.avg > Myeloid.exp.avg) & (Endo.exp.avg > MyoEpi.exp.avg) & 
           (Endo.exp.avg > PVL.exp.avg) & (Endo.exp.avg > Fibro.exp.avg) &
           (Endo.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "Endothelial Cells"
    
    
  }
  
  else if ((Epi.exp.avg > NK.exp.avg) & (Epi.exp.avg > T.exp.avg) &
           (Epi.exp.avg > B.exp.avg) & (Epi.exp.avg > Plasma.exp.avg) &
           (Epi.exp.avg > Myeloid.exp.avg) & (Epi.exp.avg > MyoEpi.exp.avg) & 
           (Epi.exp.avg > PVL.exp.avg) & (Epi.exp.avg > Fibro.exp.avg) &
           (Epi.exp.avg > Endo.exp.avg)) {
    
    cell.result[[cell]] <- "Epithelial Cells"
    
    
  }
  else 
    cell.result[[cell]] <- "Unknown"
  
  
}


#for T cells and Myeloid cells, split into subsets 
head(cell.result)
UnknownExp.Results <- t(as.data.frame(cell.result))
head(UnknownExp.Results)
table(UnknownExp.Results) 

Extreme.Final <- Extreme.Results[Extreme.Results != "Not extreme Case"]
NotExtreme.Final <- Notextreme.Results[Notextreme.Results != "Unknown"]

MarkerBased.Call.List <- c(Extreme.Final, NotExtreme.Final, zero.for.everything, cell.result)
head(MarkerBased.Call.List)

MarkerBased.dataframe <- as.data.frame(t(as.data.frame(MarkerBased.Call.List)))
head(MarkerBased.dataframe)
table(MarkerBased.dataframe$V1)


cell.result <- list()
#test cell: Aziziimmune_BC3_48068
for (cell in names(markersum.list)) {
  
  extract <- markersum.list[[cell]]
  
  #extract the genes that were expressed
  test.NK <- unlist(strsplit(extract$NK.genes,split=', ',fixed=TRUE))
  test.T <- unlist(strsplit(extract$T.genes,split=', ',fixed=TRUE))
  test.B <- unlist(strsplit(extract$B.genes,split=', ',fixed=TRUE))
  test.Plasma <- unlist(strsplit(extract$Plasma.genes,split=', ',fixed=TRUE))
  test.Mye <- unlist(strsplit(extract$Mye.genes,split=', ',fixed=TRUE))
  test.MyoEpi <- unlist(strsplit(extract$MyoEpi.genes,split=', ',fixed=TRUE))
  test.PVL <- unlist(strsplit(extract$PVL.genes,split=', ',fixed=TRUE))
  test.Fibro <- unlist(strsplit(extract$Fibro.genes,split=', ',fixed=TRUE))
  test.Endo <- unlist(strsplit(extract$Endo.genes,split=', ',fixed=TRUE))
  test.Epi <- unlist(strsplit(extract$Epi.genes,split=', ',fixed=TRUE))
  test.StrongEpi <- unlist(strsplit(extract$StrongEpi.genes,split=', ',fixed=TRUE))
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #get the expression matrix for this cell and those genes ++++++++++++++
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  #NK Expression Sum ++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  NK.test <- as.data.frame(Subset.Mat[c(test.NK, "PTPRC"), cell])
  colnames(NK.test) <- cell
  rownames(NK.test) <- c(test.NK, "PTPRC")
  NK.exp.sum <- colSums(NK.test)
  NK.exp.avg <- NK.exp.sum/69059  #average expression (denominator is #rows of dataset)
  #NK.exp.avg <- colMeans(NK.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #T Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Tcell.test <- as.data.frame(Subset.Mat[c(test.T, "PTPRC"), cell])
  colnames(Tcell.test) <- cell
  rownames(Tcell.test) <- c(test.T, "PTPRC")
  T.exp.sum <- colSums(Tcell.test)
  T.exp.avg <- T.exp.sum/69059  #average expression (denominator is #rows of dataset)
  #T.exp.avg <- colMeans(Tcell.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #B Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Bcell.test <- as.data.frame(Subset.Mat[c(test.B, "PTPRC"), cell])
  colnames(Bcell.test) <- cell
  rownames(Bcell.test) <- c(test.B, "PTPRC")
  B.exp.sum <- colSums(Bcell.test)
  B.exp.avg <- B.exp.sum/69059  #average expression (denominator is #rows of dataset)
  #B.exp.avg <- colMeans(Bcell.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #Plasma Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  Plasmacell.test <- as.data.frame(Subset.Mat[c(test.Plasma, "PTPRC"), cell])
  colnames(Plasmacell.test) <- cell
  rownames(Plasmacell.test) <- c(test.Plasma, "PTPRC")
  Plasma.exp.sum <- colSums(Plasmacell.test)
  Plasma.exp.avg <- Plasma.exp.sum/69059  #average expression (denominator is #rows of dataset)
  #Plasma.exp.avg <- colMeans(Plasmacell.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #Myeloid Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Myeloid.test <- as.data.frame(Subset.Mat[c(test.Mye, "PTPRC"), cell])
  colnames(Myeloid.test) <- cell
  rownames(Myeloid.test) <- c(test.Mye, "PTPRC")
  Myeloid.exp.sum <- colSums(Myeloid.test)
  Myeloid.exp.avg <- Myeloid.exp.sum/69059 #average expression (denominator is #rows of dataset)
  #Myeloid.exp.avg <- colMeans(Myeloid.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #Myoepi Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  MyoEpi.test <- as.data.frame(Subset.Mat[c(test.MyoEpi), cell])
  colnames(MyoEpi.test) <- cell
  rownames(MyoEpi.test) <- c(test.MyoEpi)
  MyoEpi.exp.sum <- colSums(MyoEpi.test)
  MyoEpi.exp.avg <- MyoEpi.exp.sum/69059 #average expression (denominator is #rows of dataset)
  #MyoEpi.exp.avg <- colMeans(MyoEpi.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #PVL Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  PVL.test <- as.data.frame(Subset.Mat[c(test.PVL), cell])
  colnames(PVL.test) <- cell
  rownames(PVL.test) <- c(test.PVL)
  PVL.exp.sum <- colSums(PVL.test)
  PVL.exp.avg <- PVL.exp.sum/69059 #average expression (denominator is #rows of dataset)
  #PVL.exp.avg <- colMeans(PVL.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  #Fibro Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Fibro.test <- as.data.frame(Subset.Mat[c(test.Fibro), cell])
  colnames(Fibro.test) <- cell
  rownames(Fibro.test) <- c(test.Fibro)
  Fibro.exp.sum <- colSums(Fibro.test)
  Fibro.exp.avg <- Fibro.exp.sum/69059 #average expression (denominator is #rows of dataset)
  #Fibro.exp.avg <- colMeans(Fibro.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  #Endothelial Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Endo.test <- as.data.frame(Subset.Mat[c(test.Endo), cell])
  colnames(Endo.test) <- cell
  rownames(Endo.test) <- c(test.Endo)
  Endo.exp.sum <- colSums(Endo.test)
  Endo.exp.avg <- Endo.exp.sum/69059 #average expression (denominator is #rows of dataset)
  #Endo.exp.avg <- colMeans(Endo.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  #Epithelial Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Epi.test <- as.data.frame(Subset.Mat[c(test.Epi, test.StrongEpi), cell])
  colnames(Epi.test) <- cell
  rownames(Epi.test) <- c(test.Epi, test.StrongEpi)
  Epi.exp.sum <- colSums(Epi.test)
  Epi.exp.avg <- Epi.exp.sum/69059 #average expression (denominator is #rows of dataset)
  #Epi.exp.avg <- colMeans(Epi.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Final Call ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  if ((NK.exp.avg > T.exp.avg) & (NK.exp.avg > B.exp.avg) &
      (NK.exp.avg > Plasma.exp.avg) & (NK.exp.avg > Myeloid.exp.avg) &
      (NK.exp.avg > MyoEpi.exp.avg) & (NK.exp.avg > PVL.exp.avg) & 
      (NK.exp.avg > Fibro.exp.avg) & (NK.exp.avg > Endo.exp.avg) &
      (NK.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "NK Cells"
    
    
  }
  
  else if ((T.exp.avg > NK.exp.avg) & (T.exp.avg > B.exp.avg) &
           (T.exp.avg > Plasma.exp.avg) & (T.exp.avg > Myeloid.exp.avg) &
           (T.exp.avg > MyoEpi.exp.avg) & (T.exp.avg > PVL.exp.avg) & 
           (T.exp.avg > Fibro.exp.avg) & (T.exp.avg > Endo.exp.avg) &
           (T.exp.avg > Epi.exp.avg)) {
    
    
    if(("CD8A" %in% rownames(Subset.Mat)) | ("CD8B" %in% rownames(Subset.Mat)))
    {
      CD8T.test <- as.data.frame(Subset.Mat[c("CD8A","CD8B"), cell])
      colnames(CD8T.test) <- cell
      rownames(CD8T.test) <-c("CD8A","CD8B")
      CD8.exp.sum <- colSums(CD8T.test)
      CD8.exp.avg <- CD8.exp.sum/69059 
      #CD8.exp.avg <- colMeans(CD8T.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      CD8.exp.avg = 0
    }
    
    if("CD4" %in% rownames(Subset.Mat))
    {
      CD4.test <- as.data.frame(Subset.Mat[c("CD4"), cell])
      colnames(CD4.test) <- cell
      rownames(CD4.test) <-c("CD4")
      CD4.exp.sum <- colSums(CD4.test)
      CD4.exp.avg <- CD4.exp.sum/69059 
      #CD4.exp.avg <- colMeans(CD4.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }  
    else {
      CD4.exp.avg = 0
    }
    if("FOXP3" %in% rownames(Subset.Mat))
    {
      Treg.test <- as.data.frame(Subset.Mat[c("FOXP3"), cell])
      colnames(Treg.test) <- cell
      rownames(Treg.test) <-c("FOXP3")
      Treg.exp.sum <- colSums(Treg.test)
      Treg.exp.avg <- Treg.exp.sum/69059 
      #Treg.exp.avg <- colMeans(Treg.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }  
    else {
      Treg.exp.avg = 0
    }
    
    if ((CD8.exp.avg > CD4.exp.avg) & (CD8.exp.avg > Treg.exp.avg))
    {
      cell.result[[cell]] <- "CD8+ T Cells"
    }
    
    else if ((CD4.exp.avg > CD8.exp.avg) & (CD4.exp.avg > Treg.exp.avg))
    {
      cell.result[[cell]] <- "CD4+ T Cells"
    }
    
    else if ((Treg.exp.avg > CD8.exp.avg) & (Treg.exp.avg > CD4.exp.avg))
    {
      cell.result[[cell]] <- "Regulatory T Cells"
    }
    
    else 
    {
      cell.result[[cell]] <- "Unspecific T Cells"
    }
    
  }
  
  else if ((B.exp.avg > NK.exp.avg) & (B.exp.avg > T.exp.avg) &
           (B.exp.avg > Plasma.exp.avg) & (B.exp.avg > Myeloid.exp.avg) &
           (B.exp.avg > MyoEpi.exp.avg) & (B.exp.avg > PVL.exp.avg) & 
           (B.exp.avg > Fibro.exp.avg) & (B.exp.avg > Endo.exp.avg) &
           (B.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "B Cells"
    
    
  }
  
  else if ((Plasma.exp.avg > NK.exp.avg) & (Plasma.exp.avg > T.exp.avg) &
           (Plasma.exp.avg > B.exp.avg) & (Plasma.exp.avg > Myeloid.exp.avg) &
           (Plasma.exp.avg > MyoEpi.exp.avg) & (Plasma.exp.avg > PVL.exp.avg) & 
           (Plasma.exp.avg > Fibro.exp.avg) & (Plasma.exp.avg > Endo.exp.avg) &
           (Plasma.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "Plasma Cells"
    
    
  }
  
  else if ((Myeloid.exp.avg > NK.exp.avg) & (Myeloid.exp.avg > T.exp.avg) &
           (Myeloid.exp.avg > B.exp.avg) & (Myeloid.exp.avg > Plasma.exp.avg) &
           (Myeloid.exp.avg > MyoEpi.exp.avg) & (Myeloid.exp.avg > PVL.exp.avg) & 
           (Myeloid.exp.avg > Fibro.exp.avg) & (Myeloid.exp.avg > Endo.exp.avg) &
           (Myeloid.exp.avg > Epi.exp.avg)) {
    
    
    if(("CD33" %in% rownames(Subset.Mat)) | ("ITGAM" %in% rownames(Subset.Mat)) | 
       ("FUT4" %in% rownames(Subset.Mat)) | ("CEACAM8" %in% rownames(Subset.Mat)) |
       ("IL4R" %in% rownames(Subset.Mat)) | ("HLA-DR" %in% rownames(Subset.Mat)))
    {
      MDSc.test <- as.data.frame(Subset.Mat[c("CD33", #https://www.nature.com/articles/ncomms12150; table 1
                                              "ITGAM", 
                                              "FUT4", #CD15
                                              "CEACAM8",
                                              "IL4R",
                                              "HLA-DRA"), cell])
      colnames(MDSc.test) <- cell
      rownames(MDSc.test) <-c("CD33", #https://www.nature.com/articles/ncomms12150; table 1
                              "ITGAM", 
                              "FUT4", #CD15
                              "CEACAM8",
                              "IL4R",
                              "HLA-DRA")
      MDSc.exp.sum <- colSums(MDSc.test)
      MDSc.exp.avg <- MDSc.exp.sum/69059 
      #MDSc.exp.avg <- colMeans(MDSc.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      MDSc.exp.avg = 0
    }
    if(("INHBA" %in% rownames(Subset.Mat)) | ("IL1RN" %in% rownames(Subset.Mat)) | 
       ("CCL4" %in% rownames(Subset.Mat)) | ("NLRP3" %in% rownames(Subset.Mat)) |
       ("EREG" %in% rownames(Subset.Mat)) | ("IL1B" %in% rownames(Subset.Mat)) |
       ("LYVE1" %in% rownames(Subset.Mat)) |  ("PLTP" %in% rownames(Subset.Mat)) |
       ("SELENOP" %in% rownames(Subset.Mat)) |  ("C1QC" %in% rownames(Subset.Mat)) |
       ("C1QA" %in% rownames(Subset.Mat)) |  ("APOE" %in% rownames(Subset.Mat)))
    {
      Macro.test <- as.data.frame(Subset.Mat[c("INHBA", "IL1RN", "CCL4", "NLRP3",
                                               "EREG", "IL1B", "LYVE1", "PLTP",
                                               "SELENOP", "C1QC", "C1QA", "APOE"), cell])
      colnames(Macro.test) <- cell
      rownames(Macro.test) <-c("INHBA", "IL1RN", "CCL4", "NLRP3",
                               "EREG", "IL1B", "LYVE1", "PLTP",
                               "SELENOP", "C1QC", "C1QA", "APOE")
      Macro.exp.sum <- colSums(Macro.test)
      Macro.exp.avg <- Macro.exp.sum/69059 
      #Macro.exp.avg <- colMeans(Macro.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      Macro.exp.avg = 0
    }
    
    if(("FCN1" %in% rownames(Subset.Mat)) | ("S100A9" %in% rownames(Subset.Mat)) | 
       ("S100A8" %in% rownames(Subset.Mat)) | ("FCGR3A" %in% rownames(Subset.Mat)) |
       ("LST1" %in% rownames(Subset.Mat)) | ("LILRB2" %in% rownames(Subset.Mat)))
    {
      Mono.test <- as.data.frame(Subset.Mat[c("FCN1", "S100A9", "S100A8", "FCGR3A", "LST1", "LILRB2"), cell])
      colnames(Mono.test) <- cell
      rownames(Mono.test) <-c("FCN1", "S100A9", "S100A8", "FCGR3A", "LST1", "LILRB2")
      Mono.exp.sum <- colSums(Mono.test)
      Mono.exp.avg <- Mono.exp.sum/69059 
      #Mono.exp.avg <- colMeans(Mono.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      Mono.exp.avg = 0
    }
    if(("CXCR1" %in% rownames(Subset.Mat)) | ("FUT4" %in% rownames(Subset.Mat)) | 
       ("FCGR3A" %in% rownames(Subset.Mat)) | ("CSF3R" %in% rownames(Subset.Mat)) |
       ("S100A9" %in% rownames(Subset.Mat)) | ("TNF" %in% rownames(Subset.Mat)) |
       ("CD274" %in% rownames(Subset.Mat)))
    {
      Neut.test <- as.data.frame(Subset.Mat[c("CXCR1", #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
                                              "FUT4", #https://www.rndsystems.com/resources/cell-markers/immune-cells/granulocytes/neutrophil-cell-markers
                                              "FCGR3A", #CD16
                                              "CSF3R", #TCGA, Chung
                                              "S100A9","TNF","CD274"), cell])
      colnames(Neut.test) <- cell
      rownames(Neut.test) <-c("CXCR1", #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
                              "FUT4", #https://www.rndsystems.com/resources/cell-markers/immune-cells/granulocytes/neutrophil-cell-markers
                              "FCGR3A", #CD16
                              "CSF3R", #TCGA, Chung
                              "S100A9","TNF","CD274")
      Neut.exp.sum <- colSums(Neut.test)
      Neut.exp.avg <- Neut.exp.sum/69059 
      #Neut.exp.avg <- colMeans(Neut.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      Neut.exp.avg = 0
    }
    if(("LILRA4" %in% rownames(Subset.Mat)) | ("GZMB" %in% rownames(Subset.Mat)) | 
       ("IL3RA" %in% rownames(Subset.Mat)) | ("CLEC9A" %in% rownames(Subset.Mat)) |
       ("FLT3" %in% rownames(Subset.Mat)) | ("IDO1" %in% rownames(Subset.Mat)) |
       ("CD1C" %in% rownames(Subset.Mat)) | ("FCER1A" %in% rownames(Subset.Mat)) |
       ("HLA-DQA1" %in% rownames(Subset.Mat)) | ("LAMP3" %in% rownames(Subset.Mat)) |
       ("CCR7" %in% rownames(Subset.Mat)) | ("FSCN1" %in% rownames(Subset.Mat)))
    {
      DC.test <- as.data.frame(Subset.Mat[c("LILRA4","GZMB","IL3RA",
                                            "CLEC9A", "FLT3", "IDO1", 
                                            "CD1C", "FCER1A","HLA-DQA1", 
                                            "LAMP3", "CCR7", "FSCN1"), cell])
      colnames(DC.test) <- cell
      rownames(DC.test) <-c("LILRA4","GZMB","IL3RA",
                            "CLEC9A", "FLT3", "IDO1", 
                            "CD1C", "FCER1A","HLA-DQA1", 
                            "LAMP3", "CCR7", "FSCN1")
      DC.exp.sum <- colSums(DC.test)
      DC.exp.avg <- DC.exp.sum/69059 
      #DC.exp.avg <- colMeans(DC.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      DC.exp.avg = 0
    }
    if(("KIT" %in% rownames(Subset.Mat)) | ("TPSAB1" %in% rownames(Subset.Mat)) | 
       ("CPA4" %in% rownames(Subset.Mat)))
    {
      Mast.test <- as.data.frame(Subset.Mat[c("KIT","TPSAB1","CPA4"), cell])
      colnames(Mast.test) <- cell
      rownames(Mast.test) <-c("KIT","TPSAB1","CPA4")
      Mast.exp.sum <- colSums(Mast.test)
      Mast.exp.avg <- Mast.exp.sum/69059 
      #Mast.exp.avg <- colMeans(Mast.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
    }
    else {
      Mast.exp.avg = 0
    }
    
    
    if ((MDSc.exp.avg > Macro.exp.avg) & (MDSc.exp.avg > Mono.exp.avg) & 
        (MDSc.exp.avg > Neut.exp.avg) & (MDSc.exp.avg > DC.exp.avg) &
        (MDSc.exp.avg > Mast.exp.avg))
    {
      cell.result[[cell]] <- "MDSCs"
    }
    
    else if ((Macro.exp.avg > MDSc.exp.avg) & (Macro.exp.avg > Mono.exp.avg) & 
             (Macro.exp.avg > Neut.exp.avg) & (Macro.exp.avg > DC.exp.avg) &
             (Macro.exp.avg > Mast.exp.avg))
    {
      cell.result[[cell]] <- "Macrophages"
    }
    
    else if ((Mono.exp.avg > MDSc.exp.avg) & (Mono.exp.avg > Macro.exp.avg) & 
             (Mono.exp.avg > Neut.exp.avg) & (Mono.exp.avg > DC.exp.avg) &
             (Mono.exp.avg > Mast.exp.avg))
    {
      cell.result[[cell]] <- "Monocytes"
    }
    
    else if ((Neut.exp.avg > MDSc.exp.avg) & (Neut.exp.avg > Macro.exp.avg) & 
             (Neut.exp.avg > Mono.exp.avg) & (Neut.exp.avg > DC.exp.avg) &
             (Neut.exp.avg > Mast.exp.avg))
    {
      cell.result[[cell]] <- "Neutrophils"
    }
    
    else if ((DC.exp.avg > MDSc.exp.avg) & (DC.exp.avg > Macro.exp.avg) & 
             (DC.exp.avg > Mono.exp.avg) & (DC.exp.avg > Neut.exp.avg) &
             (DC.exp.avg > Mast.exp.avg))
    {
      cell.result[[cell]] <- "Dendritic Cells"
    }
    
    else if ((Mast.exp.avg > MDSc.exp.avg) & (Mast.exp.avg > Macro.exp.avg) & 
             (Mast.exp.avg > Mono.exp.avg) & (Mast.exp.avg > Neut.exp.avg) &
             (Mast.exp.avg > DC.exp.avg))
    {
      cell.result[[cell]] <- "Mast Cells"
    }
    
    else 
    {
      cell.result[[cell]] <- "Unspecific Myeloid Cells"
    }
    
    
  }
  
  else if ((MyoEpi.exp.avg > NK.exp.avg) & (MyoEpi.exp.avg > T.exp.avg) &
           (MyoEpi.exp.avg > B.exp.avg) & (MyoEpi.exp.avg > Plasma.exp.avg) &
           (MyoEpi.exp.avg > Myeloid.exp.avg) & (MyoEpi.exp.avg > PVL.exp.avg) & 
           (MyoEpi.exp.avg > Fibro.exp.avg) & (MyoEpi.exp.avg > Endo.exp.avg) &
           (MyoEpi.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "Myoepithelial Cells"
    
    
  }
  
  else if ((PVL.exp.avg > NK.exp.avg) & (PVL.exp.avg > T.exp.avg) &
           (PVL.exp.avg > B.exp.avg) & (PVL.exp.avg > Plasma.exp.avg) &
           (PVL.exp.avg > Myeloid.exp.avg) & (PVL.exp.avg > MyoEpi.exp.avg) & 
           (PVL.exp.avg > Fibro.exp.avg) & (PVL.exp.avg > Endo.exp.avg) &
           (PVL.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "Perivascular-like (PVL) Cells"
    
    
  }
  
  else if ((Fibro.exp.avg > NK.exp.avg) & (Fibro.exp.avg > T.exp.avg) &
           (Fibro.exp.avg > B.exp.avg) & (Fibro.exp.avg > Plasma.exp.avg) &
           (Fibro.exp.avg > Myeloid.exp.avg) & (Fibro.exp.avg > MyoEpi.exp.avg) & 
           (Fibro.exp.avg > PVL.exp.avg) & (Fibro.exp.avg > Endo.exp.avg) &
           (Fibro.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "Fibroblasts"
    
    
  }
  
  else if ((Endo.exp.avg > NK.exp.avg) & (Endo.exp.avg > T.exp.avg) &
           (Endo.exp.avg > B.exp.avg) & (Endo.exp.avg > Plasma.exp.avg) &
           (Endo.exp.avg > Myeloid.exp.avg) & (Endo.exp.avg > MyoEpi.exp.avg) & 
           (Endo.exp.avg > PVL.exp.avg) & (Endo.exp.avg > Fibro.exp.avg) &
           (Endo.exp.avg > Epi.exp.avg)) {
    
    cell.result[[cell]] <- "Endothelial Cells"
    
    
  }
  
  else if ((Epi.exp.avg > NK.exp.avg) & (Epi.exp.avg > T.exp.avg) &
           (Epi.exp.avg > B.exp.avg) & (Epi.exp.avg > Plasma.exp.avg) &
           (Epi.exp.avg > Myeloid.exp.avg) & (Epi.exp.avg > MyoEpi.exp.avg) & 
           (Epi.exp.avg > PVL.exp.avg) & (Epi.exp.avg > Fibro.exp.avg) &
           (Epi.exp.avg > Endo.exp.avg)) {
    
    cell.result[[cell]] <- "Epithelial Cells"
    
    
  }
  else 
    cell.result[[cell]] <- "Unknown"
  
  
}

head(cell.result)
ExpOnly.Results <- t(as.data.frame(cell.result))
head(ExpOnly.Results)
table(ExpOnly.Results)

head(MarkerBased.dataframe)
table(MarkerBased.dataframe)


Uscore.result <- list()
for (cell in names(markersum.list)) {
  
  cell.name <- markersum.list[[cell]]
  
  NK.UScore <- sobj@meta.data[cell,"signature_1NK.mark2"]
  T.UScore <- sobj@meta.data[cell,"signature_1T_karacd4"]
  B.UScore <- sobj@meta.data[cell,"signature_1B_KaraBrech"]
  Plasma.UScore <- sobj@meta.data[cell,"signature_1plasma_colo"]
  Myeloid.UScore <- sobj@meta.data[cell,"signature_1myeloid_wuElliot"]
  Epi.UScore <- sobj@meta.data[cell,"signature_1genepi_Kara"]
  MyoEpi.UScore <- sobj@meta.data[cell,"signature_1Myoepi_wuNguyen"]
  PVL.UScore <- sobj@meta.data[cell,"signature_1newWuPVL"]
  Fibro.UScore <- sobj@meta.data[cell,"signature_1fibro_wumelan"]
  Endo.UScore <- sobj@meta.data[cell,"signature_1endokara"]
  
  CD8.UScore <- sobj@meta.data[cell,"signature_1CD8sig"]
  CD4.UScore <- sobj@meta.data[cell,"signature_1CD4sig"]
  Treg.UScore <- sobj@meta.data[cell,"signature_1Treg_sig"]
  MDSc.UScore <- sobj@meta.data[cell,"signature_1MDSc_sig"]
  Macro.UScore <- sobj@meta.data[cell,"signature_1Macro_Cheng"]
  Mono.UScore <- sobj@meta.data[cell,"signature_1Mono_sig"]
  Neut.UScore <- sobj@meta.data[cell,"signature_1Neutro_sig"]
  DC.UScore <- sobj@meta.data[cell,"signature_1DC_sig"]
  Mast.UScore <- sobj@meta.data[cell,"signature_1mast_Cheng"]
  
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Final Call ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  if ((NK.UScore > T.UScore) & (NK.UScore > B.UScore) &
      (NK.UScore > Plasma.UScore) & (NK.UScore > Myeloid.UScore) &
      (NK.UScore > MyoEpi.UScore) & (NK.UScore > PVL.UScore) & 
      (NK.UScore > Fibro.UScore) & (NK.UScore > Endo.UScore) &
      (NK.UScore > Epi.UScore)) {
    
    Uscore.result[[cell]] <- "NK Cells"
    
    
  }
  
  else if ((T.UScore > NK.UScore) & (T.UScore > B.UScore) &
           (T.UScore > Plasma.UScore) & (T.UScore > Myeloid.UScore) &
           (T.UScore > MyoEpi.UScore) & (T.UScore > PVL.UScore) & 
           (T.UScore > Fibro.UScore) & (T.UScore > Endo.UScore) &
           (T.UScore > Epi.UScore)) {
    
    
    
    if ((CD8.UScore > CD4.UScore) & (CD8.UScore > Treg.UScore))
    {
      Uscore.result[[cell]] <- "CD8+ T Cells"
    }
    
    else if ((CD4.UScore > CD8.UScore) & (CD4.UScore > Treg.UScore))
    {
      Uscore.result[[cell]] <- "CD4+ T Cells"
    }
    
    else if ((Treg.UScore > CD8.UScore) & (Treg.UScore > CD4.UScore))
    {
      Uscore.result[[cell]] <- "Regulatory T Cells"
    }
    
    else 
    {
      Uscore.result[[cell]] <- "Unspecific T Cells"
    }
    
  }
  
  else if ((B.UScore > NK.UScore) & (B.UScore > T.UScore) &
           (B.UScore > Plasma.UScore) & (B.UScore > Myeloid.UScore) &
           (B.UScore > MyoEpi.UScore) & (B.UScore > PVL.UScore) & 
           (B.UScore > Fibro.UScore) & (B.UScore > Endo.UScore) &
           (B.UScore > Epi.UScore)) {
    
    Uscore.result[[cell]] <- "B Cells"
    
    
  }
  
  else if ((Plasma.UScore > NK.UScore) & (Plasma.UScore > T.UScore) &
           (Plasma.UScore > B.UScore) & (Plasma.UScore > Myeloid.UScore) &
           (Plasma.UScore > MyoEpi.UScore) & (Plasma.UScore > PVL.UScore) & 
           (Plasma.UScore > Fibro.UScore) & (Plasma.UScore > Endo.UScore) &
           (Plasma.UScore > Epi.UScore)) {
    
    Uscore.result[[cell]] <- "Plasma Cells"
    
    
  }
  
  else if ((Myeloid.UScore > NK.UScore) & (Myeloid.UScore > T.UScore) &
           (Myeloid.UScore > B.UScore) & (Myeloid.UScore > Plasma.UScore) &
           (Myeloid.UScore > MyoEpi.UScore) & (Myeloid.UScore > PVL.UScore) & 
           (Myeloid.UScore > Fibro.UScore) & (Myeloid.UScore > Endo.UScore) &
           (Myeloid.UScore > Epi.UScore)) {
    
    
    if ((MDSc.UScore > Macro.UScore) & (MDSc.UScore > Mono.UScore) & 
        (MDSc.UScore > Neut.UScore) & (MDSc.UScore > DC.UScore) &
        (MDSc.UScore > Mast.UScore))
    {
      Uscore.result[[cell]] <- "MDSCs"
    }
    
    else if ((Macro.UScore > MDSc.UScore) & (Macro.UScore > Mono.UScore) & 
             (Macro.UScore > Neut.UScore) & (Macro.UScore > DC.UScore) &
             (Macro.UScore > Mast.UScore))
    {
      Uscore.result[[cell]] <- "Macrophages"
    }
    
    else if ((Mono.UScore > MDSc.UScore) & (Mono.UScore > Macro.UScore) & 
             (Mono.UScore > Neut.UScore) & (Mono.UScore > DC.UScore) &
             (Mono.UScore > Mast.UScore))
    {
      Uscore.result[[cell]] <- "Monocytes"
    }
    
    else if ((Neut.UScore > MDSc.UScore) & (Neut.UScore > Macro.UScore) & 
             (Neut.UScore > Mono.UScore) & (Neut.UScore > DC.UScore) &
             (Neut.UScore > Mast.UScore))
    {
      Uscore.result[[cell]] <- "Neutrophils"
    }
    
    else if ((DC.UScore > MDSc.UScore) & (DC.UScore > Macro.UScore) & 
             (DC.UScore > Mono.UScore) & (DC.UScore > Neut.UScore) &
             (DC.UScore > Mast.UScore))
    {
      Uscore.result[[cell]] <- "Dendritic Cells"
    }
    
    else if ((Mast.UScore > MDSc.UScore) & (Mast.UScore > Macro.UScore) & 
             (Mast.UScore > Mono.UScore) & (Mast.UScore > Neut.UScore) &
             (Mast.UScore > DC.UScore))
    {
      Uscore.result[[cell]] <- "Mast Cells"
    }
    
    else 
    {
      Uscore.result[[cell]] <- "Unspecific Myeloid Cells"
    }
    
    
  }
  
  else if ((MyoEpi.UScore > NK.UScore) & (MyoEpi.UScore > T.UScore) &
           (MyoEpi.UScore > B.UScore) & (MyoEpi.UScore > Plasma.UScore) &
           (MyoEpi.UScore > Myeloid.UScore) & (MyoEpi.UScore > PVL.UScore) & 
           (MyoEpi.UScore > Fibro.UScore) & (MyoEpi.UScore > Endo.UScore) &
           (MyoEpi.UScore > Epi.UScore)) {
    
    Uscore.result[[cell]] <- "Myoepithelial Cells"
    
    
  }
  
  else if ((PVL.UScore > NK.UScore) & (PVL.UScore > T.UScore) &
           (PVL.UScore > B.UScore) & (PVL.UScore > Plasma.UScore) &
           (PVL.UScore > Myeloid.UScore) & (PVL.UScore > MyoEpi.UScore) & 
           (PVL.UScore > Fibro.UScore) & (PVL.UScore > Endo.UScore) &
           (PVL.UScore > Epi.UScore)) {
    
    Uscore.result[[cell]] <- "Perivascular-like (PVL) Cells"
    
    
  }
  
  else if ((Fibro.UScore > NK.UScore) & (Fibro.UScore > T.UScore) &
           (Fibro.UScore > B.UScore) & (Fibro.UScore > Plasma.UScore) &
           (Fibro.UScore > Myeloid.UScore) & (Fibro.UScore > MyoEpi.UScore) & 
           (Fibro.UScore > PVL.UScore) & (Fibro.UScore > Endo.UScore) &
           (Fibro.UScore > Epi.UScore)) {
    
    Uscore.result[[cell]] <- "Fibroblasts"
    
    
  }
  
  else if ((Endo.UScore > NK.UScore) & (Endo.UScore > T.UScore) &
           (Endo.UScore > B.UScore) & (Endo.UScore > Plasma.UScore) &
           (Endo.UScore > Myeloid.UScore) & (Endo.UScore > MyoEpi.UScore) & 
           (Endo.UScore > PVL.UScore) & (Endo.UScore > Fibro.UScore) &
           (Endo.UScore > Epi.UScore)) {
    
    Uscore.result[[cell]] <- "Endothelial Cells"
    
    
  }
  
  else if ((Epi.UScore > NK.UScore) & (Epi.UScore > T.UScore) &
           (Epi.UScore > B.UScore) & (Epi.UScore > Plasma.UScore) &
           (Epi.UScore > Myeloid.UScore) & (Epi.UScore > MyoEpi.UScore) & 
           (Epi.UScore > PVL.UScore) & (Epi.UScore > Fibro.UScore) &
           (Epi.UScore > Endo.UScore)) {
    
    Uscore.result[[cell]] <- "Epithelial Cells"
    
    
  }
  else 
    Uscore.result[[cell]] <- "Unknown"
  
  
}

head(Uscore.result)
UScore.dataframe <- t(as.data.frame(Uscore.result))
head(UScore.dataframe)
table(UScore.dataframe) 


rownames(MarkerBased.dataframe) <- gsub(".", "-", rownames(MarkerBased.dataframe), fixed=TRUE)
rownames(ExpOnly.Results) <- gsub(".", "-", rownames(ExpOnly.Results), fixed=TRUE)
#rownames(UScore.dataframe) <- gsub(".", "-", rownames(UScore.dataframe), fixed=TRUE)


#test NA cell: Palprim_AACCGCGTCTCTGTCG-1_0135_TNBCTot #<- dash got turned into "." 
#test non-NA cell in both: Wu_CID44041_CCTACCACAATAGAGT
combo.reference <- AddMetaData(combo.reference, MarkerBased.dataframe, "MarkerBasedKara.NK")
combo.reference <- AddMetaData(combo.reference, ExpOnly.Results, "ExpBasedOnly.NK")

combo.reference <- AddMetaData(combo.reference, MarkerBased.dataframe, "MarkerBasedKara.Epi")
combo.reference <- AddMetaData(combo.reference, ExpOnly.Results, "ExpBasedOnly.Epi")

combo.reference <- AddMetaData(combo.reference, MarkerBased.dataframe, "MarkerBasedKara.Stroma")
combo.reference <- AddMetaData(combo.reference, ExpOnly.Results, "ExpBasedOnly.Stroma")

combo.reference <- AddMetaData(combo.reference, MarkerBased.dataframe, "MarkerBasedKara.B")
combo.reference <- AddMetaData(combo.reference, ExpOnly.Results, "ExpBasedOnly.B")

combo.reference <- AddMetaData(combo.reference, MarkerBased.dataframe, "MarkerBasedKara.Mye")
combo.reference <- AddMetaData(combo.reference, ExpOnly.Results, "ExpBasedOnly.Mye")

combo.reference <- AddMetaData(combo.reference, MarkerBased.dataframe, "MarkerBasedKara.T")
combo.reference <- AddMetaData(combo.reference, ExpOnly.Results, "ExpBasedOnly.T")

colnames(combo.reference@meta.data)

combo.reference$MarkerBasedKara <- apply(combo.reference@meta.data[,c(77, 79, 81, 83, 85, 87)], 1, function(x) x[!is.na(x)][1])
combo.reference@meta.data <- combo.reference@meta.data[,-c(77, 79, 81, 83, 85, 87)]

combo.reference$ExpBasedOnly <- apply(combo.reference@meta.data[,c(77:82)], 1, function(x) x[!is.na(x)][1])
combo.reference@meta.data <- combo.reference@meta.data[,-c(77:82)]

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects")
#saveRDS(combo.reference, "PrimObject_6422.rds")
combo.reference <- readRDS("PrimObject_6422.rds")

##visualize ===========================

pdf("test.pdf", width = 14.22, height = 13.01)
DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE) + ggtitle("")
DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE, group.by = "MarkerBasedKara") + ggtitle("")
DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE, group.by = "ExpBasedOnly") + ggtitle("")
dev.off()


#final is where all agree ============================

combo.reference@meta.data[1:20,c(75, 77:78)]
#75 = UCell
#77 = markerbased
#78 = expression based only
table(combo.reference$MarkerBasedKara)
#problem child row: 78936
combo.reference$RefineMarkerBased <- "hi"
for (i in 1:nrow(combo.reference@meta.data)) {
  combo.reference@meta.data[i,79] <-ifelse(is.na(combo.reference@meta.data[i,77]),
                                           as.character(combo.reference@meta.data[i,75]), as.character(combo.reference@meta.data[i,77]))
}

for (i in 1:nrow(combo.reference@meta.data)) {

  if (combo.reference@meta.data[i,79] == "Unspecific T Cells")
  {
    if ((combo.reference@meta.data[i,75] == "CD8+ T Cells") |
        (combo.reference@meta.data[i,75] == "CD4+ T Cells") |
        (combo.reference@meta.data[i,75] == "Regulatory T Cells")) 
    {
      combo.reference@meta.data[i,79] <- as.character(combo.reference@meta.data[i,75])
      
    }
    else {
      combo.reference@meta.data[i,79] <- "not Tcell by UScore"
    }  
    }
    
  else if (combo.reference@meta.data[i,79] == "T Cells")
  {
    if ((combo.reference@meta.data[i,75] == "CD8+ T Cells") |
        (combo.reference@meta.data[i,75] == "CD4+ T Cells") |
        (combo.reference@meta.data[i,75] == "Regulatory T Cells")) 
      {
          combo.reference@meta.data[i,79] <- as.character(combo.reference@meta.data[i,75])

      }
    else {
      combo.reference@meta.data[i,79] <- "not Tcell by UScore"
    }
    
  }
  else if (combo.reference@meta.data[i,79] == "Unspecific Myeloid Cells")
  {
    if ((combo.reference@meta.data[i,75] == "Macrophages") |
        (combo.reference@meta.data[i,75] == "Mast Cells") |
        (combo.reference@meta.data[i,75] == "Monocytes") |
        (combo.reference@meta.data[i,75] == "Dendritic Cells") |
        (combo.reference@meta.data[i,75] == "Neutrophils") |
        (combo.reference@meta.data[i,75] == "MDSCs")) 
    {
      combo.reference@meta.data[i,79] <- as.character(combo.reference@meta.data[i,75])
      
    }
    else {
      combo.reference@meta.data[i,79] <- "not Myeloid cell by UScore"
    }  
  }
  else if (combo.reference@meta.data[i,79] == "Myeloid Cells")
  {
    if ((combo.reference@meta.data[i,75] == "Macrophages") |
        (combo.reference@meta.data[i,75] == "Mast Cells") |
        (combo.reference@meta.data[i,75] == "Monocytes") |
        (combo.reference@meta.data[i,75] == "Dendritic Cells") |
        (combo.reference@meta.data[i,75] == "Neutrophils") |
        (combo.reference@meta.data[i,75] == "MDSCs")) 
    {
      combo.reference@meta.data[i,79] <- as.character(combo.reference@meta.data[i,75])
      
    }
    else {
      combo.reference@meta.data[i,79] <- "not Myeloid cell by UScore"
    }  
  }
  else if (combo.reference@meta.data[i,79] == "Unknown")
  {
    combo.reference@meta.data[i,79] <- as.character(combo.reference@meta.data[i,75])
  }
  else if (combo.reference@meta.data[i,79] == "Zero Expression for All Markers")
  {
    combo.reference@meta.data[i,79] <- as.character(combo.reference@meta.data[i,75])
  }
  else if (combo.reference@meta.data[i,79] == "Non-Specific Immune Cells")
  {
    if ((combo.reference@meta.data[i,75] == "CD8+ T Cells") |
        (combo.reference@meta.data[i,75] == "CD4+ T Cells") |
        (combo.reference@meta.data[i,75] == "Regulatory T Cells") | 
        (combo.reference@meta.data[i,75] == "NK Cells") |
        (combo.reference@meta.data[i,75] == "B Cells") |
        (combo.reference@meta.data[i,75] == "Plasma Cells") |
        (combo.reference@meta.data[i,75] == "MDSCs") |
        (combo.reference@meta.data[i,75] == "Macrophages") |
        (combo.reference@meta.data[i,75] == "Monocytes") |
        (combo.reference@meta.data[i,75] == "Neutrophils") |
        (combo.reference@meta.data[i,75] == "Dendritic Cells") |
        (combo.reference@meta.data[i,75] == "Mast Cells")) 
    {
      combo.reference@meta.data[i,79] <- as.character(combo.reference@meta.data[i,75])
      
    }
    else {
      combo.reference@meta.data[i,79] <- "not immune cell by UScore"
    }  }
  
  else {
    combo.reference@meta.data[i,79] <- as.character(combo.reference@meta.data[i,79])
    
  }
  }

combo.reference@meta.data[78936,c(75, 77:79)]

table(combo.reference$RefineMarkerBased)



combo.reference$RefineExpBased <- "hi"
for (i in 1:nrow(combo.reference@meta.data)) {
  combo.reference@meta.data[i,80] <-ifelse(is.na(combo.reference@meta.data[i,78]),
                                           as.character(combo.reference@meta.data[i,75]), as.character(combo.reference@meta.data[i,78]))
}

for (i in 1:nrow(combo.reference@meta.data)) {
  
  if (combo.reference@meta.data[i,80] == "Unspecific T Cells")
  {
    if ((combo.reference@meta.data[i,75] == "CD8+ T Cells") |
        (combo.reference@meta.data[i,75] == "CD4+ T Cells") |
        (combo.reference@meta.data[i,75] == "Regulatory T Cells")) 
    {
      combo.reference@meta.data[i,80] <- as.character(combo.reference@meta.data[i,75])
      
    }
    else {
      combo.reference@meta.data[i,80] <- "not Tcell by UScore"
    }  
  }
  else if (combo.reference@meta.data[i,80] == "Unspecific Myeloid Cells")
  {
    if ((combo.reference@meta.data[i,75] == "Macrophages") |
        (combo.reference@meta.data[i,75] == "Mast Cells") |
        (combo.reference@meta.data[i,75] == "Monocytes") |
        (combo.reference@meta.data[i,75] == "Dendritic Cells") |
        (combo.reference@meta.data[i,75] == "Neutrophils") |
        (combo.reference@meta.data[i,75] == "MDSCs")) 
    {
      combo.reference@meta.data[i,80] <- as.character(combo.reference@meta.data[i,75])
      
    }
    else {
      combo.reference@meta.data[i,80] <- "not Myeloid cell by UScore"
    }  
  }
  else if (combo.reference@meta.data[i,80] == "Myeloid Cells")
  {
    if ((combo.reference@meta.data[i,75] == "Macrophages") |
        (combo.reference@meta.data[i,75] == "Mast Cells") |
        (combo.reference@meta.data[i,75] == "Monocytes") |
        (combo.reference@meta.data[i,75] == "Dendritic Cells") |
        (combo.reference@meta.data[i,75] == "Neutrophils") |
        (combo.reference@meta.data[i,75] == "MDSCs")) 
    {
      combo.reference@meta.data[i,80] <- as.character(combo.reference@meta.data[i,75])
      
    }
    else {
      combo.reference@meta.data[i,80] <- "not Myeloid cell by UScore"
    }  
  }
  else if (combo.reference@meta.data[i,80] == "T Cells")
  {
    if ((combo.reference@meta.data[i,75] == "CD8+ T Cells") |
        (combo.reference@meta.data[i,75] == "CD4+ T Cells") |
        (combo.reference@meta.data[i,75] == "Regulatory T Cells")) 
    {
      combo.reference@meta.data[i,80] <- as.character(combo.reference@meta.data[i,75])
      
    }
    else {
      combo.reference@meta.data[i,80] <- "not Tcell by UScore"
    }
    
  }
  
  else if (combo.reference@meta.data[i,80] == "Unknown")
  {
    combo.reference@meta.data[i,80] <- as.character(combo.reference@meta.data[i,75])
  }
  else if (combo.reference@meta.data[i,80] == "Zero Expression for All Markers")
  {
    combo.reference@meta.data[i,80] <- as.character(combo.reference@meta.data[i,75])
  }
  else if (combo.reference@meta.data[i,80] == "Non-Specific Immune Cells")
  {
    if ((combo.reference@meta.data[i,75] == "CD8+ T Cells") |
        (combo.reference@meta.data[i,75] == "CD4+ T Cells") |
        (combo.reference@meta.data[i,75] == "Regulatory T Cells") | 
        (combo.reference@meta.data[i,75] == "NK Cells") |
        (combo.reference@meta.data[i,75] == "B Cells") |
        (combo.reference@meta.data[i,75] == "Plasma Cells") |
        (combo.reference@meta.data[i,75] == "MDSCs") |
        (combo.reference@meta.data[i,75] == "Macrophages") |
        (combo.reference@meta.data[i,75] == "Monocytes") |
        (combo.reference@meta.data[i,75] == "Neutrophils") |
        (combo.reference@meta.data[i,75] == "Dendritic Cells") |
        (combo.reference@meta.data[i,75] == "Mast Cells")) 
    {
      combo.reference@meta.data[i,80] <- as.character(combo.reference@meta.data[i,75])
      
    }
    else {
      combo.reference@meta.data[i,80] <- "not immune cell by UScore"
    }  }
  
  else {
    combo.reference@meta.data[i,80] <- as.character(combo.reference@meta.data[i,80])
    
  }
}

combo.reference@meta.data[78936,c(75, 77:80)]

table(combo.reference$RefineExpBased)


pdf("test.pdf", width = 14.22, height = 13.01)
DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE) + ggtitle("")
DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE, group.by = "RefineMarkerBased") + ggtitle("")
DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE, group.by = "RefineExpBased") + ggtitle("")
dev.off()


# agreement test ===============


#markerExpIdent.dat <- combo.reference@meta.data[,c(75,79:80)]
#head(markerExpIdent.dat)

library(DescTools)
colnames(combo.reference@meta.data)

combo.reference$celltype_final <- "hi"
for (i in 1:nrow(combo.reference@meta.data)){
  test <- c(as.character(combo.reference@meta.data[i,75]), 
          as.character(combo.reference@meta.data[i,79]), 
          as.character(combo.reference@meta.data[i,80]))
  
  test.entry <- as.character(Mode(test))
  
  if (!is.na(test.entry) == T) {
    if ((as.character(combo.reference@meta.data[i,1]) == "Savas") && (
      (test.entry != "Regulatory T Cells") && (test.entry != "CD8+ T Cells") &&
      (test.entry != "CD4+ T Cells"))){
      combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,75])
    }
    
    else if ((as.character(combo.reference@meta.data[i,1]) == "AziziT") && (
      (test.entry != "Regulatory T Cells") && (test.entry != "CD8+ T Cells") &&
      (test.entry != "CD4+ T Cells"))){
      combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,75])
    }
    
    else if ((as.character(combo.reference@meta.data[i,1]) == "Aziziimmune") && (
      (test.entry != "Mast Cells") | (test.entry != "Dendritic Cells") | (test.entry != "Neutrophils") |
      (test.entry != "Monocytes") | (test.entry != "Macrophages") | (test.entry != "MDSCs") |
      (test.entry != "Regulatory T Cells") | (test.entry != "CD4+ T Cells") | (test.entry != "CD8+ T Cells") |
      (test.entry != "Plasma Cells") | (test.entry != "B Cells") | (test.entry != "NK Cells"))){
      #if UCell has it as some immune cell
      if ((as.character(combo.reference@meta.data[i,75]) == "Mast Cells") |
          (as.character(combo.reference@meta.data[i,75]) == "Dendritic Cells") |
          (as.character(combo.reference@meta.data[i,75]) == "Neutrophils") | 
          (as.character(combo.reference@meta.data[i,75]) == "Monocytes") |
          (as.character(combo.reference@meta.data[i,75]) == "Macrophages") |
          (as.character(combo.reference@meta.data[i,75]) == "MDSCs") | 
          (as.character(combo.reference@meta.data[i,75]) == "Regulatory T Cells") |
          (as.character(combo.reference@meta.data[i,75]) == "CD4+ T Cells") |
          (as.character(combo.reference@meta.data[i,75]) == "CD8+ T Cells") |
          (as.character(combo.reference@meta.data[i,75]) == "Plasma Cells") | 
          (as.character(combo.reference@meta.data[i,75]) == "B Cells") | 
          (as.character(combo.reference@meta.data[i,75]) == "NK Cells")) 
      {
        combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,75])
      }
      #if karaavaz has it as some immune cell
      else if ((as.character(combo.reference@meta.data[i,79]) == "Mast Cells") |
          (as.character(combo.reference@meta.data[i,79]) == "Dendritic Cells") |
          (as.character(combo.reference@meta.data[i,79]) == "Neutrophils") | 
          (as.character(combo.reference@meta.data[i,79]) == "Monocytes") |
          (as.character(combo.reference@meta.data[i,79]) == "Macrophages") |
          (as.character(combo.reference@meta.data[i,79]) == "MDSCs") | 
          (as.character(combo.reference@meta.data[i,79]) == "Regulatory T Cells") |
          (as.character(combo.reference@meta.data[i,79]) == "CD4+ T Cells") |
          (as.character(combo.reference@meta.data[i,79]) == "CD8+ T Cells") |
          (as.character(combo.reference@meta.data[i,79]) == "Plasma Cells") | 
          (as.character(combo.reference@meta.data[i,79]) == "B Cells") | 
          (as.character(combo.reference@meta.data[i,79]) == "NK Cells")) {
        combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,79])
      }
      #if marker expression has it as some kind of immune cell
      else if ((as.character(combo.reference@meta.data[i,80]) == "Mast Cells") |
          (as.character(combo.reference@meta.data[i,80]) == "Dendritic Cells") |
          (as.character(combo.reference@meta.data[i,80]) == "Neutrophils") | 
          (as.character(combo.reference@meta.data[i,80]) == "Monocytes") |
          (as.character(combo.reference@meta.data[i,80]) == "Macrophages") |
          (as.character(combo.reference@meta.data[i,80]) == "MDSCs") | 
          (as.character(combo.reference@meta.data[i,80]) == "Regulatory T Cells") |
          (as.character(combo.reference@meta.data[i,80]) == "CD4+ T Cells") |
          (as.character(combo.reference@meta.data[i,80]) == "CD8+ T Cells") |
          (as.character(combo.reference@meta.data[i,80]) == "Plasma Cells") | 
          (as.character(combo.reference@meta.data[i,80]) == "B Cells") | 
          (as.character(combo.reference@meta.data[i,80]) == "NK Cells")) {
        combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,80])
      }
      else {
        combo.reference@meta.data[i,81] <- "noone called it immune"
        #combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,75])
        
      }
  }
    
    else if (test.entry == "not Tcell by UScore") {
      combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,75])
    }
    else if (test.entry == "not Myeloid cell by UScore") {
      combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,75])
      
    }
    else {
          combo.reference@meta.data[i,81] <- test.entry
    }
  }
  else {
    if ((as.character(combo.reference@meta.data[i,1]) == "Savas")){
      combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,75])
    }
    
    else if ((as.character(combo.reference@meta.data[i,1]) == "AziziT")){
      combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,75])
    }
    
    else if ((as.character(combo.reference@meta.data[i,1]) == "Aziziimmune")){
      #if UCell has it as some immune cell
      if ((as.character(combo.reference@meta.data[i,75]) == "Mast Cells") |
          (as.character(combo.reference@meta.data[i,75]) == "Dendritic Cells") |
          (as.character(combo.reference@meta.data[i,75]) == "Neutrophils") | 
          (as.character(combo.reference@meta.data[i,75]) == "Monocytes") |
          (as.character(combo.reference@meta.data[i,75]) == "Macrophages") |
          (as.character(combo.reference@meta.data[i,75]) == "MDSCs") | 
          (as.character(combo.reference@meta.data[i,75]) == "Regulatory T Cells") |
          (as.character(combo.reference@meta.data[i,75]) == "CD4+ T Cells") |
          (as.character(combo.reference@meta.data[i,75]) == "CD8+ T Cells") |
          (as.character(combo.reference@meta.data[i,75]) == "Plasma Cells") | 
          (as.character(combo.reference@meta.data[i,75]) == "B Cells") | 
          (as.character(combo.reference@meta.data[i,75]) == "NK Cells")) 
      {
        combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,75])
      }
      #if karaavaz has it as some immune cell
      else if ((as.character(combo.reference@meta.data[i,79]) == "Mast Cells") |
               (as.character(combo.reference@meta.data[i,79]) == "Dendritic Cells") |
               (as.character(combo.reference@meta.data[i,79]) == "Neutrophils") | 
               (as.character(combo.reference@meta.data[i,79]) == "Monocytes") |
               (as.character(combo.reference@meta.data[i,79]) == "Macrophages") |
               (as.character(combo.reference@meta.data[i,79]) == "MDSCs") | 
               (as.character(combo.reference@meta.data[i,79]) == "Regulatory T Cells") |
               (as.character(combo.reference@meta.data[i,79]) == "CD4+ T Cells") |
               (as.character(combo.reference@meta.data[i,79]) == "CD8+ T Cells") |
               (as.character(combo.reference@meta.data[i,79]) == "Plasma Cells") | 
               (as.character(combo.reference@meta.data[i,79]) == "B Cells") | 
               (as.character(combo.reference@meta.data[i,79]) == "NK Cells")) {
        combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,79])
      }
      #if marker expression has it as some kind of immune cell
      else if ((as.character(combo.reference@meta.data[i,80]) == "Mast Cells") |
               (as.character(combo.reference@meta.data[i,80]) == "Dendritic Cells") |
               (as.character(combo.reference@meta.data[i,80]) == "Neutrophils") | 
               (as.character(combo.reference@meta.data[i,80]) == "Monocytes") |
               (as.character(combo.reference@meta.data[i,80]) == "Macrophages") |
               (as.character(combo.reference@meta.data[i,80]) == "MDSCs") | 
               (as.character(combo.reference@meta.data[i,80]) == "Regulatory T Cells") |
               (as.character(combo.reference@meta.data[i,80]) == "CD4+ T Cells") |
               (as.character(combo.reference@meta.data[i,80]) == "CD8+ T Cells") |
               (as.character(combo.reference@meta.data[i,80]) == "Plasma Cells") | 
               (as.character(combo.reference@meta.data[i,80]) == "B Cells") | 
               (as.character(combo.reference@meta.data[i,80]) == "NK Cells")) {
        combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,80])
      }
      else {
        combo.reference@meta.data[i,81] <- "noone called it immune"
        #combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,75])
      }
    }
    
    else {
          combo.reference@meta.data[i,81] <- as.character(combo.reference@meta.data[i,75])

    }
  }
}

head(combo.reference@meta.data)
table(combo.reference$celltype_final)

#problem child : Aziziimmune_BC1_55819


pdf("test.pdf", width = 14.22, height = 13.01)
DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE) + ggtitle("")
DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE, group.by = "RefineMarkerBased") + ggtitle("")
DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE, group.by = "RefineExpBased") + ggtitle("")
DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE, group.by = "celltype_final") + ggtitle("")
DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE, group.by = "celltype_minor") + ggtitle("") +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1)) 

# DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE, group.by = "celltype_final",
#         order = c("B Cells", "Plasma Cells", "Epithelial Cells", 
#                   "CD4+ T Cells", "CD8+ T Cells", "Regulatory T Cells", "NK Cells",
#                   "Neutrophils", "Monocytes", "Macrophages", "MDSCs", "Dendritic Cells","Mast Cells",
#                   "Fibroblasts", "Myoepithelial Cells", "Perivascular-like (PVL) Cells", "Endothelial Cells"),
#         cols = c("#c082ba",
#                  "#69ba3d",
#                  "#6d3ac1",
#                  "#b4a645",
#                  "#ca4abf",
#                  "#61ae6e",
#                  "#6d6ecc",
#                  "#d24932",
#                  "#6baca2",
#                  "#cc4772",
#                  "#526531",
#                  "#4e2a61",
#                  "#c57d3c",
#                  "#7091b9",
#                  "#6f322e",
#                  "#c09086",
#                  "#353c3c")) + ggtitle("")
dev.off()

Idents(combo.reference) <- combo.reference$celltype_final

table(Idents(combo.reference))

Savas <- subset(combo.reference, subset = orig.ident == "Savas")
head(Savas@meta.data)
table(Idents(Savas))
table(Savas$RefineMarkerBased)
table(Savas$RefineExpBased)
table(Savas$celltype_UCell)


AziziT <- subset(combo.reference, subset = orig.ident == "AziziT")
head(AziziT@meta.data)
table(Idents(AziziT))

Aziziprim <- subset(combo.reference, subset = orig.ident == "Aziziimmune")
head(Aziziprim@meta.data)
table(Idents(Aziziprim))
table(Aziziprim$celltype_UCell)

aziziprim.test <- subset(Aziziprim, idents = "Epithelial Cells")

head(aziziprim.test@meta.data[,c(75, 79:81)])
table(aziziprim.test$RefineExpBased)
table(aziziprim.test$RefineMarkerBased)
table(aziziprim.test$celltype_UCell)
table(aziziprim.test$celltype_final)

#problem child : Aziziimmune_BC1_55472


combo.reference <- subset(combo.reference, idents = "noone called it immune", invert = T)




