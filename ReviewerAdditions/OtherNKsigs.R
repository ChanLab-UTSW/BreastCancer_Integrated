


NKsubDir <- "/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only"
ReviewerDir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/ReviewerAdditions"
 

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







#load nk subset =====

setwd(NKsubDir)
NK.all.combo <- readRDS("NK_withNKpathways_noZallgenedem_71422.rds") #newnewnew
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")


#looking for other NK subtypes ================== =========================================

# NK activating/inhibiting  ------------

#https://static-content.springer.com/esm/art%3A10.1038%2Fni1581/MediaObjects/41590_2008_BFni1581_MOESM203_ESM.pdf

#httr::GET("https://www.genenames.org/help/rest/", config = httr::config(connecttimeout = 120))
GeneSymbolThesarus("CD244")

DefaultAssay(NK.all.combo) <- "RNA"

##BOTH
CD244 <- list(c("CD244"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = CD244, name = "CD244", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1CD244", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KLRD1 <- list(c("KLRD1"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KLRD1, name = "KLRD1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KLRD1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

##inhibitory

KLRG1 <- list(c("KLRG1"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KLRG1, name = "KLRG1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KLRG1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KLRC1 <- list(c("KLRC1"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KLRC1, name = "KLRC1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KLRC1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KLRB1 <- list(c("KLRB1"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KLRB1, name = "KLRB1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KLRB1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

LAIR1 <- list(c("LAIR1"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = LAIR1, name = "LAIR1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1LAIR1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

LILRB1 <- list(c("LILRB1"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = LILRB1, name = "LILRB1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1LILRB1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KIR2DL3 <- list(c("KIR2DL3"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KIR2DL3, name = "KIR2DL3", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KIR2DL3", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KIR2DL2 <- list(c("KIR2DL2")) ##NA in counts assay ***
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KIR2DL2, name = "KIR2DL2", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KIR2DL2", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KIR2DL1 <- list(c("KIR2DL1"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KIR2DL1, name = "KIR2DL1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KIR2DL1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KIR3DL1 <- list(c("KIR3DL1"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KIR3DL1, name = "KIR3DL1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KIR3DL1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KIR2DL5A <- list(c("KIR2DL5A")) #**NA in counts assay
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KIR2DL5A, name = "KIR2DL5A", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KIR2DL5A", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KIR2DL5B <- list(c("KIR2DL5B")) #**NA in counts assay
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KIR2DL5B, name = "KIR2DL5B", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KIR2DL5B", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KIR3DL2 <- list(c("KIR3DL2"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KIR3DL2, name = "KIR3DL2", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KIR3DL2", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

CEACAM1 <- list(c("CEACAM1"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = CEACAM1, name = "CEACAM1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1CEACAM1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

SIGLEC7 <- list(c("SIGLEC7"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = SIGLEC7, name = "SIGLEC7", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1SIGLEC7", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)


##Activating
KLRC2 <- list(c("KLRC2"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KLRC2, name = "KLRC2", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KLRC2", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KLRK1 <- list(c("KLRK1"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KLRK1, name = "KLRK1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KLRK1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KLRF1 <- list(c("KLRF1"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KLRF1, name = "KLRF1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KLRF1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KIR2DS1 <- list(c("KIR2DS1")) ##****NA in genome
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KIR2DS1, name = "KIR2DS1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KIR2DS1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KIR2DS2 <- list(c("KIR2DS2")) ##NA in genome
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KIR2DS2, name = "KIR2DS2", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KIR2DS2", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KIR2DS4 <- list(c("KIR2DS4"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KIR2DS4, name = "KIR2DS4", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KIR2DS4", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KIR3DS1 <- list(c("KIR3DS1")) #NA in genome
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KIR3DS1, name = "KIR3DS1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KIR3DS1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

KIR2DL4 <- list(c("KIR2DL4"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = KIR2DL4, name = "KIR2DL4", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1KIR2DL4", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

NCR1 <- list(c("NCR1"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = NCR1, name = "NCR1", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1NCR1", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

NCR2 <- list(c("NCR2"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = NCR2, name = "NCR2", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1NCR2", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

NCR3 <- list(c("NCR3"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = NCR3, name = "NCR3", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1NCR3", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

FCGR3A <- list(c("FCGR3A"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = FCGR3A, name = "FCGR3A", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1FCGR3A", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

CD226 <- list(c("CD226"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = CD226, name = "CD226", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1CD226", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

SLAMF7 <- list(c("SLAMF7"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = SLAMF7, name = "SLAMF7", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1SLAMF7", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

SLAMF6 <- list(c("SLAMF6"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = SLAMF6, name = "SLAMF6", assay = "RNA")  
FeaturePlot(object = NK.all.combo, features = "signature_1SLAMF6", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)

#Moreno-Nieves ==============================

DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")

setwd(ReviewerDir)
MorenoNieves <- read.csv(file = "MorenoNievesTop20.csv")
head(MorenoNieves)


Peripheral <- MorenoNieves[,1]
Peripheral <- Peripheral[!is.na(Peripheral)]
Peripheral <- as.list(Peripheral)
Peripheral <- list(Peripheral)
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = Peripheral, name = "PeripheralNK_MorenoNieves", assay = "RNA")  


NK_1 <- MorenoNieves[,2]
NK_1 <- NK_1[!is.na(NK_1)]
NK_1 <- as.list(NK_1)
NK_1 <- list(NK_1)
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = NK_1, name = "NK_1_MorenoNieves", assay = "RNA")  


NK_2 <- MorenoNieves[,3]
NK_2 <- NK_2[!is.na(NK_2)]
NK_2 <- as.list(NK_2)
NK_2 <- list(NK_2)
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = NK_2, name = "NK_2_MorenoNieves", assay = "RNA")  


NK_ILC1_intermediate <- MorenoNieves[,4]
NK_ILC1_intermediate <- NK_ILC1_intermediate[!is.na(NK_ILC1_intermediate)]
NK_ILC1_intermediate <- as.list(NK_ILC1_intermediate)
NK_ILC1_intermediate <- list(NK_ILC1_intermediate)
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = NK_ILC1_intermediate, name = "NK_ILC1_interme_MorenoNieves", assay = "RNA")  

ILC1 <- MorenoNieves[,5]
ILC1 <- ILC1[!is.na(ILC1)]
ILC1 <- as.list(ILC1)
ILC1 <- list(ILC1)
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = ILC1, name = "ILC1_MorenoNieves", assay = "RNA")  

ILC3 <- MorenoNieves[,6]
ILC3 <- ILC3[!is.na(ILC3)]
ILC3 <- as.list(ILC3)
ILC3 <- list(ILC3)
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = ILC3, name = "ILC3_MorenoNieves", assay = "RNA")  

ieILC1 <- MorenoNieves[,7]
ieILC1 <- ieILC1[!is.na(ieILC1)]
ieILC1 <- as.list(ieILC1)
ieILC1 <- list(ieILC1)
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = ieILC1, name = "ieILC1_MorenoNieves", assay = "RNA")  

ieILC1.cycling <- MorenoNieves[,8]
ieILC1.cycling <- ieILC1.cycling[!is.na(ieILC1.cycling)]
ieILC1.cycling <- as.list(ieILC1.cycling)
ieILC1.cycling <- list(ieILC1.cycling)
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = ieILC1.cycling, name = "ieILC1.cycling_MorenoNieves", assay = "RNA")  


setwd(ReviewerDir)
pdf("MorenoNievesOverlay_11222.pdf", width = 5.11, height = 4.5)
FeaturePlot(object = NK.all.combo, features = "signature_1PeripheralNK_MorenoNieves", order = T, label = F, raster = FALSE)#min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)
FeaturePlot(object = NK.all.combo, features = "signature_1NK_1_MorenoNieves", order = TRUE, label = F, raster = FALSE)#min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)
FeaturePlot(object = NK.all.combo, features = "signature_1NK_2_MorenoNieves", order = TRUE, label = F, raster = FALSE)#min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)
FeaturePlot(object = NK.all.combo, features = "signature_1NK_ILC1_interme_MorenoNieves", order = T, label = F, raster = FALSE)#min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)
FeaturePlot(object = NK.all.combo, features = "signature_1ILC1_MorenoNieves", order = TRUE, label = F, raster = FALSE)#min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)
FeaturePlot(object = NK.all.combo, features = "signature_1ILC3_MorenoNieves", order = TRUE, label = F, raster = FALSE)# min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)
FeaturePlot(object = NK.all.combo, features = "signature_1ieILC1_MorenoNieves", order = TRUE, label = F, raster = FALSE)# min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)
FeaturePlot(object = NK.all.combo, features = "signature_1ieILC1.cycling_MorenoNieves", order = TRUE, label = F, raster = FALSE)# min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)
FeaturePlot(object = NK.all.combo, features = "signature_1CD56", order = TRUE, label = F, raster = FALSE)# min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)
FeaturePlot(object = NK.all.combo, features = "signature_1ILCIL3short", order = TRUE, label = F, raster = FALSE)# min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)
FeaturePlot(object = NK.all.combo, features = "signature_1NK2short", order = TRUE, label = F, raster = FALSE)# min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)
dev.off()





#Dogra (fig 4B) ================================

#https://www.sciencedirect.com/science/article/pii/S0092867420301033?via%3Dihub#mmc2

Dogralist <- c("ENC1", "TTC38", "PRF1", "CMKLR1",
               "BOK", "CX3CR1", "DGKK", "SLC1A7", "MYOM2",
               "MTSS1", "FGFBP2", "FCGR3A", "LGR6",
               "KIR2DL1", "KIR2DL3", "TGFBR3", "SPON2", "ADGRG1",
               "COL9A2", "GNAI1", "SRGAP3", "IL7R",
               "ITM2C", "XCL1", "SPRY1", "SYPL1",
               "SETD7", "RCAN3", "KIT", "GPR183",
               "AHR", "RASSF8", "FBXL16", "SIRPG",
               "TNFSF8", "UNC93B1", "PATJ", "SLC4A10",
               "IL23R", "GZMK", "JAML", "RUNX2",
               "INPP4B", "CCDC141", "CXCR6", "SPRY2",
               "ITGA1", "CAMK4", "MVB12B", "STYK1")
for (gene in Dogralist) {
  GeneSymbolThesarus(gene)
}

CD56briCD16min <- list(c("COL9A2", "GNAI1", "SRGAP3", "IL7R",
                         "ITM2C", "XCL1", "SPRY1", "SYPL1",
                         "SETD7", "RCAN3", "KIT", "GPR183",
                         "AHR", "RASSF8", "FBXL16", "SIRPG",
                         "TNFSF8", "UNC93B1", "PATJ", "SLC4A10",
                         "IL23R", "GZMK", "JAML", "RUNX2",
                         "INPP4B", "CCDC141", "CXCR6", "SPRY2",
                         "ITGA1", "CAMK4", "MVB12B", "STYK1", #up
                         "ENC1-", "TTC38-", "PRF1-", "CMKLR1-", #down
                         "BOK-", "CX3CR1-", "DGKK-", "SLC1A7-", "MYOM2-",
                         "MTSS1-", "FGFBP2-", "FCGR3A-", "LGR6-",
                         "KIR2DL1-", "KIR2DL3-", "TGFBR3-", "SPON2-", "ADGRG1-"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = CD56briCD16min, name = "CD56briCD16min", assay = "RNA")  


CD56dimCD16plus <- list(c("ENC1", "TTC38", "PRF1", "CMKLR1",
                          "BOK", "CX3CR1", "DGKK", "SLC1A7", "MYOM2",
                          "MTSS1", "FGFBP2", "FCGR3A", "LGR6",
                          "KIR2DL1", "KIR2DL3", "TGFBR3", "SPON2", "ADGRG1", #up
                          "COL9A2-", "GNAI1-", "SRGAP3-", "IL7R-", #down
                          "ITM2C-", "XCL1-", "SPRY1-", "SYPL1-",
                          "SETD7-", "RCAN3-", "KIT-", "GPR183-",
                          "AHR-", "RASSF8-", "FBXL16-", "SIRPG-",
                          "TNFSF8-", "UNC93B1-", "PATJ-", "SLC4A10-",
                          "IL23R-", "GZMK-", "JAML-", "RUNX2-",
                          "INPP4B-", "CCDC141-", "CXCR6-", "SPRY2-",
                          "ITGA1-", "CAMK4-", "MVB12B-", "STYK1-"))
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = CD56dimCD16plus, name = "CD56dimCD16plus", assay = "RNA")  


pdf("Dogra_CD56dimbright_11222.pdf", width = 5.11, height = 4.5)
FeaturePlot(object = NK.all.combo, features = "signature_1CD56briCD16min", order = TRUE, label = F, repel = TRUE, min.cutoff = 0, max.cutoff = 0.1, raster = FALSE)
FeaturePlot(object = NK.all.combo, features = "signature_1CD56dimCD16plus", order = TRUE, label = F, repel = TRUE, min.cutoff = 0, max.cutoff = 0.1, raster = FALSE)
dev.off()



#NK pathways ========================
#table S2 from https://bmccancer.biomedcentral.com/articles/10.1186/s12885-022-09230-y
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")


setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/GeneralAnalyses")

NKPathways <- read.csv("TableS2ChengLongGlioma.csv")
head(NKPathways)
colnames(NKPathways)

for (i in 1:ncol(NKPathways))
{
  pathway.name <- colnames(NKPathways)[i]
  
  pathway <- NKPathways[,i]
  pathway <- pathway[!is.na(pathway)]
  # for (gene in pathway) {
  #   GeneSymbolThesarus(gene)
  # }
  pathway <- as.list(pathway)
  pathway <- list(pathway)
  NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = pathway, name = pathway.name, assay = "RNA")  
  
}

colnames(NK.all.combo@meta.data)

NKgenes <- read.csv("TableS1ChengLongGlioma.csv")
head(NKgenes)

NKgenes <- NKgenes[,1]
NKgenes.list <- NKgenes[!is.na(NKgenes)]


for (gene in NKgenes.list)
{
  NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = gene, name = gene, assay = "RNA")  
}

colnames(NK.all.combo@meta.data)



