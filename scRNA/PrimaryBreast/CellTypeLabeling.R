

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# You will need the whole integrated scRNA-seq object
# This is SCTransformed and uses reference-based integration
# so after getting the data slot from the RNA assay, 
# make sure to run NormalizeData()

# The marker method is inspired by Karaayvaz et. al (2018)
# Paper link: https://www.nature.com/articles/s41467-018-06052-0 
# code link: https://github.com/Michorlab/tnbc_scrnaseq/blob/master/code/funcs_markers.R

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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





#Load in Integrated Primary Object =================

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects")
combo.reference <- readRDS("SCTRefwithUMAP_pc70res5_6222.rds")


#Cluster Labeling (UScore) =========================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This labels the clusters as a whole using the UCell package
# Paper: https://www.sciencedirect.com/science/article/pii/S2001037021002"celltype_final"6?via%3Dihub 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


DefaultAssay(combo.reference) <- "RNA"
#https://www.biostars.org/p/395951/

#https://github.com/satijalab/seurat/issues/5738
#https://github.com/satijalab/seurat/issues/5847

combo.reference <- NormalizeData(combo.reference, assay = "RNA")

DimPlot(combo.reference, reduction = "umap", label = TRUE, raster = FALSE)+ ggtitle(label = "Cell Type (PC = 40)") #+ NoLegend()


#NK __________________________________________________________________

NK.mark2 <- list(c("PTPRC","FGFBP2", "SPON2", "KLRF1", "FCGR3A", "KLRD1", "TRDC", "CD3E-",
                   "CD3G-", "CD33-", "EPCAM-",
                   "NCAM1")) #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
combo.reference <- AddModuleScore_UCell(combo.reference, features = NK.mark2, name = "NK.mark2", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = combo.reference, features = "signature_1NK.mark2", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "UCell's NK Signature (PC = 40)")
dev.off()


combo.reference <- RenameIdents(combo.reference, `175` = "NK Cells", `43` = "NK Cells", 
                                `128` = "NK Cells", `163` = "NK Cells", 
                                `98` = "NK Cells", `18` = "NK Cells", `64` = "NK Cells",
                                `178` = "NK Cells",`102` = "NK Cells")
NK <- subset(combo.reference, idents = "NK Cells")
NK <- AddModuleScore_UCell(NK, features = "CD3", name = "CD3", assay = "RNA")
NK <- AddModuleScore_UCell(NK, features = "CD3D", name = "CD3D", assay = "RNA")
NK <- AddModuleScore_UCell(NK, features = "CD8A", name = "CD8A", assay = "RNA")
NK <- AddModuleScore_UCell(NK, features = "CD4", name = "CD4", assay = "RNA")
NK <- AddModuleScore_UCell(NK, features = "CD33", name = "CD33", assay = "RNA")
NK <- AddModuleScore_UCell(NK, features = "EPCAM", name = "EPCAM", assay = "RNA")
NK <- AddModuleScore_UCell(NK, features = "PTPRC", name = "PTPRC", assay = "RNA")

summary(NK$signature_1NK.mark2)
summary(NK$signature_1CD3)
summary(NK$signature_1CD3D)
summary(NK$signature_1CD8A)
summary(NK$signature_1CD4)
summary(NK$signature_1CD33)
summary(NK$signature_1EPCAM)

noT8 <- WhichCells(NK, expression = signature_1CD8A > 0.4)
Idents(combo.reference, cells = noT8) <- "T Cells"
noTcd3d <- WhichCells(NK, expression = signature_1CD3D > 0.4)
Idents(combo.reference, cells = noTcd3d) <- "T Cells"
noTcd4 <- WhichCells(NK, expression = signature_1CD4 > 0.4)
Idents(combo.reference, cells = noTcd4) <- "T Cells"

notenoughNK <- WhichCells(NK, expression = signature_1NK.mark2 < 0.05)
Idents(combo.reference, cells = notenoughNK) <- "T Cells"


#T Cells __________________________________________________________________

T_karacd4 <- list(c("PTPRC", "CD2",
                    "CD3D",
                    "CD3E",
                    "CD3G",
                    "CD8A",
                    "CD8B",
                    "CD4",
                    "CD33-", "EPCAM-")) #myeloid progenitor
combo.reference <- AddModuleScore_UCell(combo.reference, features = T_karacd4, name = "T_karacd4", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = combo.reference, features = "signature_1T_karacd4", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Chung's T Cell Signature (PC = 40)")
dev.off()

combo.reference <- RenameIdents(combo.reference, `103` = "T Cells", `75` = "T Cells",
                                `12` = "T Cells", `209` = "T Cells", `208` = "T Cells",
                                `0` = "T Cells", `245` = "T Cells", `143` = "T Cells",
                                `190` = "T Cells", `249` = "T Cells", `21` = "T Cells",
                                `182` = "T Cells", `126` = "T Cells", `62` = "T Cells",
                                `4` = "T Cells", `66` = "T Cells", `5` = "T Cells",
                                `188` = "T Cells", `156` = "T Cells", `8` = "T Cells",
                                `146` = "T Cells", `59` = "T Cells", `250` = "T Cells",
                                `38` = "T Cells", `1` = "T Cells", `150` = "T Cells",
                                `28` = "T Cells", `213` = "T Cells", `38` = "T Cells",
                                `215` = "T Cells", `11` = "T Cells", `7` = "T Cells",
                                `48` = "T Cells", `37` = "T Cells", 
                                `20` = "T Cells", `176` = "T Cells", `9` = "T Cells",
                                `45` = "T Cells", `95` = "T Cells", `22` = "T Cells",
                                `118` = "T Cells", `67` = "T Cells", `238` = "T Cells",
                                `32` = "T Cells", `33` = "T Cells", 
                                `151` = "T Cells", `6` = "T Cells", `196` = "T Cells",
                                `15` = "T Cells", `100` = "T Cells", `232` = "T Cells",
                                `13` = "T Cells", `88` = "T Cells", `31` = "T Cells",
                                `51` = "T Cells", `123` = "T Cells", `198` = "T Cells",
                                `206` = "T Cells", `241` = "T Cells", `228` = "T Cells",
                                `60` = "T Cells", `40` = "T Cells", `24` = "T Cells",
                                `30` = "T Cells", `119` = "T Cells",
                                `99` = "T Cells", `220` = "T Cells")

#B Cells __________________________________________________________________

B_KaraBrech <- list(c("PTPRC", "CD79A",
                      "CD79B",
                      "BLNK",
                      "CD19",
                      "MS4A1",
                      "CD33-", "EPCAM-")) #alias for CD20

combo.reference <- AddModuleScore_UCell(combo.reference, features = B_KaraBrech, name = "B_KaraBrech", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = combo.reference, features = "signature_1B_KaraBrech", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) #+ ggtitle(label = "Chung's B Cell Signature (PC = 40)")
dev.off()

combo.reference <- RenameIdents(combo.reference, `219` = "B Cells", `159` = "B Cells", 
                                `216` = "B Cells", `74` = "B Cells", 
                                `172` = "B Cells", `246` = "B Cells", `2` = "B Cells",
                                `124` = "B Cells", `236` = "B Cells")


#plasma cells __________________________________________________________________

plasma_colo <- list(c("PTPRC","CD27",
                      "IGHA1",
                      "SDC1",
                      "TNFRSF17",
                      "JCHAIN",
                      "MZB1",
                      "DERL3",
                      "CD38",
                      "IGHG1",
                      "IGHG3",
                      "IGHG4",
                      "CD33-", "EPCAM-"))
combo.reference <- AddModuleScore_UCell(combo.reference, features = plasma_colo, name = "plasma_colo", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = combo.reference, features = "signature_1plasma_colo", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Ramachandran's Plasma Cell Signature (PC = 40)")
dev.off()

combo.reference <- RenameIdents(combo.reference, `135` = "Plasma Cells", `214` = "Plasma Cells",
                                `235` = "Plasma Cells", `52` = "Plasma Cells", `168` = "Plasma Cells",
                                `189` = "Plasma Cells", `210` = "Plasma Cells", `144` = "Plasma Cells")


#myeloid __________________________________________________________________


myeloid_wuElliot <- list(c("PTPRC",
                           "ITGAM",
                           "HLA-DRA",
                           "ITGAX",
                           "CD14",
                           "FCGR3A", #CD16
                           "CD1C",
                           "CD1A",
                           "CD68",
                           "CD33",
                           "EPCAM-")) #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
combo.reference <- AddModuleScore_UCell(combo.reference, features = myeloid_wuElliot, name = "myeloid_wuElliot", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = combo.reference, features = "signature_1myeloid_wuElliot", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) #+ ggtitle(label = "Wu's Myeloid Cell Signature (PC = 40)")
dev.off()

combo.reference <- RenameIdents(combo.reference,`231` = "Myeloid Cells", `109` = "Myeloid Cells",#mast cells
                                `225` = "Myeloid Cells", `116` = "Myeloid Cells", `78` = "Myeloid Cells",
                                `239` = "Myeloid Cells", `173` = "Myeloid Cells", `139` = "Myeloid Cells",
                                `193` = "Myeloid Cells", `136` = "Myeloid Cells", `205` = "Myeloid Cells",
                                `134` = "Myeloid Cells",`177` = "Myeloid Cells", `19` = "Myeloid Cells",
                                `131` = "Myeloid Cells", `141` = "Myeloid Cells", `106` = "Myeloid Cells",
                                `77` = "Myeloid Cells", `251` = "Myeloid Cells", `10` = "Myeloid Cells",
                                `130` = "Myeloid Cells", `138` = "Myeloid Cells", `107` = "Myeloid Cells", 
                                `97` = "Myeloid Cells") 


#Non-Immune _________________________________________________________


#gen epi __________________________________________________________________

genepi_baslumgen <- list(c("EGFR",
                           "FZR1",
                           "KRT14",
                           "ITGA6",
                           "KRT5",
                           "TP63",
                           "KRT17",
                           "MME",
                           "FOXA1",
                           "GATA3",
                           "MUC1",
                           "CD24",
                           "GABRP", 
                           "EPCAM",
                           "KRT8",
                           "KRT18",
                           "KRT19",
                           "PTPRC-"))

combo.reference <- AddModuleScore_UCell(combo.reference, features = genepi_baslumgen, name = "genepi_baslumgen", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = combo.reference, features = "signature_1genepi_baslumgen", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Karaayvaz's General Epithelial Signature (PC = 40)")
dev.off()

combo.reference <- RenameIdents(combo.reference, `3` = "Epithelial Cells", `17` = "Epithelial Cells",
                                `23` = "Epithelial Cells", `26` = "Epithelial Cells", `35` = "Epithelial Cells",
                                `36` = "Epithelial Cells", `39` = "Epithelial Cells", `42` = "Epithelial Cells",
                                `44` = "Epithelial Cells", `46` = "Epithelial Cells", `47` = "Epithelial Cells",
                                `49` = "Epithelial Cells", `50` = "Epithelial Cells",
                                `53` = "Epithelial Cells", `55` = "Epithelial Cells", `56` = "Epithelial Cells",
                                `58` = "Epithelial Cells", `61` = "Epithelial Cells", `65` = "Epithelial Cells",
                                `69` = "Epithelial Cells", `70` = "Epithelial Cells", `71` = "Epithelial Cells",
                                `72` = "Epithelial Cells", `73` = "Epithelial Cells", `79` = "Epithelial Cells",
                                `80` = "Epithelial Cells", `81` = "Epithelial Cells", `84` = "Epithelial Cells",
                                `85` = "Epithelial Cells", `87` = "Epithelial Cells", `89` = "Epithelial Cells",
                                `90` = "Epithelial Cells", `91` = "Epithelial Cells", `92` = "Epithelial Cells",
                                `94` = "Epithelial Cells", `101` = "Epithelial Cells", `105` = "Epithelial Cells",
                                `108` = "Epithelial Cells", `110` = "Epithelial Cells", `111` = "Epithelial Cells",
                                `112` = "Epithelial Cells", `115` = "Epithelial Cells", `117` = "Epithelial Cells",
                                `121` = "Epithelial Cells", `122` = "Epithelial Cells", `125` = "Epithelial Cells",
                                `129` = "Epithelial Cells", `133` = "Epithelial Cells", `137` = "Epithelial Cells",
                                `140` = "Epithelial Cells", `145` = "Epithelial Cells", `147` = "Epithelial Cells",
                                `149` = "Epithelial Cells", `152` = "Epithelial Cells", `153` = "Epithelial Cells",
                                `157` = "Epithelial Cells", `158` = "Epithelial Cells", `160` = "Epithelial Cells",
                                `161` = "Epithelial Cells", `162` = "Epithelial Cells", `164` = "Epithelial Cells",
                                `165` = "Epithelial Cells", `166`= "Epithelial Cells", `167`= "Epithelial Cells", `169`= "Epithelial Cells",
                                `170`= "Epithelial Cells", `171`= "Epithelial Cells", `174`= "Epithelial Cells", `183`= "Epithelial Cells",
                                `184`= "Epithelial Cells", `186`= "Epithelial Cells",`187`= "Epithelial Cells", `191`= "Epithelial Cells",
                                `192`= "Epithelial Cells", `194`= "Epithelial Cells", `195`= "Epithelial Cells",
                                `197`= "Epithelial Cells", `199`= "Epithelial Cells", `200`= "Epithelial Cells",
                                `201`= "Epithelial Cells", `203`= "Epithelial Cells", `211`= "Epithelial Cells",
                                `212`= "Epithelial Cells", `217`= "Epithelial Cells", `221`= "Epithelial Cells",
                                `222`= "Epithelial Cells", `223`= "Epithelial Cells", `226`= "Epithelial Cells",
                                `227`= "Epithelial Cells", `229`= "Epithelial Cells", `230`= "Epithelial Cells",
                                `234`= "Epithelial Cells", `240`= "Epithelial Cells", `243`= "Epithelial Cells",
                                `244`= "Epithelial Cells", `247`= "Epithelial Cells")

#myoepi __________________________________________________________________

Myoepi_wuNguyen <- list(c("KRT5",
                          "KRT14",
                          "ACTA2",
                          "TAGLN",
                          "EPCAM",
                          "PTPRC-"))
combo.reference <- AddModuleScore_UCell(combo.reference, features = Myoepi_wuNguyen, name = "Myoepi_wuNguyen", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = combo.reference, features = "signature_1Myoepi_wuNguyen", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Wu's Myoepithelial Signature (PC = 40)")
dev.off()
FeaturePlot(object = combo.reference, features = "EPCAM", order = TRUE, label = TRUE, min.cutoff = 0,  raster = FALSE) + ggtitle(label = "Wu's Myoepithelial Signature (PC = 40)")

combo.reference <- RenameIdents(combo.reference, `114` = "Myoepithelial Cells")

#PVL _________________________

PVL.mark <- list(c("MCAM", "ACTA2", "PDGFRB", "PTPRC-", "EPCAM-")) #from "Stromal subclasses resemble diverse..." (new Wu 2021)
combo.reference <- AddModuleScore_UCell(combo.reference, features = PVL.mark, name = "newWuPVL", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = combo.reference, features = "signature_1newWuPVL", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "New Wu's PVL Signature (PC = 40)")
dev.off()
combo.reference <- RenameIdents(combo.reference, `248` = "Perivascular-like (PVL) Cells",
                                `179` = "Perivascular-like (PVL) Cells", `224` = "Perivascular-like (PVL) Cells",
                                `185` = "Perivascular-like (PVL) Cells", `41` = "Perivascular-like (PVL) Cells",
                                `54` = "Perivascular-like (PVL) Cells", `132` = "Perivascular-like (PVL) Cells",
                                `93` = "Perivascular-like (PVL) Cells")

#Fibroblasts ________________________-

fibro_wumelan <- list(c("FAP",
                        "THY1",
                        "DCN",
                        "COL1A1",
                        "COL1A2",
                        "COL6A1",
                        "COL6A2",
                        "COL6A3",
                        "EPCAM-",
                        "PTPRC-"))
combo.reference <- AddModuleScore_UCell(combo.reference, features = fibro_wumelan, name = "fibro_wumelan", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = combo.reference, features = "signature_1fibro_wumelan", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Karaayvaz's Stroma Signature (PC = 40)")
dev.off()

combo.reference <- RenameIdents(combo.reference, `68` = "Fibroblasts",`204` = "Fibroblasts", 
                                `63` = "Fibroblasts", `27` = "Fibroblasts", `113` = "Fibroblasts",
                                `154` = "Fibroblasts", `34` = "Fibroblasts", `233` = "Fibroblasts",
                                `25` = "Fibroblasts", `16` = "Fibroblasts", `207` = "Fibroblasts",
                                `218` = "Fibroblasts", `202` = "Fibroblasts", `76` = "Fibroblasts",
                                `83` = "Fibroblasts", `86` = "Fibroblasts", `181` = "Fibroblasts",
                                `14` = "Fibroblasts", `155` = "Fibroblasts", `237` = "Fibroblasts",
                                `120` = "Fibroblasts", `127` = "Fibroblasts")

#emdo __________________________________________________________________

endo_kara <- list(c("PECAM1",
                    "VWF",
                    "CDH5",
                    "SELE",
                    "PTPRC-",
                    "EPCAM-"))
combo.reference <- AddModuleScore_UCell(combo.reference, features = endo_kara, name = "endokara", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = combo.reference, features = "signature_1endokara", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Karaayvaz's Endothelial Signature (PC = 40)")
dev.off()
combo.reference <- RenameIdents(combo.reference, `29` = "Endothelial Cells", `57` = "Endothelial Cells",
                                `82` = "Endothelial Cells", `148` = "Endothelial Cells", `180` = "Endothelial Cells",
                                `242` = "Endothelial Cells", `104` = "Endothelial Cells",
                                `142` = "Endothelial Cells", `96` = "Endothelial Cells")

DimPlot(combo.reference, reduction = "umap", label = FALSE, raster = FALSE, 
        order = c("NK Cells", "Epithelial Cells", "T Cells", "Stroma", "B Cells",
                  "Myoepithelial Cells", "Endothelial Cells", "Plasma Cells", "Myeloid Cells")) #+ ggtitle(label = "Integrated Primary Breast Cancer Sets (PC = 60)")


Savas <- subset(combo.reference, subset = orig.ident == "Savas")
Idents(Savas) <- "T Cells"
Savas.cells <- WhichCells(Savas, idents = "T Cells")
AziziT <- subset(combo.reference, subset = orig.ident == "AziziT")
Idents(AziziT) <- "T Cells"
AziziT.cells <- WhichCells(AziziT, idents = "T Cells")

Idents(combo.reference, Savas.cells) <- "T Cells"
Idents(combo.reference, AziziT.cells) <- "T Cells"


#T cell subsets ========================

DefaultAssay(combo.reference) <- "RNA"
combo.reference <- NormalizeData(combo.reference, assay = "RNA")

Tcellsub <- subset(combo.reference, idents = "T Cells")

table(Tcellsub$Capture.Method)
table(Tcellsub$orig.ident)
DefaultAssay(Tcellsub) <- "RNA"
Tcellsub <- NormalizeData(Tcellsub, assay = "RNA")

Tcellsub.list <- SplitObject(Tcellsub, split.by = "Capture.Method")
for (i in 1:length(Tcellsub.list)) {
  Tcellsub.list[[i]] <- SCTransform(Tcellsub.list[[i]], verbose = T, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                                                         "percent.platelet", "percent.heatshock"))
}

T.features <- SelectIntegrationFeatures(object.list = Tcellsub.list, nfeatures = 3000)
Tcellsub.list <- PrepSCTIntegration(object.list = Tcellsub.list, anchor.features = T.features)


reference.1 <- which(names(Tcellsub.list) == c("10X Genomics Chromium"))
reference.3 <- which(names(Tcellsub.list) == c("10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"))
reference.5 <- which(names(Tcellsub.list) == c("10X Genomics Chromium v2 5'"))
reference.6 <- which(names(Tcellsub.list) == c("10X Genomics Single Cell 3' v2"))

reference.list <- c(reference.1,reference.3, 
                    reference.5, reference.6)

T.anchors <- FindIntegrationAnchors(object.list = Tcellsub.list, normalization.method = "SCT",
                                    anchor.features = T.features, reference = reference.list)

Tcellsub.combo <- IntegrateData(anchorset = T.anchors, normalization.method = "SCT")


DefaultAssay(Tcellsub.combo) <- "integrated"

Tcellsub.combo <- RunPCA(Tcellsub.combo, npcs = 100, verbose = FALSE)

pdf("BatchTestTsub_62922.pdf", width = 14.22, height = 13.01)
DimPlot(Tcellsub.combo, reduction = "pca", raster = F, group.by = "Capture.Method")
DimPlot(Tcellsub.combo, reduction = "pca", raster = F, group.by = "orig.ident")

DimPlot(Tcellsub.combo, reduction = "pca", raster = F, group.by = "Capture.Method",
        cols = c("Transparent", "Transparent", "Transparent", "blue", "yellow", "green", "purple"))
dev.off()

pdf("test.pdf")
DimHeatmap(Tcellsub.combo, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(Tcellsub.combo, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(Tcellsub.combo, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(Tcellsub.combo, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(Tcellsub.combo, dims = 45:55, cells = 500, balanced = TRUE)
DimHeatmap(Tcellsub.combo, dims = 55:65, cells = 500, balanced = TRUE)
dev.off()

DefaultAssay(Tcellsub.combo) <- "integrated"

Tcellsub.combo <- FindNeighbors(Tcellsub.combo, reduction = "pca", dims = 1:60)
Tcellsub.combo <- FindClusters(Tcellsub.combo, resolution = 9)
Tcellsub.combo <- RunUMAP(Tcellsub.combo, reduction = "pca", dims = 1:60, verbose = TRUE, seed.use=123)

pdf("test.pdf", width = 14.22, height = 13.01)
DimPlot(Tcellsub.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE) #+ ggtitle(label = "NK Cells from Primary Integrated")
DimPlot(Tcellsub.combo, reduction = "umap", label = TRUE, repel = TRUE, group.by = "integrated_snn_res.5", raster = FALSE) +theme(legend.position = "None")#+ ggtitle(label = "NK Cells from Primary Integrated")
dev.off()
DimPlot(Tcellsub.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "celltype_main") + ggtitle(label = "NK Cells from lumA samples (PC = 15)")
DimPlot(Tcellsub.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "orig.ident")


DefaultAssay(Tcellsub.combo) <- "RNA"
Tcellsub.combo <- NormalizeData(Tcellsub.combo, assay = "RNA")

#CD8 __________________________________________________________________

CD8sig <- list(c("CD8A",#Wu
                 "CD8B"))
Tcellsub.combo <- AddModuleScore_UCell(Tcellsub.combo, features = CD8sig, name = "CD8sig", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = Tcellsub.combo, features = "signature_1CD8sig", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Karaayvaz's Endothelial Signature (PC = 40)")
dev.off()

Tcellsub.combo <- RenameIdents(Tcellsub.combo, `94` = "CD8+ T Cells",
                               `44` = "CD8+ T Cells", `28` = "CD8+ T Cells",
                               `36` = "CD8+ T Cells", `103` = "CD8+ T Cells",
                               `89` = "CD8+ T Cells", `64` = "CD8+ T Cells",
                               `75` = "CD8+ T Cells", `26` = "CD8+ T Cells",
                               `85` = "CD8+ T Cells", `117` = "CD8+ T Cells",
                               `116` = "CD8+ T Cells", `119` = "CD8+ T Cells",
                               `96` = "CD8+ T Cells", `8` = "CD8+ T Cells",
                               `97` = "CD8+ T Cells", `61` = "CD8+ T Cells",
                               `65` = "CD8+ T Cells", `58` = "CD8+ T Cells",
                               `101` = "CD8+ T Cells", `118` = "CD8+ T Cells",
                               `34` = "CD8+ T Cells", `77` = "CD8+ T Cells",
                               `104` = "CD8+ T Cells", `27` = "CD8+ T Cells",
                               `6` = "CD8+ T Cells", `2` = "CD8+ T Cells",
                               `105` = "CD8+ T Cells", `13` = "CD8+ T Cells",
                               `73` = "CD8+ T Cells", `21` = "CD8+ T Cells",
                               `60` = "CD8+ T Cells", `83` = "CD8+ T Cells",
                               `55` = "CD8+ T Cells", `76` = "CD8+ T Cells",
                               `95` = "CD8+ T Cells", `15` = "CD8+ T Cells",
                               `76` = "CD8+ T Cells", `95` = "CD8+ T Cells",
                               `112` = "CD8+ T Cells", `57` = "CD8+ T Cells",
                               `30` = "CD8+ T Cells", `106` = "CD8+ T Cells",
                               `78` = "CD8+ T Cells", `40` = "CD8+ T Cells",
                               `72` = "CD8+ T Cells", `17` = "CD8+ T Cells",
                               `109` = "CD8+ T Cells", `81` = "CD8+ T Cells",
                               `37` = "CD8+ T Cells", `20` = "CD8+ T Cells",
                               `31` = "CD8+ T Cells", `81` = "CD8+ T Cells",
                               `50` = "CD8+ T Cells", `53` = "CD8+ T Cells")

#CD4 __________________________________________________________________

CD4sig <- list(c("CD4"))
Tcellsub.combo <- AddModuleScore_UCell(Tcellsub.combo, features = CD4sig, name = "CD4sig", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = Tcellsub.combo, features = "signature_1CD4sig", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Karaayvaz's Endothelial Signature (PC = 40)")
dev.off()
Tcellsub.combo <- RenameIdents(Tcellsub.combo, `0` = "CD4+ T Cells",
                               `1` = "CD4+ T Cells", `3` = "CD4+ T Cells",
                               `5` = "CD4+ T Cells", `7` = "CD4+ T Cells",
                               `9` = "CD4+ T Cells", `10` = "CD4+ T Cells",
                               `11` = "CD4+ T Cells", `12` = "CD4+ T Cells",
                               `14` = "CD4+ T Cells", `16` = "CD4+ T Cells",
                               `18` = "CD4+ T Cells", `23` = "CD4+ T Cells",
                               `24` = "CD4+ T Cells", `29` = "CD4+ T Cells",
                               `32` = "CD4+ T Cells", `33` = "CD4+ T Cells",
                               `35` = "CD4+ T Cells", `38` = "CD4+ T Cells",
                               `39` = "CD4+ T Cells", `41` = "CD4+ T Cells",
                               `43` = "CD4+ T Cells", `45` = "CD4+ T Cells",
                               `46` = "CD4+ T Cells", `47` = "CD4+ T Cells",
                               `49` = "CD4+ T Cells", `51` = "CD4+ T Cells",
                               `52` = "CD4+ T Cells", `54` = "CD4+ T Cells",
                               `59` = "CD4+ T Cells", `62` = "CD4+ T Cells",
                               `63` = "CD4+ T Cells", `67` = "CD4+ T Cells",
                               `68` = "CD4+ T Cells", `69` = "CD4+ T Cells",
                               `70` = "CD4+ T Cells", `74` = "CD4+ T Cells",
                               `79` = "CD4+ T Cells", `80` = "CD4+ T Cells",
                               `82` = "CD4+ T Cells", `84` = "CD4+ T Cells",
                               `87` = "CD4+ T Cells", `88` = "CD4+ T Cells",
                               `90` = "CD4+ T Cells", `91` = "CD4+ T Cells",
                               `92` = "CD4+ T Cells", `93` = "CD4+ T Cells",
                               `98` = "CD4+ T Cells", `99` = "CD4+ T Cells",
                               `102` = "CD4+ T Cells", `107` = "CD4+ T Cells",
                               `108` = "CD4+ T Cells", `114` = "CD4+ T Cells",
                               `115` = "CD4+ T Cells")



#Treg __________________________________________________________________

Treg_sig <- list(c("FOXP3"))
Tcellsub.combo <- AddModuleScore_UCell(Tcellsub.combo, features = Treg_sig, name = "Treg_sig", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = Tcellsub.combo, features = "signature_1Treg_sig", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Karaayvaz's Endothelial Signature (PC = 40)")
dev.off()
Tcellsub.combo <- RenameIdents(Tcellsub.combo, `4` = "Regulatory T Cells", `110` = "Regulatory T Cells",
                               `22` = "Regulatory T Cells", `71` = "Regulatory T Cells",
                               `48` = "Regulatory T Cells", `100` = "Regulatory T Cells",
                               `66` = "Regulatory T Cells", `86` = "Regulatory T Cells",
                               `42` = "Regulatory T Cells", `56` = "Regulatory T Cells",
                               `25` = "Regulatory T Cells", `111` = "Regulatory T Cells",
                               `113` = "Regulatory T Cells", `19` = "Regulatory T Cells")



CD8 <- WhichCells(Tcellsub.combo, idents = "CD8+ T Cells")
CD4 <- WhichCells(Tcellsub.combo, idents = "CD4+ T Cells")
Treg <- WhichCells(Tcellsub.combo, idents = "Regulatory T Cells")

Idents(combo.reference, cells = CD8) <- "CD8+ T Cells"
Idents(combo.reference, cells = CD4) <- "CD4+ T Cells"
Idents(combo.reference, cells = Treg) <- "Regulatory T Cells"

#myeloid cell subset ==============================

DefaultAssay(combo.reference) <- "RNA"
combo.reference <- NormalizeData(combo.reference, assay = "RNA")

Myecellsub <- subset(combo.reference, idents = "Myeloid Cells")

table(Myecellsub$Capture.Method)
table(Myecellsub$orig.ident)
DefaultAssay(Myecellsub) <- "RNA"
Myecellsub <- NormalizeData(Myecellsub, assay = "RNA")

Myecellsub.list <- SplitObject(Myecellsub, split.by = "Capture.Method")
for (i in names(Myecellsub.list)) {
  Myecellsub.list[[i]] <- SCTransform(Myecellsub.list[[i]], verbose = T, 
                                      vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                          "percent.platelet", "percent.heatshock"))
}


Mye.features <- SelectIntegrationFeatures(object.list = Myecellsub.list, nfeatures = 3000)
Myecellsub.list <- PrepSCTIntegration(object.list = Myecellsub.list, anchor.features = Mye.features)


reference.1 <- which(names(Myecellsub.list) == c("10X Genomics Chromium"))
reference.3 <- which(names(Myecellsub.list) == c("10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"))
reference.5 <- which(names(Myecellsub.list) == c("10X Genomics Chromium v2 5'"))
reference.6 <- which(names(Myecellsub.list) == c("10X Genomics Single Cell 3' v2"))

reference.list <- c(reference.1,reference.3, 
                    reference.5, reference.6)

Mye.anchors <- FindIntegrationAnchors(object.list = Myecellsub.list, normalization.method = "SCT",
                                      anchor.features = Mye.features, reference = reference.list)

Myecellsub.combo <- IntegrateData(anchorset = Mye.anchors, normalization.method = "SCT", k.weight = 62)

DefaultAssay(Myecellsub.combo) <- "integrated"

Myecellsub.combo <- RunPCA(Myecellsub.combo, npcs = 100, verbose = FALSE)

pdf("BatchTestMyesub_62822.pdf", width = 14.22, height = 13.01)
DimPlot(Myecellsub.combo, reduction = "pca", raster = F, group.by = "Capture.Method")
DimPlot(Myecellsub.combo, reduction = "pca", raster = F, group.by = "orig.ident")

DimPlot(Myecellsub.combo, reduction = "pca", raster = F, group.by = "Capture.Method",
        cols = c("Transparent", "Transparent", "Transparent", "blue", "yellow", "green", "purple"))
dev.off()

pdf("test.pdf")
DimHeatmap(Myecellsub.combo, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(Myecellsub.combo, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(Myecellsub.combo, dims = 25:35, cells = 500, balanced = TRUE)
DimHeatmap(Myecellsub.combo, dims = 35:45, cells = 500, balanced = TRUE)
DimHeatmap(Myecellsub.combo, dims = 45:55, cells = 500, balanced = TRUE)
DimHeatmap(Myecellsub.combo, dims = 55:65, cells = 500, balanced = TRUE)
dev.off()


DefaultAssay(Myecellsub.combo) <- "integrated"
Myecellsub.combo <- FindNeighbors(Myecellsub.combo, reduction = "pca", dims = 1:50)
Myecellsub.combo <- FindClusters(Myecellsub.combo, resolution = 3)
Myecellsub.combo <- RunUMAP(Myecellsub.combo, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use=123)

pdf("test.pdf")
DimPlot(Myecellsub.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE) #+ ggtitle(label = "NK Cells from Primary Integrated")
dev.off()

DefaultAssay(Myecellsub.combo) <- "RNA"
Myecellsub.combo <- NormalizeData(Myecellsub.combo, assay = "RNA")

#MDSc __________________________________________________________________

#https://link.springer.com/article/10.1007/s00262-011-1161-9/tables/1 <- other

MDSc_sig <- list(c("CD33", #https://www.nature.com/articles/ncomms12150; table 1
                   "ITGAM", #CD11b, https://www.rndsystems.com/resources/cell-markers/immune-cells/myeloid--derived-suppressor-cells/monocytic-myeloid--derived-suppressor-cell--mdsc--markers
                   "LY6G",
                   "LY6C",
                   "FUT4", #CD15
                   "CEACAM8",
                   "IL4R",
                   "HLA-DRA",
                   "CD3D-",
                   "CD14-",
                   "CD19-",
                   "NCAM1-", #CD56
                   "EPCAM-")) #noepi
Myecellsub.combo <- AddModuleScore_UCell(Myecellsub.combo, features = MDSc_sig, name = "MDSc_sig", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = Myecellsub.combo, features = "signature_1MDSc_sig", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Karaayvaz's Endothelial Signature (PC = 40)")
dev.off()

Myecellsub.combo <-RenameIdents(Myecellsub.combo, `32` = "MDSCs", `19` = "MDSCs", `8` = "MDSCs",
                                `9` = "MDSCs", `24` = "MDSCs")

#Macrophage __________________________________________________________________

#https://link.springer.com/article/10.1007/s00262-011-1161-9/tables/1 <- other

Macro_Cheng <- list(c("INHBA", "IL1RN", "CCL4", "NLRP3",
                      "EREG", "IL1B", "LYVE1", "PLTP",
                      "SELENOP", "C1QC", "C1QA", "APOE"))

Myecellsub.combo <- AddModuleScore_UCell(Myecellsub.combo, features = Macro_Cheng, name = "Macro_Cheng", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = Myecellsub.combo, features = "signature_1Macro_Cheng", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) 
dev.off()

Myecellsub.combo <- RenameIdents(Myecellsub.combo, `0` = "Macrophages", `1` = "Macrophages", 
                                 `2` = "Macrophages", `3` = "Macrophages", 
                                 `4` = "Macrophages", `5` = "Macrophages", `6` = "Macrophages",
                                 `7` = "Macrophages", `10` = "Macrophages", `11` = "Macrophages",
                                 `12` = "Macrophages", `13` = "Macrophages", `16` = "Macrophages",
                                 `20` = "Macrophages", `21` = "Macrophages", `22` = "Macrophages",
                                 `23` = "Macrophages", `25` = "Macrophages", `26` = "Macrophages",
                                 `27` = "Macrophages", `28` = "Macrophages", `30` = "Macrophages",
                                 `36` = "Macrophages", `37` = "Macrophages", `38` = "Macrophages",
                                 `39` = "Macrophages", `41` = "Macrophages", `42` = "Macrophages",
                                 `43` = "Macrophages")
#Monocytes __________________________________________________________________

#https://link.springer.com/article/10.1007/s00262-011-1161-9/tables/1 <- other

Mono_sig <- list(c("FCN1", "S100A9", "S100A8", "FCGR3A", "LST1", "LILRB2", 
                   "HLA-DRA-")) #https://aacrjournals.org/cancerimmunolres/article/5/1/3/468526/Myeloid-Derived-Suppressor-CellsMasters-of

Myecellsub.combo <- AddModuleScore_UCell(Myecellsub.combo, features = Mono_sig, name = "Mono_sig", assay = "RNA")
pdf("test.pdf")
FeaturePlot(object = Myecellsub.combo, features = "signature_1Mono_sig", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE)
dev.off()

Myecellsub.combo <-RenameIdents(Myecellsub.combo, `18` = "Monocytes", `17` = "Monocytes", `14` = "Monocytes")


#Neutrophil __________________________________________________________________

#https://link.springer.com/article/10.1007/s00262-011-1161-9/tables/1 <- other


Neutro_sig <- list(c("CXCR1", #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
                     "CD15", #https://www.rndsystems.com/resources/cell-markers/immune-cells/granulocytes/neutrophil-cell-markers
                     "FCGR3A", #CD16
                     "CD14-",
                     "CSF3R", #TCGA, Chung
                     "S100A9","CD24A","TNF","CD274")) #scrna seq class
Myecellsub.combo <- AddModuleScore_UCell(Myecellsub.combo, features = Neutro_sig, name = "Neutro_sig", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = Myecellsub.combo, features = "signature_1Neutro_sig", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) 
dev.off()

Myecellsub.combo <- RenameIdents(Myecellsub.combo, `34` = "Neutrophils")

#Dendritic Cells __________________________________________________________________

#https://link.springer.com/article/10.1007/s00262-011-1161-9/tables/1 <- other


DC_sig <- list(c("LILRA4","GZMB","IL3RA",
                 "CLEC9A", "FLT3", "IDO1", 
                 "CD1C", "FCER1A","HLA-DQA1", 
                 "LAMP3", "CCR7", "FSCN1")) #Cheng

Myecellsub.combo <- AddModuleScore_UCell(Myecellsub.combo, features = DC_sig, name = "DC_sig", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = Myecellsub.combo, features = "signature_1DC_sig", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) 
dev.off()

Myecellsub.combo <- RenameIdents(Myecellsub.combo, `29` = "Dendritic Cells",
                                 `35` = "Dendritic Cells")

#mast cells __________________________________________________________________

mast_Cheng <- list(c("KIT","TPSAB1","CPA4"))

Myecellsub.combo <- AddModuleScore_UCell(Myecellsub.combo, features = mast_Cheng, name = "mast_Cheng", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = Myecellsub.combo, features = "signature_1mast_Cheng", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Chung's Mast Cell Signature (PC = 40)")
dev.off()

Myecellsub.combo <- RenameIdents(Myecellsub.combo, `33` = "Mast Cells", `15` = "Mast Cells",
                                 `31` = "Mast Cells", `40` = "Mast Cells")



MDSC <- WhichCells(Myecellsub.combo, idents = "MDSCs")
Macro <- WhichCells(Myecellsub.combo, idents = "Macrophages")
Mono <- WhichCells(Myecellsub.combo, idents = "Monocytes")
Neut <- WhichCells(Myecellsub.combo, idents = "Neutrophils")
DC <- WhichCells(Myecellsub.combo, idents = "Dendritic Cells")
Mast <- WhichCells(Myecellsub.combo, idents = "Mast Cells")

Idents(combo.reference, cells = MDSC) <- "MDSCs"
Idents(combo.reference, cells = Macro) <- "Macrophages"
Idents(combo.reference, cells = Mono) <- "Monocytes"
Idents(combo.reference, cells = Neut) <- "Neutrophils"
Idents(combo.reference, cells = DC) <- "Dendritic Cells"
Idents(combo.reference, cells = Mast) <- "Mast Cells"


combo.reference$celltype_UCell <- Idents(combo.reference)
#combo.reference$celltype_final <- Idents(combo.reference)

saveRDS(combo.reference, "finalfixed_primobj_withmarkersigs_62822_savecopysafe.rds")

#Marker and Expression Tests ========================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This will get the expression matrix of each celltype 
# as determined by UCell's cluster assignment

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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
  NK.exp.avg <- NK.exp.sum/nrow(combo.reference)  #average expression (denominator is #rows of dataset)
  #NK.exp.avg <- colMeans(NK.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #T Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Tcell.test <- as.data.frame(Subset.Mat[c(test.T, "PTPRC"), cell])
  colnames(Tcell.test) <- cell
  rownames(Tcell.test) <- c(test.T, "PTPRC")
  T.exp.sum <- colSums(Tcell.test)
  T.exp.avg <- T.exp.sum/nrow(combo.reference)  #average expression (denominator is #rows of dataset)
  #T.exp.avg <- colMeans(Tcell.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #B Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Bcell.test <- as.data.frame(Subset.Mat[c(test.B, "PTPRC"), cell])
  colnames(Bcell.test) <- cell
  rownames(Bcell.test) <- c(test.B, "PTPRC")
  B.exp.sum <- colSums(Bcell.test)
  B.exp.avg <- B.exp.sum/nrow(combo.reference)  #average expression (denominator is #rows of dataset)
  #B.exp.avg <- colMeans(Bcell.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #Plasma Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  Plasmacell.test <- as.data.frame(Subset.Mat[c(test.Plasma, "PTPRC"), cell])
  colnames(Plasmacell.test) <- cell
  rownames(Plasmacell.test) <- c(test.Plasma, "PTPRC")
  Plasma.exp.sum <- colSums(Plasmacell.test)
  Plasma.exp.avg <- Plasma.exp.sum/nrow(combo.reference)  #average expression (denominator is #rows of dataset)
  #Plasma.exp.avg <- colMeans(Plasmacell.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #Myeloid Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Myeloid.test <- as.data.frame(Subset.Mat[c(test.Mye, "PTPRC"), cell])
  colnames(Myeloid.test) <- cell
  rownames(Myeloid.test) <- c(test.Mye, "PTPRC")
  Myeloid.exp.sum <- colSums(Myeloid.test)
  Myeloid.exp.avg <- Myeloid.exp.sum/nrow(combo.reference) #average expression (denominator is #rows of dataset)
  #Myeloid.exp.avg <- colMeans(Myeloid.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #Myoepi Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  MyoEpi.test <- as.data.frame(Subset.Mat[c(test.MyoEpi), cell])
  colnames(MyoEpi.test) <- cell
  rownames(MyoEpi.test) <- c(test.MyoEpi)
  MyoEpi.exp.sum <- colSums(MyoEpi.test)
  MyoEpi.exp.avg <- MyoEpi.exp.sum/nrow(combo.reference) #average expression (denominator is #rows of dataset)
  #MyoEpi.exp.avg <- colMeans(MyoEpi.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #PVL Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  PVL.test <- as.data.frame(Subset.Mat[c(test.PVL), cell])
  colnames(PVL.test) <- cell
  rownames(PVL.test) <- c(test.PVL)
  PVL.exp.sum <- colSums(PVL.test)
  PVL.exp.avg <- PVL.exp.sum/nrow(combo.reference) #average expression (denominator is #rows of dataset)
  #PVL.exp.avg <- colMeans(PVL.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  #Fibro Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Fibro.test <- as.data.frame(Subset.Mat[c(test.Fibro), cell])
  colnames(Fibro.test) <- cell
  rownames(Fibro.test) <- c(test.Fibro)
  Fibro.exp.sum <- colSums(Fibro.test)
  Fibro.exp.avg <- Fibro.exp.sum/nrow(combo.reference) #average expression (denominator is #rows of dataset)
  #Fibro.exp.avg <- colMeans(Fibro.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  #Endothelial Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Endo.test <- as.data.frame(Subset.Mat[c(test.Endo), cell])
  colnames(Endo.test) <- cell
  rownames(Endo.test) <- c(test.Endo)
  Endo.exp.sum <- colSums(Endo.test)
  Endo.exp.avg <- Endo.exp.sum/nrow(combo.reference) #average expression (denominator is #rows of dataset)
  #Endo.exp.avg <- colMeans(Endo.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  #Epithelial Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Epi.test <- as.data.frame(Subset.Mat[c(test.Epi, test.StrongEpi), cell])
  colnames(Epi.test) <- cell
  rownames(Epi.test) <- c(test.Epi, test.StrongEpi)
  Epi.exp.sum <- colSums(Epi.test)
  Epi.exp.avg <- Epi.exp.sum/nrow(combo.reference) #average expression (denominator is #rows of dataset)
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
      CD8.exp.avg <- CD8.exp.sum/nrow(combo.reference) 
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
      CD4.exp.avg <- CD4.exp.sum/nrow(combo.reference) 
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
      Treg.exp.avg <- Treg.exp.sum/nrow(combo.reference) 
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
      MDSc.exp.avg <- MDSc.exp.sum/nrow(combo.reference) 
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
      Macro.exp.avg <- Macro.exp.sum/nrow(combo.reference) 
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
      Mono.exp.avg <- Mono.exp.sum/nrow(combo.reference) 
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
      Neut.exp.avg <- Neut.exp.sum/nrow(combo.reference) 
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
      DC.exp.avg <- DC.exp.sum/nrow(combo.reference) 
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
      Mast.exp.avg <- Mast.exp.sum/nrow(combo.reference) 
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
#test cell: Aziziimmune_BC3_4"RefineExpBased"68
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
  NK.exp.avg <- NK.exp.sum/nrow(combo.reference)  #average expression (denominator is #rows of dataset)
  #NK.exp.avg <- colMeans(NK.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #T Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Tcell.test <- as.data.frame(Subset.Mat[c(test.T, "PTPRC"), cell])
  colnames(Tcell.test) <- cell
  rownames(Tcell.test) <- c(test.T, "PTPRC")
  T.exp.sum <- colSums(Tcell.test)
  T.exp.avg <- T.exp.sum/nrow(combo.reference)  #average expression (denominator is #rows of dataset)
  #T.exp.avg <- colMeans(Tcell.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #B Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Bcell.test <- as.data.frame(Subset.Mat[c(test.B, "PTPRC"), cell])
  colnames(Bcell.test) <- cell
  rownames(Bcell.test) <- c(test.B, "PTPRC")
  B.exp.sum <- colSums(Bcell.test)
  B.exp.avg <- B.exp.sum/nrow(combo.reference)  #average expression (denominator is #rows of dataset)
  #B.exp.avg <- colMeans(Bcell.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #Plasma Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  Plasmacell.test <- as.data.frame(Subset.Mat[c(test.Plasma, "PTPRC"), cell])
  colnames(Plasmacell.test) <- cell
  rownames(Plasmacell.test) <- c(test.Plasma, "PTPRC")
  Plasma.exp.sum <- colSums(Plasmacell.test)
  Plasma.exp.avg <- Plasma.exp.sum/nrow(combo.reference)  #average expression (denominator is #rows of dataset)
  #Plasma.exp.avg <- colMeans(Plasmacell.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #Myeloid Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Myeloid.test <- as.data.frame(Subset.Mat[c(test.Mye, "PTPRC"), cell])
  colnames(Myeloid.test) <- cell
  rownames(Myeloid.test) <- c(test.Mye, "PTPRC")
  Myeloid.exp.sum <- colSums(Myeloid.test)
  Myeloid.exp.avg <- Myeloid.exp.sum/nrow(combo.reference) #average expression (denominator is #rows of dataset)
  #Myeloid.exp.avg <- colMeans(Myeloid.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #Myoepi Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  MyoEpi.test <- as.data.frame(Subset.Mat[c(test.MyoEpi), cell])
  colnames(MyoEpi.test) <- cell
  rownames(MyoEpi.test) <- c(test.MyoEpi)
  MyoEpi.exp.sum <- colSums(MyoEpi.test)
  MyoEpi.exp.avg <- MyoEpi.exp.sum/nrow(combo.reference) #average expression (denominator is #rows of dataset)
  #MyoEpi.exp.avg <- colMeans(MyoEpi.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  #PVL Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  PVL.test <- as.data.frame(Subset.Mat[c(test.PVL), cell])
  colnames(PVL.test) <- cell
  rownames(PVL.test) <- c(test.PVL)
  PVL.exp.sum <- colSums(PVL.test)
  PVL.exp.avg <- PVL.exp.sum/nrow(combo.reference) #average expression (denominator is #rows of dataset)
  #PVL.exp.avg <- colMeans(PVL.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  #Fibro Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Fibro.test <- as.data.frame(Subset.Mat[c(test.Fibro), cell])
  colnames(Fibro.test) <- cell
  rownames(Fibro.test) <- c(test.Fibro)
  Fibro.exp.sum <- colSums(Fibro.test)
  Fibro.exp.avg <- Fibro.exp.sum/nrow(combo.reference) #average expression (denominator is #rows of dataset)
  #Fibro.exp.avg <- colMeans(Fibro.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  #Endothelial Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Endo.test <- as.data.frame(Subset.Mat[c(test.Endo), cell])
  colnames(Endo.test) <- cell
  rownames(Endo.test) <- c(test.Endo)
  Endo.exp.sum <- colSums(Endo.test)
  Endo.exp.avg <- Endo.exp.sum/nrow(combo.reference) #average expression (denominator is #rows of dataset)
  #Endo.exp.avg <- colMeans(Endo.test) #average expression (taken from Cheng) #https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Signature_gene_sets_analysis.R
  
  
  
  #Epithelial Expression Sum +++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  Epi.test <- as.data.frame(Subset.Mat[c(test.Epi, test.StrongEpi), cell])
  colnames(Epi.test) <- cell
  rownames(Epi.test) <- c(test.Epi, test.StrongEpi)
  Epi.exp.sum <- colSums(Epi.test)
  Epi.exp.avg <- Epi.exp.sum/nrow(combo.reference) #average expression (denominator is #rows of dataset)
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
      CD8.exp.avg <- CD8.exp.sum/nrow(combo.reference) 
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
      CD4.exp.avg <- CD4.exp.sum/nrow(combo.reference) 
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
      Treg.exp.avg <- Treg.exp.sum/nrow(combo.reference) 
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
      MDSc.exp.avg <- MDSc.exp.sum/nrow(combo.reference) 
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
      Macro.exp.avg <- Macro.exp.sum/nrow(combo.reference) 
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
      Mono.exp.avg <- Mono.exp.sum/nrow(combo.reference) 
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
      Neut.exp.avg <- Neut.exp.sum/nrow(combo.reference) 
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
      DC.exp.avg <- DC.exp.sum/nrow(combo.reference) 
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
      Mast.exp.avg <- Mast.exp.sum/nrow(combo.reference) 
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


rownames(MarkerBased.dataframe) <- gsub(".", "-", rownames(MarkerBased.dataframe), fixed=TRUE)
rownames(ExpOnly.Results) <- gsub(".", "-", rownames(ExpOnly.Results), fixed=TRUE)


#test NA cell: Palprim_AACCGCGTCTCTGTCG-1_0135_TNBCTot #<- dash got turned into "." 
#test non-NA cell in both: Wu_CID44041_CCTACCACAATAGAGT
combo.reference <- AddMetaData(combo.reference, MarkerBased.dataframe, "MarkerBasedKara.NK")
combo.reference <- AddMetaData(combo.reference, ExpOnly.Results, "ExpBasedOnly.NK")

rm(MarkerBased.dataframe)
rm(ExpOnly.Results)

combo.reference <- AddMetaData(combo.reference, MarkerBased.dataframe, "MarkerBasedKara.Epi")
combo.reference <- AddMetaData(combo.reference, ExpOnly.Results, "ExpBasedOnly.Epi")

rm(MarkerBased.dataframe)
rm(ExpOnly.Results)

combo.reference <- AddMetaData(combo.reference, MarkerBased.dataframe, "MarkerBasedKara.Stroma")
combo.reference <- AddMetaData(combo.reference, ExpOnly.Results, "ExpBasedOnly.Stroma")

rm(MarkerBased.dataframe)
rm(ExpOnly.Results)

combo.reference <- AddMetaData(combo.reference, MarkerBased.dataframe, "MarkerBasedKara.B")
combo.reference <- AddMetaData(combo.reference, ExpOnly.Results, "ExpBasedOnly.B")

rm(MarkerBased.dataframe)
rm(ExpOnly.Results)

combo.reference <- AddMetaData(combo.reference, MarkerBased.dataframe, "MarkerBasedKara.Mye")
combo.reference <- AddMetaData(combo.reference, ExpOnly.Results, "ExpBasedOnly.Mye")

rm(MarkerBased.dataframe)
rm(ExpOnly.Results)

combo.reference <- AddMetaData(combo.reference, MarkerBased.dataframe, "MarkerBasedKara.T")
combo.reference <- AddMetaData(combo.reference, ExpOnly.Results, "ExpBasedOnly.T")

rm(MarkerBased.dataframe)
rm(ExpOnly.Results)

colnames(combo.reference@meta.data)

combo.reference$MarkerBasedKara <- apply(combo.reference@meta.data[,c(78, 80, 82, 84, 86, 88)], 1, function(x) x[!is.na(x)][1])
combo.reference@meta.data <- combo.reference@meta.data[,-c(78, 80, 82, 84, 86, 88)]

combo.reference$ExpBasedOnly <- apply(combo.reference@meta.data[,c(78:83)], 1, function(x) x[!is.na(x)][1])
combo.reference@meta.data <- combo.reference@meta.data[,-c(78:83)]

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/")
#saveRDS(combo.reference, "PrimObject_6422.rds")
combo.reference <- readRDS("PrimObject_62922.rds")

#Refine Expression and Marker Methods ============================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This will assign NA values to the larger cluster label, and 
# try to further categorize T and myeloid subsets by again looking
# at the larger cluster

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


combo.reference@meta.data[1:20,c(75, 77:80)]
#75 = UCell
#77 = markerbased
#78 = expression based only
table(combo.reference$MarkerBasedKara)

combo.reference$RefineMarkerBased <- "hi"
for (i in 1:nrow(combo.reference@meta.data)) {
  combo.reference@meta.data[i,"RefineMarkerBased"] <-ifelse(is.na(combo.reference@meta.data[i,"MarkerBasedKara"]),
                                           as.character(combo.reference@meta.data[i,"celltype_UCell"]), as.character(combo.reference@meta.data[i,"MarkerBasedKara"]))
}

for (i in 1:nrow(combo.reference@meta.data)) {
  
  if (combo.reference@meta.data[i,"RefineMarkerBased"] == "Unspecific T Cells")
  {
    if ((combo.reference@meta.data[i,"celltype_UCell"] == "CD8+ T Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "CD4+ T Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Regulatory T Cells")) 
    {
      combo.reference@meta.data[i,"RefineMarkerBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      
    }
    else {
      combo.reference@meta.data[i,"RefineMarkerBased"] <- "not Tcell by UScore"
    }  
  }
  
  else if (combo.reference@meta.data[i,"RefineMarkerBased"] == "T Cells")
  {
    if ((combo.reference@meta.data[i,"celltype_UCell"] == "CD8+ T Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "CD4+ T Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Regulatory T Cells")) 
    {
      combo.reference@meta.data[i,"RefineMarkerBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      
    }
    else {
      combo.reference@meta.data[i,"RefineMarkerBased"] <- "not Tcell by UScore"
    }
    
  }
  else if (combo.reference@meta.data[i,"RefineMarkerBased"] == "Unspecific Myeloid Cells")
  {
    if ((combo.reference@meta.data[i,"celltype_UCell"] == "Macrophages") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Mast Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Monocytes") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Dendritic Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Neutrophils") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "MDSCs")) 
    {
      combo.reference@meta.data[i,"RefineMarkerBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      
    }
    else {
      combo.reference@meta.data[i,"RefineMarkerBased"] <- "not Myeloid cell by UScore"
    }  
  }
  else if (combo.reference@meta.data[i,"RefineMarkerBased"] == "Myeloid Cells")
  {
    if ((combo.reference@meta.data[i,"celltype_UCell"] == "Macrophages") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Mast Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Monocytes") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Dendritic Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Neutrophils") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "MDSCs")) 
    {
      combo.reference@meta.data[i,"RefineMarkerBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      
    }
    else {
      combo.reference@meta.data[i,"RefineMarkerBased"] <- "not Myeloid cell by UScore"
    }  
  }
  else if (combo.reference@meta.data[i,"RefineMarkerBased"] == "Unknown")
  {
    combo.reference@meta.data[i,"RefineMarkerBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
  }
  else if (combo.reference@meta.data[i,"RefineMarkerBased"] == "Zero Expression for All Markers")
  {
    combo.reference@meta.data[i,"RefineMarkerBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
  }
  else if (combo.reference@meta.data[i,"RefineMarkerBased"] == "Non-Specific Immune Cells")
  {
    if ((combo.reference@meta.data[i,"celltype_UCell"] == "CD8+ T Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "CD4+ T Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Regulatory T Cells") | 
        (combo.reference@meta.data[i,"celltype_UCell"] == "NK Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "B Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Plasma Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "MDSCs") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Macrophages") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Monocytes") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Neutrophils") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Dendritic Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Mast Cells")) 
    {
      combo.reference@meta.data[i,"RefineMarkerBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      
    }
    else {
      combo.reference@meta.data[i,"RefineMarkerBased"] <- "not immune cell by UScore"
    }  }
  
  else {
    combo.reference@meta.data[i,"RefineMarkerBased"] <- as.character(combo.reference@meta.data[i,"RefineMarkerBased"])
    
  }
}

combo.reference@meta.data[78936,c(75, 77:79)]

table(combo.reference$RefineMarkerBased)



combo.reference$RefineExpBased <- "hi"
for (i in 1:nrow(combo.reference@meta.data)) {
  combo.reference@meta.data[i,"RefineExpBased"] <-ifelse(is.na(combo.reference@meta.data[i,"ExpBasedOnly"]),
                                           as.character(combo.reference@meta.data[i,"celltype_UCell"]), as.character(combo.reference@meta.data[i,"ExpBasedOnly"]))
}

for (i in 1:nrow(combo.reference@meta.data)) {
  
  if (combo.reference@meta.data[i,"RefineExpBased"] == "Unspecific T Cells")
  {
    if ((combo.reference@meta.data[i,"celltype_UCell"] == "CD8+ T Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "CD4+ T Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Regulatory T Cells")) 
    {
      combo.reference@meta.data[i,"RefineExpBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      
    }
    else {
      combo.reference@meta.data[i,"RefineExpBased"] <- "not Tcell by UScore"
    }  
  }
  else if (combo.reference@meta.data[i,"RefineExpBased"] == "Unspecific Myeloid Cells")
  {
    if ((combo.reference@meta.data[i,"celltype_UCell"] == "Macrophages") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Mast Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Monocytes") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Dendritic Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Neutrophils") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "MDSCs")) 
    {
      combo.reference@meta.data[i,"RefineExpBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      
    }
    else {
      combo.reference@meta.data[i,"RefineExpBased"] <- "not Myeloid cell by UScore"
    }  
  }
  else if (combo.reference@meta.data[i,"RefineExpBased"] == "Myeloid Cells")
  {
    if ((combo.reference@meta.data[i,"celltype_UCell"] == "Macrophages") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Mast Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Monocytes") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Dendritic Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Neutrophils") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "MDSCs")) 
    {
      combo.reference@meta.data[i,"RefineExpBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      
    }
    else {
      combo.reference@meta.data[i,"RefineExpBased"] <- "not Myeloid cell by UScore"
    }  
  }
  else if (combo.reference@meta.data[i,"RefineExpBased"] == "T Cells")
  {
    if ((combo.reference@meta.data[i,"celltype_UCell"] == "CD8+ T Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "CD4+ T Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Regulatory T Cells")) 
    {
      combo.reference@meta.data[i,"RefineExpBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      
    }
    else {
      combo.reference@meta.data[i,"RefineExpBased"] <- "not Tcell by UScore"
    }
    
  }
  
  else if (combo.reference@meta.data[i,"RefineExpBased"] == "Unknown")
  {
    combo.reference@meta.data[i,"RefineExpBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
  }
  else if (combo.reference@meta.data[i,"RefineExpBased"] == "Zero Expression for All Markers")
  {
    combo.reference@meta.data[i,"RefineExpBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
  }
  else if (combo.reference@meta.data[i,"RefineExpBased"] == "Non-Specific Immune Cells")
  {
    if ((combo.reference@meta.data[i,"celltype_UCell"] == "CD8+ T Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "CD4+ T Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Regulatory T Cells") | 
        (combo.reference@meta.data[i,"celltype_UCell"] == "NK Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "B Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Plasma Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "MDSCs") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Macrophages") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Monocytes") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Neutrophils") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Dendritic Cells") |
        (combo.reference@meta.data[i,"celltype_UCell"] == "Mast Cells")) 
    {
      combo.reference@meta.data[i,"RefineExpBased"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      
    }
    else {
      combo.reference@meta.data[i,"RefineExpBased"] <- "not immune cell by UScore"
    }  }
  
  else {
    combo.reference@meta.data[i,"RefineExpBased"] <- as.character(combo.reference@meta.data[i,"RefineExpBased"])
    
  }
}


table(combo.reference$RefineExpBased)



#agreement test ===============

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This will create the final labels by basically taking a vote between 
# the three methods

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(DescTools)
colnames(combo.reference@meta.data)

combo.reference$celltype_final <- "hi"
for (i in 1:nrow(combo.reference@meta.data)){
  test <- c(as.character(combo.reference@meta.data[i,"celltype_UCell"]), 
            as.character(combo.reference@meta.data[i,"RefineMarkerBased"]), 
            as.character(combo.reference@meta.data[i,"RefineExpBased"]))
  
  test.entry <- as.character(Mode(test))
  
  if (!is.na(test.entry) == T) {
    if ((as.character(combo.reference@meta.data[i,"orig.ident"]) == "Savas") && (
      (test.entry != "Regulatory T Cells") && (test.entry != "CD8+ T Cells") &&
      (test.entry != "CD4+ T Cells"))){
      combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
    }
    
    else if ((as.character(combo.reference@meta.data[i,"orig.ident"]) == "AziziT") && (
      (test.entry != "Regulatory T Cells") && (test.entry != "CD8+ T Cells") &&
      (test.entry != "CD4+ T Cells"))){
      combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
    }
    
    else if ((as.character(combo.reference@meta.data[i,"orig.ident"]) == "Aziziimmune") && (
      (test.entry != "Mast Cells") | (test.entry != "Dendritic Cells") | (test.entry != "Neutrophils") |
      (test.entry != "Monocytes") | (test.entry != "Macrophages") | (test.entry != "MDSCs") |
      (test.entry != "Regulatory T Cells") | (test.entry != "CD4+ T Cells") | (test.entry != "CD8+ T Cells") |
      (test.entry != "Plasma Cells") | (test.entry != "B Cells") | (test.entry != "NK Cells"))){
      #if UCell has it as some immune cell
      if ((as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Mast Cells") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Dendritic Cells") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Neutrophils") | 
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Monocytes") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Macrophages") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "MDSCs") | 
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Regulatory T Cells") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "CD4+ T Cells") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "CD8+ T Cells") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Plasma Cells") | 
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "B Cells") | 
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "NK Cells")) 
      {
        combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      }
      #if karaavaz has it as some immune cell
      else if ((as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Mast Cells") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Dendritic Cells") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Neutrophils") | 
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Monocytes") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Macrophages") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "MDSCs") | 
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Regulatory T Cells") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "CD4+ T Cells") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "CD8+ T Cells") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Plasma Cells") | 
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "B Cells") | 
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "NK Cells")) {
        combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"RefineMarkerBased"])
      }
      #if marker expression has it as some kind of immune cell
      else if ((as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Mast Cells") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Dendritic Cells") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Neutrophils") | 
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Monocytes") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Macrophages") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "MDSCs") | 
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Regulatory T Cells") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "CD4+ T Cells") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "CD8+ T Cells") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Plasma Cells") | 
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "B Cells") | 
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "NK Cells")) {
        combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"RefineExpBased"])
      }
      else {
        combo.reference@meta.data[i,"celltype_final"] <- "noone called it immune"

      }
    }
    
    else if (test.entry == "not Tcell by UScore") {
      combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
    }
    else if (test.entry == "not Myeloid cell by UScore") {
      combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      
    }
    else {
      combo.reference@meta.data[i,"celltype_final"] <- test.entry
    }
  }
  else {
    if ((as.character(combo.reference@meta.data[i,"orig.ident"]) == "Savas")){
      combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
    }
    
    else if ((as.character(combo.reference@meta.data[i,"orig.ident"]) == "AziziT")){
      combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
    }
    
    else if ((as.character(combo.reference@meta.data[i,"orig.ident"]) == "Aziziimmune")){
      #if UCell has it as some immune cell
      if ((as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Mast Cells") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Dendritic Cells") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Neutrophils") | 
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Monocytes") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Macrophages") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "MDSCs") | 
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Regulatory T Cells") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "CD4+ T Cells") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "CD8+ T Cells") |
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "Plasma Cells") | 
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "B Cells") | 
          (as.character(combo.reference@meta.data[i,"celltype_UCell"]) == "NK Cells")) 
      {
        combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      }
      #if karaavaz has it as some immune cell
      else if ((as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Mast Cells") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Dendritic Cells") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Neutrophils") | 
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Monocytes") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Macrophages") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "MDSCs") | 
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Regulatory T Cells") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "CD4+ T Cells") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "CD8+ T Cells") |
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "Plasma Cells") | 
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "B Cells") | 
               (as.character(combo.reference@meta.data[i,"RefineMarkerBased"]) == "NK Cells")) {
        combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"RefineMarkerBased"])
      }
      #if marker expression has it as some kind of immune cell
      else if ((as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Mast Cells") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Dendritic Cells") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Neutrophils") | 
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Monocytes") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Macrophages") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "MDSCs") | 
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Regulatory T Cells") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "CD4+ T Cells") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "CD8+ T Cells") |
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "Plasma Cells") | 
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "B Cells") | 
               (as.character(combo.reference@meta.data[i,"RefineExpBased"]) == "NK Cells")) {
        combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"RefineExpBased"])
      }
      else {
        combo.reference@meta.data[i,"celltype_final"] <- "noone called it immune"
      }
    }
    
    else {
      combo.reference@meta.data[i,"celltype_final"] <- as.character(combo.reference@meta.data[i,"celltype_UCell"])
      
    }
  }
}

head(combo.reference@meta.data)
table(combo.reference$celltype_final)


Idents(combo.reference) <- combo.reference$celltype_final

table(Idents(combo.reference))

#double check immune only datasets are in fact, immune only: 

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



combo.reference <- subset(combo.reference, idents = "noone called it immune", invert = T)

saveRDS(combo.reference, "FinalPrimObject_62922.rds")


