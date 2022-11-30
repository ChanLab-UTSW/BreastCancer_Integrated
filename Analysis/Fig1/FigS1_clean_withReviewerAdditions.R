


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This will generate all of the sub-figures for Figure S1,
# and is the code for QC and batch correction testing.


##note: S1A and S1B were made in prism 
##(datasetpie_72122.pzfx in "Project - single cell almanac -> Supplemental_Data -> Figure1)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



finalQC_Dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/PrimaryBreast_Processing/Individual_Dataset_Objects/finalQC"
Integration_Results_Dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/PrimaryBreast_Processing/Integration_Objects"


# Libraries -------------------------------------------------------------
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




# ================================================================== ======
#S1C age barplot ========================
# ================================================================== ======
# load in Prim object ------

setwd(PrimDir)
combo.reference <- readRDS("PrimObject_withreprog_noZallgenedem_71322.rds") ##preprint
setwd("/project/InternalMedicine/Chan_lab/shared/")
combo.reference <- readRDS("PrimObj_withGEmeta_with_correctassays_withfinalreprog_111622.rds") ##newnew

DefaultAssay(combo.reference) <- "RNA"
combo.reference <- NormalizeData(combo.reference, assay = "RNA")


# generate barplot -------------------------------------------------------


sobjlists <- FetchData(object = combo.reference, vars = c("samples", "Patient", "Age"))
table(sobjlists$Age)
sobjlists <- sobjlists %>% dplyr::group_by(Patient, Age)  %>%
  dplyr::summarise(Nb = n()) %>%
  dplyr::mutate(C = sum(Nb))  

sobjlists <- sobjlists[!grepl("NA", sobjlists$Age),]
sobjlists <- sobjlists[!grepl("NA", sobjlists$Age),]
table(sobjlists$Age)
class(sobjlists$Age)
sobjlists$Age[which(sobjlists$Age == "40-45")] <- "43"
sobjlists$Age[which(sobjlists$Age == "46-50")] <- "48"
sobjlists$Age[which(sobjlists$Age == "50-55")] <- "53"
sobjlists$Age[which(sobjlists$Age == "60-65")] <- "63"
sobjlists$Age[which(sobjlists$Age == "70-75")] <- "73"
sobjlists$Age[which(sobjlists$Age == "76-80")] <- "78"
sobjlists$Age[which(sobjlists$Age == "80-85")] <- "83"
sobjlists$Age[which(sobjlists$Age == "86-90")] <- "88"

sobjlists$Age <- as.numeric(sobjlists$Age)

p <- ggplot(sobjlists, aes(x = Age)) +
  geom_bar(fill = "#ef899f", color = "black") + 
  scale_x_binned() +
  #scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12), limits = c(0, 12)) +
  theme_classic() +
  xlab("Age Group") + ylab("# Patients") + 
  theme(axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 16), 
        axis.title.x = element_text(margin=margin(t = 10)), 
        axis.title.y = element_text(margin=margin(r = 10))) +
  RotatedAxis()
ggsave("patientnumbar_byage_71022.pdf", p, width = 4, height = 6.2)


# ================================================================== ======
#S1D-S1E PC plots ========================
# ================================================================== ======

# Merge Datasets (fixed) =============================================================

setwd(finalQC_Dir)
AziziPrim.fixed <- readRDS("AziziPrim_fixed_62722.rds")
AziziT.fixed <- readRDS("AziziT_fixed_62722.rds")
Karaayvaz.fixed <- readRDS("Karaayvaz_fixed_62722.rds")
Pal.fixed <- readRDS("Pal_fixed_62722.rds")
Qian.fixed <- readRDS("Qian_fixed_62722.rds")
Savas.fixed <- readRDS("Savas_fixed_62722.rds")
OldWu.fixed <- readRDS("OldWu_fixed_62722.rds")
NewWu.fixed <- readRDS("NewWu_fixed_62722.rds")
Xu.fixed <- readRDS("Xu_fixed_62722.rds")


combo <- merge(Savas.fixed, y = c(AziziT.fixed, AziziPrim.fixed, OldWu.fixed, 
                                  NewWu.fixed, Qian.fixed, Xu.fixed, Karaayvaz.fixed, Pal.fixed), 
               add.cell.ids = c("Savas", "AziziT", "Aziziimmune", "Wu", 
                                "Wu2021prim", "Qian", "Xu", "Karaayvaz",  "Palprim"))
head(combo@meta.data)

table(combo$Treatment.Status)
table(combo$Capture.Method)

combo$Treatment.Status[which(combo$Treatment.Status == "Niave")] <- "Naive"
combo$Capture.Method[which(combo$Capture.Method == "10X Genomics Chromium v2")] <- "10X Genomics Single Cell 3' v2"
combo$Capture.Method[which(combo$Capture.Method == "10X Genomics Chromium Single Cell 5' Library and Bead Kit")] <- "10X Genomics Chromium v2 5'"
DefaultAssay(combo) <- "RNA"



# Integration and PC plot generation ==============================================================

table(combo$Doublet.Call)
combo <- subset(combo, subset = Doublet.Call == "Doublet", invert = T)

combo.list <- SplitObject(combo, split.by = "Capture.Method")

#https://github.com/satijalab/sctransform/issues/94 <- sctransform for read counts (smart-seq2)
for (i in 1:length(combo.list)) {
  combo.list[[i]] <- SCTransform(combo.list[[i]], verbose = T, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                                                   "percent.platelet", "percent.heatshock"))
}


features <- SelectIntegrationFeatures(object.list = combo.list, nfeatures = 3000)
combo.list <- PrepSCTIntegration(object.list = combo.list, anchor.features = features)


reference.1 <-  which(names(combo.list) == c("10X Genomics Chromium"))
reference.2 <-  which(names(combo.list) == c("10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"))
reference.3 <-  which(names(combo.list) == c("10X Genomics Chromium v2 5'"))
reference.4 <-  which(names(combo.list) == c("10X Genomics Single Cell 3' v2"))

reference.list <- c(reference.1, reference.2, reference.3, reference.4)
sort(table(combo.reference$Capture.Method))
combo.anchors <- FindIntegrationAnchors(object.list = combo.list, normalization.method = "SCT",
                                        anchor.features = features, reference = reference.list)
#saveRDS(combo.anchors, file = "primOnlyYESSSREFSCTAnchors_BiggestRef_doubletRemoved_62722.rds")


combo.reference <- IntegrateData(anchorset = combo.anchors, normalization.method = "SCT")
setwd(Integration_Results_Dir)
#saveRDS(combo.reference, file = "primOnlyYESSSREFsctINTEGRATED_BiggestRef_doubletRemoved_62722.rds")

combo.reference <- readRDS(file = "primOnlyYESSSREFsctINTEGRATED_BiggestRef_doubletRemoved_62722.rds")
DefaultAssay(combo.reference) <- "integrated"

combo.reference <- RunPCA(combo.reference, npcs = 100, verbose = FALSE)

combo.reference$orig.ident[which(combo.reference$orig.ident == "AziziT")] <- "Azizi"
combo.reference$orig.ident[which(combo.reference$orig.ident == "Aziziimmune")] <- "Azizi"

sort(table(combo.reference$Capture.Method))
sort(table(combo.reference$orig.ident))
#pdf("BatchTestSCTRef_71922_long.pdf", width = 16, height = 9)
pdf("test.pdf", width = 16, height = 9)
p1 <- DimPlot(combo.reference, reduction = "pca", raster = F, group.by = "Capture.Method",
              order = c("Smart-Seq2", "inDrop v2", "Singleron GEXSCOPE Single Cell RNAseq Library Kit",
                        "10X Genomics Single Cell 3' v2", "10X Genomics Chromium v2 5'",
                        "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library",
                        "10X Genomics Chromium"))
p1
p2 <- DimPlot(combo.reference, reduction = "pca", raster = F, group.by = "orig.ident",
              order =c("Karaayvaz", "Savas", "Wu", "Xu",
                       "Azizi", "Qian", "Wu2021prim", "Pal_Prim"))
p2
dev.off()

# ================================================================== ======
#S1F BC separated UMAP ========================
# ================================================================== ======
# load in objects (REVIEWER ADDITION) ++++++++++++++++++++++++++++++++++++++++ ====

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/ReviewerAdditions/Pre_post_integration_compare")
combo.unint <- readRDS("PrimObject_nointegration_102822.rds") ##unintegrated

setwd("/project/InternalMedicine/Chan_lab/shared/")
combo.reference <- readRDS("PrimObj_withGEmeta_with_correctassays_withfinalreprog_111622.rds") ##integrated


# UMAPs (BCsub, orig.ident, samples) -------

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/ReviewerAdditions/Pre_post_integration_compare")


pdf("preInt_BCsub_primUMAP_112522.pdf", width = 16.22, height = 14.22)
DimPlot(combo.unint, reduction = "umap", group.by = "BC.Subtype", raster = FALSE,
        cols = c("#D8197D", "#56A008", "#063970")) + ggtitle(label = " ") +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()

pdf("postInt_BCsub_primUMAP_112522.pdf", width = 16.22, height = 14.22)
DimPlot(combo.reference, reduction = "umap", group.by = "BC.Subtype", raster = FALSE,
        cols = c("#D8197D", "#56A008", "#063970")) + ggtitle(label = " ") +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()



combo.unint$orig.ident[which(combo.unint$orig.ident == "AziziT")] <- "Azizi"
combo.unint$orig.ident[which(combo.unint$orig.ident == "Aziziimmune")] <- "Azizi"


pdf("preInt_OrigStudy_primUMAP_112522.pdf", width = 16.22, height = 14.22)
DimPlot(combo.unint, reduction = "umap", group.by = "orig.ident", raster = FALSE,
        cols = c("#1432EF", "#E74124", "#54BE38",
                 "#A941DB", "#EC9336", "#000000",
                 "#997938", "#071887")) + ggtitle(label = " ") +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()


combo.reference$orig.ident[which(combo.reference$orig.ident == "AziziT")] <- "Azizi"
combo.reference$orig.ident[which(combo.reference$orig.ident == "Aziziimmune")] <- "Azizi"

pdf("postInt_OrigStudy_primUMAP_112522.pdf", width = 16.22, height = 14.22)
DimPlot(combo.reference, reduction = "umap", group.by = "orig.ident", raster = FALSE,
        cols = c("#1432EF", "#E74124", "#54BE38",
                 "#A941DB", "#EC9336", "#000000",
                 "#997938", "#071887")) + ggtitle(label = " ") +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()


pdf("preInt_samples_primUMAP_112522.pdf", width = 16.22, height = 14.22)
DimPlot(combo.unint, reduction = "umap", group.by = "samples", raster = FALSE,
        cols = c(turbo(119))) + ggtitle(label = " ") +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()

pdf("postInt_samples_primUMAP_112522.pdf", width = 16.22, height = 14.22)
DimPlot(combo.reference, reduction = "umap", group.by = "samples", raster = FALSE,
        cols = c(turbo(119))) + ggtitle(label = " ") +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()

#PC plots (BCsub, orig.ident, samples) ====

pdf("preInt_PCAplots_bcsubTechorigIdent_112522.pdf", width = 16, height = 9)
p0 <- DimPlot(combo.unint, reduction = "pca", raster = F, group.by = "BC.Subtype", 
              order = c("HER2+", "TNBC", "HR+"),
              cols = c("#56A008", "#063970", "#D8197D"))
p0
p1 <- DimPlot(combo.unint, reduction = "pca", raster = F, group.by = "Capture.Method",
              order = c("Smart-Seq2", "inDrop v2", "Singleron GEXSCOPE Single Cell RNAseq Library Kit",
                        "10X Genomics Single Cell 3' v2", "10X Genomics Chromium v2 5'",
                        "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library",
                        "10X Genomics Chromium"), cols = c(turbo(7)))
p1
p2 <- DimPlot(combo.unint, reduction = "pca", raster = F, group.by = "orig.ident",
              order =c("Karaayvaz", "Savas", "Wu", "Xu",
                       "Azizi", "Qian", "Wu2021prim", "Pal_Prim"),
              cols = c("#E74124", "#EC9336", "#000000", "#071887", "#1432EF",
                       "#A941DB", "#997938", "#54BE38"))
p2
dev.off()


pdf("postInt_PCAplots_bcsubTechorigIdent_112522.pdf", width = 16, height = 9)
p0 <- DimPlot(combo.reference, reduction = "pca", raster = F, group.by = "BC.Subtype", 
              order = c("HER2+", "TNBC", "HR+"),
              cols = c("#56A008", "#063970", "#D8197D"))
p0
p1 <- DimPlot(combo.reference, reduction = "pca", raster = F, group.by = "Capture.Method",
              order = c("Smart-Seq2", "inDrop v2", "Singleron GEXSCOPE Single Cell RNAseq Library Kit",
                        "10X Genomics Single Cell 3' v2", "10X Genomics Chromium v2 5'",
                        "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library",
                        "10X Genomics Chromium"), cols = c(turbo(7)))
p1
p2 <- DimPlot(combo.reference, reduction = "pca", raster = F, group.by = "orig.ident",
              order =c("Karaayvaz", "Savas", "Wu", "Xu",
                       "Azizi", "Qian", "Wu2021prim", "Pal_Prim"),
              cols = c("#E74124", "#EC9336", "#000000", "#071887", "#1432EF",
                              "#A941DB", "#997938", "#54BE38"))
p2
dev.off()


