

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only")


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




#SCT =======================================================

#just filtered objects =======================================================


setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Individual_Datasets/DoubletFinder/DoubletFiles_62622/")
AziziPrim <- readRDS("AziziPrim_withDoublet_62622.rds")
AziziT <- readRDS("AziziT_withDoublet_62622.rds")
Karaayvaz <- readRDS("Karaayvaz_withDoublet_62622.rds")
Pal <- readRDS("Pal_withDoublet_62622.rds")
Qian <- readRDS("Qian_withDouble_62622t.rds")
Savas <- readRDS("Savas_withDoublet_62622.rds")
OldWu <- readRDS("OldWu_withDoublet_62622.rds")
Wu2021 <- readRDS("Wu2021_withDoublet_62622.rds")
Xu <- readRDS("Xu_withDoublet_62622.rds")
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/df.converted")


#Check For Duplicates (- and .) : AziziPrim ============================================

AziziPrim.df <- GetAssayData(AziziPrim, slot = "counts", assay = "RNA")
AziziPrim.df[1:5,1:5]
dim(AziziPrim.df)

old.gene.names <- rownames(AziziPrim.df)
new.gene.names <- gsub("[.]", "-", old.gene.names)
rownames(AziziPrim.df) <- new.gene.names

new.gene.names.dupl <- new.gene.names[duplicated(new.gene.names)]
new.gene.names.dupl #is zero


which(rownames(AziziPrim.df) == "CMP21-97G8.1") 
which(rownames(AziziPrim.df) == "CMP21-97G8-1") 
which(rownames(AziziPrim.df) == "RP11-34P13.7") 
which(rownames(AziziPrim.df) == "HLA-A") 
which(rownames(AziziPrim.df) == "HLA.A") 

AziziPrim.meta <- AziziPrim@meta.data
AziziPrim.fixed <- CreateSeuratObject(AziziPrim.df, meta.data = AziziPrim.meta)
head(AziziPrim.fixed@meta.data)
identical(AziziPrim.fixed@meta.data, AziziPrim@meta.data)

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/df.converted")
saveRDS(AziziPrim.fixed, "AziziPrim_fixed_62722.rds")


#Check For Duplicates (- and .) : AziziT ============================================

AziziT.df <- GetAssayData(AziziT, slot = "counts", assay = "RNA")
AziziT.df[1:5,1:5]
dim(AziziT.df)

old.gene.names <- rownames(AziziT.df)
new.gene.names <- gsub("[.]", "-", old.gene.names)
rownames(AziziT.df) <- new.gene.names
new.gene.names.dupl <- new.gene.names[duplicated(new.gene.names)]
new.gene.names.dupl #isn't zero

which(rownames(AziziT.df) == "TRAF3IP2-AS1") 
which(rownames(AziziT.df) == "KRTAP10-11") 

AziziT.df.duplicated.1 <- as.data.frame(t(AziziT.df[c(11303,11305),]))
AziziT.df.duplicated.1$TRAF3IP2.AS1.new <- apply(AziziT.df.duplicated.1[,c(1:2)], 1, function(x) x[!is.na(x)][1])
AziziT.df.duplicated.1[1:2,1:3]
AziziT.df.duplicated.1 <- AziziT.df.duplicated.1[,-c(1:2), drop = F]
colnames(AziziT.df.duplicated.1) <- "TRAF3IP2-AS1"
AziziT.df.duplicated.1 <- t(AziziT.df.duplicated.1)
AziziT.df.duplicated.1 <- Matrix(AziziT.df.duplicated.1, sparse = TRUE)
AziziT.df.duplicated.1[1,1:5]
rownames(AziziT.df.duplicated.1)

AziziT.df.duplicated.2 <- as.data.frame(t(AziziT.df[c(33505,33507),]))
AziziT.df.duplicated.2$KRTAP10.11.new <- apply(AziziT.df.duplicated.2[,c(1:2)], 1, function(x) x[!is.na(x)][1])
AziziT.df.duplicated.2[1:2,1:3]
AziziT.df.duplicated.2 <- AziziT.df.duplicated.2[,-c(1:2), drop = F]
colnames(AziziT.df.duplicated.2) <- "KRTAP10-11"
AziziT.df.duplicated.2 <- t(AziziT.df.duplicated.2)
AziziT.df.duplicated.2 <- Matrix(AziziT.df.duplicated.2, sparse = TRUE)
AziziT.df.duplicated.2[1,1:5]
rownames(AziziT.df.duplicated.2)

AziziT.df.nodupl <- AziziT.df[-c(11303,11305, 33505, 33507),]
AziziT.df.nodupl <- rbind(AziziT.df.nodupl, AziziT.df.duplicated.1, AziziT.df.duplicated.2)
final.check <- rownames(AziziT.df.nodupl)

final.check <- final.check[duplicated(final.check)]
final.check

which(rownames(AziziT.df.nodupl) == "TRAF3IP2-AS1") 
which(rownames(AziziT.df.nodupl) == "KRTAP10-11") 
dim(AziziT.df.nodupl)


which(rownames(AziziT.df.nodupl) == "CMP21-97G8.1") 
which(rownames(AziziT.df.nodupl) == "CMP21-97G8-1") 
which(rownames(AziziT.df.nodupl) == "RP11-34P13.7") 
which(rownames(AziziT.df.nodupl) == "HLA-A") 
which(rownames(AziziT.df.nodupl) == "HLA.A") 


AziziT.meta <- AziziT@meta.data

AziziT.fixed <- CreateSeuratObject(AziziT.df.nodupl, meta.data = AziziT.meta)
head(AziziT.fixed@meta.data)
identical(AziziT.fixed@meta.data, AziziT@meta.data)
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/df.converted")
saveRDS(AziziT.fixed, "AziziT_fixed_62722.rds")

#Check For Duplicates (- and .) : Karaayvaz ============================================

Karaayvaz.df <- GetAssayData(Karaayvaz, slot = "counts", assay = "RNA")
Karaayvaz.df[1:5,1:5]
dim(Karaayvaz.df)

old.gene.names <- rownames(Karaayvaz.df)
new.gene.names <- gsub("[.]", "-", old.gene.names)
rownames(Karaayvaz.df) <- new.gene.names
new.gene.names.dupl <- new.gene.names[duplicated(new.gene.names)]
new.gene.names.dupl #is zero

which(rownames(Karaayvaz.df) == "CMP21-97G8.1") 
which(rownames(Karaayvaz.df) == "CMP21-97G8-1") 
which(rownames(Karaayvaz.df) == "RP11-34P13.7") 
which(rownames(Karaayvaz.df) == "HLA-A") 
which(rownames(Karaayvaz.df) == "HLA.A") 


Karaayvaz.meta <- Karaayvaz@meta.data

Karaayvaz.fixed <- CreateSeuratObject(Karaayvaz.df, meta.data = Karaayvaz.meta)
head(Karaayvaz.fixed@meta.data)
identical(Karaayvaz.fixed@meta.data, Karaayvaz@meta.data)
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/df.converted")
saveRDS(Karaayvaz.fixed, "Karaayvaz_fixed_62722.rds")


#Check For Duplicates (- and .) : Pal ============================================

Pal.df <- GetAssayData(Pal, slot = "counts", assay = "RNA")
Pal.df[1:5,1:5]
dim(Pal.df)

old.gene.names <- rownames(Pal.df)
new.gene.names <- gsub("[.]", "-", old.gene.names)
rownames(Pal.df) <- new.gene.names
new.gene.names.dupl <- new.gene.names[duplicated(new.gene.names)]
new.gene.names.dupl #isn't zero

which(rownames(Pal.df) == "KRTAP10-11") 

Pal.df.duplicated <- as.data.frame(t(Pal.df[c(33396,33398),]))
Pal.df.duplicated$KRTAP10.11.new <- apply(Pal.df.duplicated[,c(1:2)], 1, function(x) x[!is.na(x)][1])
Pal.df.duplicated[1:2,1:3]
Pal.df.duplicated <- Pal.df.duplicated[,-c(1:2), drop = F]
colnames(Pal.df.duplicated) <- "KRTAP10-11"
Pal.df.duplicated <- t(Pal.df.duplicated)
Pal.df.duplicated <- Matrix(Pal.df.duplicated, sparse = TRUE)
Pal.df.duplicated[1,1:5]
rownames(Pal.df.duplicated)

Pal.df.nodupl <- Pal.df[-c(33396, 33398),]
Pal.df.nodupl <- rbind(Pal.df.nodupl, Pal.df.duplicated)
final.check <- rownames(Pal.df.nodupl)

final.check <- final.check[duplicated(final.check)]
final.check

which(rownames(Pal.df.nodupl) == "KRTAP10-11") 
which(rownames(Pal.df.nodupl) == "CMP21-97G8.1") 
which(rownames(Pal.df.nodupl) == "CMP21-97G8-1") 
which(rownames(Pal.df.nodupl) == "RP11-34P13.7") 
which(rownames(Pal.df.nodupl) == "HLA-A") 
which(rownames(Pal.df.nodupl) == "HLA.A") 
dim(Pal.df.nodupl)


Pal.meta <- Pal@meta.data

Pal.fixed <- CreateSeuratObject(Pal.df.nodupl, meta.data = Pal.meta)
head(Pal.fixed@meta.data)
identical(Pal.fixed@meta.data, Pal@meta.data)
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/df.converted")
saveRDS(Pal.fixed, "Pal_fixed_62722.rds")

#Check For Duplicates (- and .) : Qian ============================================

Qian.df <- GetAssayData(Qian, slot = "counts", assay = "RNA")
Qian.df[1:5,1:5]
dim(Qian.df)

old.gene.names <- rownames(Qian.df)
new.gene.names <- gsub("[.]", "-", old.gene.names)
rownames(Qian.df) <- new.gene.names
new.gene.names.dupl <- new.gene.names[duplicated(new.gene.names)]
new.gene.names.dupl #is zero

which(rownames(Qian.df) == "KRTAP10-11") 
which(rownames(Qian.df) == "CMP21-97G8.1") 
which(rownames(Qian.df) == "CMP21-97G8-1") 
which(rownames(Qian.df) == "RP11-34P13.7") 
which(rownames(Qian.df) == "HLA-A") 
which(rownames(Qian.df) == "HLA.A") 

Qian.meta <- Qian@meta.data

Qian.fixed <- CreateSeuratObject(Qian.df, meta.data = Qian.meta)
head(Qian.fixed@meta.data)
identical(Qian.fixed@meta.data, Qian@meta.data)
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/df.converted")
saveRDS(Qian.fixed, "Qian_fixed_62722.rds")


#Check For Duplicates (- and .) : Savas ============================================

Savas.df <- GetAssayData(Savas, slot = "counts", assay = "RNA")
Savas.df[1:5,1:5]
dim(Savas.df)

old.gene.names <- rownames(Savas.df)
new.gene.names <- gsub("[.]", "-", old.gene.names)
rownames(Savas.df) <- new.gene.names
new.gene.names.dupl <- new.gene.names[duplicated(new.gene.names)]
new.gene.names.dupl #isn't zero

which(rownames(Savas.df) == "TRAF3IP2-AS1") 
which(rownames(Savas.df) == "KRTAP10-11") 

Savas.df.duplicated.1 <- as.data.frame(t(Savas.df[c(11303,11305),]))
Savas.df.duplicated.1$TRAF3IP2.AS1.new <- apply(Savas.df.duplicated.1[,c(1:2)], 1, function(x) x[!is.na(x)][1])
Savas.df.duplicated.1[1:2,1:3]
Savas.df.duplicated.1 <- Savas.df.duplicated.1[,-c(1:2), drop = F]
colnames(Savas.df.duplicated.1) <- "TRAF3IP2-AS1"
Savas.df.duplicated.1 <- t(Savas.df.duplicated.1)
Savas.df.duplicated.1 <- Matrix(Savas.df.duplicated.1, sparse = TRUE)
Savas.df.duplicated.1[1,1:5]
rownames(Savas.df.duplicated.1)

Savas.df.duplicated.2 <- as.data.frame(t(Savas.df[c(33505,33507),]))
Savas.df.duplicated.2$KRTAP10.11.new <- apply(Savas.df.duplicated.2[,c(1:2)], 1, function(x) x[!is.na(x)][1])
Savas.df.duplicated.2[1:2,1:3]
Savas.df.duplicated.2 <- Savas.df.duplicated.2[,-c(1:2), drop = F]
colnames(Savas.df.duplicated.2) <- "KRTAP10-11"
Savas.df.duplicated.2 <- t(Savas.df.duplicated.2)
Savas.df.duplicated.2 <- Matrix(Savas.df.duplicated.2, sparse = TRUE)
Savas.df.duplicated.2[1,1:5]
rownames(Savas.df.duplicated.2)

Savas.df.nodupl <- Savas.df[-c(11303,11305, 33505, 33507),]
Savas.df.nodupl <- rbind(Savas.df.nodupl, Savas.df.duplicated.1, Savas.df.duplicated.2)
final.check <- rownames(Savas.df.nodupl)

final.check <- final.check[duplicated(final.check)]
final.check

which(rownames(Savas.df.nodupl) == "TRAF3IP2-AS1") 
which(rownames(Savas.df.nodupl) == "KRTAP10-11") 

which(rownames(Savas.df.nodupl) == "CMP21-97G8.1") 
which(rownames(Savas.df.nodupl) == "CMP21-97G8-1") 
which(rownames(Savas.df.nodupl) == "RP11-34P13.7") 
which(rownames(Savas.df.nodupl) == "HLA-A") 
which(rownames(Savas.df.nodupl) == "HLA.A") 
dim(Savas.df.nodupl)



Savas.meta <- Savas@meta.data

Savas.fixed <- CreateSeuratObject(Savas.df.nodupl, meta.data = Savas.meta)
head(Savas.fixed@meta.data)
identical(Savas.fixed@meta.data, Savas@meta.data)
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/df.converted")
saveRDS(Savas.fixed, "Savas_fixed_62722.rds")


#Check For Duplicates (- and .) : OldWu ============================================

OldWu.df <- GetAssayData(OldWu, slot = "counts", assay = "RNA")
OldWu.df[1:5,1:5]
dim(OldWu.df)

old.gene.names <- rownames(OldWu.df)
new.gene.names <- gsub("[.]", "-", old.gene.names)
rownames(OldWu.df) <- new.gene.names
new.gene.names.dupl <- new.gene.names[duplicated(new.gene.names)]
new.gene.names.dupl #is zero


which(rownames(OldWu.df) == "CMP21-97G8.1") 
which(rownames(OldWu.df) == "CMP21-97G8-1") 
which(rownames(OldWu.df) == "RP11-34P13.7") 
which(rownames(OldWu.df) == "HLA-A") 
which(rownames(OldWu.df) == "HLA.A") 


OldWu.meta <- OldWu@meta.data

OldWu.fixed <- CreateSeuratObject(OldWu.df, meta.data = OldWu.meta)
head(OldWu.fixed@meta.data)
identical(OldWu.fixed@meta.data, OldWu@meta.data)
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/df.converted")
saveRDS(OldWu.fixed, "OldWu_fixed_62722.rds")


#Check For Duplicates (- and .) : Wu2021 ============================================

Wu2021.df <- GetAssayData(Wu2021, slot = "counts", assay = "RNA")
Wu2021.df[1:5,1:5]
dim(Wu2021.df)

old.gene.names <- rownames(Wu2021.df)
new.gene.names <- gsub("[.]", "-", old.gene.names)
rownames(Wu2021.df) <- new.gene.names
new.gene.names.dupl <- new.gene.names[duplicated(new.gene.names)]
new.gene.names.dupl #is zero


which(rownames(Wu2021.df) == "CMP21-97G8.1") 
which(rownames(Wu2021.df) == "CMP21-97G8-1") 
which(rownames(Wu2021.df) == "RP11-34P13.7") 
which(rownames(Wu2021.df) == "HLA-A") 
which(rownames(Wu2021.df) == "HLA.A") 


Wu2021.meta <- Wu2021@meta.data

Wu2021.fixed <- CreateSeuratObject(Wu2021.df, meta.data = Wu2021.meta)
head(Wu2021.fixed@meta.data)
identical(Wu2021.fixed@meta.data, Wu2021@meta.data)

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/df.converted")
saveRDS(Wu2021.fixed, "NewWu_fixed_62722.rds")


#Check For Duplicates (- and .) : Xu ============================================

Xu.df <- GetAssayData(Xu, slot = "counts", assay = "RNA")
Xu.df[1:5,1:5]
dim(Xu.df)

old.gene.names <- rownames(Xu.df)
new.gene.names <- gsub("[.]", "-", old.gene.names)
rownames(Xu.df) <- new.gene.names
new.gene.names.dupl <- new.gene.names[duplicated(new.gene.names)]
new.gene.names.dupl #is zero


which(rownames(Xu.df) == "CMP21-97G8.1") 
which(rownames(Xu.df) == "CMP21-97G8-1") 
which(rownames(Xu.df) == "RP11-34P13.7") 
which(rownames(Xu.df) == "HLA-A") 
which(rownames(Xu.df) == "HLA.A") 

Xu.meta <- Xu@meta.data

Xu.fixed <- CreateSeuratObject(Xu.df, meta.data = Xu.meta)
head(Xu.fixed@meta.data)
identical(Xu.fixed@meta.data, Xu@meta.data)

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/df.converted")
saveRDS(Xu.fixed, "Xu_fixed_62722.rds")


#Merge Datasets (fixed) =============================================================

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/df.converted")
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



# ==============================================================

table(combo$Doublet.Call)
combo <- subset(combo, subset = Doublet.Call == "Doublet", invert = T)

combo.list <- SplitObject(combo, split.by = "Capture.Method")

#combo.list <- lapply(X = combo.list, FUN = SCTransform)
#https://github.com/satijalab/sctransform/issues/94 <- sctransform for read counts (smart-seq2)
for (i in 1:length(combo.list)) {
  combo.list[[i]] <- SCTransform(combo.list[[i]], verbose = T, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                                                   "percent.platelet", "percent.heatshock"))
}
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/")

saveRDS(combo.list, file = "SCTobjectlist_nodoublets_62722.rds")

combo.list <- readRDS("SCTobjectlist_nodoublets_62722.rds")
features <- SelectIntegrationFeatures(object.list = combo.list, nfeatures = 3000)
combo.list <- PrepSCTIntegration(object.list = combo.list, anchor.features = features)


reference.1 <-  which(names(combo.list) == c("10X Genomics Chromium"))
reference.2 <-  which(names(combo.list) == c("10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"))
reference.3 <-  which(names(combo.list) == c("10X Genomics Chromium v2 5'"))
reference.4 <-  which(names(combo.list) == c("10X Genomics Single Cell 3' v2"))

reference.list <- c(reference.1, reference.2, reference.3, reference.4)
# reference.list <-  which(names(combo.list) == c("10X Genomics Chromium",
#                                                 "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library", 
#                                                 "10X Genomics Chromium v2 5'", 
#                                                 "10X Genomics Single Cell 3' v2"))
sort(table(combo.reference$Capture.Method))
combo.anchors <- FindIntegrationAnchors(object.list = combo.list, normalization.method = "SCT",
                                        anchor.features = features, reference = reference.list)
saveRDS(combo.anchors, file = "primOnlyYESSSREFSCTAnchors_BiggestRef_doubletRemoved_62722.rds")

combo.reference <- IntegrateData(anchorset = combo.anchors, normalization.method = "SCT")
saveRDS(combo.reference, file = "primOnlyYESSSREFsctINTEGRATED_BiggestRef_doubletRemoved_62722.rds")
#saveRDS(combo.reference, file = "primaryREFERENCEintegratedPalQianXUNoChung_subsetsnoBRCA2REPLACEDGENES.rds")

combo.reference <- readRDS(file = "primOnlyYESSSREFsctINTEGRATED_BiggestRef_doubletRemoved_62722.rds")
DefaultAssay(combo.reference) <- "integrated"

#combo.reference <- ScaleData(combo.reference, assay = "RNA", verbose = TRUE)
combo.reference <- RunPCA(combo.reference, npcs = 100, verbose = FALSE)

pdf("BatchTestSCTRef_62822.pdf", width = 14.22, height = 13.01)
DimPlot(combo.reference, reduction = "pca", raster = F, group.by = "Capture.Method")
DimPlot(combo.reference, reduction = "pca", raster = F, group.by = "orig.ident")

DimPlot(combo.reference, reduction = "pca", raster = F, group.by = "Capture.Method",
        cols = c("Transparent", "Transparent", "Transparent", "blue", "yellow", "green", "purple"))
DimPlot(combo.reference, reduction = "pca", raster = F, group.by = "orig.ident",
        cols = c("red", "blue", "green", "Transparent", "yellow", "pink", "purple", "brown", "grey"))
dev.off()

pdf("test.pdf", width = 14.22, height = 13.01)
ElbowPlot(combo.reference, ndims = 100)
DimHeatmap(combo.reference, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(combo.reference, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(combo.reference, dims = 30:45, cells = 500, balanced = TRUE)
DimHeatmap(combo.reference, dims = 45:55, cells = 500, balanced = TRUE)
DimHeatmap(combo.reference, dims = 55:65, cells = 500, balanced = TRUE)
DimHeatmap(combo.reference, dims = 65:75, cells = 500, balanced = TRUE)
dev.off()

DefaultAssay(combo.reference) <- "integrated"

combo.reference <- FindNeighbors(combo.reference, reduction = "pca", dims = 1:60, nn.method = "rann")
combo.reference <- FindClusters(combo.reference, resolution = 12)
combo.reference <- RunUMAP(combo.reference, reduction = "pca", dims = 1:60, verbose = TRUE, seed.use = 123)
#saveRDS(combo.reference, "SCTRefwithUMAP_pc70res5_6222.rds")

combo.reference$BC.Subtype[which(combo.reference$BC.Subtype == "ER+")] <- "HR+"

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects")
combo.old <- readRDS("PrimObject_SCTnormRNA_withreprog_61422.rds")
combo.reference <- AddMetaData(combo.reference, metadata = combo.old@meta.data[,"celltype_final", drop = F], col.name = "old.labels")

pdf("test.pdf", width = 24.22, height = 23.01)
p <- DimPlot(combo.reference, reduction = "umap", label = F, repel = T, raster = FALSE) + ggtitle("")
DimPlot(combo.reference, reduction = "umap", label = TRUE, repel = T, raster = FALSE, group.by = "old.labels") + ggtitle("")
dev.off()
pdf("PrimUMAP_withTMyeloidSubs_withlegend.pdf", width = 14.22, height = 13.01)
p <- DimPlot(combo.reference, reduction = "umap", label = F, repel = T, raster = FALSE, 
             order = c("NK Cells", "Epithelial Cells", "Myoepithelial Cells",
                       "MDSCs","Dendritic Cells","Regulatory T Cells",
                       "Perivascular-like (PVL) Cells","B Cells", "CD4+ T Cells", "Neutrophils",
                       "Monocytes",  "Mast Cells","CD8+ T Cells", "Fibroblasts",
                       "Endothelial Cells", "Plasma Cells", "Macrophages")) + ggtitle(label = " ")
#https://stackoverflow.com/questions/36474643/graphical-parameters-in-ggplot2-how-to-change-axis-tick-thickness
pdf("PRIMUMAPfig1a_62922.pdf", width = 16, height = 9)
pdf("FINALFINALPRIMUMAPfig1a_62922.pdf", width = 16, height = 9)
pdf("test.pdf", width = 16, height = 9)
p + theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()
DimPlot(combo.reference, reduction = "umap", group.by = "BC.Subtype", order = c("HER2+", "HR+", "TNBC"), raster = FALSE) + ggtitle(label = "Integrated Primary by Breast Cancer Subtype (PC = 40)")
DimPlot(combo.reference, reduction = "umap", label = TRUE, group.by = "celltype_minor", raster = FALSE)+ ggtitle(label = "Cell Type (PC = 40)")
DimPlot(combo.reference, reduction = "umap", label = TRUE, group.by = "celltype_main", raster = FALSE)+ ggtitle(label = "Cell Type (PC = 40)")


setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Individual_Datasets/DoubletFinder")
combo.reference <- readRDS("comboSCTRef_withsigs_6222.rds")


combo.reference <- AddMetaData(combo.reference, combo.reference.old@meta.data[,"celltype_final", drop = F], "old.labels")
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects")
combo.reference <- readRDS("SCTRefwithUMAP_pc70res5_6222.rds")

#batch effect testing ===========================================

testAzizi <- subset(combo.reference, subset = orig.ident == "Aziziimmune")
DimPlot(testAzizi, reduction = "umap", label = TRUE, raster = FALSE)+ ggtitle(label = "Azizi_immune batch test (PC = 40)")

testAziziT <- subset(combo.reference, subset = orig.ident == "AziziT")
DimPlot(testAziziT, reduction = "umap", label = TRUE, raster = FALSE)+ ggtitle(label = "AziziT batch test (PC = 40)")

testKara <- subset(combo.reference, subset = orig.ident == "Karaayvaz")
DimPlot(testKara, reduction = "umap", label = TRUE, raster = FALSE)+ ggtitle(label = "Karaayvaz batch test (PC = 40)")

testPal <- subset(combo.reference, subset = orig.ident == "Pal_Prim")
DimPlot(testPal, reduction = "umap", label = TRUE, raster = FALSE)+ ggtitle(label = "Pal batch test (PC = 40)")

testQian <- subset(combo.reference, subset = orig.ident == "Qian")
DimPlot(testQian, reduction = "umap", label = TRUE, raster = FALSE)+ ggtitle(label = "Qian batch test (PC = 40)")

testSavas <- subset(combo.reference, subset = orig.ident == "Savas")
DimPlot(testSavas, reduction = "umap", label = TRUE, raster = FALSE)+ ggtitle(label = "Savas batch test (PC = 40)")

testoldwu <- subset(combo.reference, subset = orig.ident == "Wu")
DimPlot(testoldwu, reduction = "umap", label = TRUE, raster = FALSE)+ ggtitle(label = "Old Wu batch test (PC = 40)")

testnewwu <- subset(combo.reference, subset = orig.ident == "Wu2021prim")
DimPlot(testnewwu, reduction = "umap", label = TRUE, raster = FALSE)+ ggtitle(label = "New Wu batch test (PC = 40)")

testxu <- subset(combo.reference, subset = orig.ident == "Xu")
DimPlot(testxu, reduction = "umap", label = TRUE, raster = FALSE)+ ggtitle(label = "Xu batch test (PC = 40)")



#better markers =========================

DefaultAssay(combo.reference) <- "RNA"
#https://www.biostars.org/p/395951/

#https://github.com/satijalab/seurat/issues/5738
#https://github.com/satijalab/seurat/issues/5847

combo.reference <- NormalizeData(combo.reference, assay = "RNA")
#combo.reference <- ScaleData(combo.reference, assay = "RNA", vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
#                                                                                 "percent.platelet", "percent.heatshock"))

DimPlot(combo.reference, reduction = "umap", label = TRUE, raster = FALSE)+ ggtitle(label = "Cell Type (PC = 40)") #+ NoLegend()


#NK __________________________________________________________________

NK.mark2 <- list(c("PTPRC","FGFBP2", "SPON2", "KLRF1", "FCGR3A", "KLRD1", "TRDC", "CD3E-",
                   "CD3G-", "CD33-", "EPCAM-",
                   "NCAM1")) #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
combo.reference <- AddModuleScore_UCell(combo.reference, features = NK.mark2, name = "NK.mark2", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = combo.reference, features = "signature_1NK.mark2", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "UCell's NK Signature (PC = 40)")
dev.off()

Idents(combo.reference) <- combo.reference$integrated_snn_res.12
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
# Epi <- WhichCells(NK, expression = signature_1EPCAM > 0.5 & signature_1PTPRC == 0)
# Idents(combo.reference, cells = Epi) <- "Epithelial Cells_NK"


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

Idents(combo.reference) <- combo.reference$integrated_snn_res.12
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
# Tc <- subset(combo.reference, idents = "T Cells")
# # Tc <- AddModuleScore_UCell(Tc, features = "CD33", name = "CD33", assay = "RNA")
# Tc <- AddModuleScore_UCell(Tc, features = "EPCAM", name = "EPCAM", assay = "RNA")
# Tc <- AddModuleScore_UCell(Tc, features = "PTPRC", name = "PTPRC", assay = "RNA")
# # 
# # summary(Tc$signature_1T_karacd4)
# # summary(Tc$signature_1CD33)
# # summary(Tc$signature_1EPCAM)
# # 
# # notenougT <- WhichCells(Tc, expression = signature_1T_karacd4 < 0.05)
# # Idents(combo.reference, cells = notenougT) <- "Unknown_T"
# # 
# Epi <- WhichCells(Tc, expression = signature_1EPCAM > 0.5 & signature_1PTPRC == 0)
# # Idents(combo.reference, cells = Epi) <- "Epithelial Cells_T"


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
# Bc <- subset(combo.reference, idents = "B Cells")
# Bc <- AddModuleScore_UCell(Bc, features = "CD33", name = "CD33", assay = "RNA")
# Bc <- AddModuleScore_UCell(Bc, features = "EPCAM", name = "EPCAM", assay = "RNA")
# Bc <- AddModuleScore_UCell(Bc, features = "PTPRC", name = "PTPRC", assay = "RNA")
# 
# summary(Bc$signature_1B_KaraBrech)
# summary(Bc$signature_1CD33)
# summary(Bc$signature_1EPCAM)
# 
# notenougB <- WhichCells(Bc, expression = signature_1B_KaraBrech < 0.05)
# Idents(combo.reference, cells = notenougB) <- "Unknown_B"
# Epi <- WhichCells(Bc, expression = signature_1EPCAM > 0.5 & signature_1PTPRC == 0)
# Idents(combo.reference, cells = Epi) <- "Epithelial Cells_B"
# 

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
# Plasma <- subset(combo.reference, idents = "Plasma Cells")
# Plasma <- AddModuleScore_UCell(Plasma, features = "CD33", name = "CD33", assay = "RNA")
# Plasma <- AddModuleScore_UCell(Plasma, features = "EPCAM", name = "EPCAM", assay = "RNA")
# Plasma <- AddModuleScore_UCell(Plasma, features = "PTPRC", name = "PTPRC", assay = "RNA")
# 
# summary(Plasma$signature_1plasma_colo)
# summary(Plasma$signature_1CD33)
# summary(Plasma$signature_1EPCAM)
# 
# notenougPlasma <- WhichCells(Plasma, expression = signature_1plasma_colo < 0.05)
# Idents(combo.reference, cells = notenougPlasma) <- "Unknown_Plasma"
# Epi <- WhichCells(Plasma, expression = signature_1EPCAM > 0.5 & signature_1PTPRC == 0)
# Idents(combo.reference, cells = Epi) <- "Epithelial Cells_Plasma"
# 
# 

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

#mast cells __________________________________________________________________

mast_azizitcga <- list(c("PTPRC","KIT",
                         "ENPP3",
                         "ITGAM",
                         "CD33",
                         "IL7R",
                         "MS4A2",
                         "TPSAB1",
                         "CPA4",
                         "HECA",
                         "TPSB2",
                         "EPCAM-"))
combo.reference <- AddModuleScore_UCell(combo.reference, features = mast_azizitcga, name = "mast_azizitcga", assay = "RNA")  
pdf("test.pdf")
FeaturePlot(object = combo.reference, features = "signature_1mast_azizitcga", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Chung's Mast Cell Signature (PC = 40)")
dev.off()


#Non-Immune _________________________________________________________


#gen epi __________________________________________________________________

genepi_baslumgen <- list(c("EPCAM", "PTPRC-"))
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

#stroma __________________________________________________________________

strom_karawu <- list(c("FAP",
                       "COL1A1",
                       "COL3A1",
                       "COL5A1",
                       "ACTA2",
                       "TAGLN",
                       "LUM",
                       "FBLN1",
                       "COL6A3",
                       "COL1A2",
                       "COL6A1",
                       "COL6A2",
                       "PDGFRB",
                       "FGFR1", #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
                       "FGFR2", #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
                       "FGFR3", #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
                       "FGFR4", #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
                       "JAG1",
                       "PTPRC-",
                       "EPCAM-")) #https://docs.abcam.com/pdf/immunology/Guide-to-human-CD-antigens.pdf
combo.reference <- AddModuleScore_UCell(combo.reference, features = strom_karawu, name = "strom_karawu", assay = "RNA")  
FeaturePlot(object = combo.reference, features = "signature_1strom_karawu", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "Karaayvaz's Stroma Signature (PC = 40)")
# 
# 
# combo.reference <- RenameIdents(combo.reference,`22` = "Stroma",`61` = "Stroma", 
#                                 `122` = "Stroma", `30` = "Stroma", `114` = "Stroma",
#                                 `55` = "Stroma", `24` = "Stroma", `40` = "Stroma",
#                                 `39` = "Stroma", `90` = "Stroma", `104` = "Stroma",
#                                 `97` = "Stroma", `26` = "Stroma", `50` = "Stroma")
# 

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

saveRDS(combo.reference, "finalfixed_primobj_withmarkersigs_62822_savecopysafe.rds")

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
saveRDS(Tcellsub.combo, file = "TsubSCTREF_62822.rds")

Tcellsub.combo <-readRDS(file = "TsubSCTREF_62822.rds")
DefaultAssay(Tcellsub.combo) <- "integrated"

#combo.reference <- ScaleData(combo.reference, assay = "RNA", verbose = TRUE)
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

#Epi.all.combo <- readRDS(file = "Epiprim_withUMAPrenamednewmarkes.rds") #start here to get actual clustering

Tcellsub.combo <- FindNeighbors(Tcellsub.combo, reduction = "pca", dims = 1:60)

resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
Tcellsub.combo <- FindClusters(Tcellsub.combo, resolution = resolution.range)
pdf("test.pdf", width = 14.22, height = 13.01)
clustree(Tcellsub.combo, prefix = "integrated_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")
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

pdf("test.pdf")
dev.off()
saveRDS(Tcellsub.combo, file = "TcellsubcomboSCT_6422.rds")
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only")
Tcellsub.combo <-readRDS(file = "Tcellsubcombo.rds")


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
saveRDS(Myecellsub.combo, file = "MyesubSCTREF_62822.rds")

DefaultAssay(Myecellsub.combo) <- "integrated"

#combo.reference <- ScaleData(combo.reference, assay = "RNA", verbose = TRUE)
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

#Epi.all.combo <- readRDS(file = "Epiprim_withUMAPrenamednewmarkes.rds") #start here to get actual clustering

Myecellsub.combo <- FindNeighbors(Myecellsub.combo, reduction = "pca", dims = 1:60)

resolution.range <- c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6)
Myecellsub.combo <- FindClusters(Myecellsub.combo, resolution = resolution.range)
pdf("test.pdf", width = 14.22, height = 13.01)
clustree(Myecellsub.combo, prefix = "integrated_snn_res.", node_colour = "sc3_stability", layout = "sugiyama")
dev.off()

DefaultAssay(Myecellsub.combo) <- "integrated"
Myecellsub.combo <- FindNeighbors(Myecellsub.combo, reduction = "pca", dims = 1:50)
Myecellsub.combo <- FindClusters(Myecellsub.combo, resolution = 3)
Myecellsub.combo <- RunUMAP(Myecellsub.combo, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use=123)

pdf("test.pdf")
DimPlot(Myecellsub.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE) #+ ggtitle(label = "NK Cells from Primary Integrated")
DimPlot(Myecellsub.combo, reduction = "umap", label = TRUE, repel = TRUE, group.by = "old.labels", raster = FALSE) +theme(legend.position = "None")#+ ggtitle(label = "NK Cells from Primary Integrated")
dev.off()
DimPlot(Myecellsub.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "celltype_main") + ggtitle(label = "NK Cells from lumA samples (PC = 15)")
DimPlot(Myecellsub.combo, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE, group.by = "orig.ident")


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

pdf("test.pdf")
dev.off()



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

saveRDS(combo.reference, "PrimObject_6422.rds")

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects")
combo.reference <- readRDS("PrimObject_6422.rds")

#subtype call per cell ===================================================

DefaultAssay(combo.reference) <- "RNA"
#DefaultAssay(combo.reference) <- "SCT"

#https://github.com/satijalab/seurat/issues/5847
Mydata <- subset(combo.reference, idents = "Epithelial Cells")

#Mydata <- subset(combo.reference, idents = c("Unspecified Epithelial Cells"))
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects")
sigdat <- read.csv("NatGen_Supplementary_table_S4.csv")
temp_allgenes <- c(as.vector(sigdat[,"Basal_SC"]),
                   as.vector(sigdat[,"Her2E_SC"]),
                   as.vector(sigdat[,"LumA_SC"]),
                   as.vector(sigdat[,"LumB_SC"]))
temp_allgenes <- unique(temp_allgenes[!temp_allgenes == ""])

#Read in the single cell RDS object as 'Mydata'
# Mydata <- readRDS("path_to_seurat_object.Rdata")
Mydata <- NormalizeData(Mydata, assay = "RNA")
Mydata <- ScaleData(Mydata, features=temp_allgenes, assay = "RNA")

tocalc<-as.data.frame(Mydata@assays$RNA@scale.data)
#tocalc<-as.data.frame(Mydata@assays$SCT@scale.data)

#calculate mean scsubtype scores
outdat <- matrix(0,
                 nrow=ncol(sigdat),
                 ncol=ncol(tocalc),
                 dimnames=list(colnames(sigdat),
                               colnames(tocalc)))
for(i in 1:ncol(sigdat)){
  
  # sigdat[i,!is.na(sigdat[i,])]->module
  # row <- as.character(unlist(module))
  row <- as.character(sigdat[,i])
  row<-unique(row[row != ""])
  genes<-which(rownames(tocalc) %in% row)
  
  temp<-apply(tocalc[genes,],2,function(x){mean(as.numeric(x),na.rm=TRUE)})
  
  outdat[i,]<-as.numeric(temp)
  
}

final<-outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
final<-as.data.frame(final)
is.num <- sapply(final, is.numeric);final[is.num] <- lapply(final[is.num], round, 4)
finalm<-as.matrix(final)

##Scaling scores function before calling the highest Call
center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
  get_average <- function(v) sum(v * row.w)/sum(row.w)
  average <- apply(x, 2, get_average)
  sweep(x, 2, average)
}

##Obtaining the highest call
finalmt<-as.data.frame(t(finalm))
finalm.sweep.t<-center_sweep(finalmt)
Finalnames<-colnames(finalm.sweep.t)[max.col(finalm.sweep.t,ties.method="first")]
finalm.sweep.t$SCSubtypeCall <- Finalnames

##Writing out output files (rownames remain the same for both)

write.csv(finalm.sweep.t, "Mydata_Scores_62922.csv")
#write.csv(finalm.sweep.t, "Mydata_Scores_6622_includingAzzi880.csv")

#write.table(finalm.sweep.t, "Mydata_Scores.txt", sep="\t")
#write.table(Finalnames, "Mydata_CALLS.txt", sep="\t")

Mydata <- AddMetaData(Mydata, finalm.sweep.t[,5], col.name = "sc50.Pred")
combo.reference <- AddMetaData(combo.reference, finalm.sweep.t[,5, drop = F], col.name = "sc50.Pred")

cancer.epi <- Mydata
cancer.epi$samples <- paste(cancer.epi$samples, cancer.epi$BC.Subtype, sep = "_")
sobjlists <- FetchData(object = cancer.epi, vars = c("samples","Patient", "BC.Subtype", "sc50.Pred"))
library(reshape2)

sobjlists = melt(sobjlists, id.vars = c("samples","Patient", "BC.Subtype", "sc50.Pred"))
sobjlists$cells <- rownames(sobjlists)
head(sobjlists)

counts_call <- table(sobjlists$samples, sobjlists$sc50.Pred)
counts_call <- as.data.frame.matrix(counts_call) 
counts_call$max <- colnames(counts_call)[apply(counts_call,1,which.max)]
head(counts_call)
#write.csv(counts_call, file = "call_count_persample_6622_ucelllabel.csv")
#write.csv(counts_call, file = "call_count_persample_6622_old_label.csv")
write.csv(counts_call, file = "call_count_persample_62922_RNAnormscale.csv")
#write.csv(counts_call, file = "call_count_persample_6622_withAzizi880.csv")


#labels ==================================================================

#no shortcut really for this one -> ran loop, did ablebits to merge .csv files into one workbook
#then foundmost common call for cells in sample, and that's the call for the sample

combo.reference$pam50.pred <- combo.reference$samples
combo.reference$pam50.pred <- paste(combo.reference$pam50.pred, combo.reference$BC.Subtype, sep = "_")

#Pal
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0001_HR_HR+")] <- "LumA"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0025_HR_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0029_7C_HR_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0029_9C_HR_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0031_Her2_HER2+")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0032_HR_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0040_HR_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0042_HR_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0043_HR_HR+")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0056_HR_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0064_HR_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0069_Her2_HER2+")] <- "Her2" 
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0106_TNBC_TNBC")] <- "Basal" 
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0114_T3_HR_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0114_TNBC_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0125_HR_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0126_TNBC_TNBC")] <- "Basal" 
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0135_TNBC_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0151_HR_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0161_Her2_HER2+")] <- "LumB" 
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0163_HR_HR+")] <- "LumA"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0167_HR_HR+")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0173_HR_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0176_Her2_HER2+")] <- "LumB" #*******
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0308_Her2_HER2+")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0319_HR_HR+")] <- "Basal" #*****
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0337_Her2_HER2+")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0360_HR_HR+")] <- "LumB"

#Qian
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_1_HER2+")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_10_TNBC")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_11_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_12_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_13_HR+")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_14_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_2_TNBC")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_4_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_5_TNBC")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_6_HER2+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_7_TNBC")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_8_HER2+")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC_9_HER2+")] <- "Her2"

combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID3586_HER2+")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID3838_HER2+")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID3921_HER2+")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID3941_HR+")] <- "LumA"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID3946_TNBC")] <- "TNBC"#"Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID3948_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4040_HR+")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4067_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4290A_HR+")] <- "LumA"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID44041_TNBC")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4461_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4463_HR+")] <- "LumA"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4465_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4471_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4495_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID44971_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4515_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID45171_HER2+")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4530N_HR+")] <- "LumA"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4535_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4398_HR+")] <- "HR+"


combo.reference$pam50.pred[which(combo.reference$pam50.pred == "P1prim_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "P2prim_HR+")] <- "LumA"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "P3prim_HER2+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "P4prim_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "P5prim_HER2+")] <- "Her2"

combo.reference$pam50.pred[which(combo.reference$pam50.pred == "Patient 1_TNBC")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "Patient 2_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "Patient 3_TNBC")] <- "Her2"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "Patient 5_TNBC")] <- "Basal"

combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT039_P1_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT039_P10_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT039_P11_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT039_P121_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT058_P1_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT081_P1_TNBC")] <- "Basal"#"Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT081_P3_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT081_P5_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT089_P1_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT089_P3_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT089_P7_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT089_P9_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT126_P3_TNBC")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "PT126_P7_TNBC")] <- "Basal"


combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC10_T1_1_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC11_T1_1_HER2+")] <- "HER2+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC11_T2_1_HER2+")] <- "HER2+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC11_T2_HER2+")] <- "HER2+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC1_TUMOR_1_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC1_TUMOR_3_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC1_TUMOR_4_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC2_TUMOR_1_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC2_TUMOR_2_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC2_TUMOR_3_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC2_TUMOR_4_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC3_TUMOR_1_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC3_TUMOR_3_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC3_TUMOR_5_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC4_TUMOR_1_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC4_TUMOR_2_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC4_TUMOR_3_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC4_TUMOR_4_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC4_TUMOR_5_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC4_TUMOR_6_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC5_TUMOR_1_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC5_TUMOR_2_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC5_TUMOR_3_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC5_TUMOR_4_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC6_TUMOR_1_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC6_TUMOR_2_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC6_TUMOR_3_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC7_TUMOR_2_HER2+")] <- "HER2+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC7_TUMOR_3_HER2+")] <- "HER2+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC8_TUMOR_1_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC8_TUMOR_2_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC8_TUMOR_3_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC10_T1_1_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC9_T1_HR+")] <- "HR+"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "BC9_T2_2_HR+")] <- "HR+"


combo.reference$pam50.pred[which(combo.reference$pam50.pred == "TIL20_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "TIL32_TNBC")] <- "TNBC"

table(combo.reference$pam50.pred)

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/")
#saveRDS(combo.reference, "PrimObject_6722_noambigAzizi_SCTnormRNA.rds")
saveRDS(combo.reference, "FINALprimObj_withsubtypes.rds")
#combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4535_ER+")] <- "LumB"
#combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4530N_ER+")] <- "LumA"

#combo.reference$pam50.pred[which(combo.reference$pam50.pred == "lumA")] <- "LumA"
#combo.reference$pam50.pred[which(combo.reference$pam50.pred == "lumB")] <- "LumB"
#combo.reference$pam50.pred[which(combo.reference$pam50.pred == "Normal-like")] <- "TNBC"
