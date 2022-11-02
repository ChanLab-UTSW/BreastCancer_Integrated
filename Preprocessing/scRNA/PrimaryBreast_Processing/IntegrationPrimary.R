

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This is the third step in processing the scRNA-seq  
# primary breast datasets, and requires all filtering done
# from "IndividualDatasetLoading_Primary.R" and "DoubletFinder_Primary.R".

# This is SCTransformed and uses reference-based integration.
# The technology used for cell capture is the batch being corrected for.

# The output of this R script will be the input 
# of "CellTypeLabeling.R"

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DoubletResult_Dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/PrimaryBreast_Processing/Individual_Dataset_Objects/withDoubletMeta"
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




#just filtered objects =======================================================

setwd(DoubletResult_Dir)
AziziPrim <- readRDS("AziziPrim_withDoublet_62622.rds")
AziziT <- readRDS("AziziT_withDoublet_62622.rds")
Karaayvaz <- readRDS("Karaayvaz_withDoublet_62622.rds")
Pal <- readRDS("Pal_withDoublet_62622.rds")
Qian <- readRDS("Qian_withDouble_62622t.rds")
Savas <- readRDS("Savas_withDoublet_62622.rds")
OldWu <- readRDS("OldWu_withDoublet_62622.rds")
Wu2021 <- readRDS("Wu2021_withDoublet_62622.rds")
Xu <- readRDS("Xu_withDoublet_62622.rds")


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

setwd(finalQC_Dir)
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

setwd(finalQC_Dir)
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

setwd(finalQC_Dir)
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

setwd(finalQC_Dir)
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

setwd(finalQC_Dir)
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

setwd(finalQC_Dir)
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

setwd(finalQC_Dir)
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


setwd(finalQC_Dir)
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


setwd(finalQC_Dir)
saveRDS(Xu.fixed, "Xu_fixed_62722.rds")


#Merge Datasets (fixed) =============================================================

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



# ==============================================================

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
setwd(Integration_Results_Dir)
saveRDS(combo.anchors, file = "primOnlyANCHORS_4ref_fromdoubletfixed627_102822.rds")
#saveRDS(combo.anchors, file = "primOnlyYESSSREFSCTAnchors_BiggestRef_doubletRemoved_62722.rds")


combo.reference <- IntegrateData(anchorset = combo.anchors, normalization.method = "SCT")
saveRDS(combo.reference, file = "primOnlyINTEGRATED_4ref_fromdoubletfixed627_102822.rds")
#saveRDS(combo.reference, file = "primOnlyYESSSREFsctINTEGRATED_BiggestRef_doubletRemoved_62722.rds")

#combo.reference <- readRDS(file = "primOnlyYESSSREFsctINTEGRATED_BiggestRef_doubletRemoved_62722.rds")
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

#ggsave("BatchTestSCTRef_origident_71922_long.pdf", plot = p2, width = 16, height = 9)



ElbowPlot(combo.reference, ndims = 100)
DimHeatmap(combo.reference, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(combo.reference, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(combo.reference, dims = 30:45, cells = 500, balanced = TRUE)
DimHeatmap(combo.reference, dims = 45:55, cells = 500, balanced = TRUE)
DimHeatmap(combo.reference, dims = 55:65, cells = 500, balanced = TRUE)
DimHeatmap(combo.reference, dims = 65:75, cells = 500, balanced = TRUE)


DefaultAssay(combo.reference) <- "integrated"

combo.reference <- FindNeighbors(combo.reference, reduction = "pca", dims = 1:70, nn.method = "rann")#, k.param= 500)
combo.reference <- FindClusters(combo.reference, resolution = 12)
combo.reference <- RunUMAP(combo.reference, reduction = "pca", dims = 1:70, verbose = TRUE, seed.use = 123)

combo.reference$BC.Subtype[which(combo.reference$BC.Subtype == "ER+")] <- "HR+"
combo.reference$Patient[which(combo.reference$Patient == "BC11_1")] <- "BC11"
combo.reference$Patient[which(combo.reference$Patient == "BC11_2")] <- "BC11"
Idents(combo.reference) <- combo.reference$celltype_final



#pdf("PrimUMAP_withTMyeloidSubs_withlegend.pdf", width = 14.22, height = 13.01)
p <- DimPlot(combo.reference, reduction = "umap", label = F, repel = T, raster = FALSE) + ggtitle(label = " ") #,
cols = rev(c("#2C1439","#FF94FA","#7488E0",
             "#000000","#6BD2AA","#FFD61F",
             "#00536A","#FF2AB4","#2C1439","#00BFB7",
             "#D47100","#9300C9","#AA4300","#DD7D72",
             "#3F9BFF","#C07FFF","#8CBF00"))) + ggtitle(label = " ")

p <- DimPlot(combo.reference, reduction = "umap", label = F, repel = T, raster = FALSE, 
             order = c("NK Cells",  "Myoepithelial Cells",
                       "MDSCs","Dendritic Cells","Regulatory T Cells",
                       "Perivascular-like (PVL) Cells","B Cells", "CD4+ T Cells", "Neutrophils",
                       "Monocytes",  "Macrophages", "Mast Cells","CD8+ T Cells", "Fibroblasts",
                       "Endothelial Cells", "Plasma Cells", "Epithelial Cells"),
             cols = c("#D47100","#AA4300","#00BFB7",
                      "#DD7D72","#C07FFF","#FF94FA","#9300C9","#6BD2AA",
                      "#E90000","#3F9BFF","#FF2AB4","#2C1439",
                      "#8CBF00","#FFD61F","#000000",
                      "#7488E0","#00536A")) + ggtitle(label = " ")
#https://stackoverflow.com/questions/36474643/graphical-parameters-in-ggplot2-how-to-change-axis-tick-thickness
pdf("PRIMUMAPfig1a_72522.pdf", width = 16.22, height = 14.22)
pdf("test.pdf", width = 16.22, height = 14.22)
p + theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()

pdf("BCsub_wholeUMAP_72822.pdf", width = 16.22, height = 14.22)
DimPlot(combo.reference, reduction = "umap", group.by = "BC.Subtype", raster = FALSE) + ggtitle(label = " ") +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()
DimPlot(combo.reference, reduction = "umap", label = TRUE, group.by = "celltype_minor", raster = FALSE)+ ggtitle(label = "Cell Type (PC = 40)")
DimPlot(combo.reference, reduction = "umap", label = TRUE, group.by = "celltype_main", raster = FALSE)+ ggtitle(label = "Cell Type (PC = 40)")

pdf("PrimUMAPbigger_63022.pdf", width = 19, height = 15.3)
pdf("test.pdf", width = 19, height = 15.3)
p + #guides(color = guide_legend(override.aes = list(size=10), ncol=1) )
  theme(axis.line = element_line(colour = 'black', size = 2)) + 
  theme(axis.ticks = element_line(colour = "black", size = 2)) +
  theme(text = element_text(size = 35)) +
  theme(axis.text = element_text(size = 20))
dev.off()

setwd(Integration_Results_Dir)
saveRDS(combo.reference, "FinalPrim_noZgeneDem_71222.rds")


# Figure S1.C_ BAR patients per age group -------------------------------------------------------


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

