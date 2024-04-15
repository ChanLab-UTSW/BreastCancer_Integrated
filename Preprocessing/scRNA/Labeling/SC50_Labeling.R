
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This is the fifth and final step in processing the scRNA-seq  
# primary breast datasets, and requires the epithelial cells
# to be labeled. This applies Wu (2021)'s sc50 Subtype
# classifier.

# This same script is also used to label the subtypes of
# Bassez's antiPD1 epithelial cells.

# Wu (2021)
# Paper link: https://www.nature.com/articles/s41588-021-00911-1
# code link: https://github.com/Swarbricklab-code/BrCa_cell_atlas/blob/main/scSubtype/Highest_calls.R

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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





# subtype call per cell (Primary) ===================================================

DefaultAssay(combo.reference) <- "RNA"

#https://github.com/satijalab/seurat/issues/5847
Mydata <- subset(combo.reference, idents = "Epithelial Cells")

#Mydata <- subset(combo.reference, idents = c("Unspecified Epithelial Cells"))
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


write.csv(finalm.sweep.t, "PrimObject_Mydata_Scores_71222.csv")

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

write.csv(counts_call, file = "call_count_persample_71222_RNAnormscale.csv")


# labels (Primary) ==================================================================


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
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "0043_HR_HR+")] <- "LumB"
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
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID3946_TNBC")] <- "TNBC"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID3948_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4040_HR+")] <- "Basal"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4067_HR+")] <- "LumB"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID4290A_HR+")] <- "LumA"
combo.reference$pam50.pred[which(combo.reference$pam50.pred == "CID44041_TNBC")] <- "Basal"
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

saveRDS(combo.reference, "PrimObject_FINALnoZgenedem_71222.rds")

# pam50 labeling (antiPD1) =================================


DefaultAssay(sobj1) <- "RNA"
Mydata <- subset(sobj1, idents = "Epithelial Cells")
DefaultAssay(sobj2) <- "RNA"
Mydata <- subset(sobj2, idents = "Epithelial Cells")

sigdat <- read.csv("NatGen_Supplementary_table_S4.csv")

temp_allgenes <- c(as.vector(sigdat[,"Basal_SC"]),
                   as.vector(sigdat[,"Her2E_SC"]),
                   as.vector(sigdat[,"LumA_SC"]),
                   as.vector(sigdat[,"LumB_SC"]))
temp_allgenes <- unique(temp_allgenes[!temp_allgenes == ""])

#Read in the single cell RDS object as 'Mydata'
Mydata <- NormalizeData(Mydata, assay = "RNA")
Mydata <- ScaleData(Mydata, features=temp_allgenes, assay = "RNA")

tocalc<-as.data.frame(Mydata@assays$RNA@scale.data)

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

write.csv(finalm.sweep.t, "Mydata_Scores_antiPD1nochemo_71322.csv")
write.csv(finalm.sweep.t, "Mydata_Scores_antiPD1yeschemo_71322.csv")

Mydata <- AddMetaData(Mydata, finalm.sweep.t[,5], col.name = "sc50.Pred")
sobj1 <- AddMetaData(sobj1, finalm.sweep.t[,5, drop = F], col.name = "sc50.Pred")
sobj2 <- AddMetaData(sobj2, finalm.sweep.t[,5, drop = F], col.name = "sc50.Pred")

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
write.csv(counts_call, file = "call_count_persample_antiPD1nochemo_71322.csv")
write.csv(counts_call, file = "call_count_persample_antiPD1yeschemo_71322.csv")

#label subtypes (antiPD1) ===================================================

sobj1$pam50.pred <- sobj1$Patient
sobj1$pam50.pred <- paste(sobj1$pam50.pred, sobj1$BC.Subtype, sep = "_")

sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_1_TNBC")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_10_TNBC")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_11_TNBC")] <- "Her2"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_12_HR+")] <- "LumB"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_13_HER2+")] <- "Her2"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_14_TNBC")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_15_TNBC")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_16_TNBC")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_17_HR+")] <- "LumB"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_18_HR+")] <- "LumA"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_19_TNBC")] <- "LumB"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_2_TNBC")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_20_HR+")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_21_HR+")] <- "LumA"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_22_HR+")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_23_HER2+")] <- "LumA"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_24_HR+")] <- "LumA"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_25_TNBC")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_26_TNBC")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_27_HR+")] <- "LumB"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_28_HER2+")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_29_HR+")] <- "LumA"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_3_HR+")] <- "Her2"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_30_HR+")] <- "LumB"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_31_TNBC")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_4_HR+")] <- "LumB"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_5_HR+")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_6_HR+")] <- "LumB"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_7_HR+")] <- "LumA"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_8_TNBC")] <- "Basal"
sobj1$pam50.pred[which(sobj1$pam50.pred == "BIOKEY_9_TNBC")] <- "Basal"

table(sobj1$pam50.pred)

sobj2$pam50.pred <- sobj2$Patient
sobj2$pam50.pred <- paste(sobj2$pam50.pred, sobj2$BC.Subtype, sep = "_")

sobj2$pam50.pred[which(sobj2$pam50.pred == "BIOKEY_32_HR+")] <- "LumB"
sobj2$pam50.pred[which(sobj2$pam50.pred == "BIOKEY_33_TNBC")] <- "Basal"
sobj2$pam50.pred[which(sobj2$pam50.pred == "BIOKEY_34_TNBC")] <- "LumB"
sobj2$pam50.pred[which(sobj2$pam50.pred == "BIOKEY_35_TNBC")] <- "LumA"
sobj2$pam50.pred[which(sobj2$pam50.pred == "BIOKEY_36_TNBC")] <- "Basal"
sobj2$pam50.pred[which(sobj2$pam50.pred == "BIOKEY_37_HR+")] <- "LumB"
sobj2$pam50.pred[which(sobj2$pam50.pred == "BIOKEY_38_HER2+")] <- "LumB"
sobj2$pam50.pred[which(sobj2$pam50.pred == "BIOKEY_39_TNBC")] <- "Her2"
sobj2$pam50.pred[which(sobj2$pam50.pred == "BIOKEY_40_HR+")] <- "LumB"
sobj2$pam50.pred[which(sobj2$pam50.pred == "BIOKEY_41_TNBC")] <- "Basal"
sobj2$pam50.pred[which(sobj2$pam50.pred == "BIOKEY_42_HR+")] <- "LumB"

table(sobj2$pam50.pred)



saveRDS(sobj1, file = "antiPD1nochemo_withpam50_nozallgenedem_71322.rds")
saveRDS(sobj2, file = "antiPD1YESchemo_withpam50_nozallgenedem_71322.rds")
