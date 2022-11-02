


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This will generate all of the sub-figures for Figure S3,
# and contains additional analyses for the NK subsets and rNK population.


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


PrimDir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/PrimaryBreast_Processing/Main_Prim_Object"
NKsubDir <- "/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only"
DEGdir <- "/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only"

rNKdir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig2"


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






# ================================================================== ======

# FUNCTIONS ----- ------------------------------

DEG_Remove_mito <- function (df){
  df_rm_mito <- df[!grepl("^MT-|^MT.",rownames(df)),]
  return(df_rm_mito)
}

# ================================================================== ======
#S3A NKsub clustree  ========================
# ================================================================== ======
# load in Prim object ------

setwd(PrimDir)
combo.reference <- readRDS("PrimObject_withreprog_noZallgenedem_71322.rds")
DefaultAssay(combo.reference) <- "RNA"
combo.reference <- NormalizeData(combo.reference, assay = "RNA")


# NK clustering up until clustree ------------------

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

# ================================================================== ======
#S3B NKsub UMAPs  ========================
# ================================================================== ======

# load in NKsub ----------

setwd(NKsubDir)
NK.all.combo <- readRDS("NKall_PC17manhattan0.4res_71322.rds")
DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")

Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub



# NK-0 ------

setwd(DEGdir)
c0.mark <- readRDS("NK.c0markers_71322.rds")

library(RColorBrewer)
NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = c0.mark, name = "c0.mark.final2", assay = "RNA")
#pdf("NKc0_7522.pdf", width = 5.11, height = 4.5)
pdf("NKc0_71922.pdf", width = 5.11, height = 4.5)
p <- FeaturePlot(object = NK.all.combo, features = "signature_1c0.mark", order = TRUE, label = FALSE, 
                 repel = TRUE, min.cutoff = 0, raster = FALSE, pt.size = 1.5) +
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_color_viridis(option = "D")+
  ggtitle(label = " ")
p + theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 25)) +
  theme(axis.text = element_text(size = 20))
dev.off()

# NK-1 ----

setwd(DEGdir)
c1.mark <- readRDS("NK.c1markers_71322.rds")

NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = c1.mark, name = "c1.mark", assay = "RNA")
#pdf("NKc1_7422.pdf", width = 7.11, height = 6.5)
pdf("NKc1_71922.pdf", width = 5.11, height = 4.5)
p <- FeaturePlot(object = NK.all.combo, features = "signature_1c1.mark", order = TRUE, 
                 label = FALSE, repel = TRUE, min.cutoff = 0, raster = FALSE, pt.size = 1.5) +
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_color_viridis(option = "D")+
  ggtitle(label = " ")
p + theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 25)) +
  theme(axis.text = element_text(size = 20))
dev.off()

# NK-2 ----

setwd(DEGdir)
c2.mark <- readRDS("NK.c2markers_71322.rds")

NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = c2.mark, name = "c2.mark", assay = "RNA")
pdf("NKc2_71922.pdf",width = 5.11, height = 4.5)
p <- FeaturePlot(object = NK.all.combo, features = "signature_1c2.mark", order = TRUE, 
                 label = FALSE, repel = TRUE, min.cutoff = 0, raster = FALSE,
                 pt.size = 1.5) +
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_color_viridis(option = "D")+
  ggtitle(label = " ")
p + theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 25)) +
  theme(axis.text = element_text(size = 20))
dev.off()


# NK-3 ----

setwd(DEGdir)
c3.mark <- readRDS("NK.c3markers_71322.rds")

NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = c3.mark, name = "c3.mark", assay = "RNA")
pdf("NKc3_71922.pdf", width = 5.11, height = 4.5)
p <- FeaturePlot(object = NK.all.combo, features = "signature_1c3.mark", 
                 order = TRUE, label = FALSE, repel = TRUE, min.cutoff = 0, 
                 raster = FALSE, pt.size = 1.5) +
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_color_viridis(option = "D")+
  ggtitle(label = " ")
p + theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 25)) +
  theme(axis.text = element_text(size = 20))
dev.off()


# NK-4 ----

setwd(DEGdir)
c4.mark <- readRDS("NK.c4markers_71322.rds")

NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = c4.mark, name = "c4.mark", assay = "RNA")
pdf("NKc4_71922.pdf", width = 5.11, height = 4.5)
p <- FeaturePlot(object = NK.all.combo, features = "signature_1c4.mark", 
                 order = TRUE, label = FALSE, repel = TRUE, 
                 min.cutoff = 0, raster = FALSE, pt.size = 1.5) +
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_color_viridis(option = "D")+
  ggtitle(label = " ")
p + theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 25)) +
  theme(axis.text = element_text(size = 20))
dev.off()


# NK-5 ----

setwd(DEGdir)
c5.mark <- readRDS("NK.c5markers_71322.rds")

NK.all.combo <- AddModuleScore_UCell(NK.all.combo, features = c5.mark, name = "c5.mark", assay = "RNA")
pdf("NKc5_71922.pdf", width = 5.11, height = 4.5)
p <- FeaturePlot(object = NK.all.combo, features = "signature_1c5.mark", 
                 order = TRUE, label = FALSE, repel = TRUE, 
                 min.cutoff = 0, raster = FALSE, pt.size = 1.5) +
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_color_viridis(option = "D")+
  ggtitle(label = " ")
p + theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 25)) +
  theme(axis.text = element_text(size = 20))
dev.off()

# ================================================================== ======
#S3C NKsub MA plots  ========================
# ================================================================== ======


# DEGs -------

Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub


DefaultAssay(NK.all.combo) <- "RNA"
NK.all.combo <- NormalizeData(NK.all.combo, assay = "RNA")


NKsub.mark <- FindAllMarkers(NK.all.combo,
                              test.use = "MAST")

setwd(DEGdir)
write.csv(NKsub.mark, "NKsubDEG_ABSOLUTELYNOTHRESH.csv")

setwd(DEGdir)
j.markers_DGE_filtered <- read.csv("NKsubDEG_ABSOLUTELYNOTHRESH.csv")
j.markers_DGE_filtered <- DEG_Remove_mito(j.markers_DGE_filtered)
write.csv(j.markers_DGE_filtered, "NKsub_NOTHRESHnomito.csv")


# MA plot NK-0 -------

Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub
# dfsample <- read.csv("NKclustDEGs_NOTHRESH.csv", header = T, row.names = 1)
dfsample <- read.csv("rNKvnon_DEGs_NOTHRESH_71422.csv", header = T, row.names = 1)
dfsample <- DEG_Remove_mito(dfsample)

avgexp <- AverageExpression(object = NK.all.combo, features = unique(dfsample$gene), assays = c("RNA", "SCT"))
avgexp <- as.data.frame(avgexp[['RNA']])
avgexp <- avgexp[which(rownames(avgexp) %in% unique(dfsample$gene)),]

head(avgexp)

##repeat following for each NKsubset
c0DEGs <- dfsample[dfsample$cluster=="0",]
avgexpc0 <- avgexp[rownames(avgexp) %in% c0DEGs$gene,"0", drop = F]


##repeat following for each NKsubset
data <- cbind(c0DEGs, avgexpc0)
data_input <- data[,c(9,2, 5, 7)]


head(data)
colnames(data)

which(data_input$gene == "KLRG1")
which(data_input$gene == "TIGIT")
which(data_input$gene == "NR4A3")

colnames(data_input) <- c("baseMean", "log2FoldChange", "padj", "gene")

library(ggpubr)
options(ggrepel.max.overlaps = 15)
#options(ggrepel.max.overlaps = 300)

p <- ggmaplot(data_input, 
              fdr = 0.05, fc = 1.5, 
              #fdr = 0.05, fc = 1.7, 
              size = 1,
              palette = c("#B31B21", "#1465AC", "darkgray"),
              top = 0,              
              #top = 100,
              # label.select = c("NR4A1", "NR4A2", "NR4A3", "RHOB", #up
              #                  "HSPA1A", "FOS", "DUSP1", "JUN", "DNAJB1", "HSPA1B",
              #                  "FOSB", "TNFAIP3",
              #                  "MUCL1", "CYBA"), #down
              #label.select = c("NR4A1", "NR4A2", "NR4A3"),
              genenames = as.vector(data_input$gene),
              
              font.label = c("plain", 11),
              select.top.method = "fc",
              ggtheme = ggplot2::theme_classic()) +
  xlim(0,7) + ylim(-2.5,4) + NoLegend()

ggsave("c0_MAplot_72822.pdf", plot = p, width = 4.3, height = 4)


# MA plot NK-1 -------

Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub
# dfsample <- read.csv("NKclustDEGs_NOTHRESH.csv", header = T, row.names = 1)
dfsample <- read.csv("rNKvnon_DEGs_NOTHRESH_71422.csv", header = T, row.names = 1)
dfsample <- DEG_Remove_mito(dfsample)

avgexp <- AverageExpression(object = NK.all.combo, features = unique(dfsample$gene), assays = c("RNA", "SCT"))
avgexp <- as.data.frame(avgexp[['RNA']])
avgexp <- avgexp[which(rownames(avgexp) %in% unique(dfsample$gene)),]

head(avgexp)


c1DEGs <- dfsample[dfsample$cluster=="1",]
avgexpc1 <- avgexp[rownames(avgexp) %in% c1DEGs$gene,"1", drop = F]


data <- cbind(c1DEGs, avgexpc1)
data_input <- data[,c(9,2, 5, 7)]


head(data)
colnames(data)

which(data_input$gene == "KLRG1")
which(data_input$gene == "TIGIT")
which(data_input$gene == "NR4A3")

colnames(data_input) <- c("baseMean", "log2FoldChange", "padj", "gene")

library(ggpubr)
options(ggrepel.max.overlaps = 15)
#options(ggrepel.max.overlaps = 300)

p <- ggmaplot(data_input, 
              fdr = 0.05, fc = 1.5, 
              #fdr = 0.05, fc = 1.7, 
              size = 1,
              palette = c("#B31B21", "#1465AC", "darkgray"),
              top = 0,              
              #top = 100,
              # label.select = c("NR4A1", "NR4A2", "NR4A3", "RHOB", #up
              #                  "HSPA1A", "FOS", "DUSP1", "JUN", "DNAJB1", "HSPA1B",
              #                  "FOSB", "TNFAIP3",
              #                  "MUCL1", "CYBA"), #down
              #label.select = c("NR4A1", "NR4A2", "NR4A3"),
              genenames = as.vector(data_input$gene),
              
              font.label = c("plain", 11),
              select.top.method = "fc",
              ggtheme = ggplot2::theme_classic()) +
  xlim(0,7) + ylim(-2.5,4) + NoLegend()

ggsave("c1_MAplot_72822.pdf", plot = p, width = 4.3, height = 4)


# MA plot NK-2 -------

Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub
# dfsample <- read.csv("NKclustDEGs_NOTHRESH.csv", header = T, row.names = 1)
dfsample <- read.csv("rNKvnon_DEGs_NOTHRESH_71422.csv", header = T, row.names = 1)
dfsample <- DEG_Remove_mito(dfsample)

avgexp <- AverageExpression(object = NK.all.combo, features = unique(dfsample$gene), assays = c("RNA", "SCT"))
avgexp <- as.data.frame(avgexp[['RNA']])
avgexp <- avgexp[which(rownames(avgexp) %in% unique(dfsample$gene)),]

head(avgexp)

c2DEGs <- dfsample[dfsample$cluster=="2",]
avgexpc2 <- avgexp[rownames(avgexp) %in% c2DEGs$gene,"2", drop = F]

data <- cbind(c2DEGs, avgexpc2)
data_input <- data[,c(9,2, 5, 7)]


head(data)
colnames(data)

which(data_input$gene == "KLRG1")
which(data_input$gene == "TIGIT")
which(data_input$gene == "NR4A3")

colnames(data_input) <- c("baseMean", "log2FoldChange", "padj", "gene")

library(ggpubr)
options(ggrepel.max.overlaps = 15)
#options(ggrepel.max.overlaps = 300)

p <- ggmaplot(data_input, 
              fdr = 0.05, fc = 1.5, 
              #fdr = 0.05, fc = 1.7, 
              size = 1,
              palette = c("#B31B21", "#1465AC", "darkgray"),
              top = 0,              
              #top = 100,
              # label.select = c("NR4A1", "NR4A2", "NR4A3", "RHOB", #up
              #                  "HSPA1A", "FOS", "DUSP1", "JUN", "DNAJB1", "HSPA1B",
              #                  "FOSB", "TNFAIP3",
              #                  "MUCL1", "CYBA"), #down
              #label.select = c("NR4A1", "NR4A2", "NR4A3"),
              genenames = as.vector(data_input$gene),
              
              font.label = c("plain", 11),
              select.top.method = "fc",
              ggtheme = ggplot2::theme_classic()) +
  xlim(0,7) + ylim(-2.5,4) + NoLegend()

ggsave("c2_MAplot_72822.pdf", plot = p, width = 4.3, height = 4)


# MA plot NK-3 -------

Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub
# dfsample <- read.csv("NKclustDEGs_NOTHRESH.csv", header = T, row.names = 1)
dfsample <- read.csv("rNKvnon_DEGs_NOTHRESH_71422.csv", header = T, row.names = 1)
dfsample <- DEG_Remove_mito(dfsample)

avgexp <- AverageExpression(object = NK.all.combo, features = unique(dfsample$gene), assays = c("RNA", "SCT"))
avgexp <- as.data.frame(avgexp[['RNA']])
avgexp <- avgexp[which(rownames(avgexp) %in% unique(dfsample$gene)),]

head(avgexp)

c3DEGs <- dfsample[dfsample$cluster=="3",]
avgexpc3 <- avgexp[rownames(avgexp) %in% c3DEGs$gene,"3", drop = F]

data <- cbind(c3DEGs, avgexpc3)
data_input <- data[,c(9,2, 5, 7)]


head(data)
colnames(data)

which(data_input$gene == "KLRG1")
which(data_input$gene == "TIGIT")
which(data_input$gene == "NR4A3")

colnames(data_input) <- c("baseMean", "log2FoldChange", "padj", "gene")

library(ggpubr)
options(ggrepel.max.overlaps = 15)
#options(ggrepel.max.overlaps = 300)

p <- ggmaplot(data_input, 
              fdr = 0.05, fc = 1.5, 
              #fdr = 0.05, fc = 1.7, 
              size = 1,
              palette = c("#B31B21", "#1465AC", "darkgray"),
              top = 0,              
              #top = 100,
              # label.select = c("NR4A1", "NR4A2", "NR4A3", "RHOB", #up
              #                  "HSPA1A", "FOS", "DUSP1", "JUN", "DNAJB1", "HSPA1B",
              #                  "FOSB", "TNFAIP3",
              #                  "MUCL1", "CYBA"), #down
              #label.select = c("NR4A1", "NR4A2", "NR4A3"),
              genenames = as.vector(data_input$gene),
              
              font.label = c("plain", 11),
              select.top.method = "fc",
              ggtheme = ggplot2::theme_classic()) +
  xlim(0,7) + ylim(-2.5,4) + NoLegend()

ggsave("c3_MAplot_72822.pdf", plot = p, width = 4.3, height = 4)


# MA plot NK-4 -------

Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub
# dfsample <- read.csv("NKclustDEGs_NOTHRESH.csv", header = T, row.names = 1)
dfsample <- read.csv("rNKvnon_DEGs_NOTHRESH_71422.csv", header = T, row.names = 1)
dfsample <- DEG_Remove_mito(dfsample)

avgexp <- AverageExpression(object = NK.all.combo, features = unique(dfsample$gene), assays = c("RNA", "SCT"))
avgexp <- as.data.frame(avgexp[['RNA']])
avgexp <- avgexp[which(rownames(avgexp) %in% unique(dfsample$gene)),]

head(avgexp)

c4DEGs <- dfsample[dfsample$cluster=="4",]
avgexpc4 <- avgexp[rownames(avgexp) %in% c4DEGs$gene,"4", drop = F]

data <- cbind(c4DEGs, avgexpc4)
data_input <- data[,c(9,2, 5, 7)]


head(data)
colnames(data)

which(data_input$gene == "KLRG1")
which(data_input$gene == "TIGIT")
which(data_input$gene == "NR4A3")

colnames(data_input) <- c("baseMean", "log2FoldChange", "padj", "gene")

library(ggpubr)
options(ggrepel.max.overlaps = 15)
#options(ggrepel.max.overlaps = 300)

p <- ggmaplot(data_input, 
              fdr = 0.05, fc = 1.5, 
              #fdr = 0.05, fc = 1.7, 
              size = 1,
              palette = c("#B31B21", "#1465AC", "darkgray"),
              top = 0,              
              #top = 100,
              # label.select = c("NR4A1", "NR4A2", "NR4A3", "RHOB", #up
              #                  "HSPA1A", "FOS", "DUSP1", "JUN", "DNAJB1", "HSPA1B",
              #                  "FOSB", "TNFAIP3",
              #                  "MUCL1", "CYBA"), #down
              #label.select = c("NR4A1", "NR4A2", "NR4A3"),
              genenames = as.vector(data_input$gene),
              
              font.label = c("plain", 11),
              select.top.method = "fc",
              ggtheme = ggplot2::theme_classic()) +
  xlim(0,7) + ylim(-2.5,4) + NoLegend()

ggsave("c4_MAplot_72822.pdf", plot = p, width = 4.3, height = 4)


# MA plot NK-5 -------

Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub
# dfsample <- read.csv("NKclustDEGs_NOTHRESH.csv", header = T, row.names = 1)
dfsample <- read.csv("rNKvnon_DEGs_NOTHRESH_71422.csv", header = T, row.names = 1)
dfsample <- DEG_Remove_mito(dfsample)

avgexp <- AverageExpression(object = NK.all.combo, features = unique(dfsample$gene), assays = c("RNA", "SCT"))
avgexp <- as.data.frame(avgexp[['RNA']])
avgexp <- avgexp[which(rownames(avgexp) %in% unique(dfsample$gene)),]

head(avgexp)

c5DEGs <- dfsample[dfsample$cluster=="5",]
avgexpc5 <- avgexp[rownames(avgexp) %in% c5DEGs$gene,"5", drop = F]

data <- cbind(c5DEGs, avgexpc5)
data_input <- data[,c(9,2, 5, 7)]


head(data)
colnames(data)

which(data_input$gene == "KLRG1")
which(data_input$gene == "TIGIT")
which(data_input$gene == "NR4A3")

colnames(data_input) <- c("baseMean", "log2FoldChange", "padj", "gene")

library(ggpubr)
options(ggrepel.max.overlaps = 15)
#options(ggrepel.max.overlaps = 300)

p <- ggmaplot(data_input, 
              fdr = 0.05, fc = 1.5, 
              #fdr = 0.05, fc = 1.7, 
              size = 1,
              palette = c("#B31B21", "#1465AC", "darkgray"),
              top = 0,              
              #top = 100,
              # label.select = c("NR4A1", "NR4A2", "NR4A3", "RHOB", #up
              #                  "HSPA1A", "FOS", "DUSP1", "JUN", "DNAJB1", "HSPA1B",
              #                  "FOSB", "TNFAIP3",
              #                  "MUCL1", "CYBA"), #down
              #label.select = c("NR4A1", "NR4A2", "NR4A3"),
              genenames = as.vector(data_input$gene),
              
              font.label = c("plain", 11),
              select.top.method = "fc",
              ggtheme = ggplot2::theme_classic()) +
  xlim(0,7) + ylim(-2.5,4) + NoLegend()

ggsave("c5_MAplot_72822.pdf", plot = p, width = 4.3, height = 4)


# ================================================================== ======
#S3D NKsub GSEA  ========================
# ================================================================== ======



# GSEA LILY ===============================

library(AnnotationDbi)
library(msigdbr)
library(clusterProfiler)
library(dplyr)
# library(EnsDb.Hsapiens.v79)
library(org.Hs.eg.db)
# orgdb = "EnsDb.Hsapiens.v79"
org.db = "org.Hs.eg.db"

setwd(DEGdir)
dfsample <- read.csv("NKclustDEGs_NOTHRESH_71322.csv", header = T, row.names = 1)
dfsample$min.pct.diff = abs(dfsample$pct.1 - dfsample$pct.2)
dfsample <- DEG_Remove_mito(dfsample)
dfsample <- dfsample[which(dfsample$avg_log2FC > 0),]
#dfsample <- dfsample[which(abs(dfsample$min.pct.diff) >= 0.1),]

dfsample <- split(dfsample$gene, dfsample$cluster)

dfsample$`0` <- bitr(dfsample$`0`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`1` <- bitr(dfsample$`1`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`2` <- bitr(dfsample$`2`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`3` <- bitr(dfsample$`3`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`4` <- bitr(dfsample$`4`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`5` <- bitr(dfsample$`5`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)

genelist <- list("NK0" = dfsample$`0`$ENTREZID,
                 "NK1" = dfsample$`1`$ENTREZID,
                 "NK2" = dfsample$`2`$ENTREZID,
                 "NK3" = dfsample$`3`$ENTREZID,
                 "NK4" = dfsample$`4`$ENTREZID,
                 "NK5" = dfsample$`5`$ENTREZID)

m_df <- msigdbr(species = "Homo sapiens", category ="H") %>%#, subcategory = "CGP") %>%
  dplyr::select(gs_name, entrez_gene)

m_df <- m_df[-which(grepl("_UP", m_df$gs_name, fixed = TRUE)),]

GOclusterplot <- compareCluster(genelist,
                                fun = enricher,
                                TERM2GENE = m_df,
                                pvalueCutoff = 0.05#,
                                #qvalueCutoff=0.2,
                                # minGSSize = 10,
                                # maxGSSize = 500
)

p <- dotplot(GOclusterplot, includeAll = T)
p <- p + scale_y_discrete(label = function(x) {
  x %>% sub("HALLMARK_", "",.) %>%
    sub("BIOCARTA_", "",.) %>%
    sub("PID_", "",.) %>%
    sub("REACTOME_", "",.) %>%
    gsub("_", " ", .) %>%
    stringr::str_trunc(30, "right")
})

p <- p + scale_x_discrete(label = function(x) {
  x %>% stringr::str_sub(., 1,4)
})

ggsave("NKsubGSEA_orghsdb_72522.pdf", plot = p, width = 15, height = 10, units = "cm")
#ggsave("NKsubGSEA_71322.pdf", plot = p, width = 15, height = 10, units = "cm")

# ================================================================== ======
#S3F rNK genes per NKsub  ========================
# ================================================================== ======



# dotplot ==============================================
#https://github.com/satijalab/seurat/issues/3174
origup <- c("CALD1", "CLU", "ALOX12", "LTBP1", "CAVIN2",
            "PARVB","GP6", "SCD", "ITGAX", "NR4A3", "CCL4",
            "CR2", "HEATR9", "XDH", "RASGRP2", "MID1", "JUN",
            "CMKLR1", "DUSP1", "FOS", "ABCA1", "TNFAIP3",
            "NR4A1", "KLRG1", "DTX1", "NHSL2", "GFRA2",
            "FAM81A", "CX3CR1", "RHPN1", "HES1", "F5", 
            "GAS2L1", "THBS1", "MYLK", "TMTC1", "FOSB",
            "NR4A2", "MPIG6B", "SLC6A4", "PLXNA4", "VWF",
            "TUBB1", "SLC7A5")
origdown <- c("BCAT1", "ALDH1L2", "COX6A2", "PYCR1", #down
              "LHFPL2", "AHRR", "EXTL1", "ASNS", "CHAC1",
              "MTHFD2", "NEK6", "SLC6A9", "FMNL2", "ASB2",
              "SLC7A3", "AVIL", "CDH1", "CISH", "LGALS3",
              "GPT2", "CXCR6", "TRIB3", "CDKN1A", "ATF5",
              "SLC1A4", "PMEPA1", "CEMIP2", "OSGIN1",
              "ZNF503", "ITGA1", "ISG20", "PACSIN1",
              "TBC1D16", "RN7SL1", "SH3PXD2B", "SCN3B",
              "OSBPL1A", "ME1", "HPGDS", "PPP2R2C",
              "CLBA1", "HMOX1", "NQO1", "CARS1", "SSTR2",
              "SNORA23")
# NKup <- c("GNLY", #NK specific
#           "PRF1",
#           "KLRD1",
#           "NKG7")
# NKdown <- c("HSPB1",
#             "COX6C",
#             "ADIRF",
#             "CST3",
#             "EPCAM", #additions
#             "CD33",
#             "CD3D",
#             "MS4A1")
# Cancerup <- c("XCL1",
#               "MUCL1",
#               "CCL4L2",
#               "CCL3")
# Cancerdown <- c("CXCL13",
#                 "LAIR2",
#                 "CD8A",
#                 "MARCKSL1")
# 


#Idents(NK.all.combo) <- NK.all.combo$BC.Subtype
Idents(NK.all.combo) <- NK.all.combo$celltype_NKsub
features <- list("Original Mouse Reprog Up" = origup,
                 "Original Mouse Reprog Down" = origdown)#,

a <- DotPlot(object = NK.all.combo, features=features, #cluster.idents=T,
             dot.scale = 10) + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_colour_gradient2(low="steelblue", mid="lightgrey", high="red")
ggsave("NKdotplotsep.pdf", a, width = 20, height = 5)
ggsave("NKdotplotsep_updownonly.pdf", a, width = 11.8, height = 4.3)

pdf("ReprogDotplot_71922.pdf", width = 26.3, height = 4.9)

a+ theme(axis.line = element_line(colour = 'black', size = 1.5)) +
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(size = 18)) +
  theme(strip.text = element_text(size=20))
dev.off()
