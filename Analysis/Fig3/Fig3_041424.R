
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# This will generate all of the sub-figures for Figure 3, -------
# and is the analysis of cancer epithelial cell transcriptional
# heterogeneity defined by 10 GEs.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# SETUP: -----------------------------------------------------------------------
# Libraries -------------------------------------------------------------------

library(BiocManager)
library(GEOquery) 
library(plyr)
library(dplyr) 
library(Matrix)
library(Seurat)
library(ggplot2)
library(cowplot) 
library(multtest)
library(msigdbr)
library(fgsea)
library(monocle3)
library(velocyto.R)
library(loomR)
library(clustree)
library(tibble)
library(SeuratData)
library(matrixStats)
library(sparseMatrixStats)
library(DESeq2)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(viridis)
library(gridExtra)
library(ggplotify)
library(multtest)
library(metap)
library(writexl)
library(Rcpp)
library(RcppZiggurat)
library(Rfast)
library(ggh4x)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(cola)
library(msigdbr)
library(UCell)
library(RColorBrewer)

# Functions -------------------------------------------------------------------

# DEG_remove_mito 
DEG_Remove_mito <- function (df){ø
  df_rm_mito <- df[!grepl("^MT-|^MT.",df$gene),]
  return(df_rm_mito)
}

# DEG_remove_heat 
DEG_Remove_heat <- function (df){
  df_rm_hsp <- df[!grepl("^HSP",rownames(df)),]
  return(df_rm_hsp)
}

# plotSimilarityMatrix 
# create a similarity matrix plot from a dataframe
# adapted from klic package (adjusted heatmap parameters)
plotSimilarityMatrix = function(X, y = NULL, 
                                min.val = 0, 
                                max.val = 1,
                                clusLabels = NULL, 
                                colX = NULL, colY = NULL, 
                                clr = FALSE, clc = FALSE, 
                                annotation_col = NULL, 
                                annotation_row = NULL, 
                                annotation_colors = NULL, 
                                myLegend = NULL, 
                                fileName = "posteriorSimilarityMatrix", 
                                savePNG = FALSE, 
                                semiSupervised = FALSE, 
                                showObsNames = FALSE) {
  
  if (!is.null(y)) {
    # Check if the rownames correspond to the ones in the similarity matrix
    check <- sum(1 - rownames(X) %in% row.names(y))
    if (check == 1)
      stop("X and y must have the same row names.")
  }
  
  if (!is.null(clusLabels)) {
    if (!is.integer(clusLabels))
      stop("Cluster labels must be integers.")
    
    n_clusters <- length(table(clusLabels))
    riordina <- NULL
    for (i in 1:n_clusters) {
      riordina <- c(riordina, which(clusLabels == i))
    }
    
    X <- X[riordina, riordina]
    y <- y[riordina, ]
    y <- as.data.frame(y)
  }
  
  if (savePNG)
    grDevices::png(paste(fileName, ".png", sep = ""))
  
  if (!is.null(y)) {
    ht <- ComplexHeatmap::pheatmap(X, legend = TRUE,  
                                   color = rev(brewer.pal(11, "RdBu")), 
                                   breaks = seq(min.val, max.val, length.out = 11),
                                   cluster_rows = clr, 
                                   cluster_cols = clc, 
                                   #annotation_col = y,
                                   show_rownames = showObsNames, 
                                   show_colnames = showObsNames, 
                                   drop_levels = TRUE, 
                                   treeheight_row = -1,
                                   treeheight_col = -1,
                                   #annotation_row = annotation_row, 
                                   #annotation_col = annotation_col, 
                                   annotation_colors = annotation_colors)
  } else {
    ht <- ComplexHeatmap::pheatmap(X, legend = TRUE,
                                   color = rev(brewer.pal(11, "RdBu")),  
                                   breaks = seq(min.val, max.val, length.out = 11),
                                   cluster_rows = clr, 
                                   cluster_cols = clc,
                                   show_rownames = showObsNames, 
                                   show_colnames = showObsNames,
                                   treeheight_row = -1,
                                   treeheight_col = -1,
                                   drop_levels = TRUE, 
                                   #annotation_row = annotation_row, 
                                   #annotation_col = annotation_col, 
                                   annotation_colors = annotation_colors)
  }
  
  if (savePNG)
    grDevices::dev.off()
  
  return(ht)
}

# simil function 
# calculate similarity matrix 
# param df is dataframe of all cells for similarity matrix computation
# param drop is list of genes to drop from analysis
# param file is file name for output
# param method ("jaccard" or "corr") for method used
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

# simil_plot 
# plot interactive similarity heatmap 
# param a is an .rds file generated by simil_calc
# param annot is a vector used for heatmap annotation (metadata column from Seurat object)

simil_plot <- function(a, min.val, max.val, annot) {
  jc <- readRDS(a) #read in similarity matrix
  
  # add annotation
  if (length(annot) > 0) {
    row_annot <- annot[rownames(jc), , drop = FALSE] # select rows to use for heatmap annotation
    col_annot <- annot[colnames(jc), , drop = FALSE] # select columns to use for heatmap annotation
    colors <- mako(n_distinct(annot)) # create color vector to use for annotation
    names(colors) <- base::unique(annot)[[1]]
    colors <- list(colors, colors)
    names(colors) <- c(as.name(names(annot)), as.name(names(annot)))
  } 
  else {
    row_annot <- NULL
    col_annot <- NULL
    colors <- NULL
  }
  
  # generate simialrity matrix plot (using plotSimilarity Matrix)
  hm <- plotSimilarityMatrix(jc, clr = TRUE, clc = TRUE, 
                             min.val = min.val, max.val = max.val,
                             annotation_row = row_annot, 
                             annotation_col = col_annot, 
                             annotation_colors = colors,
                             showObsNames = T) # plot full matrix
  return(hm)
  #hm <- draw(hm)
  #htShiny(hm) # generate interactive heatmap 
}

# simil_sc50 
setwd("/work/InternalMedicine/s437775/simil")

# read in sc50 genes and convert to base::unique vector
sc50 <- read.csv("NatGen_Supplementary_table_S4.csv")
sc50 <- base::unique(unlist(sc50))
sc50 <- gsub("\\.", "-", sc50)
sc50 <- which(rownames(cancer.epi) %in% sc50) # get indices of sc50 genes

# calculate Jaccard similarity matrix 
# param df is dataframe for comparison
# param file is file name for output
# param method ("jaccard" or "corr") for similarity method used

simil_sc50 <- function(df, file, method) {
  jc <- df[sc50, , drop = FALSE] # subset to sc50 genes
  # jc[jc > 0.5] <- 1
  # jc[jc < 0.5 & jc > -0.5] <- 0
  # jc[jc < -0.5] <- -1
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

# simil_GE function 
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")

GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- which(rownames(cancer.epi) %in% GElist$gene) # get indices of GE genes

# calculate Jaccard similarity matrix 
# param df is dataframe for comparison
# param file is file name for output
# param method ("jaccard" or "corr") for similarity method used
simil_GE <- function(df, file, method) {
  jc <- df[GElist, , drop = FALSE] # subset to sc50 genes
  # jc[jc > 0.5] <- 1
  # jc[jc < 0.5 & jc > -0.5] <- 0
  # jc[jc < -0.5] <- -1
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


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ***Load in cancer epithelial cell object --------
# Load in Seurat object -------------------------------------------------------

# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/")
# combo.reference <- readRDS("PrimObject_FINALnoZgenedem_71222.rds") ##newnewnew

# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/")
# cancer.epi <- readRDS("Episubset_FINALnoZallgenedem_71322.rds") 
# cancer.epi <- NormalizeData(cancer.epi, assay = "RNA")

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
cancer.epi <- readRDS("cancerepi_withGEmetadata_111422.rds")
DefaultAssay(cancer.epi) <- "RNA"

# Scale data ------------------------------------------------------------------

DefaultAssay(cancer.epi) <- "RNA"
cancer.epi <- NormalizeData(cancer.epi, assay = "RNA")
cancer.epi <- FindVariableFeatures(cancer.epi, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)
all.genes <- rownames(cancer.epi)
cancer.epi <- ScaleData(cancer.epi, features = all.genes)

# Drop samples with too few cells ---------------------------------------------

# create dataframe of Seurat object metadata
sobjlists <- FetchData(object = cancer.epi, 
                       vars = c("samples",
                                "Patient", 
                                "BC.Subtype", 
                                "sc50.Pred"))

# determine number of cancer epi cells per sample
drop_samples <- sobjlists %>% dplyr::count(sobjlists$samples) # generate dataframe with counts of # cells per sample
drop_samples <- drop_samples[drop_samples$n < 50, ] # get list of samples with n < 10 cells
nrow(drop_samples)

# drop samples with too few cells 
sobjlists <- sobjlists[!(sobjlists$samples %in% drop_samples$`sobjlists$samples`), ] 

# create dataframe for bar graph analysis
sobjlists <- sobjlists %>% dplyr::group_by(samples, 
                                           Patient, 
                                           BC.Subtype, 
                                           sc50.Pred) %>% # group cells by samples and BC.Subtype
  dplyr::summarise(Nb = n()) %>% # add column with number of cells by sc50.Pred
  dplyr::mutate(C = sum(Nb)) %>% # add column with total number of cells in the sample
  dplyr::mutate(percent = Nb/C*100) # add column with % of cells by sc50.Pred out of all cells in sample

# reorder sobjlists
sobjlists <- sobjlists[order(sobjlists$BC.Subtype),] # order by BC Subtype

# get list of tumor samples with enough cells for analysis
samples <- unique(sobjlists$samples)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ***GE generation -------------------------------------------------
# All cancer epi clustering ---------------------------------------------------

# prepare cancer epithelial cell object
Epi.all.combo <- cancer.epi
DefaultAssay(Epi.all.combo) <- "RNA"
Epi.all.combo <- FindNeighbors(Epi.all.combo, reduction = "pca", dims = 1:50)

# cluster cells at various resolutions 
resolution.range <- c(0.01, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6, 1.8, 2.0)
Epi.all.combo <- FindClusters(Epi.all.combo, 
                              resolution = resolution.range)

# generate cancer epithelial cell UMAP
Epi.all.combo <- FindNeighbors(Epi.all.combo, reduction = "pca", dims = 1:50)
Epi.all.combo <- FindClusters(Epi.all.combo, resolution = resolution.range)
Epi.all.combo <- RunUMAP(Epi.all.combo, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use=123)

# plot UMAP by clustering
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE) + SeuratAxes() + NoLegend()
ggsave("cancerepi_UMAP.pdf", plot = as.ggplot(p), width = 5.8, height = 5.5)

# plot UMAP by original dataset
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "orig.ident") + SeuratAxes()
ggsave("cancerepi_UMAP_origdataset.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

# plot UMAP by clinical subtype
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "BC.Subtype") + SeuratAxes()
ggsave("cancerepi_UMAP_BCsubtype.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

# plot UMAP by sc50 subtype
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "sc50.Pred") + SeuratAxes()
ggsave("cancerepi_UMAP_sc50subtype.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

# update cancer epithelial cell object
cancer.epi <- Epi.all.combo

# list of clinically actionable targets
targets <- list("ESR1",
                "ERBB2", #HER2
                "PIK3CA",
                c("NTRK1", "NTRK2", "NTRK3"),
                "CD274", #PD-L1
                "ERBB3", #HER3
                "EGFR",
                c("FGFR1", "FGFR2", "FGFR3", "FGFR4"),
                "TACSTD2", #TROP2
                c("CDK4", "CDK6"), 
                "AR",
                "NECTIN2", 
                "LAG3")

# Add UCell score for clinically actionable targets
cancer.epi <- AddModuleScore_UCell(cancer.epi,
                                   features = targets,
                                   name = names(targets),
                                   assay = "RNA")

cancer.epi@meta.data <- cancer.epi@meta.data[,-c(99:4942,4944:15186)]

colnames(cancer.epi@meta.data)[102:114] <- c("ESR1",
                                             "ERBB2", 
                                             "PIK3CA",
                                             "NTRK", 
                                             "CD274", 
                                             "ERBB3", 
                                             "EGFR",
                                             "FGFR",
                                             "TACSTD2", 
                                             "CDK",
                                             "AR",
                                             "NECTIN2", 
                                             "LAG3")

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
saveRDS(cancer.epi, "cancerepi_withtargets_110922.rds")

# Unsupervised DGE generation -------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
j <- readRDS("cancerepi_withtargets_110922.rds")

# initialize objects for unsupervised DGE gene signature generation
all_DGE <- data.frame()

# perform DGE analysis at each of the various clustering resolutions
for (i in colnames(j@meta.data)[94:98]) {
  Idents(j) <- j@meta.data[,i]
  
  # DGE analysis (cluster biomarkers)
  j.markers_DGE <- FindAllMarkers(j, only.pos = T, 
                                  min.cells.group = 50, 
                                  min.diff.pct = 0.25, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST")
  j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
  j.markers_DGE$res <- strsplit(i, "res.")[[1]][2]
  j.markers_DGE$cluster_res <- paste0(j.markers_DGE$cluster, "_", j.markers_DGE$res)
  
  all_DGE <- rbind(all_DGE, j.markers_DGE)
}

# save full output of all unsupervised DGE signatures
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(all_DGE,
           path = "cancerepi_DGE_unsupervised_3.xlsx", 
           col_names = TRUE, 
           format_headers = TRUE)

# Unsupervised sample-level DGE generation -----------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
j <- readRDS("cancerepi_withtargets_110922.rds")
cancer.epi <- j

# initialize objects for unsupervised DGE gene signature generation
all_DGE <- data.frame()

for (i in samples[49:64]) {
  subset <- subset(j, subset = samples == i)
  DefaultAssay(subset) <- "RNA"
  subset <- NormalizeData(subset, assay = "RNA")
  subset <- FindVariableFeatures(subset, 
                                 selection.method = "vst", 
                                 nfeaøtures = 2000)
  all.genes <- rownames(subset)
  subset <- ScaleData(subset, features = all.genes)
  
  subset <- FindNeighbors(subset, reduction = "pca", dims = 1:30)
  resolution.range <- c(0.01, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6, 1.8, 2.0)
  subset <- FindClusters(subset, 
                         graph.name = "integrated_snn", 
                         resolution = resolution.range)
  subset <- RunUMAP(subset, reduction = "pca", dims = 1:30, verbose = TRUE, seed.use=123)
  
  # perform DGE analysis at each of the various clustering resolutions
  for (k in colnames(subset@meta.data)[(grep("integrated_snn", colnames(subset@meta.data)[85:115]) + 84)]) {
    Idents(subset) <- subset@meta.data[,k]
    
    # DGE analysis (cluster biomarkers)
    j.markers_DGE <- FindAllMarkers(subset, only.pos = T, 
                                    min.cells.group = 5, 
                                    min.diff.pct = 0.1, 
                                    #logfc.threshold = 0.25,
                                    test.use = "MAST")
    if(dim(j.markers_DGE)[1] > 0) {
      j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
      j.markers_DGE$res <- strsplit(k, "res.")[[1]][2]
      j.markers_DGE$cluster_res <- paste0(i, "_",j.markers_DGE$cluster, "_", j.markers_DGE$res)
    }
    
    all_DGE <- rbind(all_DGE, j.markers_DGE)
  }
}

# save full output of all unsupervised DGE signatures
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(all_DGE,
           path = "cancerepi_DGE_unsupervised_samplelevel_4.xlsx", 
           col_names = TRUE, 
           format_headers = TRUE)

# Supervised DGE generation based on targets ----------------------------------

# perform DGE analysis for each target 
all_DGE <- data.frame()
j <- cancer.epi 

for (i in c(102:108)) {
  # cells in top 25% by target expression are classified as "high"/"1" vs. bottom 75% as "low"/0
  cutoff <- quantile(j@meta.data[,i],0.9)
  j@meta.data[(cancer.epi@meta.data[,i] > cutoff),i] <- "high"
  j@meta.data[(cancer.epi@meta.data[,i] <= cutoff),i] <- "med"
  j@meta.data[(cancer.epi@meta.data[,i] <= 0),i] <- "low"
  
  Idents(j) <- j@meta.data[,i]
  
  # perform DGE analysis
  j.markers_DGE <- FindAllMarkers(j, 
                                  min.cells.group = 5, 
                                  min.pct = 0.2, 
                                  logfc.threshold = 0.1,
                                  test.use = "MAST")
  j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
  j.markers_DGE$res <- NA
  j.markers_DGE$cluster_res <- paste0(colnames(j@meta.data)[i], "_", j.markers_DGE$cluster)
  
  all_DGE <- rbind(all_DGE, j.markers_DGE)
  View(j.markers_DGE)
}

# save full output of all supervised DGE signatures from clinical targets
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(all_DGE,
           path = paste0("cancerepi_DGE_supervised_1.xlsx"), 
           col_names = TRUE, 
           format_headers = TRUE)

# Supervised DGE generation based on SC50 subtype -----------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
j <- readRDS("cancerepi_withtargets_110922.rds")
cancer.epi <- j

# set idents to SC50 subtype
Idents(j) <- j@meta.data$sc50.Pred

# perform DGE analysis for each SC50 subtype
j.markers_DGE <- FindAllMarkers(j, min.cells.group = 5, 
                                min.pct = 0.2, 
                                logfc.threshold = 0.1,  
                                test.use = "MAST")
j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
j.markers_DGE$res <- NA
j.markers_DGE$cluster_res <- j.markers_DGE$cluster

# save full output of all unsupervised DGE signatures from SC50 subtype
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(j.markers_DGE,
           path = "cancerepi_DGE_sc50.xlsx", 
           col_names = TRUE,
           format_headers = TRUE)

# Calculate Jaccard similarity ------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")

# create Jaccard similarity matrix for supervised and unsupervised DGE lists
files <- list.files(pattern = ".xlsx")

cancer_DGEs <- data.frame()
for (i in files) {
  cancer_DGEs <- rbind(cancer_DGEs, readxl::read_xlsx(i))
}

# filter DGE signatures
cancer_DGEs <- DEG_Remove_mito(cancer_DGEs) # remove mitochondrial genes
#cancer_DGEs <- DEG_Remove_heat(cancer_DGEs) # remove HSP gene
cancer_DGEs <- cancer_DGEs[which(cancer_DGEs$p_val_adj < 0.05),] # filter with p-value threshold
cancer_DGEs <- cancer_DGEs[which(cancer_DGEs$avg_log2FC > 0),] # filter with log2FC threshold
cancer_DGEs <- cancer_DGEs[which(cancer_DGEs$pct.1 > 0.2),] # filter with percent threshold

# select top n genes per signature by adjusted p-value (includes ties)
cancer_DGEs <- cancer_DGEs %>%
  group_by(cluster_res) %>%
  top_n(-200, p_val_adj)

# select top n genes per signature by avg_log2FC
cancer_DGEs <- cancer_DGEs %>% 
  group_by(cluster_res) %>%
  top_n(200, avg_log2FC) 

# get number of genes in each DGE signature
DGE_counts <- cancer_DGEs %>% dplyr::count(cancer_DGEs$cluster_res)
View(DGE_counts)

# prepare signature matrix to use for calculating pairwise Jaccard indices between signatures
cancer_DGEs_forjaccard <- as.data.frame(unique(cancer_DGEs$gene))
colnames(cancer_DGEs_forjaccard) <- c('gene')
rownames(cancer_DGEs_forjaccard) <- cancer_DGEs_forjaccard$gene
for (i in unique(cancer_DGEs$cluster_res)) {
  j <- as.data.frame(cancer_DGEs[which(cancer_DGEs$cluster_res == i), ]$gene)
  j$i <- 1
  colnames(j) <- c("gene", i)
  rownames(j) <- j$gene
  cancer_DGEs_forjaccard <- left_join(cancer_DGEs_forjaccard, j, by = "gene")
}

rownames(cancer_DGEs_forjaccard) <- cancer_DGEs_forjaccard$gene
cancer_DGEs_forjaccard <- cancer_DGEs_forjaccard[,-1]
cancer_DGEs_forjaccard[is.na(cancer_DGEs_forjaccard)] <- 0

# drop unsupervised DGE signatures with fewer than n genes
drop <- vector()
for (i in 1:(dim(cancer_DGEs_forjaccard)[2])) {
  if (DGE_counts[which(DGE_counts$cluster_res ==
                       colnames(cancer_DGEs_forjaccard)[i]),3] < 20) {
    drop <- append(drop, i)
  }
}
cancer_DGEs_forjaccard <- cancer_DGEs_forjaccard[, -drop]

# calculate Jaccard similarity between all DGE signatures
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
simil(cancer_DGEs_forjaccard, NULL, "cancerepi_DGEs_forjaccard_unsupervised_only.rds", "jaccard")
cancerepi_jaccard <- as.data.frame(readRDS("cancerepi_DGEs_forjaccard_unsupervised_only.rds"))

# remove redundant unsupervised DGE signatures
drop <- vector()
for (i in c(1:(length(colnames(cancerepi_jaccard))))) {
  for (j in c((i+1):(length(colnames(cancerepi_jaccard))))) {
    if (cancerepi_jaccard[i,j] > 0.9) {
      drop <- append(drop, j)
    }
  }
}
drop <- unique(drop)
cancerepi_jaccard_noredundant <- cancerepi_jaccard[-drop, -drop]
saveRDS(cancerepi_jaccard_noredundant, "cancerepi_DGEs_forjaccard_noredundant_unsupervised_only.rds")

# plot unsupervised DGE signature Jaccard similarities
p <- simil_plot("cancerepi_DGEs_forjaccard_noredundant_unsupervised_only.rds", 0, 1, NULL)
p <- as.ggplot(p)
ggsave("GE_prelim_unsupervised_only.pdf", plot = p, width = 40, height = 40)

# cola package to define GEs ----------------------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")

set.seed(123)

rh = consensus_partition(cancerepi_jaccard_noredundant, 
                         top_value_method = "ATC",
                         partition_method = "skmeans",
                         top_n = c(dim(cancerepi_jaccard_noredundant)[1]), 
                         p_sampling = 0.8,
                         max_k = 15)

k <- suggest_best_k(rh) #k <- 10

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
pdf("GE_stats_hclust.pdf", width = 7.3, height = 4)
select_partition_number(rh, mark_best = F)
dev.off()

write.csv(get_stats(rh), "GE_stats_hclust.csv")

pdf("GE_consensusplot.pdf", width = 7, height = 4)
consensus_heatmap(rh, k = k)
dev.off()

write.csv(get_classes(rh, k = k), "GE_classes_hclust.csv")

# Define GEs ------------------------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
GEs <- read.csv("GE_classes_hclust.csv", stringsAsFactors = FALSE, na.strings = "unknown")[,c(1,2)]
colnames(GEs) <- c("cluster_res", "GE")

cancer_DGEs_GE <- left_join(cancer_DGEs, GEs, by = "cluster_res")

cancer_DGEs_GE <- cancer_DGEs_GE[-which(is.na(cancer_DGEs_GE$GE)),] 
cancer_DGEs_GE$totavg_log2FC <- NA
cancer_DGEs_GE$totmin.pct.diff <- NA

for (i in unique(cancer_DGEs_GE$gene)) {
  mean <-  mean(cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$avg_log2FC)
  cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$totavg_log2FC <- mean
  
  mean <-  mean(cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$percent)
  cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$totmin.pct.diff <- mean
}

# Filter GEs ------------------------------------------------------------------

cancer_GE <- cancer_DGEs_GE %>% dplyr::group_by(GE, gene) %>%
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 

#cancer_GE <- cancer_GE[which(cancer_GE$percent > 0.1),]

cancer_GE <- cancer_GE %>%
  group_by(GE) %>%
  top_n(350, Nb)

cancer_GE <- left_join(cancer_GE, unique(cancer_DGEs_GE[,c(7,12)]))
cancer_GE <- cancer_GE %>% 
  group_by(GE) %>%
  top_n(200, totavg_log2FC) 

GE_summary <- cancer_GE %>% dplyr::group_by(GE) %>%
  dplyr::summarise(Nb = n())

View(GE_summary)

cancer_GE$old_GE <- cancer_GE$GE
cancer_GE$GE[which(cancer_GE$old_GE == 1)] <- 1
cancer_GE$GE[which(cancer_GE$old_GE == 2)] <- 3
cancer_GE$GE[which(cancer_GE$old_GE == 3)] <- 6
cancer_GE$GE[which(cancer_GE$old_GE == 4)] <- 4
cancer_GE$GE[which(cancer_GE$old_GE == 5)] <- 5
cancer_GE$GE[which(cancer_GE$old_GE == 6)] <- 9
cancer_GE$GE[which(cancer_GE$old_GE == 7)] <- 2
cancer_GE$GE[which(cancer_GE$old_GE == 8)] <- 7
cancer_GE$GE[which(cancer_GE$old_GE == 9)] <- 10
cancer_GE$GE[which(cancer_GE$old_GE == 10)] <- 8

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(cancer_GE, "cancer_GE.xlsx")

# Add GE scores for cancer cells ----------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- split(GElist$gene,GElist$GE)
length(GElist) # get number of GEs

cancer.epi <- AddModuleScore_UCell(cancer.epi,
                                   features = GElist,
                                   assay = "RNA")

# setwd("/work/InternalMedicine/s437775/simil")
# saveRDS(cancer.epi, file = "cancerepi_withGEs_061422.rds")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Fig 3A: GE heatmap ------

set.seed(123)
expdata <- t(cancer.epi@meta.data[,c(1,4,14,83:84,115:124)])
collapse_expdata <- as.data.frame(rownames(expdata))

for (i in samples) {
  subset <- expdata[,which(expdata[2,] == i)]
  subset <- subset[,sample.int(dim(subset)[2],min(dim(subset)[2],20000))]
  collapse_expdata <- cbind(collapse_expdata, subset)
}

labels <- rownames(collapse_expdata)[6:15]
labels <- labels %>% gsub("raw_", "",.)
labels <- labels %>% gsub("X", "GE",.)

collapse_expdata <- collapse_expdata[,-1]
collapse_zscore <- collapse_expdata[-c(1:5),]
collapse_zscore <- as.matrix(sapply(collapse_zscore, as.numeric))
collapse_zscore <- t(apply(collapse_zscore, 1, function(x) (x-mean(x))/sd(x)))
sort <- rbind(apply(collapse_zscore, 2, function(x) which.max(x)),
              apply(collapse_zscore, 2, function(x) max(x)), 
              collapse_zscore)
sort <- sort[,order(sort[2,],decreasing = T)] # sort by max z-score
sort <- sort[,order(sort[1,],decreasing = F)] #sort by GE
collapse_zscore <- sort[-c(1,2),]
collapse_expdata <- collapse_expdata[,colnames(collapse_zscore)]

orig_anno <- t(as.matrix(collapse_expdata[1,]))
colnames(orig_anno) <- c("origin")
sample_anno <- t(as.matrix(collapse_expdata[2,]))
colnames(sample_anno) <- c("sample")
BC_anno <- t(as.matrix(collapse_expdata[3,]))
colnames(BC_anno) <- c("BC subtype")
pam50_anno <- t(as.matrix(collapse_expdata[5,]))
colnames(pam50_anno) <- c("PAM50 subtype")
sc50_anno <- t(as.matrix(collapse_expdata[4,]))
colnames(sc50_anno) <- c("SC50 subtype") 
anno <- HeatmapAnnotation("orig" = orig_anno,
                          "sample" = sample_anno,
                          "BC" = BC_anno,
                          "pam50" = pam50_anno,
                          "sc50" = sc50_anno, 
                          show_legend = c("sample" = FALSE), 
                          col = list("orig" = c("Pal_Prim" = "#18A900",
                                                "Qian" = "#FFC300",
                                                "Wu" = "#C70039", 
                                                "Karaayvaz" = "#eb7d34",
                                                "Wu2021prim" = "#006CA9",
                                                "Xu" = "#A300DB"),
                                     "BC" = c("HER2+" = "#700639", "HR+" = "#397006", "TNBC" = "#063970"),
                                     "pam50" = c("Basal" = "#2596be", "Her2" = "#9925be", "LumA" = "#49be25", "LumB" = "#be4d25"),
                                     "sc50" = c("Basal_SC" = "#76b5c5", "Her2E_SC" = "#ad76c5", "LumA_SC" = "#c58676", "LumB_SC" = "#8dc576"))
)

p <- Heatmap(collapse_zscore,
             heatmap_legend_param = list("title" = "UCell\nZ-score"),
             cluster_rows = F,
             cluster_columns = F,
             show_column_dend = F,
             show_row_dend = F,
             show_column_names = F,
             show_row_names = T,
             row_labels = labels,
             use_raster = T,
             raster_quality = 4, 
             top_annotation = anno)

p <- as.ggplot(p)
ggsave("GE_heatmap.pdf", plot = p, width = 12, height = 5)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Fig 3B: GE barplot -------

expdata <- t(cancer.epi@meta.data[,c(1,4,14,83:84,115:124)])

zscore <- expdata[-c(1,2,3,4,5),]
zscore <- t(as.matrix(apply(zscore, 1, as.numeric)))
zscore <- t(apply(zscore, 1, function(x) (x-mean(x))/sd(x)))
sort <- rbind(apply(zscore, 2, function(x) which.max(x)),
              apply(zscore, 2, function(x) max(x)), 
              zscore, 
              expdata[c(1,2,3,4,5),])
sort <- as.data.frame(t(sort))
colnames(sort) <- c("GE", "maxZscore", "GE1", "GE2", "GE3", "GE4", "GE5", 
                    "GE6", "GE7", "GE8", "GE9", "GE10",
                    "origin", "sample", "BC_subtype", "pam50_subtype", "sc50_subtype")

sc50_plot <- sort %>% dplyr::group_by(GE, sc50_subtype) %>% # group cells by samples and BC.Subtype
  dplyr::summarise(Nb = n()) %>% # add column with number of cells by sc50.Pred
  dplyr::mutate(C = sum(Nb)) %>% # add column with total number of cells in the sample
  dplyr::mutate(percent = Nb/C) # add column with % of cells by sc50.Pred out of all cells in sample
sc50_plot$GE <- as.numeric(sc50_plot$GE)
sc50_plot <- sc50_plot[order(sc50_plot$GE,decreasing = F),] # sort by GE
sc50_plot$GE <- factor(sc50_plot$GE, levels = unique(sc50_plot$GE))

p <- ggplot(sc50_plot, aes(fill=sc50_subtype, y=percent, x=GE)) + 
  geom_bar(width = 0.7, stat="identity") + #position = "dodge"
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), labels = (scales::percent)) + 
  scale_fill_manual(values = c("#76b5c5", "#ad76c5", "#c58676", "#8dc576")) + 
  xlab("") + ylab("Proportion of Cells") + 
  scale_x_discrete(labels=c("1" = "GE 1", 
                            "2" = "GE 2",
                            "3" = "GE 3", 
                            "4" = "GE 4",
                            "5" = "GE 5", 
                            "6" = "GE 6",
                            "7" = "GE 7", 
                            "8" = "GE 8",
                            "9" = "GE 9", 
                            "10" = "GE 10")) 
ggsave("GE_bysubtype.pdf", plot = p, width = 7, height = 2)

stacked_plot <- sort %>% dplyr::group_by(sample, BC_subtype, GE) %>% # group cells by samples and BC.Subtype
  dplyr::summarise(Nb = n()) %>% # add column with number of cells by sc50.Pred
  dplyr::mutate(C = sum(Nb)) %>% # add column with total number of cells in the sample
  dplyr::mutate(percent = Nb/C) # add column with % of cells by sc50.Pred out of all cells in sample
stacked_plot$GE <- paste0("GE ", stacked_plot$GE)
stacked_plot <- stacked_plot[order(stacked_plot$BC_subtype,decreasing = F),] # sort by GE
stacked_plot$sample <- factor(stacked_plot$sample, levels = unique(stacked_plot$sample))
stacked_plot$GE <- factor(stacked_plot$GE, levels = c("GE 1", "GE 2", "GE 3", "GE 4", 
                                                      "GE 5", "GE 6", "GE 7", "GE 8", 
                                                      "GE 9", "GE 10"))

p <- ggplot(stacked_plot, aes(fill=GE, y=percent, x=sample)) + 
  geom_bar(width = 0.7, stat="identity") + 
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), labels = (scales::percent)) + 
  scale_fill_manual(values = c("#66C2A5","#ED936B","#919FC7","#DA8EC0",
                               "#B0D867","#F9DA56","#E0C59A","#B3B3B3",
                               "#90C786", "#DDDDDD")) +
  xlab("") + ylab("Proportion of Cells")+ 
  facet_nested( ~ BC_subtype, 
                scales = "free", 
                space = "free", 
                switch = "x") + 
  theme(strip.text.x = element_text(size = 14, face = "italic"), 
        strip.background = element_rect(#colour = "white", size = 1, 
          fill = "#EEEEEE"), 
        #panel.border = element_rect(color = "black", fill = NA, size = 1), 
        panel.spacing.x = unit(0.3, "lines")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("GE_bysample.pdf", plot = p, width = 20, height = 3)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Fig 3C: Hallmark gene set analysis for GEs ---------------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
dfsample <- readxl::read_xlsx("cancer_GE.xlsx")
dfsample <- split(dfsample$gene, dfsample$GE)

library(EnsDb.Hsapiens.v79)
orgdb = "org.Hs.eg.db"

dfsample$`1` <- bitr(dfsample$`1`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)#, drop = FALSE)
dfsample$`2` <- bitr(dfsample$`2`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)#, drop = FALSE)
dfsample$`3` <- bitr(dfsample$`3`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`4` <- bitr(dfsample$`4`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`5` <- bitr(dfsample$`5`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`6` <- bitr(dfsample$`6`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`7` <- bitr(dfsample$`7`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`8` <- bitr(dfsample$`8`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`9` <- bitr(dfsample$`9`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`10` <- bitr(dfsample$`10`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)

genelist <- list("1" = dfsample$`1`$ENTREZID,
                 "2" = dfsample$`2`$ENTREZID,
                 "3" = dfsample$`3`$ENTREZID,
                 "4" = dfsample$`4`$ENTREZID,
                 "5" = dfsample$`5`$ENTREZID,
                 "6" = dfsample$`6`$ENTREZID,
                 "7" = dfsample$`7`$ENTREZID,
                 "8" = dfsample$`8`$ENTREZID, 
                 "9" = dfsample$`9`$ENTREZID,
                 "10" = dfsample$`10`$ENTREZID)

m_df <- msigdbr(species = "Homo sapiens", category ="H") %>%#, subcategory = "CGP") %>% 
  dplyr::select(gs_name, entrez_gene)

GOclusterplot <- compareCluster(genelist, 
                                fun = enricher, 
                                TERM2GENE = m_df, 
                                pvalueCutoff = .05)

p <- dotplot(GOclusterplot, includeAll = F)
p <- p + scale_y_discrete(label = function(x) {
  x %>% sub("HALLMARK_", "",.) %>% 
    sub("BIOCARTA_", "",.) %>% 
    sub("PID_", "",.) %>% 
    sub("REACTOME_", "",.) %>% 
    gsub("_", " ", .) %>% 
    stringr::str_trunc(40, "right")
}) 
p <- p + scale_x_discrete(label = function(x) {
  x %>% paste0("GE",.) %>% 
    stringr::str_sub(., 1,4)
}) 
p <- p + theme(axis.text.x = element_text(angle = 90, hjust=1))

ggsave("GE_clusterplot_hallmark.pdf", plot = p, width = 20, height = 15, units = "cm")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# *** Nichenet analysis 
# ***Nichenet analysis -------------------------------
# Create epi Nichenet object ------------

# Create epithelial cell object with GE idents
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
cancer.epi <- readRDS("cancerepi_withGEmetadata_111422.rds")
DefaultAssay(cancer.epi) <- "RNA"

GEidents <- data.frame()
for (i in c(1:10)) {
  j <- cancer.epi$maxGE
  j[which(j == i)] <- "high"
  j[which(j != "high")] <- "low"
  j <- as.data.frame(t(j))
  rownames(j) <- c(paste0("GE", i, "_idents"))
  GEidents <- rbind(GEidents, j)
}

colnames(GEidents) <- colnames(cancer.epi)
cancer.epi <- AddMetaData(cancer.epi, t(GEidents), col.name = rownames(GEidents))

# Create non-epithelial cell object
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
combo.reference <- readRDS("primobj_withCNVlabels_111022.rds")
combo.reference$celltype_withreprog <- as.character(combo.reference$celltype_withreprog)
combo.reference$celltype_withreprog[which(combo.reference$celltype_final == "Cancer Epithelial Cells")] <- "Cancer Epithelial Cells"
DefaultAssay(combo.reference) <- "RNA"

cancer.other <- subset(combo.reference, subset = celltype_withreprog != "Cancer Epithelial Cells")
cancer.other$celltype_withreprog <- Idents(cancer.other)

GEidents <- data.frame()
for (i in c(1:10)) {
  j <- as.data.frame(t(as.data.frame(cancer.other$celltype_withreprog)))
  rownames(j) <- c(paste0("GE", i, "_idents"))
  GEidents <- rbind(GEidents, j)
}

colnames(GEidents) <- colnames(cancer.other)
cancer.other <- AddMetaData(cancer.other, t(GEidents), col.name = rownames(GEidents))

nichenet_obj <- merge(cancer.epi, cancer.other)
nichenet_obj$celltype_withreprog[which(is.na(nichenet_obj$celltype_withreprog))] <- "Cancer Epithelial Cells"
nichenet_obj$celltype_final[which(nichenet_obj$celltype_withreprog == "Cancer Epithelial Cells")] <- "Cancer Epithelial Cells"

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
saveRDS(nichenet_obj, "nichenetobj_112022.rds")

# Load in NicheNet object ---------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/")
nichenet_obj <- readRDS("nichenetobj_112022.rds")

all_cell <- unique(nichenet_obj$celltype_withreprog)[-1]

# Run NicheNet -----------

setwd("/endosome/work/InternalMedicine/s437775/simil")

ligand_target_matrix = readRDS("ligand_target_matrix.rds")
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")

for (cell in all_cell[16:18]) {
  # NicheNet GE1 -------------
  
  Idents(nichenet_obj) <- nichenet_obj$GE1_idents
  temp_nichenet1 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet1, file = paste0(cell, "_nichenet1.Rdata"))
  
  # NicheNet GE2 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE2_idents
  temp_nichenet2 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet2, file = paste0(cell, "_nichenet2.Rdata"))
  
  # NicheNet GE3 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE3_idents
  temp_nichenet3 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet3, file = paste0(cell, "_nichenet3.Rdata"))
  
  # NicheNet GE4 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE4_idents
  temp_nichenet4 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet4, file = paste0(cell, "_nichenet4.Rdata"))
  
  # NicheNet GE5 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE5_idents
  temp_nichenet5 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet5, file = paste0(cell, "_nichenet5.Rdata"))
  
  # NicheNet GE6 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE6_idents
  temp_nichenet6 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet6, file = paste0(cell, "_nichenet6.Rdata"))
  
  # NicheNet GE7 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE7_idents
  temp_nichenet7 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet7, file = paste0(cell, "_nichenet7.Rdata"))
  
  # NicheNet GE8 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE8_idents
  temp_nichenet8 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet8, file = paste0(cell, "_nichenet8.Rdata"))
  
  # NicheNet GE9 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE9_idents
  temp_nichenet9 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet9, file = paste0(cell, "_nichenet9.Rdata"))
  
  # NicheNet GE10 -------------
  
  Idents(nichenet_obj) <- nichenet_obj$GE10_idents
  temp_nichenet10 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet10, file = paste0(cell, "_nichenet10.Rdata"))
  
}

# Read in Nichenet results --------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")
files <- list.files(pattern = ".Rdata")
#files <- files[-c(51, 172)]

# Nichenet object subsets ----------

all_cell <- unique(nichenet_obj$celltype_withreprog)[-1]

primsubset <- subset(nichenet_obj, subset = celltype_withreprog %in% all_cell)
Idents(primsubset) <- primsubset$celltype_withreprog

episubset <- subset(nichenet_obj, subset = celltype_final == "Cancer Epithelial Cells") 

# NicheNet circos plots -----------

all_receptors <- data.frame()
all_targets <- data.frame()

for (j in files) {
  temp_nichenet <- readRDS(j)
  cell <- strsplit(j, "_")[[1]][1]
  GE <- as.numeric(substring(strsplit(strsplit(j, "_")[[1]][2], ".R")[[1]],9)[1])
  
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
  Idents(episubset) <- episubset@meta.data[,136 + GE]
  keep_lig <- AverageExpression(subset(episubset, idents = c("high", "low")), features = gsub("\\.", "-", unique(LR_matrix$ligand)), slot = "data")$RNA
  keep_lig <- as.data.frame(keep_lig[which(keep_lig[,'high'] > keep_lig[,'low']),])
  keep_lig$ligand <- rownames(keep_lig)
  keep_lig$ligexp <- keep_lig$high / keep_lig$low
  keep_lig <- keep_lig[,-c(1:2)]
  LR_matrix <- LR_matrix[which(LR_matrix$ligand %in% keep_lig$ligand),]
  
  cellsubset <- subset(primsubset, subset = celltype_withreprog == cell)
  cellsubset <- GetAssayData(cellsubset, assay = "RNA", slot = "data")
  cellsubset <- as.data.frame(cellsubset[which(rownames(cellsubset) %in% gsub("\\.", "-", unique(LR_matrix$receptor))),])
  cellsubset[cellsubset > 0] <- 1
  cellsubset$percent <- apply(cellsubset, 1, function(x) sum(x))/dim(cellsubset)[2]
  cellsubset$receptor <- rownames(cellsubset)
  keep_rec <- as.data.frame(cbind(cellsubset$receptor, cellsubset$percent))
  colnames(keep_rec) <- c("receptor", "recpercent")
  LR_matrix <- LR_matrix[which(LR_matrix$receptor %in% keep_rec$receptor),]
  
  LR_matrix <- left_join(LR_matrix, keep_lig)
  LR_matrix <- left_join(LR_matrix, keep_rec)
  LR_matrix$LRexp <- LR_matrix$ligexp * as.numeric(LR_matrix$recpercent)
  dim(LR_matrix)
  
  lr.df.top <- LR_matrix %>% top_n(50, LRexp)
  
  if (cell %in% c("NK Cells", "Reprogrammed NK Cells")) {
    activating <- LR_matrix[which(LR_matrix$receptor %in% c("CD160", "CD226", "CD244", "CRTAM",
                                                            "FCGR3A", "KLRC2", "KLRK1", "NCR1", "NCR2", "NCR3",
                                                            "TNFRSF9")),]
    inactivating <- LR_matrix[which(LR_matrix$receptor %in% c("CD96", "KIR2DL1", "KIR2DL2", "KIR2DL3","KIR2DL4", "KLRA",
                                                              "KLRC1", "LAG3", "LILRB1", "PDCD1", "PVRIG", "TIGIT")),]
    lr.df.top <- unique(rbind(lr.df.top, activating, inactivating))
  }
  
  #lr.df.top$weight <- lr.df.top$LRexp
  lr.df.top <- lr.df.top[,c(1:3)]
  
  setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet_NK")
  
  #generate circos plot with top 10 ligands by ligexp
  pdf(paste0("NK", GE, "_", cell, "_circos.pdf"), width = 8, height = 8)
  
  circos.par(canvas.ylim=c(-1.5,1.5), # edit  canvas size
             track.margin = c(0.01, 0)) # adjust bottom and top margin
  
  chordDiagram(lr.df.top,
               directional = 1,
               link.sort = TRUE,
               link.decreasing = FALSE,
               link.visible = lr.df.top$weight > 0,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.05))
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
  }, bg.border = NA)
  circos.clear()
  
  dev.off()
  
  targets <- temp_nichenet$ligand_target_df[which(temp_nichenet$ligand_target_df$ligand %in% LR_matrix$ligand),]
  receptors <- LR_matrix
  targets <- as.data.frame(targets)
  targets$GE <- paste0("NK",GE)
  targets$cell <- cell
  all_targets <- rbind(all_targets, targets)
  
  receptors <- as.data.frame(receptors)
  receptors$GE <- paste0("NK",GE)
  receptors$cell <- cell
  all_receptors <- rbind(all_receptors, receptors)
  
}

write.csv(all_receptors, "nichenet_receptors.csv")
write.csv(all_targets, "nichenet_targets.csv")


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ***CellChat analysis --------------
# Libraries ------------

#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork) 
options(stringsAsFactors = FALSE)
library(NMF)
library(ggalluvial)
library(purrr)
library(ComplexHeatmap)
library(viridis)
library(limma)
library(tidyr)

# Load in Nichenet object ------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
nichenet_obj <- readRDS("nichenetobj_112022.rds")
DefaultAssay(nichenet_obj) <- "RNA"
nichenet_obj <- NormalizeData(nichenet_obj, verbose = FALSE)

nichenet_obj$maxGE <- paste0("GE", nichenet_obj$maxGE)
nichenet_obj$maxGE[which(nichenet_obj$maxGE == "GENA")] <- nichenet_obj$celltype_withreprog[which(nichenet_obj$maxGE == "GENA")]
data.input <- GetAssayData(nichenet_obj, assay = "RNA", slot = "data") # normalized data matrix

# Run CellChat -----------

cellchats <- list()
nets <- list()

for (i in c(82)) {
  labels <- nichenet_obj@meta.data[,i]
  meta <- data.frame(group = labels, row.names = rownames(nichenet_obj@meta.data))
  
  # Create CellChat object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

  # Load CellChat database
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  #showDatabaseCategory(CellChatDB)
  
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  cellchat@DB <- CellChatDB.use
  
  # Pre-process data
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  
  cellchats <- append(cellchats, cellchat)
  
  net <- subsetCommunication(cellchat, 
                             slot.name = "net", 
                             thresh = 0.05)
  nets <- append(nets, net)
}

saveRDS(cellchats, "cellchat_obj_113022.rds")
saveRDS(nets, "cellchat_nets_113022.rds")

# CellChat output ----------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
cellchats <- readRDS("cellchat_obj_113022.rds")

net <- data.frame(stringsAsFactors = F)

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet_NK")
all_receptors <- read.csv("nichenet_receptors.csv", row.names = 1, stringsAsFactors = F)[,-1]

net <- subsetCommunication(cellchats[[1]],
                           slot.name = "net",
                           thresh = 0.05,
                           targets.use = unique(all_receptors$cell))

write.csv(net, "cellchat_receptors.csv")

# -----
# -----
# Fig 3D: GE-immune decoder matrix from combined NicheNet and CellChat results  -------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")
all_receptors <- read.csv("nichenet_receptors.csv", row.names = 1, stringsAsFactors = F)[,-1]
net <- read.csv("cellchat_receptors.csv", row.names = 1, stringsAsFactors = F)

nichenet_receptors <- all_receptors
newgenes <- alias2SymbolTable(nichenet_receptors$receptor, species = "Hs")
newgenes <- ifelse(is.na(newgenes),nichenet_receptors$receptor,newgenes)
nichenet_receptors$receptor <- newgenes

nichenet_receptors <- nichenet_receptors  %>%
  group_by(cell) %>%
  top_n(2300, weight)

nichenet_receptors <- nichenet_receptors  %>%
  group_by(GE) %>%
  top_n(500, LRexp)

write.csv(nichenet_receptors, "cellchat_nichenet_receptors_list.csv")

nichenet_receptors <- nichenet_receptors[, c(2,6,7)]
colnames(nichenet_receptors) <- c("receptor", "GE", "cell")

cellchat_receptors <- net[which(net$prob >= 0.0),]
cellchat_receptors <- cellchat_receptors[,c(4,1,2)]
colnames(cellchat_receptors) <- c("receptor", "GE", "cell")
cellchat_receptors = separate_rows(cellchat_receptors, 1, sep = "_")
cellchat_receptors = separate_rows(cellchat_receptors, 1, sep = ":")
newgenes <- alias2SymbolTable(cellchat_receptors$receptor, species = "Hs")
newgenes <- ifelse(is.na(newgenes),cellchat_receptors$receptor,newgenes)
cellchat_receptors$receptor <- newgenes

overlap <- cellchat_receptors[which(do.call(paste0, cellchat_receptors) %in% do.call(paste0, nichenet_receptors)),]
overlap <- rbind(overlap,nichenet_receptors)
overlap <- unique(overlap)

overlap <- overlap %>% dplyr::group_by(GE, cell) %>%
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 
overlap <- split(overlap[,2:3], overlap$GE)
names <- names(overlap)
overlap <- purrr::reduce(overlap, full_join, by = "cell") %>% replace(., is.na(.), 0);
overlap <- as.data.frame(overlap)
colnames(overlap) <- c("cell", names)
rownames(overlap) <- overlap$cell
overlap <- overlap[,c(7:17)]
overlap <- overlap[,-1]

overlap <- overlap[-which(rownames(overlap) == "Epithelial Cells"),]
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")

pdf("cellchat_nichenet_heatmap.pdf")
Heatmap(overlap, 
        col = circlize::colorRamp2(seq(min(overlap),max(overlap), 
                                       length.out = 500), maGEa(500)), 
        cluster_rows = T)
dev.off()

pdf("cellchat_nichenet_heatmap_scaled.pdf")
Heatmap(overlap / rowSums(overlap), 
        col = circlize::colorRamp2(seq(min(overlap / rowSums(overlap)),
                                       #max(overlap / rowSums(overlap)), 
                                       0.2,
                                       length.out = 500), maGEa(500)), 
        cluster_rows = T)
dev.off()

write.csv(overlap, "cellchat_nichenet_receptors_sum.csv")

# -----
# -----

# Fig 3E: Cell line GE analysis -----------

# Load in files -----------------------

# from CCLE DepMap portal 
fullmat <- read.csv("CCLE_RNAseq_reads.csv")
rownames(fullmat) <- fullmat$X
fullmat <- fullmat[,-1]
colnames(fullmat) <- sub('\\.\\..*', '', colnames(fullmat))
fullmat <- t(fullmat)
fullmat[1:5,1:5]

sample.info <- read.csv("sample_info.csv")
rownames(sample.info) <- sample.info$DepMap_ID
head(sample.info)
cell.lines.wanted <- colnames(fullmat)
sample.info <- sample.info[rownames(sample.info) %in% cell.lines.wanted, ]
dim(sample.info)

cell_lines <- CreateSeuratObject(fullmat, assay = "RNA", meta.data = sample.info)
DefaultAssay(cell_lines) <- "RNA"
cell_lines <- NormalizeData(cell_lines)

# Subset cell lines and add GE module scores --------------- 

setwd("/endosome/work/InternalMedicine/s437775/simil")

# Load in supplemental data from Sheffer et al. for breast cancer cell lines 
paper_lines <- read.csv("cell_line_screens.csv", stringsAsFactors = F)

colnames(paper_lines) <- c("Cell.Lines", "Tissue",
                           "24hr_AUC", "48hr_AUC", "72hr_AUC",
                           "24hr_sensitivity", "48hr_sensitivity", "72hr_sensitivity",
                           "Mesenchymal", "Epithelial", "MSI_high",
                           "B7H6_protein_scores", "B7H6_high",
                           "HLA_protein_scores", "HLA_negative")

breast_lines <- subset(cell_lines, 
                       subset = stripped_cell_line_name %in% paper_lines$Cell.Lines)

newnames <- breast_lines@meta.data[,c(1,6)]
newnames$orig.ident <- rownames(newnames)
colnames(newnames) <- c("orig_names", "Cell.Lines")

paper_lines <- left_join(paper_lines,newnames)
paper_lines <- paper_lines[!is.na(paper_lines$orig_names),]
rownames(paper_lines) <- paper_lines$orig_names

breast_lines <- AddMetaData(breast_lines, paper_lines[,-16])

colnames(breast_lines@meta.data)

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- split(GElist$gene,GElist$GE)
length(GElist) # get number of GEs

breast_lines <- AddModuleScore_UCell(breast_lines,
                                     features = GElist,
                                     assay = "RNA")

# Z-score GE scores 
breast_lines@meta.data[41:50] <- apply(breast_lines@meta.data[41:50], 2, function(x) (x-mean(x))/sd(x))

# GE vs 24 hr sensitivity ----------------

breast_lines$X72hr_killing <- 1 - breast_lines$X72hr_AUC
breast_lines$X48hr_killing <- 1 - breast_lines$X48hr_AUC
breast_lines$X24hr_killing <- 1 - breast_lines$X24hr_AUC

p <- list() 
p.val <- vector()
corr <- vector()

breast_lines$include <- 0
breast_lines$include[which(breast_lines$X24hr_sensitivity < breast_lines$X72hr_sensitivity)] <- 1
breast_lines$include[which(breast_lines$X24hr_sensitivity < breast_lines$X48hr_sensitivity)] <- 1

for (i in c(1:10)) {
  j <- colnames(breast_lines@meta.data)[40 + i]
  subset <- breast_lines@meta.data[which(breast_lines$include == 0),]
  p[[i]] <- ggscatter(subset, x = "X24hr_killing", y = j,
                      add = "reg.line",
                      conf.int = TRUE,
                      add.params = list(color = "blue", fill = "lightgray"),
                      xlab = "Sensitivity to NK Cell Killing", 
                      ylab = paste0("GE", i, " Expression")) + 
    stat_cor(method = "spearman", 
             label.sep = "\n", label.y = 0, 
             label.x = 0.5) + 
    xlim(min(subset$X24hr_killing), max(subset$X24hr_killing)) 
  
  p.val <- append(p.val, cor.test(subset$X24hr_killing, subset[,40 + i], method = "spearman")$p.value)
  corr <- append(corr, cor.test(subset$X24hr_killing, subset[,40 + i], method = "spearman")$estimate)
}

p.val.adj <- p.adjust(p.val, method = "BH")
p.val.adj
p.val
dim(subset)

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
pdf("cell_lines_24hrAUCvsUCell_112523.pdf", width = 5, height = 5)
print(p)
dev.off()


# -----
# -----

# Fig 3F: Nichenet circos plots for NK cells vs. GE1 and GE6 -----------

# Refer to Nichenet analysis section above for code to generate circos plots
                                       # -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# This will generate all of the sub-figures for Figure 3, -------
# and is the analysis of cancer epithelial cell transcriptional
# heterogeneity defined by 10 GEs.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# SETUP: -----------------------------------------------------------------------
# Libraries -------------------------------------------------------------------

library(BiocManager)
library(GEOquery) 
library(plyr)
library(dplyr) 
library(Matrix)
library(Seurat)
library(ggplot2)
library(cowplot) 
library(multtest)
library(msigdbr)
library(fgsea)
library(monocle3)
library(velocyto.R)
library(loomR)
library(clustree)
library(tibble)
library(SeuratData)
library(matrixStats)
library(sparseMatrixStats)
library(DESeq2)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(viridis)
library(gridExtra)
library(ggplotify)
library(multtest)
library(metap)
library(writexl)
library(Rcpp)
library(RcppZiggurat)
library(Rfast)
library(ggh4x)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(cola)
library(msigdbr)
library(UCell)
library(RColorBrewer)

# Functions -------------------------------------------------------------------

# DEG_remove_mito 
DEG_Remove_mito <- function (df){ø
  df_rm_mito <- df[!grepl("^MT-|^MT.",df$gene),]
  return(df_rm_mito)
}

# DEG_remove_heat 
DEG_Remove_heat <- function (df){
  df_rm_hsp <- df[!grepl("^HSP",rownames(df)),]
  return(df_rm_hsp)
}

# plotSimilarityMatrix 
# create a similarity matrix plot from a dataframe
# adapted from klic package (adjusted heatmap parameters)
plotSimilarityMatrix = function(X, y = NULL, 
                                min.val = 0, 
                                max.val = 1,
                                clusLabels = NULL, 
                                colX = NULL, colY = NULL, 
                                clr = FALSE, clc = FALSE, 
                                annotation_col = NULL, 
                                annotation_row = NULL, 
                                annotation_colors = NULL, 
                                myLegend = NULL, 
                                fileName = "posteriorSimilarityMatrix", 
                                savePNG = FALSE, 
                                semiSupervised = FALSE, 
                                showObsNames = FALSE) {
  
  if (!is.null(y)) {
    # Check if the rownames correspond to the ones in the similarity matrix
    check <- sum(1 - rownames(X) %in% row.names(y))
    if (check == 1)
      stop("X and y must have the same row names.")
  }
  
  if (!is.null(clusLabels)) {
    if (!is.integer(clusLabels))
      stop("Cluster labels must be integers.")
    
    n_clusters <- length(table(clusLabels))
    riordina <- NULL
    for (i in 1:n_clusters) {
      riordina <- c(riordina, which(clusLabels == i))
    }
    
    X <- X[riordina, riordina]
    y <- y[riordina, ]
    y <- as.data.frame(y)
  }
  
  if (savePNG)
    grDevices::png(paste(fileName, ".png", sep = ""))
  
  if (!is.null(y)) {
    ht <- ComplexHeatmap::pheatmap(X, legend = TRUE,  
                                   color = rev(brewer.pal(11, "RdBu")), 
                                   breaks = seq(min.val, max.val, length.out = 11),
                                   cluster_rows = clr, 
                                   cluster_cols = clc, 
                                   #annotation_col = y,
                                   show_rownames = showObsNames, 
                                   show_colnames = showObsNames, 
                                   drop_levels = TRUE, 
                                   treeheight_row = -1,
                                   treeheight_col = -1,
                                   #annotation_row = annotation_row, 
                                   #annotation_col = annotation_col, 
                                   annotation_colors = annotation_colors)
  } else {
    ht <- ComplexHeatmap::pheatmap(X, legend = TRUE,
                                   color = rev(brewer.pal(11, "RdBu")),  
                                   breaks = seq(min.val, max.val, length.out = 11),
                                   cluster_rows = clr, 
                                   cluster_cols = clc,
                                   show_rownames = showObsNames, 
                                   show_colnames = showObsNames,
                                   treeheight_row = -1,
                                   treeheight_col = -1,
                                   drop_levels = TRUE, 
                                   #annotation_row = annotation_row, 
                                   #annotation_col = annotation_col, 
                                   annotation_colors = annotation_colors)
  }
  
  if (savePNG)
    grDevices::dev.off()
  
  return(ht)
}

# simil function 
# calculate similarity matrix 
# param df is dataframe of all cells for similarity matrix computation
# param drop is list of genes to drop from analysis
# param file is file name for output
# param method ("jaccard" or "corr") for method used
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

# simil_plot 
# plot interactive similarity heatmap 
# param a is an .rds file generated by simil_calc
# param annot is a vector used for heatmap annotation (metadata column from Seurat object)

simil_plot <- function(a, min.val, max.val, annot) {
  jc <- readRDS(a) #read in similarity matrix
  
  # add annotation
  if (length(annot) > 0) {
    row_annot <- annot[rownames(jc), , drop = FALSE] # select rows to use for heatmap annotation
    col_annot <- annot[colnames(jc), , drop = FALSE] # select columns to use for heatmap annotation
    colors <- mako(n_distinct(annot)) # create color vector to use for annotation
    names(colors) <- base::unique(annot)[[1]]
    colors <- list(colors, colors)
    names(colors) <- c(as.name(names(annot)), as.name(names(annot)))
  } 
  else {
    row_annot <- NULL
    col_annot <- NULL
    colors <- NULL
  }
  
  # generate simialrity matrix plot (using plotSimilarity Matrix)
  hm <- plotSimilarityMatrix(jc, clr = TRUE, clc = TRUE, 
                             min.val = min.val, max.val = max.val,
                             annotation_row = row_annot, 
                             annotation_col = col_annot, 
                             annotation_colors = colors,
                             showObsNames = T) # plot full matrix
  return(hm)
  #hm <- draw(hm)
  #htShiny(hm) # generate interactive heatmap 
}

# simil_sc50 
setwd("/work/InternalMedicine/s437775/simil")

# read in sc50 genes and convert to base::unique vector
sc50 <- read.csv("NatGen_Supplementary_table_S4.csv")
sc50 <- base::unique(unlist(sc50))
sc50 <- gsub("\\.", "-", sc50)
sc50 <- which(rownames(cancer.epi) %in% sc50) # get indices of sc50 genes

# calculate Jaccard similarity matrix 
# param df is dataframe for comparison
# param file is file name for output
# param method ("jaccard" or "corr") for similarity method used

simil_sc50 <- function(df, file, method) {
  jc <- df[sc50, , drop = FALSE] # subset to sc50 genes
  # jc[jc > 0.5] <- 1
  # jc[jc < 0.5 & jc > -0.5] <- 0
  # jc[jc < -0.5] <- -1
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

# simil_GE function 
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")

GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- which(rownames(cancer.epi) %in% GElist$gene) # get indices of GE genes

# calculate Jaccard similarity matrix 
# param df is dataframe for comparison
# param file is file name for output
# param method ("jaccard" or "corr") for similarity method used
simil_GE <- function(df, file, method) {
  jc <- df[GElist, , drop = FALSE] # subset to sc50 genes
  # jc[jc > 0.5] <- 1
  # jc[jc < 0.5 & jc > -0.5] <- 0
  # jc[jc < -0.5] <- -1
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


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ***Load in cancer epithelial cell object --------
# Load in Seurat object -------------------------------------------------------

# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/")
# combo.reference <- readRDS("PrimObject_FINALnoZgenedem_71222.rds") ##newnewnew

# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/")
# cancer.epi <- readRDS("Episubset_FINALnoZallgenedem_71322.rds") 
# cancer.epi <- NormalizeData(cancer.epi, assay = "RNA")

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
cancer.epi <- readRDS("cancerepi_withGEmetadata_111422.rds")
DefaultAssay(cancer.epi) <- "RNA"

# Scale data ------------------------------------------------------------------

DefaultAssay(cancer.epi) <- "RNA"
cancer.epi <- NormalizeData(cancer.epi, assay = "RNA")
cancer.epi <- FindVariableFeatures(cancer.epi, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)
all.genes <- rownames(cancer.epi)
cancer.epi <- ScaleData(cancer.epi, features = all.genes)

# Drop samples with too few cells ---------------------------------------------

# create dataframe of Seurat object metadata
sobjlists <- FetchData(object = cancer.epi, 
                       vars = c("samples",
                                "Patient", 
                                "BC.Subtype", 
                                "sc50.Pred"))

# determine number of cancer epi cells per sample
drop_samples <- sobjlists %>% dplyr::count(sobjlists$samples) # generate dataframe with counts of # cells per sample
drop_samples <- drop_samples[drop_samples$n < 50, ] # get list of samples with n < 10 cells
nrow(drop_samples)

# drop samples with too few cells 
sobjlists <- sobjlists[!(sobjlists$samples %in% drop_samples$`sobjlists$samples`), ] 

# create dataframe for bar graph analysis
sobjlists <- sobjlists %>% dplyr::group_by(samples, 
                                           Patient, 
                                           BC.Subtype, 
                                           sc50.Pred) %>% # group cells by samples and BC.Subtype
  dplyr::summarise(Nb = n()) %>% # add column with number of cells by sc50.Pred
  dplyr::mutate(C = sum(Nb)) %>% # add column with total number of cells in the sample
  dplyr::mutate(percent = Nb/C*100) # add column with % of cells by sc50.Pred out of all cells in sample

# reorder sobjlists
sobjlists <- sobjlists[order(sobjlists$BC.Subtype),] # order by BC Subtype

# get list of tumor samples with enough cells for analysis
samples <- unique(sobjlists$samples)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ***GE generation -------------------------------------------------
# All cancer epi clustering ---------------------------------------------------

# prepare cancer epithelial cell object
Epi.all.combo <- cancer.epi
DefaultAssay(Epi.all.combo) <- "RNA"
Epi.all.combo <- FindNeighbors(Epi.all.combo, reduction = "pca", dims = 1:50)

# cluster cells at various resolutions 
resolution.range <- c(0.01, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6, 1.8, 2.0)
Epi.all.combo <- FindClusters(Epi.all.combo, 
                              resolution = resolution.range)

# generate cancer epithelial cell UMAP
Epi.all.combo <- FindNeighbors(Epi.all.combo, reduction = "pca", dims = 1:50)
Epi.all.combo <- FindClusters(Epi.all.combo, resolution = resolution.range)
Epi.all.combo <- RunUMAP(Epi.all.combo, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use=123)

# plot UMAP by clustering
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE) + SeuratAxes() + NoLegend()
ggsave("cancerepi_UMAP.pdf", plot = as.ggplot(p), width = 5.8, height = 5.5)

# plot UMAP by original dataset
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "orig.ident") + SeuratAxes()
ggsave("cancerepi_UMAP_origdataset.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

# plot UMAP by clinical subtype
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "BC.Subtype") + SeuratAxes()
ggsave("cancerepi_UMAP_BCsubtype.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

# plot UMAP by sc50 subtype
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "sc50.Pred") + SeuratAxes()
ggsave("cancerepi_UMAP_sc50subtype.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

# update cancer epithelial cell object
cancer.epi <- Epi.all.combo

# list of clinically actionable targets
targets <- list("ESR1",
                "ERBB2", #HER2
                "PIK3CA",
                c("NTRK1", "NTRK2", "NTRK3"),
                "CD274", #PD-L1
                "ERBB3", #HER3
                "EGFR",
                c("FGFR1", "FGFR2", "FGFR3", "FGFR4"),
                "TACSTD2", #TROP2
                c("CDK4", "CDK6"), 
                "AR",
                "NECTIN2", 
                "LAG3")

# Add UCell score for clinically actionable targets
cancer.epi <- AddModuleScore_UCell(cancer.epi,
                                   features = targets,
                                   name = names(targets),
                                   assay = "RNA")

cancer.epi@meta.data <- cancer.epi@meta.data[,-c(99:4942,4944:15186)]

colnames(cancer.epi@meta.data)[102:114] <- c("ESR1",
                                             "ERBB2", 
                                             "PIK3CA",
                                             "NTRK", 
                                             "CD274", 
                                             "ERBB3", 
                                             "EGFR",
                                             "FGFR",
                                             "TACSTD2", 
                                             "CDK",
                                             "AR",
                                             "NECTIN2", 
                                             "LAG3")

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
saveRDS(cancer.epi, "cancerepi_withtargets_110922.rds")

# Unsupervised DGE generation -------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
j <- readRDS("cancerepi_withtargets_110922.rds")

# initialize objects for unsupervised DGE gene signature generation
all_DGE <- data.frame()

# perform DGE analysis at each of the various clustering resolutions
for (i in colnames(j@meta.data)[94:98]) {
  Idents(j) <- j@meta.data[,i]
  
  # DGE analysis (cluster biomarkers)
  j.markers_DGE <- FindAllMarkers(j, only.pos = T, 
                                  min.cells.group = 50, 
                                  min.diff.pct = 0.25, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST")
  j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
  j.markers_DGE$res <- strsplit(i, "res.")[[1]][2]
  j.markers_DGE$cluster_res <- paste0(j.markers_DGE$cluster, "_", j.markers_DGE$res)
  
  all_DGE <- rbind(all_DGE, j.markers_DGE)
}

# save full output of all unsupervised DGE signatures
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(all_DGE,
           path = "cancerepi_DGE_unsupervised_3.xlsx", 
           col_names = TRUE, 
           format_headers = TRUE)

# Unsupervised sample-level DGE generation -----------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
j <- readRDS("cancerepi_withtargets_110922.rds")
cancer.epi <- j

# initialize objects for unsupervised DGE gene signature generation
all_DGE <- data.frame()

for (i in samples[49:64]) {
  subset <- subset(j, subset = samples == i)
  DefaultAssay(subset) <- "RNA"
  subset <- NormalizeData(subset, assay = "RNA")
  subset <- FindVariableFeatures(subset, 
                                 selection.method = "vst", 
                                 nfeaøtures = 2000)
  all.genes <- rownames(subset)
  subset <- ScaleData(subset, features = all.genes)
  
  subset <- FindNeighbors(subset, reduction = "pca", dims = 1:30)
  resolution.range <- c(0.01, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6, 1.8, 2.0)
  subset <- FindClusters(subset, 
                         graph.name = "integrated_snn", 
                         resolution = resolution.range)
  subset <- RunUMAP(subset, reduction = "pca", dims = 1:30, verbose = TRUE, seed.use=123)
  
  # perform DGE analysis at each of the various clustering resolutions
  for (k in colnames(subset@meta.data)[(grep("integrated_snn", colnames(subset@meta.data)[85:115]) + 84)]) {
    Idents(subset) <- subset@meta.data[,k]
    
    # DGE analysis (cluster biomarkers)
    j.markers_DGE <- FindAllMarkers(subset, only.pos = T, 
                                    min.cells.group = 5, 
                                    min.diff.pct = 0.1, 
                                    #logfc.threshold = 0.25,
                                    test.use = "MAST")
    if(dim(j.markers_DGE)[1] > 0) {
      j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
      j.markers_DGE$res <- strsplit(k, "res.")[[1]][2]
      j.markers_DGE$cluster_res <- paste0(i, "_",j.markers_DGE$cluster, "_", j.markers_DGE$res)
    }
    
    all_DGE <- rbind(all_DGE, j.markers_DGE)
  }
}

# save full output of all unsupervised DGE signatures
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(all_DGE,
           path = "cancerepi_DGE_unsupervised_samplelevel_4.xlsx", 
           col_names = TRUE, 
           format_headers = TRUE)

# Supervised DGE generation based on targets ----------------------------------

# perform DGE analysis for each target 
all_DGE <- data.frame()
j <- cancer.epi 

for (i in c(102:108)) {
  # cells in top 25% by target expression are classified as "high"/"1" vs. bottom 75% as "low"/0
  cutoff <- quantile(j@meta.data[,i],0.9)
  j@meta.data[(cancer.epi@meta.data[,i] > cutoff),i] <- "high"
  j@meta.data[(cancer.epi@meta.data[,i] <= cutoff),i] <- "med"
  j@meta.data[(cancer.epi@meta.data[,i] <= 0),i] <- "low"
  
  Idents(j) <- j@meta.data[,i]
  
  # perform DGE analysis
  j.markers_DGE <- FindAllMarkers(j, 
                                  min.cells.group = 5, 
                                  min.pct = 0.2, 
                                  logfc.threshold = 0.1,
                                  test.use = "MAST")
  j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
  j.markers_DGE$res <- NA
  j.markers_DGE$cluster_res <- paste0(colnames(j@meta.data)[i], "_", j.markers_DGE$cluster)
  
  all_DGE <- rbind(all_DGE, j.markers_DGE)
  View(j.markers_DGE)
}

# save full output of all supervised DGE signatures from clinical targets
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(all_DGE,
           path = paste0("cancerepi_DGE_supervised_1.xlsx"), 
           col_names = TRUE, 
           format_headers = TRUE)

# Supervised DGE generation based on SC50 subtype -----------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
j <- readRDS("cancerepi_withtargets_110922.rds")
cancer.epi <- j

# set idents to SC50 subtype
Idents(j) <- j@meta.data$sc50.Pred

# perform DGE analysis for each SC50 subtype
j.markers_DGE <- FindAllMarkers(j, min.cells.group = 5, 
                                min.pct = 0.2, 
                                logfc.threshold = 0.1,  
                                test.use = "MAST")
j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
j.markers_DGE$res <- NA
j.markers_DGE$cluster_res <- j.markers_DGE$cluster

# save full output of all unsupervised DGE signatures from SC50 subtype
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(j.markers_DGE,
           path = "cancerepi_DGE_sc50.xlsx", 
           col_names = TRUE,
           format_headers = TRUE)

# Calculate Jaccard similarity ------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")

# create Jaccard similarity matrix for supervised and unsupervised DGE lists
files <- list.files(pattern = ".xlsx")

cancer_DGEs <- data.frame()
for (i in files) {
  cancer_DGEs <- rbind(cancer_DGEs, readxl::read_xlsx(i))
}

# filter DGE signatures
cancer_DGEs <- DEG_Remove_mito(cancer_DGEs) # remove mitochondrial genes
#cancer_DGEs <- DEG_Remove_heat(cancer_DGEs) # remove HSP gene
cancer_DGEs <- cancer_DGEs[which(cancer_DGEs$p_val_adj < 0.05),] # filter with p-value threshold
cancer_DGEs <- cancer_DGEs[which(cancer_DGEs$avg_log2FC > 0),] # filter with log2FC threshold
cancer_DGEs <- cancer_DGEs[which(cancer_DGEs$pct.1 > 0.2),] # filter with percent threshold

# select top n genes per signature by adjusted p-value (includes ties)
cancer_DGEs <- cancer_DGEs %>%
  group_by(cluster_res) %>%
  top_n(-200, p_val_adj)

# select top n genes per signature by avg_log2FC
cancer_DGEs <- cancer_DGEs %>% 
  group_by(cluster_res) %>%
  top_n(200, avg_log2FC) 

# get number of genes in each DGE signature
DGE_counts <- cancer_DGEs %>% dplyr::count(cancer_DGEs$cluster_res)
View(DGE_counts)

# prepare signature matrix to use for calculating pairwise Jaccard indices between signatures
cancer_DGEs_forjaccard <- as.data.frame(unique(cancer_DGEs$gene))
colnames(cancer_DGEs_forjaccard) <- c('gene')
rownames(cancer_DGEs_forjaccard) <- cancer_DGEs_forjaccard$gene
for (i in unique(cancer_DGEs$cluster_res)) {
  j <- as.data.frame(cancer_DGEs[which(cancer_DGEs$cluster_res == i), ]$gene)
  j$i <- 1
  colnames(j) <- c("gene", i)
  rownames(j) <- j$gene
  cancer_DGEs_forjaccard <- left_join(cancer_DGEs_forjaccard, j, by = "gene")
}

rownames(cancer_DGEs_forjaccard) <- cancer_DGEs_forjaccard$gene
cancer_DGEs_forjaccard <- cancer_DGEs_forjaccard[,-1]
cancer_DGEs_forjaccard[is.na(cancer_DGEs_forjaccard)] <- 0

# drop unsupervised DGE signatures with fewer than n genes
drop <- vector()
for (i in 1:(dim(cancer_DGEs_forjaccard)[2])) {
  if (DGE_counts[which(DGE_counts$cluster_res ==
                       colnames(cancer_DGEs_forjaccard)[i]),3] < 20) {
    drop <- append(drop, i)
  }
}
cancer_DGEs_forjaccard <- cancer_DGEs_forjaccard[, -drop]

# calculate Jaccard similarity between all DGE signatures
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
simil(cancer_DGEs_forjaccard, NULL, "cancerepi_DGEs_forjaccard_unsupervised_only.rds", "jaccard")
cancerepi_jaccard <- as.data.frame(readRDS("cancerepi_DGEs_forjaccard_unsupervised_only.rds"))

# remove redundant unsupervised DGE signatures
drop <- vector()
for (i in c(1:(length(colnames(cancerepi_jaccard))))) {
  for (j in c((i+1):(length(colnames(cancerepi_jaccard))))) {
    if (cancerepi_jaccard[i,j] > 0.9) {
      drop <- append(drop, j)
    }
  }
}
drop <- unique(drop)
cancerepi_jaccard_noredundant <- cancerepi_jaccard[-drop, -drop]
saveRDS(cancerepi_jaccard_noredundant, "cancerepi_DGEs_forjaccard_noredundant_unsupervised_only.rds")

# plot unsupervised DGE signature Jaccard similarities
p <- simil_plot("cancerepi_DGEs_forjaccard_noredundant_unsupervised_only.rds", 0, 1, NULL)
p <- as.ggplot(p)
ggsave("GE_prelim_unsupervised_only.pdf", plot = p, width = 40, height = 40)

# cola package to define GEs ----------------------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")

set.seed(123)

rh = consensus_partition(cancerepi_jaccard_noredundant, 
                         top_value_method = "ATC",
                         partition_method = "skmeans",
                         top_n = c(dim(cancerepi_jaccard_noredundant)[1]), 
                         p_sampling = 0.8,
                         max_k = 15)

k <- suggest_best_k(rh) #k <- 10

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
pdf("GE_stats_hclust.pdf", width = 7.3, height = 4)
select_partition_number(rh, mark_best = F)
dev.off()

write.csv(get_stats(rh), "GE_stats_hclust.csv")

pdf("GE_consensusplot.pdf", width = 7, height = 4)
consensus_heatmap(rh, k = k)
dev.off()

write.csv(get_classes(rh, k = k), "GE_classes_hclust.csv")

# Define GEs ------------------------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
GEs <- read.csv("GE_classes_hclust.csv", stringsAsFactors = FALSE, na.strings = "unknown")[,c(1,2)]
colnames(GEs) <- c("cluster_res", "GE")

cancer_DGEs_GE <- left_join(cancer_DGEs, GEs, by = "cluster_res")

cancer_DGEs_GE <- cancer_DGEs_GE[-which(is.na(cancer_DGEs_GE$GE)),] 
cancer_DGEs_GE$totavg_log2FC <- NA
cancer_DGEs_GE$totmin.pct.diff <- NA

for (i in unique(cancer_DGEs_GE$gene)) {
  mean <-  mean(cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$avg_log2FC)
  cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$totavg_log2FC <- mean
  
  mean <-  mean(cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$percent)
  cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$totmin.pct.diff <- mean
}

# Filter GEs ------------------------------------------------------------------

cancer_GE <- cancer_DGEs_GE %>% dplyr::group_by(GE, gene) %>%
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 

#cancer_GE <- cancer_GE[which(cancer_GE$percent > 0.1),]

cancer_GE <- cancer_GE %>%
  group_by(GE) %>%
  top_n(350, Nb)

cancer_GE <- left_join(cancer_GE, unique(cancer_DGEs_GE[,c(7,12)]))
cancer_GE <- cancer_GE %>% 
  group_by(GE) %>%
  top_n(200, totavg_log2FC) 

GE_summary <- cancer_GE %>% dplyr::group_by(GE) %>%
  dplyr::summarise(Nb = n())

View(GE_summary)

cancer_GE$old_GE <- cancer_GE$GE
cancer_GE$GE[which(cancer_GE$old_GE == 1)] <- 1
cancer_GE$GE[which(cancer_GE$old_GE == 2)] <- 3
cancer_GE$GE[which(cancer_GE$old_GE == 3)] <- 6
cancer_GE$GE[which(cancer_GE$old_GE == 4)] <- 4
cancer_GE$GE[which(cancer_GE$old_GE == 5)] <- 5
cancer_GE$GE[which(cancer_GE$old_GE == 6)] <- 9
cancer_GE$GE[which(cancer_GE$old_GE == 7)] <- 2
cancer_GE$GE[which(cancer_GE$old_GE == 8)] <- 7
cancer_GE$GE[which(cancer_GE$old_GE == 9)] <- 10
cancer_GE$GE[which(cancer_GE$old_GE == 10)] <- 8

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(cancer_GE, "cancer_GE.xlsx")

# Add GE scores for cancer cells ----------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- split(GElist$gene,GElist$GE)
length(GElist) # get number of GEs

cancer.epi <- AddModuleScore_UCell(cancer.epi,
                                   features = GElist,
                                   assay = "RNA")

# setwd("/work/InternalMedicine/s437775/simil")
# saveRDS(cancer.epi, file = "cancerepi_withGEs_061422.rds")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Fig 3A: GE heatmap ------

set.seed(123)
expdata <- t(cancer.epi@meta.data[,c(1,4,14,83:84,115:124)])
collapse_expdata <- as.data.frame(rownames(expdata))

for (i in samples) {
  subset <- expdata[,which(expdata[2,] == i)]
  subset <- subset[,sample.int(dim(subset)[2],min(dim(subset)[2],20000))]
  collapse_expdata <- cbind(collapse_expdata, subset)
}

labels <- rownames(collapse_expdata)[6:15]
labels <- labels %>% gsub("raw_", "",.)
labels <- labels %>% gsub("X", "GE",.)

collapse_expdata <- collapse_expdata[,-1]
collapse_zscore <- collapse_expdata[-c(1:5),]
collapse_zscore <- as.matrix(sapply(collapse_zscore, as.numeric))
collapse_zscore <- t(apply(collapse_zscore, 1, function(x) (x-mean(x))/sd(x)))
sort <- rbind(apply(collapse_zscore, 2, function(x) which.max(x)),
              apply(collapse_zscore, 2, function(x) max(x)), 
              collapse_zscore)
sort <- sort[,order(sort[2,],decreasing = T)] # sort by max z-score
sort <- sort[,order(sort[1,],decreasing = F)] #sort by GE
collapse_zscore <- sort[-c(1,2),]
collapse_expdata <- collapse_expdata[,colnames(collapse_zscore)]

orig_anno <- t(as.matrix(collapse_expdata[1,]))
colnames(orig_anno) <- c("origin")
sample_anno <- t(as.matrix(collapse_expdata[2,]))
colnames(sample_anno) <- c("sample")
BC_anno <- t(as.matrix(collapse_expdata[3,]))
colnames(BC_anno) <- c("BC subtype")
pam50_anno <- t(as.matrix(collapse_expdata[5,]))
colnames(pam50_anno) <- c("PAM50 subtype")
sc50_anno <- t(as.matrix(collapse_expdata[4,]))
colnames(sc50_anno) <- c("SC50 subtype") 
anno <- HeatmapAnnotation("orig" = orig_anno,
                          "sample" = sample_anno,
                          "BC" = BC_anno,
                          "pam50" = pam50_anno,
                          "sc50" = sc50_anno, 
                          show_legend = c("sample" = FALSE), 
                          col = list("orig" = c("Pal_Prim" = "#18A900",
                                                "Qian" = "#FFC300",
                                                "Wu" = "#C70039", 
                                                "Karaayvaz" = "#eb7d34",
                                                "Wu2021prim" = "#006CA9",
                                                "Xu" = "#A300DB"),
                                     "BC" = c("HER2+" = "#700639", "HR+" = "#397006", "TNBC" = "#063970"),
                                     "pam50" = c("Basal" = "#2596be", "Her2" = "#9925be", "LumA" = "#49be25", "LumB" = "#be4d25"),
                                     "sc50" = c("Basal_SC" = "#76b5c5", "Her2E_SC" = "#ad76c5", "LumA_SC" = "#c58676", "LumB_SC" = "#8dc576"))
)

p <- Heatmap(collapse_zscore,
             heatmap_legend_param = list("title" = "UCell\nZ-score"),
             cluster_rows = F,
             cluster_columns = F,
             show_column_dend = F,
             show_row_dend = F,
             show_column_names = F,
             show_row_names = T,
             row_labels = labels,
             use_raster = T,
             raster_quality = 4, 
             top_annotation = anno)

p <- as.ggplot(p)
ggsave("GE_heatmap.pdf", plot = p, width = 12, height = 5)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Fig 3B: GE barplot -------

expdata <- t(cancer.epi@meta.data[,c(1,4,14,83:84,115:124)])

zscore <- expdata[-c(1,2,3,4,5),]
zscore <- t(as.matrix(apply(zscore, 1, as.numeric)))
zscore <- t(apply(zscore, 1, function(x) (x-mean(x))/sd(x)))
sort <- rbind(apply(zscore, 2, function(x) which.max(x)),
              apply(zscore, 2, function(x) max(x)), 
              zscore, 
              expdata[c(1,2,3,4,5),])
sort <- as.data.frame(t(sort))
colnames(sort) <- c("GE", "maxZscore", "GE1", "GE2", "GE3", "GE4", "GE5", 
                    "GE6", "GE7", "GE8", "GE9", "GE10",
                    "origin", "sample", "BC_subtype", "pam50_subtype", "sc50_subtype")

sc50_plot <- sort %>% dplyr::group_by(GE, sc50_subtype) %>% # group cells by samples and BC.Subtype
  dplyr::summarise(Nb = n()) %>% # add column with number of cells by sc50.Pred
  dplyr::mutate(C = sum(Nb)) %>% # add column with total number of cells in the sample
  dplyr::mutate(percent = Nb/C) # add column with % of cells by sc50.Pred out of all cells in sample
sc50_plot$GE <- as.numeric(sc50_plot$GE)
sc50_plot <- sc50_plot[order(sc50_plot$GE,decreasing = F),] # sort by GE
sc50_plot$GE <- factor(sc50_plot$GE, levels = unique(sc50_plot$GE))

p <- ggplot(sc50_plot, aes(fill=sc50_subtype, y=percent, x=GE)) + 
  geom_bar(width = 0.7, stat="identity") + #position = "dodge"
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), labels = (scales::percent)) + 
  scale_fill_manual(values = c("#76b5c5", "#ad76c5", "#c58676", "#8dc576")) + 
  xlab("") + ylab("Proportion of Cells") + 
  scale_x_discrete(labels=c("1" = "GE 1", 
                            "2" = "GE 2",
                            "3" = "GE 3", 
                            "4" = "GE 4",
                            "5" = "GE 5", 
                            "6" = "GE 6",
                            "7" = "GE 7", 
                            "8" = "GE 8",
                            "9" = "GE 9", 
                            "10" = "GE 10")) 
ggsave("GE_bysubtype.pdf", plot = p, width = 7, height = 2)

stacked_plot <- sort %>% dplyr::group_by(sample, BC_subtype, GE) %>% # group cells by samples and BC.Subtype
  dplyr::summarise(Nb = n()) %>% # add column with number of cells by sc50.Pred
  dplyr::mutate(C = sum(Nb)) %>% # add column with total number of cells in the sample
  dplyr::mutate(percent = Nb/C) # add column with % of cells by sc50.Pred out of all cells in sample
stacked_plot$GE <- paste0("GE ", stacked_plot$GE)
stacked_plot <- stacked_plot[order(stacked_plot$BC_subtype,decreasing = F),] # sort by GE
stacked_plot$sample <- factor(stacked_plot$sample, levels = unique(stacked_plot$sample))
stacked_plot$GE <- factor(stacked_plot$GE, levels = c("GE 1", "GE 2", "GE 3", "GE 4", 
                                                      "GE 5", "GE 6", "GE 7", "GE 8", 
                                                      "GE 9", "GE 10"))

p <- ggplot(stacked_plot, aes(fill=GE, y=percent, x=sample)) + 
  geom_bar(width = 0.7, stat="identity") + 
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), labels = (scales::percent)) + 
  scale_fill_manual(values = c("#66C2A5","#ED936B","#919FC7","#DA8EC0",
                               "#B0D867","#F9DA56","#E0C59A","#B3B3B3",
                               "#90C786", "#DDDDDD")) +
  xlab("") + ylab("Proportion of Cells")+ 
  facet_nested( ~ BC_subtype, 
                scales = "free", 
                space = "free", 
                switch = "x") + 
  theme(strip.text.x = element_text(size = 14, face = "italic"), 
        strip.background = element_rect(#colour = "white", size = 1, 
          fill = "#EEEEEE"), 
        #panel.border = element_rect(color = "black", fill = NA, size = 1), 
        panel.spacing.x = unit(0.3, "lines")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("GE_bysample.pdf", plot = p, width = 20, height = 3)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Fig 3C: Hallmark gene set analysis for GEs ---------------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
dfsample <- readxl::read_xlsx("cancer_GE.xlsx")
dfsample <- split(dfsample$gene, dfsample$GE)

library(EnsDb.Hsapiens.v79)
orgdb = "org.Hs.eg.db"

dfsample$`1` <- bitr(dfsample$`1`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)#, drop = FALSE)
dfsample$`2` <- bitr(dfsample$`2`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)#, drop = FALSE)
dfsample$`3` <- bitr(dfsample$`3`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`4` <- bitr(dfsample$`4`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`5` <- bitr(dfsample$`5`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`6` <- bitr(dfsample$`6`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`7` <- bitr(dfsample$`7`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`8` <- bitr(dfsample$`8`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`9` <- bitr(dfsample$`9`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`10` <- bitr(dfsample$`10`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)

genelist <- list("1" = dfsample$`1`$ENTREZID,
                 "2" = dfsample$`2`$ENTREZID,
                 "3" = dfsample$`3`$ENTREZID,
                 "4" = dfsample$`4`$ENTREZID,
                 "5" = dfsample$`5`$ENTREZID,
                 "6" = dfsample$`6`$ENTREZID,
                 "7" = dfsample$`7`$ENTREZID,
                 "8" = dfsample$`8`$ENTREZID, 
                 "9" = dfsample$`9`$ENTREZID,
                 "10" = dfsample$`10`$ENTREZID)

m_df <- msigdbr(species = "Homo sapiens", category ="H") %>%#, subcategory = "CGP") %>% 
  dplyr::select(gs_name, entrez_gene)

GOclusterplot <- compareCluster(genelist, 
                                fun = enricher, 
                                TERM2GENE = m_df, 
                                pvalueCutoff = .05)

p <- dotplot(GOclusterplot, includeAll = F)
p <- p + scale_y_discrete(label = function(x) {
  x %>% sub("HALLMARK_", "",.) %>% 
    sub("BIOCARTA_", "",.) %>% 
    sub("PID_", "",.) %>% 
    sub("REACTOME_", "",.) %>% 
    gsub("_", " ", .) %>% 
    stringr::str_trunc(40, "right")
}) 
p <- p + scale_x_discrete(label = function(x) {
  x %>% paste0("GE",.) %>% 
    stringr::str_sub(., 1,4)
}) 
p <- p + theme(axis.text.x = element_text(angle = 90, hjust=1))

ggsave("GE_clusterplot_hallmark.pdf", plot = p, width = 20, height = 15, units = "cm")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# *** Nichenet analysis 
# ***Nichenet analysis -------------------------------
# Create epi Nichenet object ------------

# Create epithelial cell object with GE idents
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
cancer.epi <- readRDS("cancerepi_withGEmetadata_111422.rds")
DefaultAssay(cancer.epi) <- "RNA"

GEidents <- data.frame()
for (i in c(1:10)) {
  j <- cancer.epi$maxGE
  j[which(j == i)] <- "high"
  j[which(j != "high")] <- "low"
  j <- as.data.frame(t(j))
  rownames(j) <- c(paste0("GE", i, "_idents"))
  GEidents <- rbind(GEidents, j)
}

colnames(GEidents) <- colnames(cancer.epi)
cancer.epi <- AddMetaData(cancer.epi, t(GEidents), col.name = rownames(GEidents))

# Create non-epithelial cell object
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
combo.reference <- readRDS("primobj_withCNVlabels_111022.rds")
combo.reference$celltype_withreprog <- as.character(combo.reference$celltype_withreprog)
combo.reference$celltype_withreprog[which(combo.reference$celltype_final == "Cancer Epithelial Cells")] <- "Cancer Epithelial Cells"
DefaultAssay(combo.reference) <- "RNA"

cancer.other <- subset(combo.reference, subset = celltype_withreprog != "Cancer Epithelial Cells")
cancer.other$celltype_withreprog <- Idents(cancer.other)

GEidents <- data.frame()
for (i in c(1:10)) {
  j <- as.data.frame(t(as.data.frame(cancer.other$celltype_withreprog)))
  rownames(j) <- c(paste0("GE", i, "_idents"))
  GEidents <- rbind(GEidents, j)
}

colnames(GEidents) <- colnames(cancer.other)
cancer.other <- AddMetaData(cancer.other, t(GEidents), col.name = rownames(GEidents))

nichenet_obj <- merge(cancer.epi, cancer.other)
nichenet_obj$celltype_withreprog[which(is.na(nichenet_obj$celltype_withreprog))] <- "Cancer Epithelial Cells"
nichenet_obj$celltype_final[which(nichenet_obj$celltype_withreprog == "Cancer Epithelial Cells")] <- "Cancer Epithelial Cells"

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
saveRDS(nichenet_obj, "nichenetobj_112022.rds")

# Load in NicheNet object ---------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/")
nichenet_obj <- readRDS("nichenetobj_112022.rds")

all_cell <- unique(nichenet_obj$celltype_withreprog)[-1]

# Run NicheNet -----------

setwd("/endosome/work/InternalMedicine/s437775/simil")

ligand_target_matrix = readRDS("ligand_target_matrix.rds")
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")

for (cell in all_cell[16:18]) {
  # NicheNet GE1 -------------
  
  Idents(nichenet_obj) <- nichenet_obj$GE1_idents
  temp_nichenet1 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet1, file = paste0(cell, "_nichenet1.Rdata"))
  
  # NicheNet GE2 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE2_idents
  temp_nichenet2 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet2, file = paste0(cell, "_nichenet2.Rdata"))
  
  # NicheNet GE3 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE3_idents
  temp_nichenet3 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet3, file = paste0(cell, "_nichenet3.Rdata"))
  
  # NicheNet GE4 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE4_idents
  temp_nichenet4 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet4, file = paste0(cell, "_nichenet4.Rdata"))
  
  # NicheNet GE5 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE5_idents
  temp_nichenet5 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet5, file = paste0(cell, "_nichenet5.Rdata"))
  
  # NicheNet GE6 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE6_idents
  temp_nichenet6 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet6, file = paste0(cell, "_nichenet6.Rdata"))
  
  # NicheNet GE7 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE7_idents
  temp_nichenet7 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet7, file = paste0(cell, "_nichenet7.Rdata"))
  
  # NicheNet GE8 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE8_idents
  temp_nichenet8 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet8, file = paste0(cell, "_nichenet8.Rdata"))
  
  # NicheNet GE9 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE9_idents
  temp_nichenet9 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet9, file = paste0(cell, "_nichenet9.Rdata"))
  
  # NicheNet GE10 -------------
  
  Idents(nichenet_obj) <- nichenet_obj$GE10_idents
  temp_nichenet10 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet10, file = paste0(cell, "_nichenet10.Rdata"))
  
}

# Read in Nichenet results --------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")
files <- list.files(pattern = ".Rdata")
#files <- files[-c(51, 172)]

# Nichenet object subsets ----------

all_cell <- unique(nichenet_obj$celltype_withreprog)[-1]

primsubset <- subset(nichenet_obj, subset = celltype_withreprog %in% all_cell)
Idents(primsubset) <- primsubset$celltype_withreprog

episubset <- subset(nichenet_obj, subset = celltype_final == "Cancer Epithelial Cells") 

# NicheNet circos plots -----------

all_receptors <- data.frame()
all_targets <- data.frame()

for (j in files) {
  temp_nichenet <- readRDS(j)
  cell <- strsplit(j, "_")[[1]][1]
  GE <- as.numeric(substring(strsplit(strsplit(j, "_")[[1]][2], ".R")[[1]],9)[1])
  
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
  Idents(episubset) <- episubset@meta.data[,136 + GE]
  keep_lig <- AverageExpression(subset(episubset, idents = c("high", "low")), features = gsub("\\.", "-", unique(LR_matrix$ligand)), slot = "data")$RNA
  keep_lig <- as.data.frame(keep_lig[which(keep_lig[,'high'] > keep_lig[,'low']),])
  keep_lig$ligand <- rownames(keep_lig)
  keep_lig$ligexp <- keep_lig$high / keep_lig$low
  keep_lig <- keep_lig[,-c(1:2)]
  LR_matrix <- LR_matrix[which(LR_matrix$ligand %in% keep_lig$ligand),]
  
  cellsubset <- subset(primsubset, subset = celltype_withreprog == cell)
  cellsubset <- GetAssayData(cellsubset, assay = "RNA", slot = "data")
  cellsubset <- as.data.frame(cellsubset[which(rownames(cellsubset) %in% gsub("\\.", "-", unique(LR_matrix$receptor))),])
  cellsubset[cellsubset > 0] <- 1
  cellsubset$percent <- apply(cellsubset, 1, function(x) sum(x))/dim(cellsubset)[2]
  cellsubset$receptor <- rownames(cellsubset)
  keep_rec <- as.data.frame(cbind(cellsubset$receptor, cellsubset$percent))
  colnames(keep_rec) <- c("receptor", "recpercent")
  LR_matrix <- LR_matrix[which(LR_matrix$receptor %in% keep_rec$receptor),]
  
  LR_matrix <- left_join(LR_matrix, keep_lig)
  LR_matrix <- left_join(LR_matrix, keep_rec)
  LR_matrix$LRexp <- LR_matrix$ligexp * as.numeric(LR_matrix$recpercent)
  dim(LR_matrix)
  
  lr.df.top <- LR_matrix %>% top_n(50, LRexp)
  
  if (cell %in% c("NK Cells", "Reprogrammed NK Cells")) {
    activating <- LR_matrix[which(LR_matrix$receptor %in% c("CD160", "CD226", "CD244", "CRTAM",
                                                            "FCGR3A", "KLRC2", "KLRK1", "NCR1", "NCR2", "NCR3",
                                                            "TNFRSF9")),]
    inactivating <- LR_matrix[which(LR_matrix$receptor %in% c("CD96", "KIR2DL1", "KIR2DL2", "KIR2DL3","KIR2DL4", "KLRA",
                                                              "KLRC1", "LAG3", "LILRB1", "PDCD1", "PVRIG", "TIGIT")),]
    lr.df.top <- unique(rbind(lr.df.top, activating, inactivating))
  }
  
  #lr.df.top$weight <- lr.df.top$LRexp
  lr.df.top <- lr.df.top[,c(1:3)]
  
  setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet_NK")
  
  #generate circos plot with top 10 ligands by ligexp
  pdf(paste0("NK", GE, "_", cell, "_circos.pdf"), width = 8, height = 8)
  
  circos.par(canvas.ylim=c(-1.5,1.5), # edit  canvas size
             track.margin = c(0.01, 0)) # adjust bottom and top margin
  
  chordDiagram(lr.df.top,
               directional = 1,
               link.sort = TRUE,
               link.decreasing = FALSE,
               link.visible = lr.df.top$weight > 0,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.05))
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
  }, bg.border = NA)
  circos.clear()
  
  dev.off()
  
  targets <- temp_nichenet$ligand_target_df[which(temp_nichenet$ligand_target_df$ligand %in% LR_matrix$ligand),]
  receptors <- LR_matrix
  targets <- as.data.frame(targets)
  targets$GE <- paste0("NK",GE)
  targets$cell <- cell
  all_targets <- rbind(all_targets, targets)
  
  receptors <- as.data.frame(receptors)
  receptors$GE <- paste0("NK",GE)
  receptors$cell <- cell
  all_receptors <- rbind(all_receptors, receptors)
  
}

write.csv(all_receptors, "nichenet_receptors.csv")
write.csv(all_targets, "nichenet_targets.csv")


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ***CellChat analysis --------------
# Libraries ------------

#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork) 
options(stringsAsFactors = FALSE)
library(NMF)
library(ggalluvial)
library(purrr)
library(ComplexHeatmap)
library(viridis)
library(limma)
library(tidyr)

# Load in Nichenet object ------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
nichenet_obj <- readRDS("nichenetobj_112022.rds")
DefaultAssay(nichenet_obj) <- "RNA"
nichenet_obj <- NormalizeData(nichenet_obj, verbose = FALSE)

nichenet_obj$maxGE <- paste0("GE", nichenet_obj$maxGE)
nichenet_obj$maxGE[which(nichenet_obj$maxGE == "GENA")] <- nichenet_obj$celltype_withreprog[which(nichenet_obj$maxGE == "GENA")]
data.input <- GetAssayData(nichenet_obj, assay = "RNA", slot = "data") # normalized data matrix

# Run CellChat -----------

cellchats <- list()
nets <- list()

for (i in c(82)) {
  labels <- nichenet_obj@meta.data[,i]
  meta <- data.frame(group = labels, row.names = rownames(nichenet_obj@meta.data))
  
  # Create CellChat object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

  # Load CellChat database
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  #showDatabaseCategory(CellChatDB)
  
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  cellchat@DB <- CellChatDB.use
  
  # Pre-process data
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  
  cellchats <- append(cellchats, cellchat)
  
  net <- subsetCommunication(cellchat, 
                             slot.name = "net", 
                             thresh = 0.05)
  nets <- append(nets, net)
}

saveRDS(cellchats, "cellchat_obj_113022.rds")
saveRDS(nets, "cellchat_nets_113022.rds")

# CellChat output ----------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
cellchats <- readRDS("cellchat_obj_113022.rds")

net <- data.frame(stringsAsFactors = F)

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet_NK")
all_receptors <- read.csv("nichenet_receptors.csv", row.names = 1, stringsAsFactors = F)[,-1]

net <- subsetCommunication(cellchats[[1]],
                           slot.name = "net",
                           thresh = 0.05,
                           targets.use = unique(all_receptors$cell))

write.csv(net, "cellchat_receptors.csv")

# -----
# -----
# Fig 3D: GE-immune decoder matrix from combined NicheNet and CellChat results  -------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")
all_receptors <- read.csv("nichenet_receptors.csv", row.names = 1, stringsAsFactors = F)[,-1]
net <- read.csv("cellchat_receptors.csv", row.names = 1, stringsAsFactors = F)

nichenet_receptors <- all_receptors
newgenes <- alias2SymbolTable(nichenet_receptors$receptor, species = "Hs")
newgenes <- ifelse(is.na(newgenes),nichenet_receptors$receptor,newgenes)
nichenet_receptors$receptor <- newgenes

nichenet_receptors <- nichenet_receptors  %>%
  group_by(cell) %>%
  top_n(2300, weight)

nichenet_receptors <- nichenet_receptors  %>%
  group_by(GE) %>%
  top_n(500, LRexp)

write.csv(nichenet_receptors, "cellchat_nichenet_receptors_list.csv")

nichenet_receptors <- nichenet_receptors[, c(2,6,7)]
colnames(nichenet_receptors) <- c("receptor", "GE", "cell")

cellchat_receptors <- net[which(net$prob >= 0.0),]
cellchat_receptors <- cellchat_receptors[,c(4,1,2)]
colnames(cellchat_receptors) <- c("receptor", "GE", "cell")
cellchat_receptors = separate_rows(cellchat_receptors, 1, sep = "_")
cellchat_receptors = separate_rows(cellchat_receptors, 1, sep = ":")
newgenes <- alias2SymbolTable(cellchat_receptors$receptor, species = "Hs")
newgenes <- ifelse(is.na(newgenes),cellchat_receptors$receptor,newgenes)
cellchat_receptors$receptor <- newgenes

overlap <- cellchat_receptors[which(do.call(paste0, cellchat_receptors) %in% do.call(paste0, nichenet_receptors)),]
overlap <- rbind(overlap,nichenet_receptors)
overlap <- unique(overlap)

overlap <- overlap %>% dplyr::group_by(GE, cell) %>%
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 
overlap <- split(overlap[,2:3], overlap$GE)
names <- names(overlap)
overlap <- purrr::reduce(overlap, full_join, by = "cell") %>% replace(., is.na(.), 0);
overlap <- as.data.frame(overlap)
colnames(overlap) <- c("cell", names)
rownames(overlap) <- overlap$cell
overlap <- overlap[,c(7:17)]
overlap <- overlap[,-1]

overlap <- overlap[-which(rownames(overlap) == "Epithelial Cells"),]
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")

pdf("cellchat_nichenet_heatmap.pdf")
Heatmap(overlap, 
        col = circlize::colorRamp2(seq(min(overlap),max(overlap), 
                                       length.out = 500), maGEa(500)), 
        cluster_rows = T)
dev.off()

pdf("cellchat_nichenet_heatmap_scaled.pdf")
Heatmap(overlap / rowSums(overlap), 
        col = circlize::colorRamp2(seq(min(overlap / rowSums(overlap)),
                                       #max(overlap / rowSums(overlap)), 
                                       0.2,
                                       length.out = 500), maGEa(500)), 
        cluster_rows = T)
dev.off()

write.csv(overlap, "cellchat_nichenet_receptors_sum.csv")

# -----
# -----

# Fig 3E: Cell line GE analysis -----------

# Load in files -----------------------

# from CCLE DepMap portal 
fullmat <- read.csv("CCLE_RNAseq_reads.csv")
rownames(fullmat) <- fullmat$X
fullmat <- fullmat[,-1]
colnames(fullmat) <- sub('\\.\\..*', '', colnames(fullmat))
fullmat <- t(fullmat)
fullmat[1:5,1:5]

sample.info <- read.csv("sample_info.csv")
rownames(sample.info) <- sample.info$DepMap_ID
head(sample.info)
cell.lines.wanted <- colnames(fullmat)
sample.info <- sample.info[rownames(sample.info) %in% cell.lines.wanted, ]
dim(sample.info)

cell_lines <- CreateSeuratObject(fullmat, assay = "RNA", meta.data = sample.info)
DefaultAssay(cell_lines) <- "RNA"
cell_lines <- NormalizeData(cell_lines)

# Subset cell lines and add GE module scores --------------- 

setwd("/endosome/work/InternalMedicine/s437775/simil")

# Load in supplemental data from Sheffer et al. for breast cancer cell lines 
paper_lines <- read.csv("cell_line_screens.csv", stringsAsFactors = F)

colnames(paper_lines) <- c("Cell.Lines", "Tissue",
                           "24hr_AUC", "48hr_AUC", "72hr_AUC",
                           "24hr_sensitivity", "48hr_sensitivity", "72hr_sensitivity",
                           "Mesenchymal", "Epithelial", "MSI_high",
                           "B7H6_protein_scores", "B7H6_high",
                           "HLA_protein_scores", "HLA_negative")

breast_lines <- subset(cell_lines, 
                       subset = stripped_cell_line_name %in% paper_lines$Cell.Lines)

newnames <- breast_lines@meta.data[,c(1,6)]
newnames$orig.ident <- rownames(newnames)
colnames(newnames) <- c("orig_names", "Cell.Lines")

paper_lines <- left_join(paper_lines,newnames)
paper_lines <- paper_lines[!is.na(paper_lines$orig_names),]
rownames(paper_lines) <- paper_lines$orig_names

breast_lines <- AddMetaData(breast_lines, paper_lines[,-16])

colnames(breast_lines@meta.data)

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- split(GElist$gene,GElist$GE)
length(GElist) # get number of GEs

breast_lines <- AddModuleScore_UCell(breast_lines,
                                     features = GElist,
                                     assay = "RNA")

# Z-score GE scores 
breast_lines@meta.data[41:50] <- apply(breast_lines@meta.data[41:50], 2, function(x) (x-mean(x))/sd(x))

# GE vs 24 hr sensitivity ----------------

breast_lines$X72hr_killing <- 1 - breast_lines$X72hr_AUC
breast_lines$X48hr_killing <- 1 - breast_lines$X48hr_AUC
breast_lines$X24hr_killing <- 1 - breast_lines$X24hr_AUC

p <- list() 
p.val <- vector()
corr <- vector()

breast_lines$include <- 0
breast_lines$include[which(breast_lines$X24hr_sensitivity < breast_lines$X72hr_sensitivity)] <- 1
breast_lines$include[which(breast_lines$X24hr_sensitivity < breast_lines$X48hr_sensitivity)] <- 1

for (i in c(1:10)) {
  j <- colnames(breast_lines@meta.data)[40 + i]
  subset <- breast_lines@meta.data[which(breast_lines$include == 0),]
  p[[i]] <- ggscatter(subset, x = "X24hr_killing", y = j,
                      add = "reg.line",
                      conf.int = TRUE,
                      add.params = list(color = "blue", fill = "lightgray"),
                      xlab = "Sensitivity to NK Cell Killing", 
                      ylab = paste0("GE", i, " Expression")) + 
    stat_cor(method = "spearman", 
             label.sep = "\n", label.y = 0, 
             label.x = 0.5) + 
    xlim(min(subset$X24hr_killing), max(subset$X24hr_killing)) 
  
  p.val <- append(p.val, cor.test(subset$X24hr_killing, subset[,40 + i], method = "spearman")$p.value)
  corr <- append(corr, cor.test(subset$X24hr_killing, subset[,40 + i], method = "spearman")$estimate)
}

p.val.adj <- p.adjust(p.val, method = "BH")
p.val.adj
p.val
dim(subset)

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
pdf("cell_lines_24hrAUCvsUCell_112523.pdf", width = 5, height = 5)
print(p)
dev.off()


# -----
# -----

# Fig 3F: Nichenet circos plots for NK cells vs. GE1 and GE6 -----------

# Refer to Nichenet analysis section above for code to generate circos plots# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# This will generate all of the sub-figures for Figure 3, -------
# and is the analysis of cancer epithelial cell transcriptional
# heterogeneity defined by 10 GEs.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# SETUP: -----------------------------------------------------------------------
# Libraries -------------------------------------------------------------------

library(BiocManager)
library(GEOquery) 
library(plyr)
library(dplyr) 
library(Matrix)
library(Seurat)
library(ggplot2)
library(cowplot) 
library(multtest)
library(msigdbr)
library(fgsea)
library(monocle3)
library(velocyto.R)
library(loomR)
library(clustree)
library(tibble)
library(SeuratData)
library(matrixStats)
library(sparseMatrixStats)
library(DESeq2)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(viridis)
library(gridExtra)
library(ggplotify)
library(multtest)
library(metap)
library(writexl)
library(Rcpp)
library(RcppZiggurat)
library(Rfast)
library(ggh4x)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(cola)
library(msigdbr)
library(UCell)
library(RColorBrewer)

# Functions -------------------------------------------------------------------

# DEG_remove_mito 
DEG_Remove_mito <- function (df){ø
  df_rm_mito <- df[!grepl("^MT-|^MT.",df$gene),]
  return(df_rm_mito)
}

# DEG_remove_heat 
DEG_Remove_heat <- function (df){
  df_rm_hsp <- df[!grepl("^HSP",rownames(df)),]
  return(df_rm_hsp)
}

# plotSimilarityMatrix 
# create a similarity matrix plot from a dataframe
# adapted from klic package (adjusted heatmap parameters)
plotSimilarityMatrix = function(X, y = NULL, 
                                min.val = 0, 
                                max.val = 1,
                                clusLabels = NULL, 
                                colX = NULL, colY = NULL, 
                                clr = FALSE, clc = FALSE, 
                                annotation_col = NULL, 
                                annotation_row = NULL, 
                                annotation_colors = NULL, 
                                myLegend = NULL, 
                                fileName = "posteriorSimilarityMatrix", 
                                savePNG = FALSE, 
                                semiSupervised = FALSE, 
                                showObsNames = FALSE) {
  
  if (!is.null(y)) {
    # Check if the rownames correspond to the ones in the similarity matrix
    check <- sum(1 - rownames(X) %in% row.names(y))
    if (check == 1)
      stop("X and y must have the same row names.")
  }
  
  if (!is.null(clusLabels)) {
    if (!is.integer(clusLabels))
      stop("Cluster labels must be integers.")
    
    n_clusters <- length(table(clusLabels))
    riordina <- NULL
    for (i in 1:n_clusters) {
      riordina <- c(riordina, which(clusLabels == i))
    }
    
    X <- X[riordina, riordina]
    y <- y[riordina, ]
    y <- as.data.frame(y)
  }
  
  if (savePNG)
    grDevices::png(paste(fileName, ".png", sep = ""))
  
  if (!is.null(y)) {
    ht <- ComplexHeatmap::pheatmap(X, legend = TRUE,  
                                   color = rev(brewer.pal(11, "RdBu")), 
                                   breaks = seq(min.val, max.val, length.out = 11),
                                   cluster_rows = clr, 
                                   cluster_cols = clc, 
                                   #annotation_col = y,
                                   show_rownames = showObsNames, 
                                   show_colnames = showObsNames, 
                                   drop_levels = TRUE, 
                                   treeheight_row = -1,
                                   treeheight_col = -1,
                                   #annotation_row = annotation_row, 
                                   #annotation_col = annotation_col, 
                                   annotation_colors = annotation_colors)
  } else {
    ht <- ComplexHeatmap::pheatmap(X, legend = TRUE,
                                   color = rev(brewer.pal(11, "RdBu")),  
                                   breaks = seq(min.val, max.val, length.out = 11),
                                   cluster_rows = clr, 
                                   cluster_cols = clc,
                                   show_rownames = showObsNames, 
                                   show_colnames = showObsNames,
                                   treeheight_row = -1,
                                   treeheight_col = -1,
                                   drop_levels = TRUE, 
                                   #annotation_row = annotation_row, 
                                   #annotation_col = annotation_col, 
                                   annotation_colors = annotation_colors)
  }
  
  if (savePNG)
    grDevices::dev.off()
  
  return(ht)
}

# simil function 
# calculate similarity matrix 
# param df is dataframe of all cells for similarity matrix computation
# param drop is list of genes to drop from analysis
# param file is file name for output
# param method ("jaccard" or "corr") for method used
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

# simil_plot 
# plot interactive similarity heatmap 
# param a is an .rds file generated by simil_calc
# param annot is a vector used for heatmap annotation (metadata column from Seurat object)

simil_plot <- function(a, min.val, max.val, annot) {
  jc <- readRDS(a) #read in similarity matrix
  
  # add annotation
  if (length(annot) > 0) {
    row_annot <- annot[rownames(jc), , drop = FALSE] # select rows to use for heatmap annotation
    col_annot <- annot[colnames(jc), , drop = FALSE] # select columns to use for heatmap annotation
    colors <- mako(n_distinct(annot)) # create color vector to use for annotation
    names(colors) <- base::unique(annot)[[1]]
    colors <- list(colors, colors)
    names(colors) <- c(as.name(names(annot)), as.name(names(annot)))
  } 
  else {
    row_annot <- NULL
    col_annot <- NULL
    colors <- NULL
  }
  
  # generate simialrity matrix plot (using plotSimilarity Matrix)
  hm <- plotSimilarityMatrix(jc, clr = TRUE, clc = TRUE, 
                             min.val = min.val, max.val = max.val,
                             annotation_row = row_annot, 
                             annotation_col = col_annot, 
                             annotation_colors = colors,
                             showObsNames = T) # plot full matrix
  return(hm)
  #hm <- draw(hm)
  #htShiny(hm) # generate interactive heatmap 
}

# simil_sc50 
setwd("/work/InternalMedicine/s437775/simil")

# read in sc50 genes and convert to base::unique vector
sc50 <- read.csv("NatGen_Supplementary_table_S4.csv")
sc50 <- base::unique(unlist(sc50))
sc50 <- gsub("\\.", "-", sc50)
sc50 <- which(rownames(cancer.epi) %in% sc50) # get indices of sc50 genes

# calculate Jaccard similarity matrix 
# param df is dataframe for comparison
# param file is file name for output
# param method ("jaccard" or "corr") for similarity method used

simil_sc50 <- function(df, file, method) {
  jc <- df[sc50, , drop = FALSE] # subset to sc50 genes
  # jc[jc > 0.5] <- 1
  # jc[jc < 0.5 & jc > -0.5] <- 0
  # jc[jc < -0.5] <- -1
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

# simil_GE function 
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")

GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- which(rownames(cancer.epi) %in% GElist$gene) # get indices of GE genes

# calculate Jaccard similarity matrix 
# param df is dataframe for comparison
# param file is file name for output
# param method ("jaccard" or "corr") for similarity method used
simil_GE <- function(df, file, method) {
  jc <- df[GElist, , drop = FALSE] # subset to sc50 genes
  # jc[jc > 0.5] <- 1
  # jc[jc < 0.5 & jc > -0.5] <- 0
  # jc[jc < -0.5] <- -1
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


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ***Load in cancer epithelial cell object --------
# Load in Seurat object -------------------------------------------------------

# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/")
# combo.reference <- readRDS("PrimObject_FINALnoZgenedem_71222.rds") ##newnewnew

# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/")
# cancer.epi <- readRDS("Episubset_FINALnoZallgenedem_71322.rds") 
# cancer.epi <- NormalizeData(cancer.epi, assay = "RNA")

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
cancer.epi <- readRDS("cancerepi_withGEmetadata_111422.rds")
DefaultAssay(cancer.epi) <- "RNA"

# Scale data ------------------------------------------------------------------

DefaultAssay(cancer.epi) <- "RNA"
cancer.epi <- NormalizeData(cancer.epi, assay = "RNA")
cancer.epi <- FindVariableFeatures(cancer.epi, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)
all.genes <- rownames(cancer.epi)
cancer.epi <- ScaleData(cancer.epi, features = all.genes)

# Drop samples with too few cells ---------------------------------------------

# create dataframe of Seurat object metadata
sobjlists <- FetchData(object = cancer.epi, 
                       vars = c("samples",
                                "Patient", 
                                "BC.Subtype", 
                                "sc50.Pred"))

# determine number of cancer epi cells per sample
drop_samples <- sobjlists %>% dplyr::count(sobjlists$samples) # generate dataframe with counts of # cells per sample
drop_samples <- drop_samples[drop_samples$n < 50, ] # get list of samples with n < 10 cells
nrow(drop_samples)

# drop samples with too few cells 
sobjlists <- sobjlists[!(sobjlists$samples %in% drop_samples$`sobjlists$samples`), ] 

# create dataframe for bar graph analysis
sobjlists <- sobjlists %>% dplyr::group_by(samples, 
                                           Patient, 
                                           BC.Subtype, 
                                           sc50.Pred) %>% # group cells by samples and BC.Subtype
  dplyr::summarise(Nb = n()) %>% # add column with number of cells by sc50.Pred
  dplyr::mutate(C = sum(Nb)) %>% # add column with total number of cells in the sample
  dplyr::mutate(percent = Nb/C*100) # add column with % of cells by sc50.Pred out of all cells in sample

# reorder sobjlists
sobjlists <- sobjlists[order(sobjlists$BC.Subtype),] # order by BC Subtype

# get list of tumor samples with enough cells for analysis
samples <- unique(sobjlists$samples)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ***GE generation -------------------------------------------------
# All cancer epi clustering ---------------------------------------------------

# prepare cancer epithelial cell object
Epi.all.combo <- cancer.epi
DefaultAssay(Epi.all.combo) <- "RNA"
Epi.all.combo <- FindNeighbors(Epi.all.combo, reduction = "pca", dims = 1:50)

# cluster cells at various resolutions 
resolution.range <- c(0.01, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6, 1.8, 2.0)
Epi.all.combo <- FindClusters(Epi.all.combo, 
                              resolution = resolution.range)

# generate cancer epithelial cell UMAP
Epi.all.combo <- FindNeighbors(Epi.all.combo, reduction = "pca", dims = 1:50)
Epi.all.combo <- FindClusters(Epi.all.combo, resolution = resolution.range)
Epi.all.combo <- RunUMAP(Epi.all.combo, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use=123)

# plot UMAP by clustering
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE) + SeuratAxes() + NoLegend()
ggsave("cancerepi_UMAP.pdf", plot = as.ggplot(p), width = 5.8, height = 5.5)

# plot UMAP by original dataset
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "orig.ident") + SeuratAxes()
ggsave("cancerepi_UMAP_origdataset.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

# plot UMAP by clinical subtype
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "BC.Subtype") + SeuratAxes()
ggsave("cancerepi_UMAP_BCsubtype.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

# plot UMAP by sc50 subtype
p <- DimPlot(Epi.all.combo, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "sc50.Pred") + SeuratAxes()
ggsave("cancerepi_UMAP_sc50subtype.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

# update cancer epithelial cell object
cancer.epi <- Epi.all.combo

# list of clinically actionable targets
targets <- list("ESR1",
                "ERBB2", #HER2
                "PIK3CA",
                c("NTRK1", "NTRK2", "NTRK3"),
                "CD274", #PD-L1
                "ERBB3", #HER3
                "EGFR",
                c("FGFR1", "FGFR2", "FGFR3", "FGFR4"),
                "TACSTD2", #TROP2
                c("CDK4", "CDK6"), 
                "AR",
                "NECTIN2", 
                "LAG3")

# Add UCell score for clinically actionable targets
cancer.epi <- AddModuleScore_UCell(cancer.epi,
                                   features = targets,
                                   name = names(targets),
                                   assay = "RNA")

cancer.epi@meta.data <- cancer.epi@meta.data[,-c(99:4942,4944:15186)]

colnames(cancer.epi@meta.data)[102:114] <- c("ESR1",
                                             "ERBB2", 
                                             "PIK3CA",
                                             "NTRK", 
                                             "CD274", 
                                             "ERBB3", 
                                             "EGFR",
                                             "FGFR",
                                             "TACSTD2", 
                                             "CDK",
                                             "AR",
                                             "NECTIN2", 
                                             "LAG3")

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
saveRDS(cancer.epi, "cancerepi_withtargets_110922.rds")

# Unsupervised DGE generation -------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
j <- readRDS("cancerepi_withtargets_110922.rds")

# initialize objects for unsupervised DGE gene signature generation
all_DGE <- data.frame()

# perform DGE analysis at each of the various clustering resolutions
for (i in colnames(j@meta.data)[94:98]) {
  Idents(j) <- j@meta.data[,i]
  
  # DGE analysis (cluster biomarkers)
  j.markers_DGE <- FindAllMarkers(j, only.pos = T, 
                                  min.cells.group = 50, 
                                  min.diff.pct = 0.25, 
                                  logfc.threshold = 0.25,
                                  test.use = "MAST")
  j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
  j.markers_DGE$res <- strsplit(i, "res.")[[1]][2]
  j.markers_DGE$cluster_res <- paste0(j.markers_DGE$cluster, "_", j.markers_DGE$res)
  
  all_DGE <- rbind(all_DGE, j.markers_DGE)
}

# save full output of all unsupervised DGE signatures
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(all_DGE,
           path = "cancerepi_DGE_unsupervised_3.xlsx", 
           col_names = TRUE, 
           format_headers = TRUE)

# Unsupervised sample-level DGE generation -----------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
j <- readRDS("cancerepi_withtargets_110922.rds")
cancer.epi <- j

# initialize objects for unsupervised DGE gene signature generation
all_DGE <- data.frame()

for (i in samples[49:64]) {
  subset <- subset(j, subset = samples == i)
  DefaultAssay(subset) <- "RNA"
  subset <- NormalizeData(subset, assay = "RNA")
  subset <- FindVariableFeatures(subset, 
                                 selection.method = "vst", 
                                 nfeaøtures = 2000)
  all.genes <- rownames(subset)
  subset <- ScaleData(subset, features = all.genes)
  
  subset <- FindNeighbors(subset, reduction = "pca", dims = 1:30)
  resolution.range <- c(0.01, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6, 1.8, 2.0)
  subset <- FindClusters(subset, 
                         graph.name = "integrated_snn", 
                         resolution = resolution.range)
  subset <- RunUMAP(subset, reduction = "pca", dims = 1:30, verbose = TRUE, seed.use=123)
  
  # perform DGE analysis at each of the various clustering resolutions
  for (k in colnames(subset@meta.data)[(grep("integrated_snn", colnames(subset@meta.data)[85:115]) + 84)]) {
    Idents(subset) <- subset@meta.data[,k]
    
    # DGE analysis (cluster biomarkers)
    j.markers_DGE <- FindAllMarkers(subset, only.pos = T, 
                                    min.cells.group = 5, 
                                    min.diff.pct = 0.1, 
                                    #logfc.threshold = 0.25,
                                    test.use = "MAST")
    if(dim(j.markers_DGE)[1] > 0) {
      j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
      j.markers_DGE$res <- strsplit(k, "res.")[[1]][2]
      j.markers_DGE$cluster_res <- paste0(i, "_",j.markers_DGE$cluster, "_", j.markers_DGE$res)
    }
    
    all_DGE <- rbind(all_DGE, j.markers_DGE)
  }
}

# save full output of all unsupervised DGE signatures
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(all_DGE,
           path = "cancerepi_DGE_unsupervised_samplelevel_4.xlsx", 
           col_names = TRUE, 
           format_headers = TRUE)

# Supervised DGE generation based on targets ----------------------------------

# perform DGE analysis for each target 
all_DGE <- data.frame()
j <- cancer.epi 

for (i in c(102:108)) {
  # cells in top 25% by target expression are classified as "high"/"1" vs. bottom 75% as "low"/0
  cutoff <- quantile(j@meta.data[,i],0.9)
  j@meta.data[(cancer.epi@meta.data[,i] > cutoff),i] <- "high"
  j@meta.data[(cancer.epi@meta.data[,i] <= cutoff),i] <- "med"
  j@meta.data[(cancer.epi@meta.data[,i] <= 0),i] <- "low"
  
  Idents(j) <- j@meta.data[,i]
  
  # perform DGE analysis
  j.markers_DGE <- FindAllMarkers(j, 
                                  min.cells.group = 5, 
                                  min.pct = 0.2, 
                                  logfc.threshold = 0.1,
                                  test.use = "MAST")
  j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
  j.markers_DGE$res <- NA
  j.markers_DGE$cluster_res <- paste0(colnames(j@meta.data)[i], "_", j.markers_DGE$cluster)
  
  all_DGE <- rbind(all_DGE, j.markers_DGE)
  View(j.markers_DGE)
}

# save full output of all supervised DGE signatures from clinical targets
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(all_DGE,
           path = paste0("cancerepi_DGE_supervised_1.xlsx"), 
           col_names = TRUE, 
           format_headers = TRUE)

# Supervised DGE generation based on SC50 subtype -----------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
j <- readRDS("cancerepi_withtargets_110922.rds")
cancer.epi <- j

# set idents to SC50 subtype
Idents(j) <- j@meta.data$sc50.Pred

# perform DGE analysis for each SC50 subtype
j.markers_DGE <- FindAllMarkers(j, min.cells.group = 5, 
                                min.pct = 0.2, 
                                logfc.threshold = 0.1,  
                                test.use = "MAST")
j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
j.markers_DGE$res <- NA
j.markers_DGE$cluster_res <- j.markers_DGE$cluster

# save full output of all unsupervised DGE signatures from SC50 subtype
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(j.markers_DGE,
           path = "cancerepi_DGE_sc50.xlsx", 
           col_names = TRUE,
           format_headers = TRUE)

# Calculate Jaccard similarity ------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")

# create Jaccard similarity matrix for supervised and unsupervised DGE lists
files <- list.files(pattern = ".xlsx")

cancer_DGEs <- data.frame()
for (i in files) {
  cancer_DGEs <- rbind(cancer_DGEs, readxl::read_xlsx(i))
}

# filter DGE signatures
cancer_DGEs <- DEG_Remove_mito(cancer_DGEs) # remove mitochondrial genes
#cancer_DGEs <- DEG_Remove_heat(cancer_DGEs) # remove HSP gene
cancer_DGEs <- cancer_DGEs[which(cancer_DGEs$p_val_adj < 0.05),] # filter with p-value threshold
cancer_DGEs <- cancer_DGEs[which(cancer_DGEs$avg_log2FC > 0),] # filter with log2FC threshold
cancer_DGEs <- cancer_DGEs[which(cancer_DGEs$pct.1 > 0.2),] # filter with percent threshold

# select top n genes per signature by adjusted p-value (includes ties)
cancer_DGEs <- cancer_DGEs %>%
  group_by(cluster_res) %>%
  top_n(-200, p_val_adj)

# select top n genes per signature by avg_log2FC
cancer_DGEs <- cancer_DGEs %>% 
  group_by(cluster_res) %>%
  top_n(200, avg_log2FC) 

# get number of genes in each DGE signature
DGE_counts <- cancer_DGEs %>% dplyr::count(cancer_DGEs$cluster_res)
View(DGE_counts)

# prepare signature matrix to use for calculating pairwise Jaccard indices between signatures
cancer_DGEs_forjaccard <- as.data.frame(unique(cancer_DGEs$gene))
colnames(cancer_DGEs_forjaccard) <- c('gene')
rownames(cancer_DGEs_forjaccard) <- cancer_DGEs_forjaccard$gene
for (i in unique(cancer_DGEs$cluster_res)) {
  j <- as.data.frame(cancer_DGEs[which(cancer_DGEs$cluster_res == i), ]$gene)
  j$i <- 1
  colnames(j) <- c("gene", i)
  rownames(j) <- j$gene
  cancer_DGEs_forjaccard <- left_join(cancer_DGEs_forjaccard, j, by = "gene")
}

rownames(cancer_DGEs_forjaccard) <- cancer_DGEs_forjaccard$gene
cancer_DGEs_forjaccard <- cancer_DGEs_forjaccard[,-1]
cancer_DGEs_forjaccard[is.na(cancer_DGEs_forjaccard)] <- 0

# drop unsupervised DGE signatures with fewer than n genes
drop <- vector()
for (i in 1:(dim(cancer_DGEs_forjaccard)[2])) {
  if (DGE_counts[which(DGE_counts$cluster_res ==
                       colnames(cancer_DGEs_forjaccard)[i]),3] < 20) {
    drop <- append(drop, i)
  }
}
cancer_DGEs_forjaccard <- cancer_DGEs_forjaccard[, -drop]

# calculate Jaccard similarity between all DGE signatures
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
simil(cancer_DGEs_forjaccard, NULL, "cancerepi_DGEs_forjaccard_unsupervised_only.rds", "jaccard")
cancerepi_jaccard <- as.data.frame(readRDS("cancerepi_DGEs_forjaccard_unsupervised_only.rds"))

# remove redundant unsupervised DGE signatures
drop <- vector()
for (i in c(1:(length(colnames(cancerepi_jaccard))))) {
  for (j in c((i+1):(length(colnames(cancerepi_jaccard))))) {
    if (cancerepi_jaccard[i,j] > 0.9) {
      drop <- append(drop, j)
    }
  }
}
drop <- unique(drop)
cancerepi_jaccard_noredundant <- cancerepi_jaccard[-drop, -drop]
saveRDS(cancerepi_jaccard_noredundant, "cancerepi_DGEs_forjaccard_noredundant_unsupervised_only.rds")

# plot unsupervised DGE signature Jaccard similarities
p <- simil_plot("cancerepi_DGEs_forjaccard_noredundant_unsupervised_only.rds", 0, 1, NULL)
p <- as.ggplot(p)
ggsave("GE_prelim_unsupervised_only.pdf", plot = p, width = 40, height = 40)

# cola package to define GEs ----------------------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")

set.seed(123)

rh = consensus_partition(cancerepi_jaccard_noredundant, 
                         top_value_method = "ATC",
                         partition_method = "skmeans",
                         top_n = c(dim(cancerepi_jaccard_noredundant)[1]), 
                         p_sampling = 0.8,
                         max_k = 15)

k <- suggest_best_k(rh) #k <- 10

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
pdf("GE_stats_hclust.pdf", width = 7.3, height = 4)
select_partition_number(rh, mark_best = F)
dev.off()

write.csv(get_stats(rh), "GE_stats_hclust.csv")

pdf("GE_consensusplot.pdf", width = 7, height = 4)
consensus_heatmap(rh, k = k)
dev.off()

write.csv(get_classes(rh, k = k), "GE_classes_hclust.csv")

# Define GEs ------------------------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
GEs <- read.csv("GE_classes_hclust.csv", stringsAsFactors = FALSE, na.strings = "unknown")[,c(1,2)]
colnames(GEs) <- c("cluster_res", "GE")

cancer_DGEs_GE <- left_join(cancer_DGEs, GEs, by = "cluster_res")

cancer_DGEs_GE <- cancer_DGEs_GE[-which(is.na(cancer_DGEs_GE$GE)),] 
cancer_DGEs_GE$totavg_log2FC <- NA
cancer_DGEs_GE$totmin.pct.diff <- NA

for (i in unique(cancer_DGEs_GE$gene)) {
  mean <-  mean(cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$avg_log2FC)
  cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$totavg_log2FC <- mean
  
  mean <-  mean(cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$percent)
  cancer_DGEs_GE[which(cancer_DGEs_GE$gene == i),]$totmin.pct.diff <- mean
}

# Filter GEs ------------------------------------------------------------------

cancer_GE <- cancer_DGEs_GE %>% dplyr::group_by(GE, gene) %>%
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 

#cancer_GE <- cancer_GE[which(cancer_GE$percent > 0.1),]

cancer_GE <- cancer_GE %>%
  group_by(GE) %>%
  top_n(350, Nb)

cancer_GE <- left_join(cancer_GE, unique(cancer_DGEs_GE[,c(7,12)]))
cancer_GE <- cancer_GE %>% 
  group_by(GE) %>%
  top_n(200, totavg_log2FC) 

GE_summary <- cancer_GE %>% dplyr::group_by(GE) %>%
  dplyr::summarise(Nb = n())

View(GE_summary)

cancer_GE$old_GE <- cancer_GE$GE
cancer_GE$GE[which(cancer_GE$old_GE == 1)] <- 1
cancer_GE$GE[which(cancer_GE$old_GE == 2)] <- 3
cancer_GE$GE[which(cancer_GE$old_GE == 3)] <- 6
cancer_GE$GE[which(cancer_GE$old_GE == 4)] <- 4
cancer_GE$GE[which(cancer_GE$old_GE == 5)] <- 5
cancer_GE$GE[which(cancer_GE$old_GE == 6)] <- 9
cancer_GE$GE[which(cancer_GE$old_GE == 7)] <- 2
cancer_GE$GE[which(cancer_GE$old_GE == 8)] <- 7
cancer_GE$GE[which(cancer_GE$old_GE == 9)] <- 10
cancer_GE$GE[which(cancer_GE$old_GE == 10)] <- 8

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
write_xlsx(cancer_GE, "cancer_GE.xlsx")

# Add GE scores for cancer cells ----------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- split(GElist$gene,GElist$GE)
length(GElist) # get number of GEs

cancer.epi <- AddModuleScore_UCell(cancer.epi,
                                   features = GElist,
                                   assay = "RNA")

# setwd("/work/InternalMedicine/s437775/simil")
# saveRDS(cancer.epi, file = "cancerepi_withGEs_061422.rds")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Fig 3A: GE heatmap ------

set.seed(123)
expdata <- t(cancer.epi@meta.data[,c(1,4,14,83:84,115:124)])
collapse_expdata <- as.data.frame(rownames(expdata))

for (i in samples) {
  subset <- expdata[,which(expdata[2,] == i)]
  subset <- subset[,sample.int(dim(subset)[2],min(dim(subset)[2],20000))]
  collapse_expdata <- cbind(collapse_expdata, subset)
}

labels <- rownames(collapse_expdata)[6:15]
labels <- labels %>% gsub("raw_", "",.)
labels <- labels %>% gsub("X", "GE",.)

collapse_expdata <- collapse_expdata[,-1]
collapse_zscore <- collapse_expdata[-c(1:5),]
collapse_zscore <- as.matrix(sapply(collapse_zscore, as.numeric))
collapse_zscore <- t(apply(collapse_zscore, 1, function(x) (x-mean(x))/sd(x)))
sort <- rbind(apply(collapse_zscore, 2, function(x) which.max(x)),
              apply(collapse_zscore, 2, function(x) max(x)), 
              collapse_zscore)
sort <- sort[,order(sort[2,],decreasing = T)] # sort by max z-score
sort <- sort[,order(sort[1,],decreasing = F)] #sort by GE
collapse_zscore <- sort[-c(1,2),]
collapse_expdata <- collapse_expdata[,colnames(collapse_zscore)]

orig_anno <- t(as.matrix(collapse_expdata[1,]))
colnames(orig_anno) <- c("origin")
sample_anno <- t(as.matrix(collapse_expdata[2,]))
colnames(sample_anno) <- c("sample")
BC_anno <- t(as.matrix(collapse_expdata[3,]))
colnames(BC_anno) <- c("BC subtype")
pam50_anno <- t(as.matrix(collapse_expdata[5,]))
colnames(pam50_anno) <- c("PAM50 subtype")
sc50_anno <- t(as.matrix(collapse_expdata[4,]))
colnames(sc50_anno) <- c("SC50 subtype") 
anno <- HeatmapAnnotation("orig" = orig_anno,
                          "sample" = sample_anno,
                          "BC" = BC_anno,
                          "pam50" = pam50_anno,
                          "sc50" = sc50_anno, 
                          show_legend = c("sample" = FALSE), 
                          col = list("orig" = c("Pal_Prim" = "#18A900",
                                                "Qian" = "#FFC300",
                                                "Wu" = "#C70039", 
                                                "Karaayvaz" = "#eb7d34",
                                                "Wu2021prim" = "#006CA9",
                                                "Xu" = "#A300DB"),
                                     "BC" = c("HER2+" = "#700639", "HR+" = "#397006", "TNBC" = "#063970"),
                                     "pam50" = c("Basal" = "#2596be", "Her2" = "#9925be", "LumA" = "#49be25", "LumB" = "#be4d25"),
                                     "sc50" = c("Basal_SC" = "#76b5c5", "Her2E_SC" = "#ad76c5", "LumA_SC" = "#c58676", "LumB_SC" = "#8dc576"))
)

p <- Heatmap(collapse_zscore,
             heatmap_legend_param = list("title" = "UCell\nZ-score"),
             cluster_rows = F,
             cluster_columns = F,
             show_column_dend = F,
             show_row_dend = F,
             show_column_names = F,
             show_row_names = T,
             row_labels = labels,
             use_raster = T,
             raster_quality = 4, 
             top_annotation = anno)

p <- as.ggplot(p)
ggsave("GE_heatmap.pdf", plot = p, width = 12, height = 5)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Fig 3B: GE barplot -------

expdata <- t(cancer.epi@meta.data[,c(1,4,14,83:84,115:124)])

zscore <- expdata[-c(1,2,3,4,5),]
zscore <- t(as.matrix(apply(zscore, 1, as.numeric)))
zscore <- t(apply(zscore, 1, function(x) (x-mean(x))/sd(x)))
sort <- rbind(apply(zscore, 2, function(x) which.max(x)),
              apply(zscore, 2, function(x) max(x)), 
              zscore, 
              expdata[c(1,2,3,4,5),])
sort <- as.data.frame(t(sort))
colnames(sort) <- c("GE", "maxZscore", "GE1", "GE2", "GE3", "GE4", "GE5", 
                    "GE6", "GE7", "GE8", "GE9", "GE10",
                    "origin", "sample", "BC_subtype", "pam50_subtype", "sc50_subtype")

sc50_plot <- sort %>% dplyr::group_by(GE, sc50_subtype) %>% # group cells by samples and BC.Subtype
  dplyr::summarise(Nb = n()) %>% # add column with number of cells by sc50.Pred
  dplyr::mutate(C = sum(Nb)) %>% # add column with total number of cells in the sample
  dplyr::mutate(percent = Nb/C) # add column with % of cells by sc50.Pred out of all cells in sample
sc50_plot$GE <- as.numeric(sc50_plot$GE)
sc50_plot <- sc50_plot[order(sc50_plot$GE,decreasing = F),] # sort by GE
sc50_plot$GE <- factor(sc50_plot$GE, levels = unique(sc50_plot$GE))

p <- ggplot(sc50_plot, aes(fill=sc50_subtype, y=percent, x=GE)) + 
  geom_bar(width = 0.7, stat="identity") + #position = "dodge"
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), labels = (scales::percent)) + 
  scale_fill_manual(values = c("#76b5c5", "#ad76c5", "#c58676", "#8dc576")) + 
  xlab("") + ylab("Proportion of Cells") + 
  scale_x_discrete(labels=c("1" = "GE 1", 
                            "2" = "GE 2",
                            "3" = "GE 3", 
                            "4" = "GE 4",
                            "5" = "GE 5", 
                            "6" = "GE 6",
                            "7" = "GE 7", 
                            "8" = "GE 8",
                            "9" = "GE 9", 
                            "10" = "GE 10")) 
ggsave("GE_bysubtype.pdf", plot = p, width = 7, height = 2)

stacked_plot <- sort %>% dplyr::group_by(sample, BC_subtype, GE) %>% # group cells by samples and BC.Subtype
  dplyr::summarise(Nb = n()) %>% # add column with number of cells by sc50.Pred
  dplyr::mutate(C = sum(Nb)) %>% # add column with total number of cells in the sample
  dplyr::mutate(percent = Nb/C) # add column with % of cells by sc50.Pred out of all cells in sample
stacked_plot$GE <- paste0("GE ", stacked_plot$GE)
stacked_plot <- stacked_plot[order(stacked_plot$BC_subtype,decreasing = F),] # sort by GE
stacked_plot$sample <- factor(stacked_plot$sample, levels = unique(stacked_plot$sample))
stacked_plot$GE <- factor(stacked_plot$GE, levels = c("GE 1", "GE 2", "GE 3", "GE 4", 
                                                      "GE 5", "GE 6", "GE 7", "GE 8", 
                                                      "GE 9", "GE 10"))

p <- ggplot(stacked_plot, aes(fill=GE, y=percent, x=sample)) + 
  geom_bar(width = 0.7, stat="identity") + 
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), labels = (scales::percent)) + 
  scale_fill_manual(values = c("#66C2A5","#ED936B","#919FC7","#DA8EC0",
                               "#B0D867","#F9DA56","#E0C59A","#B3B3B3",
                               "#90C786", "#DDDDDD")) +
  xlab("") + ylab("Proportion of Cells")+ 
  facet_nested( ~ BC_subtype, 
                scales = "free", 
                space = "free", 
                switch = "x") + 
  theme(strip.text.x = element_text(size = 14, face = "italic"), 
        strip.background = element_rect(#colour = "white", size = 1, 
          fill = "#EEEEEE"), 
        #panel.border = element_rect(color = "black", fill = NA, size = 1), 
        panel.spacing.x = unit(0.3, "lines")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("GE_bysample.pdf", plot = p, width = 20, height = 3)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Fig 3C: Hallmark gene set analysis for GEs ---------------------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
dfsample <- readxl::read_xlsx("cancer_GE.xlsx")
dfsample <- split(dfsample$gene, dfsample$GE)

library(EnsDb.Hsapiens.v79)
orgdb = "org.Hs.eg.db"

dfsample$`1` <- bitr(dfsample$`1`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)#, drop = FALSE)
dfsample$`2` <- bitr(dfsample$`2`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)#, drop = FALSE)
dfsample$`3` <- bitr(dfsample$`3`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`4` <- bitr(dfsample$`4`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`5` <- bitr(dfsample$`5`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`6` <- bitr(dfsample$`6`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`7` <- bitr(dfsample$`7`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`8` <- bitr(dfsample$`8`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`9` <- bitr(dfsample$`9`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$`10` <- bitr(dfsample$`10`, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)

genelist <- list("1" = dfsample$`1`$ENTREZID,
                 "2" = dfsample$`2`$ENTREZID,
                 "3" = dfsample$`3`$ENTREZID,
                 "4" = dfsample$`4`$ENTREZID,
                 "5" = dfsample$`5`$ENTREZID,
                 "6" = dfsample$`6`$ENTREZID,
                 "7" = dfsample$`7`$ENTREZID,
                 "8" = dfsample$`8`$ENTREZID, 
                 "9" = dfsample$`9`$ENTREZID,
                 "10" = dfsample$`10`$ENTREZID)

m_df <- msigdbr(species = "Homo sapiens", category ="H") %>%#, subcategory = "CGP") %>% 
  dplyr::select(gs_name, entrez_gene)

GOclusterplot <- compareCluster(genelist, 
                                fun = enricher, 
                                TERM2GENE = m_df, 
                                pvalueCutoff = .05)

p <- dotplot(GOclusterplot, includeAll = F)
p <- p + scale_y_discrete(label = function(x) {
  x %>% sub("HALLMARK_", "",.) %>% 
    sub("BIOCARTA_", "",.) %>% 
    sub("PID_", "",.) %>% 
    sub("REACTOME_", "",.) %>% 
    gsub("_", " ", .) %>% 
    stringr::str_trunc(40, "right")
}) 
p <- p + scale_x_discrete(label = function(x) {
  x %>% paste0("GE",.) %>% 
    stringr::str_sub(., 1,4)
}) 
p <- p + theme(axis.text.x = element_text(angle = 90, hjust=1))

ggsave("GE_clusterplot_hallmark.pdf", plot = p, width = 20, height = 15, units = "cm")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# *** Nichenet analysis 
# ***Nichenet analysis -------------------------------
# Create epi Nichenet object ------------

# Create epithelial cell object with GE idents
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
cancer.epi <- readRDS("cancerepi_withGEmetadata_111422.rds")
DefaultAssay(cancer.epi) <- "RNA"

GEidents <- data.frame()
for (i in c(1:10)) {
  j <- cancer.epi$maxGE
  j[which(j == i)] <- "high"
  j[which(j != "high")] <- "low"
  j <- as.data.frame(t(j))
  rownames(j) <- c(paste0("GE", i, "_idents"))
  GEidents <- rbind(GEidents, j)
}

colnames(GEidents) <- colnames(cancer.epi)
cancer.epi <- AddMetaData(cancer.epi, t(GEidents), col.name = rownames(GEidents))

# Create non-epithelial cell object
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
combo.reference <- readRDS("primobj_withCNVlabels_111022.rds")
combo.reference$celltype_withreprog <- as.character(combo.reference$celltype_withreprog)
combo.reference$celltype_withreprog[which(combo.reference$celltype_final == "Cancer Epithelial Cells")] <- "Cancer Epithelial Cells"
DefaultAssay(combo.reference) <- "RNA"

cancer.other <- subset(combo.reference, subset = celltype_withreprog != "Cancer Epithelial Cells")
cancer.other$celltype_withreprog <- Idents(cancer.other)

GEidents <- data.frame()
for (i in c(1:10)) {
  j <- as.data.frame(t(as.data.frame(cancer.other$celltype_withreprog)))
  rownames(j) <- c(paste0("GE", i, "_idents"))
  GEidents <- rbind(GEidents, j)
}

colnames(GEidents) <- colnames(cancer.other)
cancer.other <- AddMetaData(cancer.other, t(GEidents), col.name = rownames(GEidents))

nichenet_obj <- merge(cancer.epi, cancer.other)
nichenet_obj$celltype_withreprog[which(is.na(nichenet_obj$celltype_withreprog))] <- "Cancer Epithelial Cells"
nichenet_obj$celltype_final[which(nichenet_obj$celltype_withreprog == "Cancer Epithelial Cells")] <- "Cancer Epithelial Cells"

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
saveRDS(nichenet_obj, "nichenetobj_112022.rds")

# Load in NicheNet object ---------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/")
nichenet_obj <- readRDS("nichenetobj_112022.rds")

all_cell <- unique(nichenet_obj$celltype_withreprog)[-1]

# Run NicheNet -----------

setwd("/endosome/work/InternalMedicine/s437775/simil")

ligand_target_matrix = readRDS("ligand_target_matrix.rds")
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")

for (cell in all_cell[16:18]) {
  # NicheNet GE1 -------------
  
  Idents(nichenet_obj) <- nichenet_obj$GE1_idents
  temp_nichenet1 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet1, file = paste0(cell, "_nichenet1.Rdata"))
  
  # NicheNet GE2 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE2_idents
  temp_nichenet2 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet2, file = paste0(cell, "_nichenet2.Rdata"))
  
  # NicheNet GE3 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE3_idents
  temp_nichenet3 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet3, file = paste0(cell, "_nichenet3.Rdata"))
  
  # NicheNet GE4 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE4_idents
  temp_nichenet4 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet4, file = paste0(cell, "_nichenet4.Rdata"))
  
  # NicheNet GE5 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE5_idents
  temp_nichenet5 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet5, file = paste0(cell, "_nichenet5.Rdata"))
  
  # NicheNet GE6 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE6_idents
  temp_nichenet6 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet6, file = paste0(cell, "_nichenet6.Rdata"))
  
  # NicheNet GE7 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE7_idents
  temp_nichenet7 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet7, file = paste0(cell, "_nichenet7.Rdata"))
  
  # NicheNet GE8 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE8_idents
  temp_nichenet8 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet8, file = paste0(cell, "_nichenet8.Rdata"))
  
  # NicheNet GE9 -------------
  Idents(nichenet_obj) <- nichenet_obj$GE9_idents
  temp_nichenet9 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet9, file = paste0(cell, "_nichenet9.Rdata"))
  
  # NicheNet GE10 -------------
  
  Idents(nichenet_obj) <- nichenet_obj$GE10_idents
  temp_nichenet10 = nichenet_seuratobj_cluster_de(
    seurat_obj = nichenet_obj,
    assay_oi = "RNA",
    receiver_affected = cell,
    receiver_reference = setdiff(all_cell, cell),
    sender = "high",
    geneset = "up", top_n_ligands = 50, top_n_targets = 200,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
  
  saveRDS(temp_nichenet10, file = paste0(cell, "_nichenet10.Rdata"))
  
}

# Read in Nichenet results --------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")
files <- list.files(pattern = ".Rdata")
#files <- files[-c(51, 172)]

# Nichenet object subsets ----------

all_cell <- unique(nichenet_obj$celltype_withreprog)[-1]

primsubset <- subset(nichenet_obj, subset = celltype_withreprog %in% all_cell)
Idents(primsubset) <- primsubset$celltype_withreprog

episubset <- subset(nichenet_obj, subset = celltype_final == "Cancer Epithelial Cells") 

# NicheNet circos plots -----------

all_receptors <- data.frame()
all_targets <- data.frame()

for (j in files) {
  temp_nichenet <- readRDS(j)
  cell <- strsplit(j, "_")[[1]][1]
  GE <- as.numeric(substring(strsplit(strsplit(j, "_")[[1]][2], ".R")[[1]],9)[1])
  
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
  Idents(episubset) <- episubset@meta.data[,136 + GE]
  keep_lig <- AverageExpression(subset(episubset, idents = c("high", "low")), features = gsub("\\.", "-", unique(LR_matrix$ligand)), slot = "data")$RNA
  keep_lig <- as.data.frame(keep_lig[which(keep_lig[,'high'] > keep_lig[,'low']),])
  keep_lig$ligand <- rownames(keep_lig)
  keep_lig$ligexp <- keep_lig$high / keep_lig$low
  keep_lig <- keep_lig[,-c(1:2)]
  LR_matrix <- LR_matrix[which(LR_matrix$ligand %in% keep_lig$ligand),]
  
  cellsubset <- subset(primsubset, subset = celltype_withreprog == cell)
  cellsubset <- GetAssayData(cellsubset, assay = "RNA", slot = "data")
  cellsubset <- as.data.frame(cellsubset[which(rownames(cellsubset) %in% gsub("\\.", "-", unique(LR_matrix$receptor))),])
  cellsubset[cellsubset > 0] <- 1
  cellsubset$percent <- apply(cellsubset, 1, function(x) sum(x))/dim(cellsubset)[2]
  cellsubset$receptor <- rownames(cellsubset)
  keep_rec <- as.data.frame(cbind(cellsubset$receptor, cellsubset$percent))
  colnames(keep_rec) <- c("receptor", "recpercent")
  LR_matrix <- LR_matrix[which(LR_matrix$receptor %in% keep_rec$receptor),]
  
  LR_matrix <- left_join(LR_matrix, keep_lig)
  LR_matrix <- left_join(LR_matrix, keep_rec)
  LR_matrix$LRexp <- LR_matrix$ligexp * as.numeric(LR_matrix$recpercent)
  dim(LR_matrix)
  
  lr.df.top <- LR_matrix %>% top_n(50, LRexp)
  
  if (cell %in% c("NK Cells", "Reprogrammed NK Cells")) {
    activating <- LR_matrix[which(LR_matrix$receptor %in% c("CD160", "CD226", "CD244", "CRTAM",
                                                            "FCGR3A", "KLRC2", "KLRK1", "NCR1", "NCR2", "NCR3",
                                                            "TNFRSF9")),]
    inactivating <- LR_matrix[which(LR_matrix$receptor %in% c("CD96", "KIR2DL1", "KIR2DL2", "KIR2DL3","KIR2DL4", "KLRA",
                                                              "KLRC1", "LAG3", "LILRB1", "PDCD1", "PVRIG", "TIGIT")),]
    lr.df.top <- unique(rbind(lr.df.top, activating, inactivating))
  }
  
  #lr.df.top$weight <- lr.df.top$LRexp
  lr.df.top <- lr.df.top[,c(1:3)]
  
  setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet_NK")
  
  #generate circos plot with top 10 ligands by ligexp
  pdf(paste0("NK", GE, "_", cell, "_circos.pdf"), width = 8, height = 8)
  
  circos.par(canvas.ylim=c(-1.5,1.5), # edit  canvas size
             track.margin = c(0.01, 0)) # adjust bottom and top margin
  
  chordDiagram(lr.df.top,
               directional = 1,
               link.sort = TRUE,
               link.decreasing = FALSE,
               link.visible = lr.df.top$weight > 0,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.05))
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
  }, bg.border = NA)
  circos.clear()
  
  dev.off()
  
  targets <- temp_nichenet$ligand_target_df[which(temp_nichenet$ligand_target_df$ligand %in% LR_matrix$ligand),]
  receptors <- LR_matrix
  targets <- as.data.frame(targets)
  targets$GE <- paste0("NK",GE)
  targets$cell <- cell
  all_targets <- rbind(all_targets, targets)
  
  receptors <- as.data.frame(receptors)
  receptors$GE <- paste0("NK",GE)
  receptors$cell <- cell
  all_receptors <- rbind(all_receptors, receptors)
  
}

write.csv(all_receptors, "nichenet_receptors.csv")
write.csv(all_targets, "nichenet_targets.csv")


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ***CellChat analysis --------------
# Libraries ------------

#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork) 
options(stringsAsFactors = FALSE)
library(NMF)
library(ggalluvial)
library(purrr)
library(ComplexHeatmap)
library(viridis)
library(limma)
library(tidyr)

# Load in Nichenet object ------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
nichenet_obj <- readRDS("nichenetobj_112022.rds")
DefaultAssay(nichenet_obj) <- "RNA"
nichenet_obj <- NormalizeData(nichenet_obj, verbose = FALSE)

nichenet_obj$maxGE <- paste0("GE", nichenet_obj$maxGE)
nichenet_obj$maxGE[which(nichenet_obj$maxGE == "GENA")] <- nichenet_obj$celltype_withreprog[which(nichenet_obj$maxGE == "GENA")]
data.input <- GetAssayData(nichenet_obj, assay = "RNA", slot = "data") # normalized data matrix

# Run CellChat -----------

cellchats <- list()
nets <- list()

for (i in c(82)) {
  labels <- nichenet_obj@meta.data[,i]
  meta <- data.frame(group = labels, row.names = rownames(nichenet_obj@meta.data))
  
  # Create CellChat object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

  # Load CellChat database
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  #showDatabaseCategory(CellChatDB)
  
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  cellchat@DB <- CellChatDB.use
  
  # Pre-process data
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  
  cellchats <- append(cellchats, cellchat)
  
  net <- subsetCommunication(cellchat, 
                             slot.name = "net", 
                             thresh = 0.05)
  nets <- append(nets, net)
}

saveRDS(cellchats, "cellchat_obj_113022.rds")
saveRDS(nets, "cellchat_nets_113022.rds")

# CellChat output ----------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
cellchats <- readRDS("cellchat_obj_113022.rds")

net <- data.frame(stringsAsFactors = F)

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet_NK")
all_receptors <- read.csv("nichenet_receptors.csv", row.names = 1, stringsAsFactors = F)[,-1]

net <- subsetCommunication(cellchats[[1]],
                           slot.name = "net",
                           thresh = 0.05,
                           targets.use = unique(all_receptors$cell))

write.csv(net, "cellchat_receptors.csv")

# -----
# -----
# Fig 3D: GE-immune decoder matrix from combined NicheNet and CellChat results  -------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")
all_receptors <- read.csv("nichenet_receptors.csv", row.names = 1, stringsAsFactors = F)[,-1]
net <- read.csv("cellchat_receptors.csv", row.names = 1, stringsAsFactors = F)

nichenet_receptors <- all_receptors
newgenes <- alias2SymbolTable(nichenet_receptors$receptor, species = "Hs")
newgenes <- ifelse(is.na(newgenes),nichenet_receptors$receptor,newgenes)
nichenet_receptors$receptor <- newgenes

nichenet_receptors <- nichenet_receptors  %>%
  group_by(cell) %>%
  top_n(2300, weight)

nichenet_receptors <- nichenet_receptors  %>%
  group_by(GE) %>%
  top_n(500, LRexp)

write.csv(nichenet_receptors, "cellchat_nichenet_receptors_list.csv")

nichenet_receptors <- nichenet_receptors[, c(2,6,7)]
colnames(nichenet_receptors) <- c("receptor", "GE", "cell")

cellchat_receptors <- net[which(net$prob >= 0.0),]
cellchat_receptors <- cellchat_receptors[,c(4,1,2)]
colnames(cellchat_receptors) <- c("receptor", "GE", "cell")
cellchat_receptors = separate_rows(cellchat_receptors, 1, sep = "_")
cellchat_receptors = separate_rows(cellchat_receptors, 1, sep = ":")
newgenes <- alias2SymbolTable(cellchat_receptors$receptor, species = "Hs")
newgenes <- ifelse(is.na(newgenes),cellchat_receptors$receptor,newgenes)
cellchat_receptors$receptor <- newgenes

overlap <- cellchat_receptors[which(do.call(paste0, cellchat_receptors) %in% do.call(paste0, nichenet_receptors)),]
overlap <- rbind(overlap,nichenet_receptors)
overlap <- unique(overlap)

overlap <- overlap %>% dplyr::group_by(GE, cell) %>%
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 
overlap <- split(overlap[,2:3], overlap$GE)
names <- names(overlap)
overlap <- purrr::reduce(overlap, full_join, by = "cell") %>% replace(., is.na(.), 0);
overlap <- as.data.frame(overlap)
colnames(overlap) <- c("cell", names)
rownames(overlap) <- overlap$cell
overlap <- overlap[,c(7:17)]
overlap <- overlap[,-1]

overlap <- overlap[-which(rownames(overlap) == "Epithelial Cells"),]
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")

pdf("cellchat_nichenet_heatmap.pdf")
Heatmap(overlap, 
        col = circlize::colorRamp2(seq(min(overlap),max(overlap), 
                                       length.out = 500), maGEa(500)), 
        cluster_rows = T)
dev.off()

pdf("cellchat_nichenet_heatmap_scaled.pdf")
Heatmap(overlap / rowSums(overlap), 
        col = circlize::colorRamp2(seq(min(overlap / rowSums(overlap)),
                                       #max(overlap / rowSums(overlap)), 
                                       0.2,
                                       length.out = 500), maGEa(500)), 
        cluster_rows = T)
dev.off()

write.csv(overlap, "cellchat_nichenet_receptors_sum.csv")

# -----
# -----

# Fig 3E: Cell line GE analysis -----------

# Load in files -----------------------

# from CCLE DepMap portal 
fullmat <- read.csv("CCLE_RNAseq_reads.csv")
rownames(fullmat) <- fullmat$X
fullmat <- fullmat[,-1]
colnames(fullmat) <- sub('\\.\\..*', '', colnames(fullmat))
fullmat <- t(fullmat)
fullmat[1:5,1:5]

sample.info <- read.csv("sample_info.csv")
rownames(sample.info) <- sample.info$DepMap_ID
head(sample.info)
cell.lines.wanted <- colnames(fullmat)
sample.info <- sample.info[rownames(sample.info) %in% cell.lines.wanted, ]
dim(sample.info)

cell_lines <- CreateSeuratObject(fullmat, assay = "RNA", meta.data = sample.info)
DefaultAssay(cell_lines) <- "RNA"
cell_lines <- NormalizeData(cell_lines)

# Subset cell lines and add GE module scores --------------- 

setwd("/endosome/work/InternalMedicine/s437775/simil")

# Load in supplemental data from Sheffer et al. for breast cancer cell lines 
paper_lines <- read.csv("cell_line_screens.csv", stringsAsFactors = F)

colnames(paper_lines) <- c("Cell.Lines", "Tissue",
                           "24hr_AUC", "48hr_AUC", "72hr_AUC",
                           "24hr_sensitivity", "48hr_sensitivity", "72hr_sensitivity",
                           "Mesenchymal", "Epithelial", "MSI_high",
                           "B7H6_protein_scores", "B7H6_high",
                           "HLA_protein_scores", "HLA_negative")

breast_lines <- subset(cell_lines, 
                       subset = stripped_cell_line_name %in% paper_lines$Cell.Lines)

newnames <- breast_lines@meta.data[,c(1,6)]
newnames$orig.ident <- rownames(newnames)
colnames(newnames) <- c("orig_names", "Cell.Lines")

paper_lines <- left_join(paper_lines,newnames)
paper_lines <- paper_lines[!is.na(paper_lines$orig_names),]
rownames(paper_lines) <- paper_lines$orig_names

breast_lines <- AddMetaData(breast_lines, paper_lines[,-16])

colnames(breast_lines@meta.data)

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- split(GElist$gene,GElist$GE)
length(GElist) # get number of GEs

breast_lines <- AddModuleScore_UCell(breast_lines,
                                     features = GElist,
                                     assay = "RNA")

# Z-score GE scores 
breast_lines@meta.data[41:50] <- apply(breast_lines@meta.data[41:50], 2, function(x) (x-mean(x))/sd(x))

# GE vs 24 hr sensitivity ----------------

breast_lines$X72hr_killing <- 1 - breast_lines$X72hr_AUC
breast_lines$X48hr_killing <- 1 - breast_lines$X48hr_AUC
breast_lines$X24hr_killing <- 1 - breast_lines$X24hr_AUC

p <- list() 
p.val <- vector()
corr <- vector()

breast_lines$include <- 0
breast_lines$include[which(breast_lines$X24hr_sensitivity < breast_lines$X72hr_sensitivity)] <- 1
breast_lines$include[which(breast_lines$X24hr_sensitivity < breast_lines$X48hr_sensitivity)] <- 1

for (i in c(1:10)) {
  j <- colnames(breast_lines@meta.data)[40 + i]
  subset <- breast_lines@meta.data[which(breast_lines$include == 0),]
  p[[i]] <- ggscatter(subset, x = "X24hr_killing", y = j,
                      add = "reg.line",
                      conf.int = TRUE,
                      add.params = list(color = "blue", fill = "lightgray"),
                      xlab = "Sensitivity to NK Cell Killing", 
                      ylab = paste0("GE", i, " Expression")) + 
    stat_cor(method = "spearman", 
             label.sep = "\n", label.y = 0, 
             label.x = 0.5) + 
    xlim(min(subset$X24hr_killing), max(subset$X24hr_killing)) 
  
  p.val <- append(p.val, cor.test(subset$X24hr_killing, subset[,40 + i], method = "spearman")$p.value)
  corr <- append(corr, cor.test(subset$X24hr_killing, subset[,40 + i], method = "spearman")$estimate)
}

p.val.adj <- p.adjust(p.val, method = "BH")
p.val.adj
p.val
dim(subset)

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
pdf("cell_lines_24hrAUCvsUCell_112523.pdf", width = 5, height = 5)
print(p)
dev.off()


# -----
# -----

# Fig 3F: Nichenet circos plots for NK cells vs. GE1 and GE6 -----------

# Refer to Nichenet analysis section above for code to generate circos plots
                                       
# -----
# -----
