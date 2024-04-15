
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# This will generate all of the sub-figures for Figure 4, -------
# and is the analysis of GE-immune interactions using spatial datasets
# and validation of InteractPrint on Bassez et al. and I-SPY2 datasets. 
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
# Fig 4A: CD8 T cell correlation with GEs across spatial samples --------


# Read in files ----------

counts <- readRDS("counts_spatial.RDS")
corr <- readRDS("corr_spatial.RDS")
meta <- readRDS("meta_spatial.RDS")
slides_all <- readRDS("all_spatial_obj.RDS")

names <- c("CID4535", "CID44971", "CID4465", "CID4290", 
           "1160920F", "1142243F",
           rep(NA, 8),
           "CytAssist_FFPE_Human_Breast_Cancer", 
           "Parent_Visium_Human_BreastCancer",
           "V1_Breast_Cancer_Block_A_Section_1", 
           "V1_Human_Invasive_Ductal_Carcinoma",
           "Visium_FFPE_Human_Breast_Cancer", 
           "Visium_Human_Breast_Cancer")

samples <- c("HR+", "TNBC", "TNBC", "HR+", "TNBC", "TNBC", 
             "HER2+", "HER2+", "HER2+", "HER2+", "HER2+", "HER2+", "HER2+", "HER2+", 
             unique(slides_all[[15]]$subtype),
             unique(slides_all[[16]]$subtype),
             unique(slides_all[[17]]$subtype),
             unique(slides_all[[18]]$subtype),
             NA, #unique(slide_all[[19]]$subtype),
             unique(slides_all[[20]]$subtype))

# Correlation analysis -------

files <- list.files(pattern = "corr.rds")

corr <- list()
corr_T <- data.frame()
for (i in files) {
  temp <- as.data.frame(readRDS(i))
  corr <- append(corr, i)
  if (i == files[1]) {
    corr_T <- temp$`CD8 T Cells`
  }
  if (i != files[1]) {
    corr_T <- cbind(corr_T, temp$`CD8+ T Cells`)
  }
}
colnames(corr_T) <- gsub("_corr.rds", "", files)
rownames(corr_T) <- rownames(temp)

corr_T <- corr_T[grep("GE", rownames(corr_T)),]

anno <- HeatmapAnnotation(subtype = samples[c(1:6,15:18,20)],
                          col = list(subtype = c("HER2+" = "#B0FC96",
                                                 "HR+" = "#96B0FC",
                                                 "TNBC" = "#FC96B0")))

pdf(paste0("CD8_corr_heatmap.pdf"), width = 6, height = 6.5)
Heatmap(corr_T, 
        top_annotation = anno,
        column_split = factor(samples[c(1:6,15:18,20)]),
        col = colorRamp2(c(-0.4, 0, 0.4), c("blue", "#F5F5F5", "red")))
dev.off()

data_cor <- data.frame()
for (j in c(1:10)) {
  for (i in c(1:6)) {
    k <- cor.test(meta[[i]][,27+j], meta[[i]]$`CD8+ T Cells`,
                  alternative = "two.sided",
                  method = "pearson",
                  exact = NULL,
                  conf.level = 0.95)
    k <- k$p.value
    names(k) <- unique(meta[[i]]$orig.ident)
    data_cor[j,i] <- k  
  }
  
  for (i in c(7:11)) {
    k <- cor.test(meta[[i]][,23+j], meta[[i]]$`CD8+ T Cells`,
                  alternative = "two.sided",
                  method = "pearson",
                  exact = NULL,
                  conf.level = 0.95)
    k <- k$p.value
    names(k) <- unique(meta[[i]]$orig.ident)
    data_cor[j,i] <- k  
  }
}

colnames(data_cor) <- c("CID4535","CID44971","CID4465","CID4290","1160920F","1142243F",
                        "CytAssist_FFPE_Human_Breast_Cancer", 
                        "Parent_Visium_Human_BreastCancer",
                        "V1_Breast_Cancer_Block_A_Section_1", 
                        "V1_Human_Invasive_Ductal_Carcinoma",
                        "Visium_Human_Breast_Cancer")
rownames(data_cor) <- c("GE1", "GE2", "GE3", "GE4", "GE5", "GE6", "GE7", "GE8", "GE9", "GE10")
write.csv(data_cor, "CD8_corr_heatmap_pvalues.csv")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Fig 4B: Representative TNBC spatial sample --------


# Spatial images -------

name <- "1160920F"
slide <- slides_all[[5]]

predict <- as.data.frame(slide@meta.data)

p <- SpatialFeaturePlot(slide, features = "GE1", 
                        images = "image", 
                        stroke = NA, 
                        alpha = c(0,0),
                        min.cutoff = 0,
                        max.cutoff = "q90",
                        image.alpha = 1) + 
  scale_fill_gradient2(mid = "red", high = "red") + NoLegend() + 
  ggtitle("H&E")
ggsave(paste0(name, "_image.pdf"), plot = as.ggplot(p), width = 3, height = 3)

# plot GE signature score across spots
for (i in c(1:10)) {
  p <- SpatialFeaturePlot(slide, features = paste0("GE",i), 
                          images = "image", 
                          stroke = NA, 
                          alpha = c(0,1),
                          min.cutoff = 0,
                          max.cutoff = "q90",
                          pt.size.factor = 2.8,
                          image.alpha = 1) + 
    scale_fill_gradient2(mid = "red", high = "red") + NoLegend() + 
    ggtitle(paste0("Areas with GE", i, "-labeled Cells"))
  ggsave(paste0(name, "_GE", i, ".pdf"), 
         plot = as.ggplot(p), width = 3, height = 3)
  
}

# plot CD8+ T cell signature score across spots
p <- SpatialFeaturePlot(slide, features = "CD8+ T Cells", 
                        images = "image", 
                        stroke = NA, 
                        alpha = c(0,1),
                        min.cutoff = 0,
                        max.cutoff = "q90",
                        pt.size.factor = 2.8,
                        image.alpha = 1) + 
  scale_fill_gradient2(mid = "red", high = "red") + NoLegend() + 
  ggtitle("Areas with CD8+ T Cells")
ggsave(paste0(name, "_CD8.pdf"), plot = as.ggplot(p), width = 3, height = 3)

min_max_norm <- function(x){ (x - min(x)) / (max(x) - min(x)) }

# plot ligand and receptor co-localization
lig <- "LTB" # or ITGB2, ALOX5AP
rec <- "TNFRSF1A" # or ITGAL, ALOX5
data <- slide[['Spatial']]@data
slide <- AddMetaData(slide, 
                     (min_max_norm(data[which(rownames(data) == lig),]) * 
                        min_max_norm(data[which(rownames(data) == rec),])), 
                     col.name = paste0(lig, "_", rec))
p <- SpatialFeaturePlot(slide, features = paste0(lig, "_", rec), 
                        images = "image", 
                        stroke = NA, 
                        alpha = c(0,1),
                        min.cutoff = "q20",
                        max.cutoff = "q90",
                        pt.size.factor = 2.8,
                        image.alpha = 1) + 
  scale_fill_gradient2(mid = "red", high = "red") + NoLegend() + 
  ggtitle(paste0(lig, "_", rec))
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/spatial")
ggsave(paste0(name,"_", lig, "_", rec, ".pdf"), plot = as.ggplot(p), width = 3, height = 3)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ***Bassez et al. dataset analysis --------
# Load in Seurat object -------------------------------------------------------

cancer.epi <- readRDS("Bassez_withGEmetadata_112022.rds")
DefaultAssay(cancer.epi) <- "RNA"

# subset to cancer epithelial cells based on inferCNV results
cancer.epi$totalCNV <- as.numeric(cancer.epi$totalCNV)
cancer.epi$corCNV <- as.numeric(cancer.epi$corCNV)
cancer.epi$cancer <- NA
cancer.epi$cancer[which((cancer.epi$totalCNV < 0.02) & (cancer.epi$corCNV < 0.4))] <- "N"
cancer.epi$cancer[which((cancer.epi$totalCNV >= 0.02) | (cancer.epi$corCNV >= 0.4))] <- "Y"
cancer.epi <- subset(cancer.epi, subset = cancer == "Y")

# Scale + UMAP ---------

DefaultAssay(cancer.epi) <- "RNA"
cancer.epi <- NormalizeData(cancer.epi, assay = "RNA")
all.genes <- rownames(cancer.epi)
cancer.epi <- ScaleData(cancer.epi, features = all.genes)

cancer.epi[["RNA"]]@meta.features <- data.frame(row.names = rownames(cancer.epi[["RNA"]]))

cancer.epi <- FindVariableFeatures(cancer.epi, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)
cancer.epi <- RunPCA(cancer.epi, npcs = 100, verbose = TRUE)

cancer.epi <- FindNeighbors(cancer.epi, reduction = "pca", dims = 1:50)

resolution.range <- c(0.01, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.6, 1.8, 2.0)
cancer.epi <- FindClusters(cancer.epi, 
                           graph.name = "integrated_snn", 
                           resolution = resolution.range)

cancer.epi <- FindNeighbors(cancer.epi, reduction = "pca", dims = 1:50)
cancer.epi <- FindClusters(cancer.epi, resolution = resolution.range)
cancer.epi <- RunUMAP(cancer.epi, reduction = "pca", dims = 1:50, verbose = TRUE, seed.use=123)

p <- DimPlot(cancer.epi, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE) + SeuratAxes() + NoLegend()
ggsave("Bassez_UMAP.pdf", plot = as.ggplot(p), width = 5.8, height = 5.5)

p <- DimPlot(cancer.epi, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "BC.Subtype") + SeuratAxes()
ggsave("Bassez_UMAP_BCsubtype.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

p <- DimPlot(cancer.epi, reduction = "umap", 
             label = F, 
             repel = TRUE, 
             raster = FALSE, 
             group.by = "sc50.Pred") + SeuratAxes()
ggsave("Bassez_UMAP_sc50.pdf", plot = as.ggplot(p), width = 7, height = 5.5)

p <- DimPlot(cancer.epi, reduction = "umap", 
             cols = c("#009933", "grey", "#cc0000"), 
             label = FALSE, 
             raster = FALSE, 
             group.by = "Expansion") + 
  scale_fill_discrete(labels = c("Expansion", "N/A", "No Expansion")) + 
  SeuratAxes()
ggsave("Bassez_UMAP_expansion.pdf", plot = as.ggplot(p), width = 6.0, height = 5.5)

# Add clonotype expansion results from Bassez et al. ---------------------------------------------------------

all <- cancer.epi

naive_clonotypes <- read.csv("1881-BIOKEY_clonotypes_combined_cohort1.csv")
naive_clonotypes <- tidyr::separate(data = naive_clonotypes, col = clonotype_id, into = c("BIOKEY", "Patient", "Tx", "Clonotype"))
naive_clonotypes$Patient <- paste(naive_clonotypes$BIOKEY, naive_clonotypes$Patient, sep="_")
naive_clonotypes <- naive_clonotypes[,-1]

chemo_clonotypes <- read.csv("1882-BIOKEY_clonotypes_combined_cohort2.csv")
chemo_clonotypes <- separate(data = chemo_clonotypes, col = clonotype_id, into = c("BIOKEY", "Patient", "Tx", "Clonotype"))
chemo_clonotypes$Patient <- paste(chemo_clonotypes$BIOKEY, chemo_clonotypes$Patient, sep="_")
chemo_clonotypes <- chemo_clonotypes[,-1]

clonotypes <- rbind(naive_clonotypes, chemo_clonotypes)
clonotype_pre <- clonotypes[which(clonotypes$Tx == "Pre"),]
clonotype_on <- clonotypes[which(clonotypes$Tx == "On"),]
clonotype_on <- clonotypes[which(clonotypes$frequency > 2),]
clonotype_expansion <- left_join(clonotype_on, clonotype_pre, by = c("Patient" = "Patient", "Clonotype" = "Clonotype"))
clonotype_expansion$freq.diff <- clonotype_expansion$frequency.x - clonotype_expansion$frequency.y
clonotype_expansion$prop.diff <- clonotype_expansion$proportion.x - clonotype_expansion$proportion.y

clonotype_expansion$expand <- "NE"
clonotype_expansion$expand[which(clonotype_expansion$freq.diff > 0)] <- "E"
clonotype_expansion$expand[which(clonotype_expansion$prop.diff > 0)] <- "E"


# create dataframe for bar graph analysis
clonotypelists <- clonotype_expansion %>% dplyr::group_by(Patient, expand) %>% # group cells by samples and BC.Subtype
  dplyr::summarise(Nb = n()) %>% # add column with number of cells by sc50.Pred
  dplyr::mutate(C = sum(Nb)) %>% # add column with total number of cells in the sample
  dplyr::mutate(percent = Nb/C*100) # add column with % of cells by sc50.Pred out of all cells in sample

clonotypelists <- clonotypelists[which(clonotypelists$expand == "E"),]

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/071322/antiPD1/pre-tx/simil_scores")

scores_pretx <- read.csv("antiPD1_pretx_corr.csv", stringsAsFactors = FALSE, na.strings = "unknown")[,-1]
scores_ontx <- read.csv("antiPD1_pretx_corr.csv", stringsAsFactors = FALSE, na.strings = "unknown")

scores <- left_join(scores_pretx, scores_ontx, by = "name")
scores <- scores[, c(1, 4, 11)]
scores[,2] <- apply(t(scores[,2]), 1, function(x)(1-(x-min(x))/(max(x)-min(x))))
scores[,3] <- apply(t(scores[,3]), 1, function(x)(1-(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))))
scores$change <- scores[, 3] - scores[, 2]
colnames(scores) <- c("Patient", "Pretx_corr", "Ontx_corr", "change")

correlation <- left_join(scores, clonotypelists, by = "Patient")
correlation$expand <- "NE"
correlation$expand[which(correlation$Nb > 30)] <- "E"

# scatterplot (with linear regression)
p <- ggscatter(correlation, x = "Nb", y = "Pretx_corr",
               add = "reg.line",
               conf.int = TRUE,
               add.params = list(color = "blue", fill = "lightgray"),
               xlab = "# clonotypes expanded (antiPD1 response)",
               ylab = "Cancer Cell ITH") +
  stat_cor(method = "pearson", label.x = 100, label.y = -0.05) #+
# scale_y_continuous(breaks = c(-0.4, -0.05),
#                    labels = c("Low", "High"),
#                    limits = c(-0.4, -0.05))

correlation_ordered <- correlation[rev(order(correlation$expand)),]
correlation_ordered[which(correlation_ordered$expand == "E"), 5] <- "Response"
correlation_ordered[which(correlation_ordered$expand == "NE"), 5] <- "No Response"
correlation_ordered <- correlation_ordered[which(correlation_ordered$Patient %in% samples),]
correlation_ordered$expand <- factor(correlation_ordered$expand,
                                     levels = c("No Response", "Response"))

# violin (with Wilcox test)
p <- ggplot(correlation_ordered, aes(x = expand, y = Pretx_corr, fill = expand)) + 
  #geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.5) +
  xlab("") + 
  ylab("") + 
  theme_bw() + 
  scale_y_continuous(breaks = c(0,1), 
                     labels = c("Low ITH", "High ITH"),
                     limits = c(0, 1.2)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual("Response",
                    values = c("#EBB8BC","#B8EBDA"),
                    labels = c("Y (expansion)", "N (no expansion)")) + 
  NoLegend() +
  # scale_y_continuous(breaks = c(-0.5, 0.05), 
  #                    labels = c("Low", "High"),
  #                    limits = c(-0.5, 0.05)) +
  stat_compare_means(comparisons = list(c("No Response", "Response")), 
                     method="wilcox.test", label="p.format", color="black")

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/071322/antiPD1/pre-tx/simil_scores")
ggsave("pretx_ITH_v_response_violin.pdf", 
       plot = p, 
       width = 3.5, height = 3)

# Add GE scores ------------

#cancer.epi@meta.data <- cancer.epi@meta.data[,-c(57:7367)]

cancer.epi@meta.data <- cancer.epi@meta.data[,-c(75:97)]

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- split(GElist$gene,GElist$GE)
length(GElist) # get number of GEs

cancer.epi <- AddModuleScore_UCell(cancer.epi,
                                   features = GElist,
                                   assay = "RNA")

colnames(cancer.epi@meta.data)[75:84] <- c("GE1", "GE2", "GE3", "GE4", "GE5", "GE6", 
                                           "GE7", "GE8", "GE9", "GE10")

cancer.epi <- AddModuleScore_UCell(cancer.epi, 
                                   features = "CD274", 
                                   assay = "RNA")
colnames(cancer.epi@meta.data)[85] <- c("CD274")

colnames(cancer.epi@meta.data)[75:85] <- paste0("raw_", colnames(cancer.epi@meta.data)[75:85])

expdata <- t(cancer.epi@meta.data[,c(1,5,8,75:85)])

zscore <- expdata[-c(1,2,3,14),]
zscore <- t(as.matrix(apply(zscore, 1, as.numeric)))
zscore <- t(apply(zscore, 1, function(x) (x-mean(x))/sd(x)))
newmetadata <- rbind(apply(zscore, 2, function(x) which.max(x)),
                     apply(zscore, 2, function(x) max(x)), 
                     zscore)
newmetadata <- as.data.frame(t(newmetadata))
colnames(newmetadata) <- c("maxGE", "maxZscore", "GE1", "GE2", "GE3", "GE4", "GE5", 
                           "GE6", "GE7", "GE8", "GE9", "GE10")
rownames(newmetadata) <- colnames(expdata)

cancer.epi <- AddMetaData(cancer.epi, newmetadata, col.name = colnames(newmetadata))

#saveRDS(cancer.epi, "Bassez_withGEmetadata_112022.rds")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ***T cell InteractPrint analysis --------------
# Read in CellChat/Nichenet output ------

set.seed(123)
setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/Nichenet")
receptor_sum <- read.csv("cellchat_nichenet_receptors_sum.csv")
nichenet_receptors <- read.csv("cellchat_nichenet_receptors_list.csv")

colnames(receptor_sum)[1] <- "cell"
receptor_sum <- receptor_sum[,c(1,2,4,5,6,7,8,9,10,11,3)]

# Incorporate negative multiplier for inhibitory GEs based on T cell checkpoints ------

T_checkpoint <- c("HAVCR2", "LAG3", "PDCD1", "CTLA4", "TIGIT", "BTLA", "KLRC1", 
                  "KLRG1")#, "TGFBR2")
inhib_GEs <- nichenet_receptors[which((nichenet_receptors$receptor %in% T_checkpoint) & (nichenet_receptors$cell == "CD8+ T Cells")),]
inhib_GEs <- unique(inhib_GEs$GE)
print(inhib_GEs)

for (i in inhib_GEs) {
  receptor_sum[3, i] <- (-1)*receptor_sum[3, i]
}

# ----------
# ----------
# Fig 4C: Bassez et al. GE and InteractPrint scores -------
# Loading in Bassez et al. dataset ---------
correlation <- cancer.epi@meta.data[,c(1,5:8,54:55,88:97,85:86)]
correlation <- correlation[which(correlation$Expansion != "n/a"),]
correlation <- correlation[which(correlation$Patient %in% samples),]

correlation_GE <- aggregate(. ~ Patient, correlation[,c(2,8:18)], sum)
correlation_num <- aggregate(. ~ Patient, correlation[,c(2,8:18)], length)
correlation_GE[,c(2:12)] <- correlation_GE[,c(2:12)] / correlation_num[,c(2:12)]

correlation_expansion <- aggregate(Expansion ~ Patient, correlation, unique)
correlation_GE <- left_join(correlation_GE, correlation_expansion, by = "Patient")

# split samples by responder vs. non-responder (based on clonotype expansion)
correlation_GE[which(correlation_GE$Expansion == "E"), 13] <- "Response"
correlation_GE[which(correlation_GE$Expansion == "NE"), 13] <- "No Response"
correlation_GE$Expansion <- factor(correlation_GE$Expansion,
                                   levels = c("No Response", "Response"))
correlation_GE <- left_join(correlation_GE, unique(correlation[,c(2,5,7)]), by = "Patient")

rownames(correlation_GE) <- correlation_GE$Patient
correlation_GE <- correlation_GE[order(correlation_GE$Expansion),]

# create annotation for heatmap
anno <- HeatmapAnnotation(expand = correlation_GE$Expansion, 
                          subtype = correlation_GE$BC.Subtype,
                          col = list(expand = c("Response" = "dark green", 
                                                "No Response" = "dark red"),
                                     subtype = c("HER2+" = "#B0FC96",
                                                 "HR+" = "#96B0FC",
                                                 "TNBC" = "#FC96B0"),
                                     subtype = c("HER2+" = "#700639", "HR+" = "#397006", "TNBC" = "#063970")#,
                          ))

# create heatmap of GEs by sample in Bassez et al. dataset
p <- Heatmap(t(correlation_GE[,c(2:11)]), 
             column_order = rownames(correlation_GE), 
             top_annotation = anno, 
             show_column_names = F,
             col = viridis(500),
             column_split = correlation_GE$BC.Subtype, 
             cluster_row_slices = F)

ggsave("Bassez_heatmap_flex_prePD1_nochemo.pdf", 
       plot = as.ggplot(p), 
       width = 10, 
       height = 3.5)

# Add on T cell InteractPrint score for Bassez et al. samples ---------

correlation <- correlation_GE
correlation_immune <- as.data.frame(as.matrix(correlation[,c(2:11)]) %*% t(as.matrix(receptor_sum[,-1])))
colnames(correlation_immune) <- receptor_sum$cell

# label samples based on anti-PD-1 therapy response (clonal expansion) 
correlation_immune$samples <- correlation$Patient
correlation_immune$Expansion <- correlation$Expansion
correlation_immune$BC.Subtype <- correlation$BC.Subtype
correlation_immune <- aggregate(. ~ BC.Subtype + Expansion + samples, correlation_immune, mean)
correlation_immune <- correlation_immune[which(correlation_immune$samples %in% samples),]
correlation_immune <- correlation_immune[order(correlation_immune$`CD8+ T Cells`),]
correlation_immune[which(correlation_immune$Expansion == "E"), 2] <- "Response"
correlation_immune[which(correlation_immune$Expansion == "NE"), 2] <- "No Response"
correlation_immune <- correlation_immune[order(correlation_immune$Expansion),]
rownames(correlation_immune) <- correlation_immune$samples

# create heatmap annotation
anno <- HeatmapAnnotation(expand = correlation_immune$Expansion, 
                          subtype = correlation_immune$BC.Subtype,
                          col = list(expand = c("Response" = "dark green", 
                                                "No Response" = "dark red"),
                                     subtype = c("HER2+" = "#B0FC96",
                                                 "HR+" = "#96B0FC",
                                                 "TNBC" = "#FC96B0")))

colnames(correlation_immune)[3] <- "Patient"
correlation_all <- left_join(correlation_GE, correlation_immune[,c(3,6)], by = "Patient")
correlation_all <- correlation_all[order(correlation_all$`CD8+ T Cells`,decreasing = F),] 
correlation_all <- correlation_all[order(correlation_all$Expansion,decreasing = F),] 
correlation_all <- correlation_all[order(correlation_all$BC.Subtype,decreasing = F),] 

# Create heatmap of GE scores for all samples in Bassez et al. ------
p <- Heatmap(t(correlation_all[,-c(1,12:16)]), 
             column_order = rownames(correlation_all), 
             top_annotation = anno, 
             col = viridis(500),
             show_column_names = F,
             column_split = correlation_all$BC.Subtype)

ggsave("Bassez_heatmap_combined1.pdf", plot = as.ggplot(p), width = 11, height = 4)

# Create heatmap of T cell InteractPrint scores for all samples in Bassez et al. -----

weighted <- as.matrix(t(correlation_all[,-c(1:15)]))
colnames(weighted) <- rownames(correlation_all)
p <- Heatmap(weighted, 
             column_order = rownames(correlation_all), 
             #top_annotation = anno, 
             col = mako(500),
             show_column_names = F,
             #row_split = 8,
             column_split = correlation_all$BC.Subtype)

ggsave("Bassez_heatmap_combined2.pdf", plot = as.ggplot(p), width = 11, height = 1)

# -------
# -------
# Fig 4D: Bassez et al. T cell InteractPrint prediction of response to anti-PD-1 therapy --------
# Correlation plot between response (expansion) and T cell InteractPrint score -------

p <- ggplot(correlation_immune, 
            aes(x = Expansion, y = `CD8+ T Cells`, fill = Expansion)) + 
  geom_boxplot(width = 0.5) +
  geom_point(position=position_dodge(width = 0.75)) +
  xlab("") + ylab("") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual("Response",
                    values = c("#EBB8BC", "#B8EBDA"),
                    labels = c("Y (expansion)", "N (no expansion)")) + 
  NoLegend() +
  stat_compare_means(comparisons = list(c("Response", "No Response")), 
                     method="wilcox.test", label="p.format", color="black")

ggsave("Bassez_correlation_weightedT.pdf", plot = as.ggplot(p), width = 3, height = 4)

# Create ROC curve for T cell InteractPrint on Bassez et al. dataset -------------

library(pROC)

pdf("Bassez_ROC_weightedT.pdf", height = 5, width = 5)
pROC_obj <- roc(correlation_all$Expansion,
                correlation_all$`CD8+ T Cells`,
                smoothed = T,
                # arguments for ci
                ci = T, ci.alpha = 0.9, stratified = T,
                # arguments for plot
                plot = TRUE, 
                auc.polygon = T, 
                max.auc.polygon = T, 
                grid = F,
                print.auc = F, 
                show.thres = F)

sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type = "shape", col = NULL)
text(0.4, 0.4, paste("AUC:", round(pROC_obj$auc, 3)))
dev.off()

pdf("Bassez_ROC_weightedTvsPDL1.pdf", width = 5, height = 5)
rocobj1 <- plot.roc(correlation_immune$Expansion,
                    correlation_immune$`CD8+ T Cells`,
                    percent=TRUE,
                    col="#1c61b6")
rocobj2 <- plot.roc(correlation_GE$Expansion,
                    correlation_GE$raw_CD274,
                    percent=TRUE, 
                    col="#008600")

testobj <- roc.test(rocobj1, rocobj2, method="bootstrap", boot.n=10000)
text(60, 30, labels=paste("p-value =", format.pval(testobj$p.value)), adj=c(0, .5))

text(85, 95, paste("AUC:", round(rocobj1$auc, 3)),
     col="#1c61b6")
text(85, 90, paste("AUC:", round(rocobj2$auc, 3)),
     col="#008600")
legend("bottomright", legend=c("Weighted CD8+ T Cell \nInteraction Score", "Average PD-L1 Expression"), col=c("#1c61b6", "#008600"), lwd=2)
dev.off()

print(testobj$p.value)
print(rocobj1$auc)

# -------
# -------
# Fig 4E: I-SPY2 GE and InteractPrint scores -----------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
data <- read.csv("ISPY2_markerbasedestimates_weighted.csv", row.names = 1)

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/validation")
# Loading in ISPY
correlation <- data[,c(16,12,17:18,22,19,65,1:10,11)]
correlation[,c(8:17)] <- t(apply(correlation[,c(8:17)], 1, function(x) (x-mean(x))/sd(x)))
colnames(correlation)[2] <- "Patient"
colnames(correlation)[1] <- "Expansion"
colnames(correlation)[5] <- "BC.Subtype"
colnames(correlation)[18] <- "raw_CD274"
correlation[which(correlation$BC.Subtype == "Basal-type"), 5] <- "TNBC"
correlation[which(correlation$BC.Subtype == "Luminal-type"), 5] <- "HR+"

correlation_GE <- aggregate(. ~ Patient, correlation[,c(2,8:18)], sum)
correlation_num <- aggregate(. ~ Patient, correlation[,c(2,8:18)], length)
correlation_GE[,c(2:12)] <- correlation_GE[,c(2:12)] / correlation_num[,c(2:12)]

correlation_expansion <- aggregate(Expansion ~ Patient, correlation, unique)
correlation_GE <- left_join(correlation_GE, correlation_expansion, by = "Patient")
correlation_GE[which(correlation_GE$Expansion == "R"), 13] <- "Response"
correlation_GE[which(correlation_GE$Expansion == "NR"), 13] <- "No Response"
correlation_GE$Expansion <- factor(correlation_GE$Expansion,
                                   levels = c("No Response", "Response"))
correlation_GE <- left_join(correlation_GE, unique(correlation[,c(2,5,7)]), by = "Patient")

p <- ggplot(correlation_GE, 
            aes(x = Expansion, y = raw_CD274, fill = Expansion)) + 
  geom_boxplot(width = 0.5) +
  geom_point(position=position_dodge(width=0.75)) +
  xlab("") + 
  ylab("") + #ylim(-0.02,0.05) +
  theme_bw() + 
  # scale_y_continuous(breaks = c(0,1), 
  #                    labels = c("Low", "High")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual("Response",
                    values = c("#EBB8BC", "#B8EBDA"),
                    labels = c("Y (expansion)", "N (no expansion)")) + 
  NoLegend() +
  stat_compare_means(comparisons = list(c("No Response", "Response")), 
                     method="wilcox.test", label="p.format", color="black")

ggsave(paste0("ISPY2_PDL1_response.pdf"),
       plot = p,
       width = 4, height = 4)

rownames(correlation_GE) <- correlation_GE$Patient
correlation_GE <- correlation_GE[order(correlation_GE$Expansion),]

anno <- HeatmapAnnotation(expand = correlation_GE$Expansion, 
                          subtype = correlation_GE$BC.Subtype,
                          #tx = correlation_immune$Treatment.Status,
                          col = list(expand = c("Response" = "dark green", 
                                                "No Response" = "dark red"),
                                     subtype = c("HER2+" = "#B0FC96",
                                                 "HR+" = "#96B0FC",
                                                 "TNBC" = "#FC96B0"),
                                     subtype = c("HER2+" = "#700639", "HR+" = "#397006", "TNBC" = "#063970")#,
                                     #tx = c("Naive (Pre antiPD1)" = "light blue"#, 
                                     # "Neoadjuvant_chemo (Pre antiPD1)" = "grey")
                          ))

p <- Heatmap(t(correlation_GE[,c(2:11)]), 
             column_order = rownames(correlation_GE), 
             top_annotation = anno, 
             show_column_names = F,
             col = viridis(500),
             column_split = correlation_GE$BC.Subtype, 
             cluster_row_slices = F)

ggsave("ISPY2_heatmap_flex_prePD1_nochemo.pdf", plot = as.ggplot(p), width = 20, height = 3.5)

# Loading in I-SPY2 dataset --------
correlation <- data[,c(12,16,22,65,17:21,1:10)]
correlation[,c(10:19)] <- t(apply(correlation[,c(10:19)], 1, function(x) (x-mean(x))/sd(x)))

# Add on T cell InteractPrint score for I-SPY2 samples ------
correlation_immune <- as.data.frame(as.matrix(correlation[,c(10:19)]) %*% t(as.matrix(receptor_sum[,-1])))
colnames(correlation_immune) <- receptor_sum$cell
correlation_immune$samples <- correlation$Patient.Identifier
correlation_immune$Expansion <- correlation$pCR
correlation_immune$BC.Subtype <- correlation$BP.subtype
rm(mean)

correlation_immune <- aggregate(. ~ BC.Subtype + Expansion + samples, correlation_immune, mean)
correlation_immune <- correlation_immune[order(correlation_immune$`CD8+ T Cells`),]
correlation_immune[which(correlation_immune$Expansion == "R"), 2] <- "Response"
correlation_immune[which(correlation_immune$Expansion == "NR"), 2] <- "No Response"
correlation_immune[which(correlation_immune$BC.Subtype == "Basal-type"), 1] <- "TNBC"
correlation_immune[which(correlation_immune$BC.Subtype == "Luminal-type"), 1] <- "HR+"

correlation_immune <- correlation_immune[order(correlation_immune$Expansion),]
tx.status <- unique(correlation[,c(2,3)])
colnames(tx.status) <- c("samples", "Treatment.Status")
correlation_immune <- left_join(correlation_immune, tx.status, by = "samples")
rownames(correlation_immune) <- correlation_immune$samples

colnames(correlation_immune)[3] <- "Patient"
correlation_all <- left_join(correlation_GE, correlation_immune[,c(3,6)], by = "Patient")
correlation_all <- correlation_all[order(correlation_all$`CD8+ T Cells`,decreasing = F),] 
correlation_all <- correlation_all[order(correlation_all$Expansion,decreasing = F),] 
correlation_all <- correlation_all[order(correlation_all$BC.Subtype,decreasing = F),] 

# create heatmap annotation
anno <- HeatmapAnnotation(expand = correlation_all$Expansion, 
                          subtype = correlation_all$BC.Subtype,
                          col = list(expand = c("Response" = "dark green", 
                                                "No Response" = "dark red"),
                                     subtype = c("HER2+" = "#B0FC96",
                                                 "HR+" = "#96B0FC",
                                                 "TNBC" = "#FC96B0")))

# Create heatmap of GE scores for all samples in I-SPY2 -------
p <- Heatmap(t(correlation_all[,-c(1,12:16)]), 
             column_order = rownames(correlation_all), 
             top_annotation = anno, 
             col = viridis(500),
             show_column_names = F,
             column_split = correlation_all$BC.Subtype)
ggsave("ISPY2_heatmap_combined1.pdf", plot = as.ggplot(p), width = 20, height = 4)

# Create heatmap of T cell InteractPrint scores for all samples in I-SPY2 -----
weighted <- as.matrix(t(correlation_all[,-c(1:15)]))
colnames(weighted) <- rownames(correlation_all)
p <- Heatmap(weighted, 
             column_order = rownames(correlation_all), 
             #top_annotation = anno, 
             col = mako(500),
             show_column_names = F,
             #row_split = 8,
             column_split = correlation_all$BC.Subtype)

ggsave("ISPY2_heatmap_combined2.pdf", plot = as.ggplot(p), width = 20, height = 1)

# -------
# -------
# Fig 4F: I-SPY2 T cell InteractPrint prediction of response to anti-PD-1 therapy --------
# Correlation plot between response (expansion) and T cell InteractPrint score -----

p <- ggplot(correlation_immune, 
            aes(x = Expansion, y = `CD8+ T Cells`, fill = Expansion)) + 
  geom_boxplot(width = 0.5) +
  geom_point(position=position_dodge(width = 0.75)) +
  xlab("") + ylab("") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual("Response",
                    values = c("#EBB8BC", "#B8EBDA"),
                    labels = c("Y (expansion)", "N (no expansion)")) + 
  NoLegend() +
  stat_compare_means(comparisons = list(c("Response", "No Response")), 
                     method="wilcox.test", label="p.format", color="black")

ggsave("ISPY2_correlation_weightedT.pdf", plot = as.ggplot(p), width = 3, height = 4)

# Create ROC curve for T cell InteractPrint on I-SPY2 dataset ------

library(pROC)

correlation_new <- correlation_all

pdf("IPSY2_IP.pdf", height = 5, width = 5)
pROC_obj <- roc(correlation_new$Expansion,
                correlation_new$`CD8+ T Cells`,
                smoothed = T,
                # arguments for ci
                ci = T, ci.alpha = 0.9, stratified = T,
                # arguments for plot
                plot = TRUE, 
                auc.polygon = T, 
                max.auc.polygon = T, 
                grid = F,
                print.auc = F, 
                show.thres = F)

sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type = "shape", col = NULL)
text(0.4, 0.4, paste("AUC:", round(pROC_obj$auc, 3)))
dev.off()

pdf("ISPY2_IPvsPDL1.pdf", width = 5, height = 5)
rocobj1 <- plot.roc(correlation_new$Expansion,
                    correlation_new$`CD8+ T Cells`,
                    percent=TRUE,
                    col="#1c61b6")
rocobj2 <- plot.roc(correlation_new$Expansion,
                    correlation_new$raw_CD274,
                    percent=TRUE, 
                    col="#008600")
testobj <- roc.test(rocobj1, rocobj2, method="bootstrap", boot.n=10000)
text(60, 30, labels=paste("p-value =", format.pval(testobj$p.value)), adj=c(0, .5))

text(85, 95, paste("AUC:", round(rocobj1$auc, 3)),
     col="#1c61b6")
text(85, 90, paste("AUC:", round(rocobj2$auc, 3)),
     col="#008600")
legend("bottomright", legend=c("Weighted CD8+ T Cell \nInteraction Score", "Average PD-L1 Expression"), col=c("#1c61b6", "#008600"), lwd=2)
dev.off()

print(testobj$p.value)
print(rocobj1$auc)

# -------
# -------
