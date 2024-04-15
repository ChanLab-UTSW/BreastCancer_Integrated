# SETUP --------
# Libraries -------------------------------------------------------------------

library(generics, lib.loc = "/home2/s437775/R/x86_64-pc-linux-gnu-library/4.1")
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(ggplotify)
library(UCell)
library(ComplexHeatmap)
library(circlize)

# ------
# ------

# VISIUM -----------------------------------------------------
# Create spatial datasets from Wu 2021 ------------------------------------------------------------

# load in GEs
GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist <- split(GElist$gene,GElist$GE)
length(GElist) # get number of GEs

# initialize objects for spatial data samples
samples <- c("CID4535", "CID44971", "CID4465", "CID4290", "1160920F", "1142243F")
slide_all <- list()

# loop through each sample to load in objects
for (i in samples) {

  # read in filtered counts matrix
  counts <- Read10X(
    paste0("/work/InternalMedicine/s437775/simil/simil_cancerepi/spatial_Wu2021/filtered_count_matrices/",
           i,
           "_filtered_count_matrix"),
    gene.column = 1,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
  )

  # read in metadata
  setwd("/work/InternalMedicine/s437775/simil/simil_cancerepi/spatial_Wu2021/metadata")
  meta <- read.csv(paste0(i, "_metadata.csv"))

  # create Seurat spatial object
  slide = CreateSeuratObject(counts = counts,
                             project = i,
                             assay = "Spatial")
  slide$slide <- 1
  slide$region <- i

  # read in image data
  image <- Read10X_Image(
    paste0("/work/InternalMedicine/s437775/simil/simil_cancerepi/spatial_Wu2021/spatial/",
           i,
           "_spatial"),
    image.name = "tissue_lowres_image.png",
    filter.matrix = TRUE)
  DefaultAssay(object = image) <- "Spatial"
  image <- image[colnames(x = slide)]
  slide[["image"]] <- image

  slide@meta.data$X <- rownames(slide@meta.data)
  slide@meta.data <- left_join(slide@meta.data, meta, by = "X")
  rownames(slide@meta.data) <- slide$X
  colnames(slide@meta.data) <- c("orig.ident", "nCount_Spatial", "nFeature_Spatial",
                                 "slide", "region", "cell", "nCount_RNA", "nFeature_RNA",
                                 "subtype", "patientid", "Classification")

  # plot spatial slide by class
  Idents(slide) <- slide$Classification
  plot2 <- SpatialDimPlot(slide, stroke = NA, image.alpha = 0.7) + theme(legend.position = "right")
  plot2 <- as.ggplot(plot2)
  ggsave(paste0(i, "_classes_all.pdf"), plot = plot2)


  idents <- setdiff(unique(slide$Classification), c(NA, "Artefact", "Uncertain", "", "Adipose tissue",
                                                    "Normal + stroma + lymphocytes", "Normal duct",
                                                    "Normal glands + lymphocytes", "Necrosis"))
  slide <- subset(slide, subset = Classification %in% idents)

  # perform SCT transform, dim reduction, clustering
  Idents(slide) <- slide$cell
  slide <- SCTransform(slide,
                       assay = "Spatial",
                       verbose = FALSE)
  slide <- RunPCA(slide, assay = "SCT", verbose = F)
  slide <- FindNeighbors(slide, reduction = "pca", dims = 1:50)
  slide <- FindClusters(slide, verbose = F)
  slide <- RunUMAP(slide, reduction = "pca", dims = 1:50)

  # plot nCount violin plot
  plot1 <- VlnPlot(slide, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot1 <- as.ggplot(plot1)
  ggsave(paste0(i, "_nCount.pdf"), plot = plot1)

  # plot spatial slide by class
  Idents(slide) <- slide$Classification
  plot2 <- SpatialDimPlot(slide, stroke = NA, image.alpha = 0.7) + theme(legend.position = "right")
  plot2 <- as.ggplot(plot2)
  ggsave(paste0(i, "_classes.pdf"), plot = plot2)

  # plot UMAP
  p <- DimPlot(slide, reduction = "umap",
               label = F,
               repel = TRUE,
               raster = FALSE,
               group.by = "Classification") + SeuratAxes()
  ggsave(paste0(i, "_UMAP.pdf"),
         plot = as.ggplot(p),
         width = 8.5, height = 5.5)
         
  # append slide to list
  slide_all <- append(slide_all, slide)
}

names(slide_all) <- samples
saveRDS(slide_all, "spatial_Wu2021_objlist_111222.rds")

# Create spatial datasets from Visium datasets ---------

# initialize objects for spatial data samples
samples <- c("CytAssist_FFPE_Human_Breast_Cancer", 
             "Parent_Visium_Human_BreastCancer",
             "V1_Breast_Cancer_Block_A_Section_1", 
             "V1_Human_Invasive_Ductal_Carcinoma",
             "Visium_FFPE_Human_Breast_Cancer", 
             "Visium_Human_Breast_Cancer")
slide_all_10x <- list()

# loop through each sample to load in objects
for (i in samples) {
  
  setwd(paste0("/work/InternalMedicine/s437775/simil/simil_cancerepi/spatial_Visium/",i))
  untar(paste0(i, "_spatial.tar.gz"))
  
  # read in image data
  image <- Read10X_Image(
    paste0("/work/InternalMedicine/s437775/simil/simil_cancerepi/spatial_Visium/", i, "/spatial"),
    image.name = "tissue_lowres_image.png",
    filter.matrix = TRUE)
  
  # read in filtered counts matrix
  slide <- Load10X_Spatial(
    paste0("/work/InternalMedicine/s437775/simil/simil_cancerepi/spatial_Visium/", i),
    filename = paste0(i, "_filtered_feature_bc_matrix.h5"),
    assay = "Spatial",
    slice = "slice1")
  
  DefaultAssay(object = image) <- "Spatial"
  image <- image[colnames(x = slide)]
  slide[["image"]] <- image
  
  slide$patientid <- i
  
  # add metadata 
  if (i == samples[1]) {
    slide$subtype <- "HR+"
    slide$stage <- "T2N1M0"
    slide$full_subtype <- "ER+/HER2+/PR-"
  }
  if (i == samples[2]) {
    slide$subtype <- "HR+"
    slide$stage <- "I"
    slide$full_subtype <- "ER+/HER2-/PR+"
  }
  if (i == samples[3]) {
    slide$subtype <- "HER2+"
    slide$stage <- "IIa"
    slide$full_subtype <- "ER+/HER2+/PR-"
  }
  if (i == samples[4]) {
    slide$subtype <- "HER2+"
    slide$stage <- "IIa"
    slide$full_subtype <- "ER+/HER2+/PR-"
  }
  if (i == samples[5]) {
    
  }
  if (i == samples[6]) {
    slide$subtype <- "HER2+" 
    slide$stage <- "T2N0M0"
    slide$full_subtype <- "ER+/HER2+/PR2"
  }
  
  # perform SCT transform, dim reduction, clustering
  slide <- SCTransform(slide,
                       assay = "Spatial",
                       verbose = FALSE)
  slide <- RunPCA(slide, assay = "SCT", verbose = F)
  slide <- FindNeighbors(slide, reduction = "pca", dims = 1:50)
  slide <- FindClusters(slide, verbose = F)
  slide <- RunUMAP(slide, reduction = "pca", dims = 1:50)
  
  # plot nCount violin plot
  plot1 <- VlnPlot(slide, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot1 <- as.ggplot(plot1)
  ggsave(paste0(i, "_nCount.pdf"), plot = plot1)
  
  # plot UMAP
  p <- DimPlot(slide, reduction = "umap",
               label = F,
               repel = TRUE,
               raster = FALSE) + SeuratAxes()
  ggsave(paste0(i, "_UMAP.pdf"),
         plot = as.ggplot(p),
         width = 8.5, height = 5.5)
  
  # append slide to list
  slide_all_10x <- append(slide_all_10x, slide)
}

names(slide_all_10x) <- samples
saveRDS(slide_all_10x, "spatial_Visium_objlist_111222.rds")

# Create spatial datasets from Andersson (did not use) ------- 

setwd("/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/spatial")
slide_Andersson <- readRDS("spatial_Andersson_objlist_111222.rds")

# Use Andersson.Rmd file 

# --------
# --------

# Load in datasets -----------

slide_Wu <- readRDS("spatial_Wu2021_objlist_111222.rds")
slide_Andersson <- readRDS("spatial_Andersson_objlist_111222.rds")
slide_Visium <- readRDS("spatial_Visium_objlist_111222.rds")

slide_all <- append(slide_Wu, slide_Andersson)
slide_all <- append(slide_all, slide_Visium)
# slide_all <- append(slide_Wu, slide_Visium)

names(slide_all) <- c(names(slide_all)[1:6],
                      unique(slide_all[[7]]$patient_id),
                      unique(slide_all[[8]]$patient_id),
                      unique(slide_all[[9]]$patient_id),
                      unique(slide_all[[10]]$patient_id),
                      unique(slide_all[[11]]$patient_id),
                      unique(slide_all[[12]]$patient_id),
                      unique(slide_all[[13]]$patient_id),
                      unique(slide_all[[14]]$patient_id),
                      names(slide_all)[15:20])

samples <- c("HR+", "TNBC", "TNBC", "HR+", "TNBC", "TNBC", 
             "HER2+", "HER2+", "HER2+", "HER2+", "HER2+", "HER2+", "HER2+", "HER2+", 
             unique(slide_all[[15]]$subtype),
             unique(slide_all[[16]]$subtype),
             unique(slide_all[[17]]$subtype),
             unique(slide_all[[18]]$subtype),
             NA, #unique(slide_all[[19]]$subtype),
             unique(slide_all[[20]]$subtype))

# Integration with scRNAseq dataset -------------------------------------------

# Load in primary object ---------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
prim_obj <- readRDS("primobj_withCNVlabels_111022.rds")

# Prepare primary object for integration with spatial slides -------

slide_name <- c("CID4535","CID44971","CID4465","CID4290","1160920F","1142243F",
                NA, NA, NA, NA, NA, NA, NA, NA, 
                "CytAssist_FFPE_Human_Breast_Cancer", 
                "Parent_Visium_Human_BreastCancer",
                "V1_Breast_Cancer_Block_A_Section_1", 
                "V1_Human_Invasive_Ductal_Carcinoma",
                NA, 
                "Visium_Human_Breast_Cancer")

prim_obj$celltype_broad <- prim_obj$celltype_final
prim_obj$celltype_broad[which(prim_obj$celltype_broad %in% c("Cancer Epithelial Cells", 
                                                             "Epithelial Cells"))] <- "Epithelial Cells"
prim_obj$celltype_broad[which(prim_obj$celltype_broad %in% c("Perivascular-like (PVL) Cells", 
                                                             "Endothelial Cells",
                                                             "Fibroblasts", 
                                                             "Myoepithelial Cells"))] <- "Stromal Cells"
prim_obj$celltype_broad[which(prim_obj$celltype_broad %in% c("CD4+ T Cells", 
                                                             "CD8+ T Cells",
                                                             "Regulatory T Cells",
                                                             "NK Cells",
                                                             "B Cells", 
                                                             "Plasma Cells"))] <- "Lymphocytes"
prim_obj$celltype_broad[which(prim_obj$celltype_broad %in% c("Dendritic Cells", 
                                                             "Macrophages",
                                                             "Mast Cells", 
                                                             "MDSCs",
                                                             "Monocytes", 
                                                             "Neutrophils"))] <- "Myeloid Cells"

prim_obj_TNBC <- subset(prim_obj, subset = BC.Subtype == "TNBC")
prim_obj_HER2 <- subset(prim_obj, subset = BC.Subtype == "HER2+")
prim_obj_HR <- subset(prim_obj, subset = BC.Subtype == "HR+")

# Initialize objects ------
new_slide_all <- list()
corrs <- list()

# Integrate Wu spatial slides with scRNAseq dataset ------------
for (i in c(1:6)) {
  # Setup slide ------- 
  
  slide <- slide_all[[i]]
  
  if (samples[i] == "HR+") {
    prim_obj_sub <- prim_obj_HR
  }
  if (samples[i] == "HER2+") {
    prim_obj_sub <- prim_obj_HER2
  }
  if (samples[i] == "TNBC") {
    prim_obj_sub <- prim_obj_TNBC
  }
  
  epicutoff <- 0.1
  
  # Perform anchor-based integration workflow ---------
  
  DefaultAssay(prim_obj_sub) <- "integrated"
  anchors <- FindTransferAnchors(reference = prim_obj_sub, 
                                 query = slide, 
                                 normalization.method = "SCT")
  
  DefaultAssay(slide) <- "Spatial"
  slide <- RunPCA(slide, assay = "SCT", verbose = F)
  predictions.assay <- TransferData(anchorset = anchors, 
                                    refdata = prim_obj_sub$celltype_broad, 
                                    prediction.assay = TRUE,
                                    weight.reduction = slide[["pca"]], 
                                    dims = 1:50)
  slide[["predictions"]] <- predictions.assay
  predict <- as.data.frame(t(GetAssayData(slide, assay = "predictions", slot = "data")))
  
  # get prediction scores for each spot and plot classes
  DefaultAssay(slide) <- "predictions"
  setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/spatial")
  
  # Add GE scores ----------------
  
  setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
  GElist <- readxl::read_xlsx("cancer_GE.xlsx")
  GElist <- split(GElist$gene,GElist$GE)
  length(GElist) # get number of GEs
  
  slide <- AddModuleScore_UCell(slide,
                                features = GElist,
                                assay = "Spatial")
  
  colnames(slide@meta.data)[16:25] <- paste0("raw_", c("GE1", "GE2", "GE3", "GE4", "GE5", "GE6", 
                                        "GE7", "GE8", "GE9", "GE10"))
  
  zscore <- apply(slide@meta.data[,c(16:25)], 2, function(x) (x-mean(x))/sd(x))
  newmetadata <- cbind(apply(zscore, 1, function(x) which.max(x)),
                       apply(zscore, 1, function(x) max(x)), 
                       zscore)
  newmetadata <- as.data.frame(newmetadata)
  colnames(newmetadata) <- c("maxGE", "maxZscore", "GE1", "GE2", "GE3", "GE4", "GE5", 
                             "GE6", "GE7", "GE8", "GE9", "GE10")
  rownames(newmetadata) <- colnames(slide)
  
  slide <- AddMetaData(slide, newmetadata, col.name = colnames(newmetadata))
  
  drop <- rownames(predict[which(predict$`Epithelial Cells` < epicutoff), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(16:37)] <- 0
  
  # Add NK cell labels --------
  
  # NK cells
  NK.mark2 <- list(c("CD3D","CD3E","CD3G","FGFBP2","KLRD1", "FCGR3A", "KLRK1", "NCAM1", 
                     "CD4-", "CD8A-", "CD8B-", "FOXP3-"))
  slide <- AddModuleScore_UCell(slide, features = NK.mark2, name = "NK.mark2", assay = "Spatial")
  
  drop <- rownames(predict[which(predict$Lymphocytes < 0.1), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(38)] <- 0
  
  # Add T cell labels -----------

  # CD8+ T Cells
  T_karacd8 <- list(c("CD2","CD3D","CD3E","CD3G","CD8A","CD8B",
                      "CD4-", "FOXP3-"))
  slide <- AddModuleScore_UCell(slide, features = T_karacd8, name = "T_karacd8", assay = "Spatial")

  # CD4+ T Cells
  T_karacd4 <- list(c("CD2","CD3D","CD3E","CD3G","CD4",
                      "CD8A-","CD8B-", "FOXP3-"))
  slide <- AddModuleScore_UCell(slide, features = T_karacd4, name = "T_karacd4", assay = "Spatial")
  
  # T reg Cells
  T_reg <- list(c("CD2","CD3D","CD3E","CD3G","FOXP3"))
  slide <- AddModuleScore_UCell(slide, features = T_reg, name = "T_reg", assay = "Spatial")
  
  drop <- rownames(predict[which(predict$Lymphocytes == 0), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(39:41)] <- 0
  
  # Add B cell labels ---------
  
  # B Cells
  B_KaraBrech <- list(c("CD79A","CD79B","BLNK","CD19","MS4A1")) 
  slide <- AddModuleScore_UCell(slide, features = B_KaraBrech, name = "B_KaraBrech", assay = "Spatial")

  # Plasma Cells
  plasma_colo <- list(c("CD27","IGHA1","SDC1","TNFRSF17","JCHAIN","MZB1","DERL3","CD38",
                        "IGHG1","IGHG3","IGHG4"))
  slide <- AddModuleScore_UCell(slide, features = plasma_colo, name = "plasma_colo", assay = "Spatial")

  drop <- rownames(predict[which(predict$Lymphocytes == 0), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(42:43)] <- 0
  
  # Add myeloid labels --------
  
  # Myeloid Cells
  myeloid_wuElliot <- list(c("PTPRC", "ITGAM","HLA-DR","ITGAX", "CD14","CD16","CD1C",
                             "CD1A","CD68","CD33",
                             "EPCAM-")) 
  slide <- AddModuleScore_UCell(slide, features = myeloid_wuElliot, name = "myeloid_wuElliot", assay = "Spatial")

  # MDSC
  MDSc_sig <- list(c("CD33", "ITGAM", "LY6G","LY6C","FUT4", "CEACAM8","IL4R","HLA-DRA",
                     "CD3D-","CD14-","CD19-","NCAM1-", "EPCAM-")) 
  slide <- AddModuleScore_UCell(slide, features = MDSc_sig, name = "MDSC_sig", assay = "Spatial")
  
  # Macrophages
  Macro_Cheng <- list(c("INHBA", "IL1RN", "CCL4", "NLRP3",
                        "EREG", "IL1B", "LYVE1", "PLTP",
                        "SELENOP", "C1QC", "C1QA", "APOE"))
  slide <- AddModuleScore_UCell(slide, features = Macro_Cheng, name = "Macrophage", assay = "Spatial")
  
  # Monocytes
  Mono_sig <- list(c("FCN1", "S100A9", "S100A8", "FCGR3A", "LST1", "LILRB2", 
                     "HLA-DRA-"))
  slide <- AddModuleScore_UCell(slide, features = Mono_sig, name = "Monocytes", assay = "Spatial")
  
  # Neutrophils 
  Neutro_sig <- list(c("CXCR1","CD15","FCGR3A", 
                       "CD14-",
                       "CSF3R", "S100A9","CD24A","TNF","CD274"))
  slide <- AddModuleScore_UCell(slide, features = Neutro_sig, name = "Neutrophils", assay = "Spatial")
  
  # Dendritic Cells
  DC_sig <- list(c("LILRA4","GZMB","IL3RA","CLEC9A", "FLT3", "IDO1", 
                   "CD1C", "FCER1A","HLA-DQA1", "LAMP3", "CCR7", "FSCN1")) 
  slide <- AddModuleScore_UCell(slide, features = DC_sig, name = "Dendritic_Cells", assay = "Spatial")
  
  # Mast Cells
  mast_Cheng <- list(c("KIT","TPSAB1","CPA4"))
  slide <- AddModuleScore_UCell(slide, features = mast_Cheng, name = "Mast_Cells", assay = "Spatial")
  
  drop <- rownames(predict[which(predict$`Myeloid Cells` == 0), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(44:50)] <- 0
  
  # Add epithelial labels ---------
  
  #Epithelial cells 
  genepi_Kara <- list(c("EPCAM", "EGFR","FZR1", "KRT14", "ITGA6",
                        "KRT5", "TP63", "KRT17", "MME", "KRT8",
                        "KRT18", "KRT19", "FOXA1", "GATA3", "MUC1",
                        "CD24", "GABRP",
                        "PTPRC-"))
  slide <- AddModuleScore_UCell(slide, features = genepi_Kara, name = "genepi_Kara", assay = "Spatial")

  drop <- rownames(predict[which(predict$`Epithelial Cells` < epicutoff), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(51)] <- 0
  
  # Add stromal labels -------
  
  # Myoepithelial Cells
  Myoepi_wuNguyen <- list(c("KRT5", "KRT14","ACTA2","TAGLN","EPCAM",
                            "PTPRC-"))
  slide <- AddModuleScore_UCell(slide, features = Myoepi_wuNguyen, name = "Myoepi_wuNguyen", assay = "Spatial")

  # Stromal 
  strom_karawu <- list(c("FAP","COL1A1","COL3A1","COL5A1","ACTA2","TAGLN","LUM",
                         "FBLN1","COL6A3","COL1A2","COL6A1","COL6A2","PDGFRB",
                         "FGFR1", "FGFR2", "FGFR3", "FGFR4", "JAG1",
                         "PTPRC-","EPCAM-")) 
  slide <- AddModuleScore_UCell(slide, features = strom_karawu, name = "strom_karawu", assay = "Spatial")

  # PVL Cells
  PVL.mark <- list(c("MCAM", "ACTA2", "PDGFRB", "PTPRC-", "EPCAM-")) 
  slide <- AddModuleScore_UCell(slide, features = PVL.mark, name = "newWuPVL", assay = "Spatial")

  # Fibroblasts
  fibro_wumelan <- list(c("FAP", "THY1", "DCN", "COL1A1","COL1A2", "COL6A1","COL6A2","COL6A3",
                          "EPCAM-","PTPRC-"))
  slide <- AddModuleScore_UCell(slide, features = fibro_wumelan, name = "fibro_wumelan", assay = "Spatial")

  # Endothelial Cells
  endo_kara <- list(c("PECAM1","VWF","CDH5","SELE","PTPRC-",
                      "EPCAM-"))
  slide <- AddModuleScore_UCell(slide, features = endo_kara, name = "endokara", assay = "Spatial")

  drop <- rownames(predict[which(predict$`Stromal Cells` == 0.1), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(52:56)] <- 0
  
  # Z-score cell type signature scores

  # zscore <- t(apply(slide@meta.data[,c(38:50)], 1, function(x) (x-mean(x))/sd(x)))
  # 
  # slide@meta.data[,c(38:50)] <- zscore
  colnames(slide@meta.data)[38:50] <- c("NK Cells", "CD8+ T Cells", "CD4+ T Cells", "B Cells", "Plasma Cells",
                                        "Myeloid Cells", "Mast Cells", "Epithelial Cells", "Myoepithelial Cells",
                                        "Stromal Cells", "PVL Cells", "Fibroblasts", "Endothelial Cells")

  # Name cell types --------
  
  colnames(slide@meta.data)[38:56] <- c("NK Cells", "CD8+ T Cells", "CD4+ T Cells", 
                                        "Regulatory T Cells", "B Cells", "Plasma Cells",
                                        "Myeloid Cells", "MDSCs", "Macrophages", 
                                        "Monocytes", "Neutrophils", "Dendritic Cells", "Mast Cells",
                                        "Epithelial Cells", "Myoepithelial Cells", "Stromal Cells", 
                                        "PVL Cells", "Fibroblasts", "Endothelial Cells")
  # Calculate correlations -------------
  
  corr <- slide@meta.data[,c(28:37, 38:56)]
  corr <- corr[which(!is.na(rowSums(corr))),]
  corr[is.na(corr)] <- 0
  
  setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/spatial")
  corr_heatmap <- simil(corr, 
                        NULL, 
                        paste0(unique(slide$orig.ident), "_corr.rds"), 
                        "corr")
  corr_heatmap <- readRDS(paste0(unique(slide$orig.ident), "_corr.rds"))
  
  GEs <- c("GE1", "GE2", "GE3", "GE4", "GE5", "GE6", "GE7", 
           "GE8", "GE9", "GE10")
  cells <- colnames(slide@meta.data)[38:56]
  
  colnames(corr_heatmap)[1:10] <- GEs
  rownames(corr_heatmap)[1:10] <- GEs
  
  corr_heatmap[is.na(corr_heatmap)] <- 0
  
  pdf(paste0(unique(slide$orig.ident), "_corrheatmap.pdf"),width = 10, height = 9)
  print(Heatmap(corr_heatmap))
  dev.off()
  
  corr_heatmap_GEs <- corr_heatmap[which(rownames(corr_heatmap) %in% GEs),
                                   which(colnames(corr_heatmap) %in% GEs)]
  pdf(paste0(unique(slide$orig.ident), "_corrheatmap_GEonly.pdf"),width = 5, height = 4.2)
  print(Heatmap(corr_heatmap_GEs))
  dev.off()
  
  corr_heatmap_Tcell <- corr_heatmap[which(rownames(corr_heatmap) %in% GEs),
                                     which(colnames(corr_heatmap) == "CD8+ T Cells")]
  pdf(paste0(unique(slide$orig.ident), "_corrheatmap_Tcell.pdf"),
      width = 5, height = 1.5)
  print(Heatmap(t(corr_heatmap_Tcell),
          col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))))
  dev.off()
  
  # Plot GEs and immune/stromal cells -----------
  
  p <- list() 
  for (k in c(1:10)) {
    j <- c(16:25, 38:56)[k]
    p[[k]] <- as.ggplot(SpatialFeaturePlot(slide, 
                                 features = colnames(slide@meta.data)[j], 
                                 images = "image", 
                                 stroke = NA, 
                                 alpha = c(0,1),
                                 min.cutoff = "q50",
                                 max.cutoff = "q90",
                                 pt.size.factor = 2.8,
                                 image.alpha = 0.3) + 
      scale_fill_gradient2(mid = "red", high = "red"))
  }
  
  for (k in c(11:29)) {
    j <- c(16:25, 38:56)[k]
    p[[k]] <- as.ggplot(SpatialFeaturePlot(slide, 
                                           features = colnames(slide@meta.data)[j], 
                                           images = "image", 
                                           stroke = NA, 
                                           alpha = c(0,1),
                                           min.cutoff = "q0",
                                           max.cutoff = "q90",
                                           pt.size.factor = 2.8,
                                           image.alpha = 0.3) + 
                          scale_fill_gradient2(mid = "red", high = "red"))
  }
  
  pdf(paste0(unique(slide$orig.ident), "_allcells.pdf"), width = 3, height = 3)
  print(p)
  dev.off()
  
  # Add slide/corrs --------- 
  
  new_slide_all[[i]] <- slide
  corrs[[i]] <- corr_heatmap
  
} 

# Integrated Visium spatial slides with scRNAseq dataset -------

for (i in c(15,16,17,18,20)) {
  # Setup slide -------

  slide <- slide_all[[i]]

  if (samples[i] == "HR+") {
    prim_obj_sub <- prim_obj_HR
  }
  if (samples[i] == "HER2+") {
    prim_obj_sub <- prim_obj_HER2
  }
  if (samples[i] == "TNBC") {
    prim_obj_sub <- prim_obj_TNBC
  }

  epicutoff <- 0.1

  # Perform anchor-based integration workflow ---------

  DefaultAssay(prim_obj_sub) <- "integrated"
  anchors <- FindTransferAnchors(reference = prim_obj_sub,
                                 query = slide,
                                 normalization.method = "SCT")

  DefaultAssay(slide) <- "Spatial"
  slide <- RunPCA(slide, assay = "SCT", verbose = F)
  predictions.assay <- TransferData(anchorset = anchors,
                                    refdata = prim_obj_sub$celltype_broad,
                                    prediction.assay = TRUE,
                                    weight.reduction = slide[["pca"]],
                                    dims = 1:50)
  slide[["predictions"]] <- predictions.assay
  predict <- as.data.frame(t(GetAssayData(slide, assay = "predictions", slot = "data")))

  # get prediction scores for each spot and plot classes
  DefaultAssay(slide) <- "predictions"
  setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/spatial")

  # Add GE scores ----------------

  setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
  GElist <- readxl::read_xlsx("cancer_GE.xlsx")
  GElist <- split(GElist$gene,GElist$GE)
  length(GElist) # get number of GEs

  slide <- AddModuleScore_UCell(slide,
                                features = GElist,
                                assay = "Spatial")
  
  slide <- AddModuleScore_UCell(slide,
                                features = GE_hub,
                                assay = "Spatial")

  slide@meta.data <- slide@meta.data[,c(1:53,55:62,54)]
  colnames(slide@meta.data)[53:62] <- paste0("raw_hub_", c("GE1", "GE2", "GE3", 
                                                           "GE4", "GE5", "GE6",
                                                           "GE7", "GE8", "GE9", "GE10"))
  zscore <- apply(slide@meta.data[,c(53:62)], 2, function(x) (x-mean(x))/sd(x))
  newmetadata <- cbind(apply(zscore, 1, function(x) which.max(x)),
                       apply(zscore, 1, function(x) max(x)),
                       zscore)
  newmetadata <- as.data.frame(newmetadata)
  colnames(newmetadata) <- c("maxGE_hub", "maxZscore_hub", 
                             "GE1_hub", "GE_hub", "GE3_hub", 
                             "GE4_hub", "GE5_hub", "GE6_hub", 
                             "GE7_hub", "GE8_hub", "GE9_hub", "GE10_hub")
  rownames(newmetadata) <- colnames(slide)

  slide <- AddMetaData(slide, newmetadata, 
                       col.name = colnames(newmetadata))

  drop <- rownames(predict[which(predict$`Epithelial Cells` < epicutoff), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(53:62)] <- 0

  # Add NK cell labels -----------
  
  # NK cells
  NK.mark2 <- list(c("CD3D","CD3E","CD3G","FGFBP2","KLRD1", "FCGR3A", "KLRK1", "NCAM1", 
                     "CD4-", "CD8A-", "CD8B-", "FOXP3-"))
  slide <- AddModuleScore_UCell(slide, features = NK.mark2, name = "NK.mark2", assay = "Spatial")
  
  drop <- rownames(predict[which(predict$Lymphocytes < 0.1), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(34)] <- 0
  
  # Add T cell labels -------
  
  # CD8+ T Cells
  T_karacd8 <- list(c("CD2","CD3D","CD3E","CD3G","CD8A","CD8B",
                      "CD4-", "FOXP3-"))
  slide <- AddModuleScore_UCell(slide, features = T_karacd8, name = "T_karacd8", assay = "Spatial")
  
  # CD4+ T Cells
  T_karacd4 <- list(c("CD2","CD3D","CD3E","CD3G","CD4",
                      "CD8A-","CD8B-", "FOXP3-"))
  slide <- AddModuleScore_UCell(slide, features = T_karacd4, name = "T_karacd4", assay = "Spatial")
  
  # T reg Cells
  T_reg <- list(c("CD2","CD3D","CD3E","CD3G","FOXP3"))
  slide <- AddModuleScore_UCell(slide, features = T_reg, name = "T_reg", assay = "Spatial")
  
  drop <- rownames(predict[which(predict$Lymphocytes == 0), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(35:37)] <- 0
  
  # Add B cell labels ------
 
  # B Cells
  B_KaraBrech <- list(c("CD79A","CD79B","BLNK","CD19","MS4A1")) 
  slide <- AddModuleScore_UCell(slide, features = B_KaraBrech, name = "B_KaraBrech", assay = "Spatial")
  
  # Plasma Cells
  plasma_colo <- list(c("CD27","IGHA1","SDC1","TNFRSF17","JCHAIN","MZB1","DERL3","CD38",
                        "IGHG1","IGHG3","IGHG4"))
  slide <- AddModuleScore_UCell(slide, features = plasma_colo, name = "plasma_colo", assay = "Spatial")
  
  drop <- rownames(predict[which(predict$Lymphocytes == 0), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(38:39)] <- 0
  
  # Add myeloid labels --------
  
  # Myeloid Cells
  myeloid_wuElliot <- list(c("PTPRC", "ITGAM","HLA-DR","ITGAX", "CD14","CD16","CD1C",
                             "CD1A","CD68","CD33",
                             "EPCAM-")) 
  slide <- AddModuleScore_UCell(slide, features = myeloid_wuElliot, name = "myeloid_wuElliot", assay = "Spatial")
  
  # MDSC
  MDSc_sig <- list(c("CD33", "ITGAM", "LY6G","LY6C","FUT4", "CEACAM8","IL4R","HLA-DRA",
                     "CD3D-","CD14-","CD19-","NCAM1-", "EPCAM-")) 
  slide <- AddModuleScore_UCell(slide, features = MDSc_sig, name = "MDSC_sig", assay = "Spatial")
  
  # Macrophages
  Macro_Cheng <- list(c("INHBA", "IL1RN", "CCL4", "NLRP3",
                        "EREG", "IL1B", "LYVE1", "PLTP",
                        "SELENOP", "C1QC", "C1QA", "APOE"))
  slide <- AddModuleScore_UCell(slide, features = Macro_Cheng, name = "Macrophage", assay = "Spatial")
  
  # Monocytes
  Mono_sig <- list(c("FCN1", "S100A9", "S100A8", "FCGR3A", "LST1", "LILRB2", 
                     "HLA-DRA-"))
  slide <- AddModuleScore_UCell(slide, features = Mono_sig, name = "Monocytes", assay = "Spatial")
  
  # Neutrophils 
  Neutro_sig <- list(c("CXCR1","CD15","FCGR3A", 
                       "CD14-",
                       "CSF3R", "S100A9","CD24A","TNF","CD274"))
  slide <- AddModuleScore_UCell(slide, features = Neutro_sig, name = "Neutrophils", assay = "Spatial")
  
  # Dendritic Cells
  DC_sig <- list(c("LILRA4","GZMB","IL3RA","CLEC9A", "FLT3", "IDO1", 
                   "CD1C", "FCER1A","HLA-DQA1", "LAMP3", "CCR7", "FSCN1")) 
  slide <- AddModuleScore_UCell(slide, features = DC_sig, name = "Dendritic_Cells", assay = "Spatial")
  
  # Mast Cells
  mast_Cheng <- list(c("KIT","TPSAB1","CPA4"))
  slide <- AddModuleScore_UCell(slide, features = mast_Cheng, name = "Mast_Cells", assay = "Spatial")
  
  drop <- rownames(predict[which(predict$`Myeloid Cells` == 0), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(40:46)] <- 0
  
  # Add epithelial labels ---------
  
  #Epithelial cells 
  genepi_Kara <- list(c("EPCAM", "EGFR","FZR1", "KRT14", "ITGA6",
                        "KRT5", "TP63", "KRT17", "MME", "KRT8",
                        "KRT18", "KRT19", "FOXA1", "GATA3", "MUC1",
                        "CD24", "GABRP",
                        "PTPRC-"))
  slide <- AddModuleScore_UCell(slide, features = genepi_Kara, name = "genepi_Kara", assay = "Spatial")
  
  drop <- rownames(predict[which(predict$`Epithelial Cells` < epicutoff), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(47)] <- 0
  
  # Add stromal labels -------
  
  # Myoepithelial Cells
  Myoepi_wuNguyen <- list(c("KRT5", "KRT14","ACTA2","TAGLN","EPCAM",
                            "PTPRC-"))
  slide <- AddModuleScore_UCell(slide, features = Myoepi_wuNguyen, name = "Myoepi_wuNguyen", assay = "Spatial")
  
  # Stromal 
  strom_karawu <- list(c("FAP","COL1A1","COL3A1","COL5A1","ACTA2","TAGLN","LUM",
                         "FBLN1","COL6A3","COL1A2","COL6A1","COL6A2","PDGFRB",
                         "FGFR1", "FGFR2", "FGFR3", "FGFR4", "JAG1",
                         "PTPRC-","EPCAM-")) 
  slide <- AddModuleScore_UCell(slide, features = strom_karawu, name = "strom_karawu", assay = "Spatial")
  
  # PVL Cells
  PVL.mark <- list(c("MCAM", "ACTA2", "PDGFRB", "PTPRC-", "EPCAM-")) 
  slide <- AddModuleScore_UCell(slide, features = PVL.mark, name = "newWuPVL", assay = "Spatial")
  
  # Fibroblasts
  fibro_wumelan <- list(c("FAP", "THY1", "DCN", "COL1A1","COL1A2", "COL6A1","COL6A2","COL6A3",
                          "EPCAM-","PTPRC-"))
  slide <- AddModuleScore_UCell(slide, features = fibro_wumelan, name = "fibro_wumelan", assay = "Spatial")
  
  # Endothelial Cells
  endo_kara <- list(c("PECAM1","VWF","CDH5","SELE","PTPRC-",
                      "EPCAM-"))
  slide <- AddModuleScore_UCell(slide, features = endo_kara, name = "endokara", assay = "Spatial")
  
  drop <- rownames(predict[which(predict$`Stromal Cells` == 0.1), ])
  slide@meta.data[which(rownames(slide@meta.data) %in% drop), c(48:52)] <- 0
  
  # Z-score cell type signature scores
  
  # zscore <- t(apply(slide@meta.data[,c(38:50)], 1, function(x) (x-mean(x))/sd(x)))
  # 
  # slide@meta.data[,c(38:50)] <- zscore
  colnames(slide@meta.data)[38:50] <- c("NK Cells", "CD8+ T Cells", "CD4+ T Cells", "B Cells", "Plasma Cells",
                                        "Myeloid Cells", "Mast Cells", "Epithelial Cells", "Myoepithelial Cells",
                                        "Stromal Cells", "PVL Cells", "Fibroblasts", "Endothelial Cells")
  
  # Name cell types --------
  
  colnames(slide@meta.data)[34:52] <- c("NK Cells", "CD8+ T Cells", "CD4+ T Cells", 
                                        "Regulatory T Cells", "B Cells", "Plasma Cells",
                                        "Myeloid Cells", "MDSCs", "Macrophages", 
                                        "Monocytes", "Neutrophils", "Dendritic Cells", "Mast Cells",
                                        "Epithelial Cells", "Myoepithelial Cells", "Stromal Cells", 
                                        "PVL Cells", "Fibroblasts", "Endothelial Cells")
  # Calculate correlations -------------

  corr <- slide@meta.data[,c(24:33, 34:52)]
  corr <- corr[which(!is.na(rowSums(corr))),]
  corr[is.na(corr)] <- 0

  corr_heatmap <- simil(corr,
                        NULL,
                        paste0(names(slide_all)[i], "_corr.rds"),
                        "corr")
  corr_heatmap <- readRDS(paste0(names(slide_all[i]), "_corr.rds"))

  GEs <- c("GE1", "GE2", "GE3", "GE4", "GE5", "GE6", "GE7",
           "GE8", "GE9", "GE10")
  cells <- colnames(slide@meta.data)[34:52]

  colnames(corr_heatmap)[1:10] <- GEs
  rownames(corr_heatmap)[1:10] <- GEs

  corr_heatmap[is.na(corr_heatmap)] <- 0

  pdf(paste0(names(slide_all)[i], "_corrheatmap.pdf"),width = 10, height = 9)
  print(Heatmap(corr_heatmap))
  dev.off()

  corr_heatmap_GEs <- corr_heatmap[which(rownames(corr_heatmap) %in% GEs),
                                   which(colnames(corr_heatmap) %in% GEs)]
  pdf(paste0(names(slide_all)[i], "_corrheatmap_GEonly.pdf"),width = 5, height = 4.2)
  print(Heatmap(corr_heatmap_GEs))
  dev.off()

  corr_heatmap_Tcell <- corr_heatmap[which(rownames(corr_heatmap) %in% GEs),
                                     which(colnames(corr_heatmap) == "CD8+ T Cells")]
  pdf(paste0(names(slide_all)[i], "_corrheatmap_Tcell.pdf"),
      width = 5, height = 1.5)
  print(Heatmap(t(corr_heatmap_Tcell),
          col = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))))
  dev.off()

  # Plot GEs and immune/stromal cells -----------
  
  slide@images[["image"]]@coordinates[["tissue"]] <- as.integer(slide@images[["image"]]@coordinates[["tissue"]])
  slide@images[["image"]]@coordinates[["row"]] <- as.integer(slide@images[["image"]]@coordinates[["row"]])
  slide@images[["image"]]@coordinates[["col"]] <- as.integer(slide@images[["image"]]@coordinates[["col"]])
  slide@images[["image"]]@coordinates[["imagerow"]] <- as.integer(slide@images[["image"]]@coordinates[["imagerow"]])
  slide@images[["image"]]@coordinates[["imagecol"]] <- as.integer(slide@images[["image"]]@coordinates[["imagecol"]])
  
  p <- list() 
  for (k in c(1:10)) {
    j <- c(53:62)[k]
    p[[k]] <- as.ggplot(SpatialFeaturePlot(slide, 
                                           features = colnames(slide@meta.data)[j], 
                                           images = "image", 
                                           stroke = NA, 
                                           alpha = c(0,1),
                                           min.cutoff = "q50",
                                           max.cutoff = "q90",
                                           pt.size.factor = 2.0,
                                           image.alpha = 0.3) + 
                          scale_fill_gradient2(mid = "red", high = "red"))
  }
  
  for (k in c(11:29)) {
    j <- c(12:21, 34:52)[k]
    p[[k]] <- as.ggplot(SpatialFeaturePlot(slide, 
                                           features = colnames(slide@meta.data)[j], 
                                           images = "image", 
                                           stroke = NA, 
                                           alpha = c(0,1),
                                           min.cutoff = "q0",
                                           max.cutoff = "q90",
                                           pt.size.factor = 2.0,
                                           image.alpha = 0.3) + 
                          scale_fill_gradient2(mid = "red", high = "red"))
  }
  
  pdf(paste0(names(slide_all)[i], "_hubGE.pdf"), width = 3, height = 3)
  print(p)
  dev.off()

  # Append slide/corrs ---------

  new_slide_all[[i]] <- slide
  corrs[[i]] <- corr_heatmap

}

# Save objects -------

saveRDS(new_slide_all, "all_spatial_obj.RDS")
saveRDS(corrs, "corr_spatial.RDS")

counts <- list()
for (slide in slide_all) {
  counts <- append(counts, list(slide[['SCT']]@data))
}
saveRDS(counts, "counts_spatial.RDS")

meta <- list()
for (i in c(1:6, 15:18,20)) {
  meta <- append(meta, list(new_slide_all[[i]]@meta.data))
}
saveRDS(meta, "meta_spatial.RDS")

# -----
# -----

# Analysis -------
# Read in files ----------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/spatial")
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

saveRDS(slides_all, "all_spatial_obj_040623.rds")


# Cell type percentages -------

pct_all <- data.frame()
for (i in c(1:11)) {
  temp <- as.data.frame(meta[[i]])
  if (i %in% c(1:6)) {
    temp <- temp[,c(38:50)]
  }
  if (i %in% c(7:11)) {
    temp <- temp[,c(34:46)]
  }
  temp[is.na(temp)] <- 0
  
  pct <- t(as.data.frame(colSums(temp)/dim(temp)[1]))
  if (i %in% c(1:6)) {
    rownames(pct) <- names(slides_all[i])
  }
  if (i %in% c(7:10)) {
    rownames(pct) <- names(slides_all[i+8])
  }
  if (i %in% c(11)) {
    rownames(pct) <- names(slides_all[20])
  }
  pct_all <- rbind(pct_all, pct)
}
colnames(pct_all) <- colnames(meta[[1]][,c(38:50)])
rownames(pct_all) <- slide_name[c(1:6, 15:18, 20)]

anno <- HeatmapAnnotation(subtype = samples[c(1:6,15:18,20)],
                          col = list(subtype = c("HER2+" = "#B0FC96",
                                                 "HR+" = "#96B0FC",
                                                 "TNBC" = "#FC96B0")))


pdf("celltypes.pdf", width = 6, height = 6)
Heatmap(t(pct_all), top_annotation = anno,
        column_split = factor(samples[c(1:6,15:18,20)]),
        col = colorRamp2(c(0, 0.2), c("#EFEFEF", "red")))
dev.off()
