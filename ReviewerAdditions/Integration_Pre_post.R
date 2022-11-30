
DoubletResult_Dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/PrimaryBreast_Processing/Individual_Dataset_Objects/withDoubletMeta"
finalQC_Dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/PrimaryBreast_Processing/Individual_Dataset_Objects/finalQC"
prepostDir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/ReviewerAdditions/Pre_post_integration_compare"


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





# FUNCTIONS ---------------------

Find_top_bot_ngenes_forPCs <- function(seurat_obj, top_or_bot = c("top", "bottom"), 
                                       PC_count = 20, gene_count = 10, save_csv = T, 
                                       project_name, output_dir)
{
  top_or_bot <- match.arg(top_or_bot, c("top", "bottom"))
  
  PC_post <- as.data.frame(seurat_obj@reductions$pca@feature.loadings)
  PC_post$gene_names <- rownames(PC_post)
  rownames(PC_post) <- 1:nrow(PC_post)
  
  if (top_or_bot == "top")
  {
    top_PC_post_final <- data.frame(matrix(ncol = as.numeric(PC_count), nrow = as.numeric(gene_count)))
    colnames(top_PC_post_final) <- colnames(PC_post)[1:PC_count]
    
    for (compon in colnames(PC_post)[1:PC_count])
    {
      ind_top <- PC_post[order(PC_post[[compon]], decreasing = TRUE)[1:gene_count], c(compon,"gene_names")]
      colnames(ind_top) <- c(compon, paste(compon, "gene_names", sep = "_"))
      top_PC_post_final <- cbind(top_PC_post_final, ind_top)
    }
    
    top_PC_post_final <- top_PC_post_final[,-c(1:PC_count)]
    rownames(top_PC_post_final) <- 1:nrow(top_PC_post_final)
    
    write.csv(top_PC_post_final, paste0(project_name,"_top",gene_count, "genes_top",PC_count, "_PCs_",Sys.Date(),".csv"))
  }
  
  else
  {
    bottom_PC_post_final <- data.frame(matrix(ncol = as.numeric(PC_count), nrow = as.numeric(gene_count)))
    colnames(bottom_PC_post_final) <- colnames(PC_post)[1:PC_count]
    
    for (compon in colnames(PC_post)[1:PC_count])
    {
      ind_bot <- PC_post[order(PC_post[[compon]], decreasing = TRUE)[(dim(PC_post)[1] - gene_count):dim(PC_post)[1]], c(compon,"gene_names")]
      colnames(ind_bot) <- c(compon, paste(compon, "gene_names", sep = "_"))
      bottom_PC_post_final <- cbind(bottom_PC_post_final, ind_bot)
    }
    
    bottom_PC_post_final <- bottom_PC_post_final[,-c(1:PC_count)]
    rownames(bottom_PC_post_final) <- 1:nrow(bottom_PC_post_final)
    
    write.csv(bottom_PC_post_final, paste0(project_name,"_bottom",gene_count, "genes_top",PC_count, "_PCs_",Sys.Date(),".csv"))
    
  }
  
  
}

#UMAP before =============================================================

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

table(combo$Doublet.Call)
combo <- subset(combo, subset = Doublet.Call == "Doublet", invert = T)

combo <- SCTransform(combo, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA", "percent.hb", 
                                                "percent.platelet", "percent.heatshock"), verbose = TRUE)


combo <- RunPCA(combo, npcs = 100, verbose = FALSE)

combo$orig.ident[which(combo$orig.ident == "AziziT")] <- "Azizi"
combo$orig.ident[which(combo$orig.ident == "Aziziimmune")] <- "Azizi"

sort(table(combo$Capture.Method))
sort(table(combo$orig.ident))

setwd(prepostDir)
#pdf("BatchTestSCTRef_71922_long.pdf", width = 16, height = 9)
pdf("BatchTest_NOintegration_102822_long.pdf", width = 16, height = 9)
p1 <- DimPlot(combo, reduction = "pca", raster = F, group.by = "Capture.Method",
              order = c("Smart-Seq2", "inDrop v2", "Singleron GEXSCOPE Single Cell RNAseq Library Kit",
                        "10X Genomics Single Cell 3' v2", "10X Genomics Chromium v2 5'",
                        "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library",
                        "10X Genomics Chromium"))
p1
p2 <- DimPlot(combo, reduction = "pca", raster = F, group.by = "orig.ident",
              order =c("Karaayvaz", "Savas", "Wu", "Xu",
                       "Azizi", "Qian", "Wu2021prim", "Pal_Prim"))
p2
dev.off()

#ggsave("BatchTestSCTRef_origident_71922_long.pdf", plot = p2, width = 16, height = 9)



ElbowPlot(combo, ndims = 100)
setwd(prepostDir)
pdf("test.pdf")
DimHeatmap(combo, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 15:25, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 30:45, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 45:55, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 55:65, cells = 500, balanced = TRUE)
DimHeatmap(combo, dims = 65:75, cells = 500, balanced = TRUE)
dev.off()


combo <- FindNeighbors(combo, reduction = "pca", dims = 1:55, nn.method = "rann")#, k.param= 500)
combo <- FindClusters(combo, resolution = 7)
combo <- RunUMAP(combo, reduction = "pca", dims = 1:55, verbose = TRUE, seed.use = 123)

combo$BC.Subtype[which(combo$BC.Subtype == "ER+")] <- "HR+"
combo$Patient[which(combo$Patient == "BC11_1")] <- "BC11"
combo$Patient[which(combo$Patient == "BC11_2")] <- "BC11"



p <- DimPlot(combo, reduction = "umap", label = F, repel = T, raster = FALSE) + ggtitle(label = " ") #,

pdf("test.pdf", width = 16.22, height = 14.22)
p + theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()

saveRDS(combo, "PrimObject_nointegration_102822.rds")


setwd(prepostDir)
pdf("UNINTEGRATED_BCsub_wholeUMAP_11322.pdf", width = 16.22, height = 14.22)
DimPlot(combo, reduction = "umap", group.by = "BC.Subtype", raster = FALSE) + ggtitle(label = " ") +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()


combo$orig.ident[which(combo.reference$orig.ident == "AziziT")] <- "Azizi"
combo$orig.ident[which(combo.reference$orig.ident == "Aziziimmune")] <- "Azizi"


#pdf("UNINTEGRATED_capturemethod_primUmAP_11322.pdf", width = 16.22, height = 14.22)
#pdf("UNINTEGRATED_libraryprep_primUmAP_11322.pdf", width = 16.22, height = 14.22)
pdf("UNINTEGRATED_dataset_primUmAP_11322.pdf", width = 16.22, height = 14.22)
DimPlot(combo, reduction = "umap", group.by = "orig.ident", raster = FALSE) + ggtitle(label = " ") +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()

#PC loadings before ------
setwd(prepostDir)
combo.unint <- readRDS("PrimObject_nointegration_102822.rds")


#top10_20PC <- Find_top_bot_ngenes_forPCs(combo.unint, "top", "PrimMain", prepostDir)

PC_post <- as.data.frame(combo.unint@reductions$pca@feature.loadings)
class(PC_post)

PC_post$gene_names <- rownames(PC_post)
rownames(PC_post) <- 1:nrow(PC_post)

top_PC_post_final <- data.frame(matrix(ncol = 20, nrow = 1))
colnames(top_PC_post_final) <- colnames(PC_post)[1:20]

for (compon in colnames(PC_post)[1:20])
{
  ind_top10 <- PC_post[order(PC_post[[compon]], decreasing = TRUE)[1:10], c(compon,"gene_names")]
  colnames(ind_top10) <- c(compon, paste(compon, "gene_names", sep = "_"))
  top_PC_post_final <- cbind(top_PC_post_final, ind_top10)
}

top_PC_post_final <- top_PC_post_final[,-c(1:20)]
rownames(top_PC_post_final) <- 1:nrow(top_PC_post_final)
write.csv(top_PC_post_final, "MainPrim_BEFOREINTEGRATION_top10genestop20PCs_andloadings_11222.csv")


dim(PC_post)
bottom_PC_post_final <- data.frame(matrix(ncol = 20, nrow = 1))
colnames(bottom_PC_post_final) <- colnames(PC_post)[1:20]

for (compon in colnames(PC_post)[1:20])
{
  ind_top10 <- PC_post[order(PC_post[[compon]], decreasing = T)[(dim(PC_post)[1] - 10):dim(PC_post)[1]], c(compon,"gene_names")]
  colnames(ind_top10) <- c(compon, paste(compon, "gene_names", sep = "_"))
  bottom_PC_post_final <- cbind(bottom_PC_post_final, ind_top10)
}

bottom_PC_post_final <- bottom_PC_post_final[-1,-c(1:20)]
rownames(bottom_PC_post_final) <- 1:nrow(bottom_PC_post_final)

write.csv(bottom_PC_post_final, "MainPrim_BEFOREINTEGRATION_BOTTOM10genestop20PCs_andloadings_11222.csv")

#UMAP after ----

pdf("BCsub_wholeUMAP_72822.pdf", width = 16.22, height = 14.22)
DimPlot(combo.reference, reduction = "umap", group.by = "BC.Subtype", raster = FALSE) + ggtitle(label = " ") +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()


combo.reference$orig.ident[which(combo.reference$orig.ident == "AziziT")] <- "Azizi"
combo.reference$orig.ident[which(combo.reference$orig.ident == "Aziziimmune")] <- "Azizi"


#pdf("capturemethod_primUmAP_102822.pdf", width = 16.22, height = 14.22)
#pdf("libraryprep_primUmAP_102822.pdf", width = 16.22, height = 14.22)
pdf("dataset_primUmAP_102822.pdf", width = 16.22, height = 14.22)
DimPlot(combo.reference, reduction = "umap", group.by = "orig.ident", raster = FALSE) + ggtitle(label = " ") +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) + 
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 15))
dev.off()




#PC loadings after -------

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/")
combo.reference <- readRDS("PrimObj_withGEmeta_with_correctassays_withfinalreprog_111622.rds") ##newnewnew


#top10_20PC <- Find_top_bot_ngenes_forPCs(combo.reference, "top", "PrimMain", prepostDir)

PC_post <- as.data.frame(combo.reference@reductions$pca@feature.loadings)
class(PC_post)

PC_post$gene_names <- rownames(PC_post)
rownames(PC_post) <- 1:nrow(PC_post)

top_PC_post_final <- data.frame(matrix(ncol = 20, nrow = 1))
colnames(top_PC_post_final) <- colnames(PC_post)[1:20]

for (compon in colnames(PC_post)[1:20])
{
  ind_top10 <- PC_post[order(PC_post[[compon]], decreasing = TRUE)[1:10], c(compon,"gene_names")]
  colnames(ind_top10) <- c(compon, paste(compon, "gene_names", sep = "_"))
  top_PC_post_final <- cbind(top_PC_post_final, ind_top10)
}

top_PC_post_final <- top_PC_post_final[,-c(1:20)]
rownames(top_PC_post_final) <- 1:nrow(top_PC_post_final)
write.csv(top_PC_post_final, "MainPrim_top10genestop20PCs_andloadings_102822.csv")


dim(PC_post)
bottom_PC_post_final <- data.frame(matrix(ncol = 1, nrow = 1))
colnames(bottom_PC_post_final) <- colnames(PC_post)[1:20]

for (compon in colnames(PC_post)[1:20])
{
  ind_top10 <- PC_post[order(PC_post[[compon]], decreasing = TRUE)[(dim(PC_post)[1] - 10):dim(PC_post)[1]], c(compon,"gene_names")]
  colnames(ind_top10) <- c(compon, paste(compon, "gene_names", sep = "_"))
  bottom_PC_post_final <- cbind(bottom_PC_post_final, ind_top10)
}

bottom_PC_post_final <- bottom_PC_post_final[,-1]
write.csv(bottom_PC_post_final, "MainPrim_BOTTOM10genestop20PCs_andloadings_102822.csv")


pdf("test.pdf")
VlnPlot(object = combo.reference, features = 'PC_1', group.by = "orig.ident", raster = F)
dev.off()

#barplot of cohort per cluster ------------------

setwd("/project/InternalMedicine/Chan_lab/shared/")
combo.reference <- readRDS("PrimObj_withGEmeta_with_correctassays_withfinalreprog_111622.rds") ##newnewnew

table(Idents(combo.reference))

Idents(combo.reference) <- combo.reference$celltype_final


write.csv(table(combo.reference$celltype_final, combo.reference$BC.Subtype), "BCsub_cellcount_percelltype_11322.csv")
write.csv(table(combo.reference$celltype_final, combo.reference$Capture.Method), "Technology_cellcount_percelltype_11322.csv")
write.csv(table(combo.reference$celltype_final, combo.reference$Library.Preparation), "LibraryPrep_cellcount_percelltype_11322.csv")
write.csv(table(combo.reference$celltype_final, combo.reference$orig.ident), "Dataset_cellcount_percelltype_11322.csv")
write.csv(table(combo.reference$celltype_final, combo.reference$Patient), "Patient_cellcount_percelltype_11322.csv")
write.csv(table(combo.reference$celltype_final, combo.reference$samples), "sample_cellcount_percelltype_11322.csv")



sobjlists <- FetchData(object = combo.reference,
                       vars = c("samples",
                                "Patient",
                                "BC.Subtype",
                                "orig.ident",
                                "Capture.Method",
                                "Library.Preparation",
                                "celltype_final"))




sobjlists <- sobjlists %>% dplyr::group_by(celltype_final,
                                           # samples,
                                           # Patient,
                                          BC.Subtype,
                                           # orig.ident,
                                           #Capture.Method,
                                           # Library.Preparation
  ) %>%
  dplyr::summarise(Nb = n()) %>%
  dplyr::mutate(C = sum(Nb)) %>%
  dplyr::mutate(percent = Nb/C*100)



library(viridis)
library(ggh4x)
bp <- ggplot(sobjlists,
             aes(x = celltype_final,
                 y = percent,
                 #group = as.factor(BC.Subtype),
                 fill = BC.Subtype)) +
  scale_fill_viridis(option = "H", discrete = T) + 
  # scale_fill_manual(name = "Dataset",
  #                   values = turbo(9),
  #                   unique(sobjlists$orig.ident)) +
  
  # scale_fill_manual(name = "Reprogrammed Status",
  #                   values = rev(c("#440154FF", "light grey")),
  #                   rev(unique(sobjlists$celltype_withreprog_simply))) +
  geom_bar(stat = "identity") +
  # geom_text(aes(label = paste0(round(sobjlists$percent_reprog), "%"),
  #               y = 120,
  #               fill = NULL),
  #           angle = 90,
  #           colour = "#666666",
  #           size = 3) +
  theme_minimal() +
  # theme(axis.text.x = element_blank(),
  #       axis.ticks.x = element_blank(),
  #       axis.title.x = element_text(size = 20),
  #       axis.text.y = element_text(size = 20),
  #       axis.title.y = element_text(size = 20)) +
  # scale_y_continuous(expand = c(0, 0),
  #                    limits = c(0, 105),#135),
  #                    breaks = seq(0,100,25)) +
  # facet_nested( ~ BC.Subtype,
  #               scales = "free",
  #               space = "free",
  #               switch = "x"
  # ) +
  theme(legend.text=element_text(size=20)) +
  theme(legend.title=element_text(size=20)) +
  ylab("% Cells") +
  xlab("Cell Types") +
  theme(strip.text.x = element_text(size = 20, face = "italic"),
        strip.background = element_rect(colour = "#FFFFFF",
                                        size = 1.5,
                                        fill = "#EEEEEE"),
        panel.spacing.x = unit(-0.1, "lines")) +
  RotatedAxis()


#pdf("proportion_dataset_perclust_112522.pdf", width = 18.4, height = 5.4)
#pdf("proportion_samples_perclust_112522.pdf", width = 25.4, height = 12.4)
#pdf("proportion_Patient_perclust_112522.pdf", width = 25.4, height = 12.4)
#pdf("proportion_CaptureMethod_perclust_112522.pdf", width = 25.4, height = 12.4)
#pdf("proportion_LibraryPrep_perclust_112522.pdf", width = 25.4, height = 12.4)
pdf("proportion_BCsub_perclust_112522.pdf", width = 25.4, height = 12.4)
bp
dev.off()


#compare original labsl -=-=======

combo.reference <- readRDS("/project/InternalMedicine/Chan_lab/shared/PrimObj_withGEmeta_with_correctassays_withfinalreprog_111622.rds")

finalQC_Dir <- "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Preprocessing/scRNA/PrimaryBreast_Processing/Individual_Dataset_Objects/finalQC"

setwd(finalQC_Dir)
AziziPrim.fixed <- readRDS("AziziPrim_fixed_62722.rds")
head(AziziPrim.fixed@meta.data)
AziziT.fixed <- readRDS("AziziT_fixed_62722.rds") #has celltype info (col: celltype_minor)
head(AziziT.fixed@meta.data)
Karaayvaz.fixed <- readRDS("Karaayvaz_fixed_62722.rds") #has celltype info (col: celltype_minor)
head(Karaayvaz.fixed@meta.data)
Pal.fixed <- readRDS("Pal_fixed_62722.rds")
head(Pal.fixed@meta.data)
Qian.fixed <- readRDS("Qian_fixed_62722.rds") #has celltype info (col: celltype_minor)
head(Qian.fixed@meta.data)
Savas.fixed <- readRDS("Savas_fixed_62722.rds") #has celltype info (col: celltype_minor)
head(Savas.fixed@meta.data)
OldWu.fixed <- readRDS("OldWu_fixed_62722.rds") #has celltype info (col: celltype_minor)
head(OldWu.fixed@meta.data)
NewWu.fixed <- readRDS("NewWu_fixed_62722.rds") #has celltype info (cols: celltype_major, celltype_minor, celltype_subset)
head(NewWu.fixed@meta.data)
Xu.fixed <- readRDS("Xu_fixed_62722.rds")
head(Xu.fixed@meta.data)


unique(combo.reference$orig.ident)

Kara.new <- subset(combo.reference, subset = orig.ident == "Karaayvaz")
unique(Kara.new$celltype_final)
unique(Kara.new$celltype_minor)
Kara.new$celltype_final[which(Kara.new$celltype_final == "Cancer Epithelial Cells")] <- "epithelial"
Kara.new$celltype_final[which(Kara.new$celltype_final == "CD4+ T Cells")] <- "Tcell"
Kara.new$celltype_final[which(Kara.new$celltype_final == "Fibroblasts")] <- "stroma"
Kara.new$celltype_final[which(Kara.new$celltype_final == "Epithelial Cells")] <- "epithelial"
Kara.new$celltype_final[which(Kara.new$celltype_final == "Perivascular-like (PVL) Cells")] <- "stroma"
Kara.new$celltype_final[which(Kara.new$celltype_final == "Plasma Cells")] <- "Bcell"



AziziT.new <- subset(combo.reference, subset = orig.ident == "AziziT")
table(Kara.new$celltype_final)
table(Kara.new$celltype_minor)

Qian.new <- subset(combo.reference, subset = orig.ident == "Qian")
table(Kara.new$celltype_final)
table(Kara.new$celltype_minor)

Savas.new <- subset(combo.reference, subset = orig.ident == "Savas")
table(Kara.new$celltype_final)
table(Kara.new$celltype_minor)

OldWu.new <- subset(combo.reference, subset = orig.ident == "Wu")
table(Kara.new$celltype_final)
table(Kara.new$celltype_minor)

NewWu.new <- subset(combo.reference, subset = orig.ident == "Wu2021prim")
table(Kara.new$celltype_final)
table(Kara.new$celltype_minor)
