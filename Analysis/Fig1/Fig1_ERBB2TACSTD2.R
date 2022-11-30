

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This will generate all of the sub-figures for Figure 1 and Figure S2,
# and is the exploratory analysis done on the cancer epithelial cells.

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# SETUP -----------------------------------------------------------------------
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

# functions ------------------------------------------------------

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

DEG_Remove_mito <- function (df){
  df_rm_mito <- df[!grepl("^MT-|^MT.",rownames(df)),]
  return(df_rm_mito)
}


# Load in object ---------------

##object for pre-print
# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only")
# cancer.epi <- readRDS("cancerEpiwithGenes_nozallgenedem_71322.rds")
# DefaultAssay(cancer.epi) <- "RNA"
# cancer.epi <- NormalizeData(cancer.epi, assay = "RNA")

##Lily's original object
#setwd("/project/InternalMedicine/Chan_lab/shared")
#cancer.epi <- readRDS("cancerepi_withGEmetadata_111122.rds")



setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig1")
cancer.epi <- readRDS("cancerepi_withinferCNV_with_Trop2Her2cats_111422.rds") ##Lily version

DefaultAssay(cancer.epi) <- "RNA"
cancer.epi <- NormalizeData(cancer.epi, assay = "RNA")
head(cancer.epi@meta.data)

##run TROP2 and HER2 labeling again
cancer.epi$TROP2_ident <- NA
cancer.epi$TROP2_ident[(cancer.epi$TACSTD2 >= quantile(cancer.epi$TACSTD2,0.9))] <- "high"
cancer.epi$TROP2_ident[(cancer.epi$TACSTD2 <= quantile(cancer.epi$TACSTD2,0.9))] <- "med"
cancer.epi$TROP2_ident[(cancer.epi$TACSTD2 <= 0)] <- "low"

cancer.epi$HER2_ident <- NA
cancer.epi$HER2_ident[(cancer.epi$ERBB2 <= quantile(cancer.epi$ERBB2, 0.975))] <- "med"
cancer.epi$HER2_ident[(cancer.epi$ERBB2 >= quantile(cancer.epi$ERBB2, 0.975))] <- "high"
cancer.epi$HER2_ident[(cancer.epi$ERBB2 <= 0)] <- "low"


# Add UScores for Clinical Targets --------------------------


HER2 <- list(c("ERBB2"))
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = HER2, name = "ERBB2", assay = "RNA")  
FeaturePlot(object = cancer.epi, features = "signature_1ERBB2", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) 


TROP2 <- list(c("TACSTD2"))
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = TROP2, name = "TACSTD2", assay = "RNA")  
FeaturePlot(object = cancer.epi, features = "signature_1TACSTD2", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) 

ER <- list(c("ESR1"))
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = ER, name = "ESR1", assay = "RNA")  
FeaturePlot(object = cancer.epi, features = "signature_1ESR1", order = TRUE, label = TRUE, repel = TRUE, raster = FALSE)# min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) 

HER3 <- list(c("ERBB3"))
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = HER3, name = "ERBB3", assay = "RNA")  
FeaturePlot(object = cancer.epi, features = "signature_1ERBB3", order = TRUE, label = TRUE, repel = TRUE, raster = FALSE)# min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) 

PIK3CA <- list(c("PIK3CA"))
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = PIK3CA, name = "PIK3CA", assay = "RNA")  
FeaturePlot(object = cancer.epi, features = "signature_1PIK3CA", order = TRUE, label = TRUE, repel = TRUE, raster = FALSE)# min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) 

CD274 <- list(c("CD274"))
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = CD274, name = "CD274", assay = "RNA")  
FeaturePlot(object = cancer.epi, features = "signature_1CD274", order = TRUE, label = TRUE, repel = TRUE,min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "CD274/PDL1 (PC = 40)")

EGFR <- list(c("EGFR"))
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = EGFR, name = "EGFR", assay = "RNA")  
FeaturePlot(object = cancer.epi, features = "signature_1EGFR", order = TRUE, label = TRUE, repel = TRUE,min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "EGFR (PC = 40)")

FGFR <- list(c("FGFR1", "FGFR2", "FGFR3"))
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = FGFR, name = "FGFR", assay = "RNA")  
FeaturePlot(object = cancer.epi, features = "signature_1FGFR", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "PIK3CA (PC = 40)")

NTRKlist <- list(c("NTRK1", "NTRK2", "NTRK3")) #NTRK1
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = NTRKlist, name = "NTRK", assay = "RNA")  
FeaturePlot(object = cancer.epi, features = "signature_1NTRK", order = TRUE, label = TRUE, repel = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) 


NECTIN2 <- list(c("NECTIN2"))
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = NECTIN2, name = "NECTIN2", assay = "RNA")  
FeaturePlot(object = cancer.epi, features = "signature_1NECTIN2", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "CD112/PVRL2/NECTIN2 (PC = 40)")

CDK <- list(c("CDK"))
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = CDK, name = "CDK", assay = "RNA")  
FeaturePlot(object = cancer.epi, features = "signature_1CDK", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "CD112/PVRL2/NECTIN2 (PC = 40)")

AR <- list(c("AKR1B1"))
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = AR, name = "AR", assay = "RNA")  
FeaturePlot(object = cancer.epi, features = "signature_1AR", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "CD112/PVRL2/NECTIN2 (PC = 40)")


LAG3 <- list(c("LAG3"))
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = LAG3, name = "LAG3")  
FeaturePlot(object = cancer.epi, features = "signature_1LAG3", order = TRUE, label = TRUE, min.cutoff = 0, max.cutoff = 0.5, raster = FALSE) + ggtitle(label = "LAG3 (PC = 40)")


colnames(cancer.epi@meta.data) <- sub("signature_1", "", colnames(cancer.epi@meta.data))
colnames(cancer.epi@meta.data)

# saveRDS(cancer.epi, "cancerepi_WHOLEwithsigs_100422.rds")

# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only")
# cancer.epi <- readRDS("cancerepi_WHOLEwithsigs_100422.rds")

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig1")
cancer.epi <- readRDS("cancerepi_withinferCNV_with_Trop2Her2cats_111422.rds")

head(cancer.epi@meta.data)

# -------------
# -------------
# TROP2 Labeling ---------------

cancer.epi$TROP2_ident <- NA
cancer.epi$TROP2_ident[(cancer.epi$TACSTD2 >= quantile(cancer.epi$TACSTD2,0.9))] <- "high"
cancer.epi$TROP2_ident[(cancer.epi$TACSTD2 <= quantile(cancer.epi$TACSTD2,0.9))] <- "med"
cancer.epi$TROP2_ident[(cancer.epi$TACSTD2 <= 0)] <- "low"

Idents(cancer.epi) <- cancer.epi$TROP2_ident


# DEG Analysis for TROP2  ---------------------------------------

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig1/Subsetting") ##for reintegrated lists

j.markers_DGE <- FindAllMarkers(cancer.epi,
                                slot = "data",
                                min.cells.group = 5,
                                min.pct = 0.2,
                                logfc.threshold = 0,
                                test.use = "MAST")

j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)

write.csv(j.markers_DGE, "reintegrated_TROP2_DGEs_logfc0mincell5minpct0.2_111822_YESmitopresent.csv")
#write.csv(j.markers_DGE, "TROP2_DGEs_logfc0mincell5minpct0.2_111522_YESmitopresent.csv") ##lily

j.markers_DGE_filtered <- DEG_Remove_mito(j.markers_DGE)
write.csv(j.markers_DGE_filtered, "reintegrated_TROP2_DGEs_logfc0mincell5minpct0.2_111822_NOmito.csv")
#write.csv(j.markers_DGE_filtered, "TROP2_DGEs_logfc0mincell5minpct0.2_111522_NOmito.csv") ##lily


# 
j.markers_DGE_nothresh <- FindAllMarkers(cancer.epi,
                                slot = "data",
                                # min.cells.group = 5,
                                # min.pct = 0.2,
                                # logfc.threshold = 0,
                                test.use = "MAST")

j.markers_DGE_nothresh$min.pct.diff <- abs(j.markers_DGE_nothresh$pct.1 - j.markers_DGE_nothresh$pct.2)

write.csv(j.markers_DGE_nothresh, "reintegrated_TROP2_DGEs_ABSOLUTELYNOTHRESH_111822.csv")
write.csv(j.markers_DGE_nothresh, "TROP2_DGEs_ABSOLUTELYNOTHRESH_111822.csv")
#write.csv(j.markers_DGE_nothresh, "TROP2_DGEs_ABSOLUTELYNOTHRESH_111522.csv") ##lily


# # write.csv(j.markers_DGE_nothresh, "TROP2_DGEs_71422_NOTHRESH.csv")
# # write.csv(j.markers_DGE_nothresh, "TROP2_DGEs_71422_logfc0mincell5minpct0.2.csv")


# TROP2 MA plot ------------

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig1")
j.markers_DGE <- read.csv("TROP2_DGEs_ABSOLUTELYNOTHRESH_111522.csv")
head(j.markers_DGE)
rownames(j.markers_DGE) <- j.markers_DGE$X
j.markers_DGE_filtered <- DEG_Remove_mito(j.markers_DGE)

Idents(cancer.epi) <- cancer.epi$TROP2_ident
avgexp <- AverageExpression(object = cancer.epi, features = unique(j.markers_DGE_filtered$gene))
avgexp <- as.data.frame(avgexp[['RNA']])

#j.markers_DGE_filtered <- j.markers_DGE_filtered[which(j.markers_DGE_filtered$cluster == "high"),]
#j.markers_DGE_filtered <- j.markers_DGE_filtered[which(j.markers_DGE_filtered$cluster == "med"),]
j.markers_DGE_filtered <- j.markers_DGE_filtered[which(j.markers_DGE_filtered$cluster == "low"),]

avgexp <- avgexp[which(rownames(avgexp) %in% unique(j.markers_DGE_filtered$gene)),]

data <- cbind(j.markers_DGE_filtered, avgexp)
colnames(data)
#data_input <- data[,c(12,3,6,8)] #high
#data_input <- data[,c(11,3,6,8)] #med
data_input <- data[,c(10,3,6,8)] #low
colnames(data_input) <- c("baseMean", "log2FoldChange", "padj", "gene")

options(ggrepel.max.overlaps = 8)
library(ggpubr)
p <- ggmaplot(data_input, 
              fdr = 0.05, fc = 1.5, 
              size = 1,
              palette = c("#B31B21", "#1465AC", "darkgray"),
              genenames = as.vector(data_input$gene),
              top = 100,
              font.label = c("plain", 11),
              select.top.method = "fc",
              ggtheme = ggplot2::theme_classic()) +
  xlim(0,7) + ylim(-2.5,4) + NoLegend()


#ggsave("TROP2_MAplot_HIGHvall_111622.pdf", plot = p, width = 4.3, height = 4)
#ggsave("TROP2_MAplot_MEDvall_111622.pdf", plot = p, width = 4.3, height = 4)
ggsave("TROP2_MAplot_LOWvall_111622.pdf", plot = p, width = 4.3, height = 4)


# ggsave("TROP2_MAplot_71422.pdf", plot = p, width = 4.3, height = 4)
# ggsave("TROP2_MAplot_LOWvall_71822.pdf", plot = p, width = 4.3, height = 4)
# ggsave("TROP2_MAplot_NEGvall_71822.pdf", plot = p, width = 4.3, height = 4)

# TROP2 GO analysis -----------------------------

library(org.Hs.eg.db)
orgdb = "org.Hs.eg.db"

dfsample <- read.csv("TROP2_DGEs_logfc0mincell5minpct0.2_111522_NOmito.csv", header = T, row.names = 1)
dfsample <- dfsample[which(abs(dfsample$min.pct.diff) >= 0.1),]
##remove mitochondrial genes
dfsample <- split(dfsample$gene, dfsample$cluster)

dfsample$high <- bitr(dfsample$high, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$med <- bitr(dfsample$med, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$low <- bitr(dfsample$low, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)

genelist <- list("high" = dfsample$high$ENTREZID,
                 "med" = dfsample$med$ENTREZID,
                 "low" = dfsample$low$ENTREZID)

m_df <- msigdbr(species = "Homo sapiens", category ="H") %>%#, subcategory = "CGP") %>% 
  dplyr::select(gs_name, entrez_gene)

m_df <- m_df[-which(grepl("_UP", m_df$gs_name, fixed = TRUE)),]

GOclusterplot <- compareCluster(genelist, 
                                fun = enricher, 
                                TERM2GENE = m_df, 
                                pvalueCutoff = 0.05
)

p <- dotplot(GOclusterplot, includeAll = F)
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


ggsave("TROP2_hallmark_orghsdb_0.1mindiffnologfc_nomito_111522.pdf", plot = p, width = 15, height = 10, units = "cm")
#ggsave("TROP2_hallmark_orghsdb_0.1mindiffnologfc_nomito_72522.pdf", plot = p, width = 15, height = 10, units = "cm")


# Barplot of TROP2 high cells by sample ------------

sobjlists <- FetchData(object = cancer.epi, 
                       vars = c("samples",
                                "Patient", 
                                "BC.Subtype", 
                                "TROP2_ident"))

sobjlists <- sobjlists %>% dplyr::group_by(samples, 
                                           Patient, 
                                           BC.Subtype, 
                                           TROP2_ident) %>% 
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 

sobjlists$percent_low <- NA
for (i in unique(sobjlists$samples)) {
  sobjlists$percent_low[which(sobjlists$samples == i)] <- sobjlists$percent[which((sobjlists$samples == i) & (sobjlists$TROP2_ident == "low"))]
}

sobjlists <- sobjlists[order(match(sobjlists$TROP2_ident, c("high", "med", "low"))),]
sobjlists <- sobjlists[order(sobjlists$percent_low, decreasing = F),] # order by HER2
sobjlists <- sobjlists[order(sobjlists$BC.Subtype),] # order by BC Subtype
sobjlists$samples <- factor(sobjlists$samples, levels = unique(sobjlists$samples))

bp <- ggplot(sobjlists, 
             aes(x = samples, 
                 y = percent, 
                 group = as.factor(BC.Subtype),
                 fill = TROP2_ident)) +
  scale_fill_manual(name = "TACSTD2 expression", 
                    values = rev(c("light grey", viridis(2))), 
                    unique(sobjlists$TROP2_ident)) +
  geom_bar(stat = "identity", width = 0.93) +
  geom_text(aes(label = paste0(round(100 - sobjlists$percent_low), "%"), 
                y = 120, 
                fill = NULL), 
            angle = 90,
            colour = "#666666",
            size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 12)) + 
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, 135), 
                     breaks = seq(0,100,25)) + 
  facet_nested( ~ BC.Subtype, 
                scales = "free", 
                space = "free", 
                switch = "x" 
  ) +
  ylab("% Cells") + xlab("Samples (IHC Subtype)") +
  theme(strip.text.x = element_text(size = 14, face = "italic"), 
        strip.background = element_rect(colour = "#FFFFFF", 
                                        size = 1.5, 
                                        fill = "#EEEEEE"), 
        panel.spacing.x = unit(-0.1, "lines"))

ggsave("TROP2_sample_barplot_112322.pdf", plot = bp, width = 14, height = 2)
#ggsave("TROP2_sample_barplot_71422.pdf", plot = bp, width = 14, height = 2)

# TROP2 correlation analysis ---------

Idents(cancer.epi) <- cancer.epi$TROP2_ident

colnames(cancer.epi@meta.data)
data_tot <- cancer.epi@meta.data[,c(4,14,136,137:149)] #"samples","BC.Subtype","TROP2_ident", gene UScores
#data_tot <- cancer.epi@meta.data[,c(4,14,137,102:114)] #"samples","BC.Subtype","TROP2_ident", gene UScores ##Lily

data <- data_tot %>% dplyr::group_by(TROP2_ident) %>% 
  dplyr::summarise(ESR1 = mean(ESR1),
                   ERBB2 = mean(ERBB2),
                   ERBB3 = mean(ERBB3),
                   PIK3CA = mean(PIK3CA),
                   NTRK = mean(NTRK),
                   CD274 = mean(CD274),
                   EGFR = mean(EGFR),
                   FGFR = mean(FGFR),
                   TACSTD2 = mean(TACSTD2),
                   CDK = mean(CDK),
                   AR = mean(AR),
                   NECTIN2 = mean(NECTIN2),
                   LAG3 = mean(LAG3)) 

data <- as.data.frame(data)
rownames(data) <- data$TROP2_ident
data <- data[,-1]

data <- t(apply(data, 2, function(x) (x-mean(x))/sd(x)))

data_pval <- vector() 
data_corval <- vector() 
colnames(data_tot)
for (i in c(4:16)) {
  j <- cor.test(data_tot$TACSTD2, data_tot[,i], 
                alternative = "two.sided",
                method = "pearson",
                exact = NULL, conf.level = 0.95)
  j.p <- j$p.value
  j.c <- j$estimate #correlation coefficient
  names(j.p) <- colnames(data_tot)[i]
  names(j.c) <- colnames(data_tot)[i]
  data_pval <- append(data_pval, j.p)
  data_corval <- append(data_corval, j.c)
  
}

data_pval
data_corval

data <- t(data)
data <- data[,-which(colnames(data) == "TACSTD2")]

p <- Heatmap(data, 
             cluster_rows = F, 
             column_split = 3, 
             row_split = c("high", "med", "low"),
             heatmap_legend_param = list(
               title = "Z-score \nAvg Expression",
               at = c(-2, 0, 2)),
             column_title_gp = gpar(col = c("red", "blue"), fontsize = 0),
             row_title_gp = gpar(fontsize = 0)
)

ggsave("reintegrated_TROP2_corr_heatmap_111822.pdf", plot = as.ggplot(p), width = 5.9, height = 3)  
#ggsave("TROP2_corr_heatmap_111622.pdf", plot = as.ggplot(p), width = 5.9, height = 3) ##Lily 
#ggsave("TROP2_corr_heatmap_71422.pdf", plot = as.ggplot(p), width = 5.9, height = 3)  ##preprint


# TROP2 clinical correlations (SAMPLE LEVEL) --------------

sobjlists <- FetchData(object = cancer.epi, 
                       vars = c("samples",
                                "Patient", 
                                "BC.Subtype", 
                                "TROP2_ident", 
                                "Age", 
                                "Stage", 
                                "Grade",
                                "Tumor.Size", 
                                "nodal_involvement", 
                                "TNM.Classification",
                                "Ki67"))

sobjlists <- sobjlists %>% dplyr::group_by(samples, 
                                           Patient, 
                                           Age,
                                           Stage,
                                           Grade,
                                           Tumor.Size,
                                           nodal_involvement,
                                           Ki67,
                                           BC.Subtype, 
                                           TNM.Classification,
                                           TROP2_ident) %>% 
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 

sobjlists$percent_low <- NA
for (i in unique(sobjlists$samples)) {
  sobjlists$percent_low[which(sobjlists$samples == i)] <- sobjlists$percent[which((sobjlists$samples == i) & (sobjlists$TROP2_ident == "low"))]
}

sobjlists <- sobjlists[order(match(sobjlists$TROP2_ident, c("high", "med", "low"))),]
sobjlists <- sobjlists[order(sobjlists$percent_low, decreasing = F),] # order by TROP2
sobjlists <- sobjlists[order(sobjlists$BC.Subtype),] # order by BC Subtype
sobjlists$samples <- factor(sobjlists$samples, levels = unique(sobjlists$samples))
sobjlists$percent_low <- 100 - sobjlists$percent_low

# Clinical variables relationship

#sobjlists <- sobjlists[which(sobjlists$BC.Subtype == "TNBC"),]
colnames(sobjlists)
age_linreg <- unique(sobjlists[,c(1,3,15)])
age_linreg$Age <- as.numeric(age_linreg$Age)
age_linreg <- age_linreg[which(!is.na(age_linreg$Age)),]

p <- ggscatter(age_linreg, x = "Age", y = "percent_low",
               add = "reg.line",
               conf.int = TRUE,
               add.params = list(color = "blue", fill = "lightgray"),
               xlab = "Age",
               ylab = "% TACSTD2+ cells") + theme_classic() +
  stat_cor(method = "pearson") + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     limits = c(0, 120))
ggsave("test.pdf", plot = p, width = 3, height = 2.5)
#ggsave("TROP2_age_scatter_allBC_71422.pdf", plot = p, width = 3, height = 2.5)

sobjlists$Age[which(sobjlists$Age < 45)] <- "<45yo"
sobjlists$Age[which((sobjlists$Age >= 45) & (sobjlists$Age <= 65))] <- "45-65yo"
sobjlists$Age[which(sobjlists$Age > 65)] <- ">65yo"
sobjlists$Age <- factor(sobjlists$Age, levels = c("<45yo", "45-65yo", ">65yo"))

p <- ggplot(unique(sobjlists[,c(1,3,15)]), aes(x = Age, y = percent_low, fill = Age)) + 
  geom_boxplot(width = 0.5) +
  ylab("% TACSTD2+ cells") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("45-65yo", ">65yo"),
                                        c("<45yo", "45-65yo"),
                                        c("<45yo", ">65yo")), 
                     method="wilcox.test", label="p.format", color="black")

ggsave("TROP2_age_corr_allBC_71422.pdf", plot = p, width = 4, height = 3.5)

sobjlists$N <- NA
sobjlists$N[which(grepl("N0", sobjlists$TNM.Classification))] <- "N0"
sobjlists$N[which(grepl("N1", sobjlists$TNM.Classification))] <- "N1"
sobjlists$N[which(grepl("N2", sobjlists$TNM.Classification))] <- "N2"
sobjlists$N[which(grepl("N3", sobjlists$TNM.Classification))] <- "N3"
sobjlists$N[which(grepl("no", sobjlists$nodal_involvement))] <- "N0"
sobjlists$N[which(grepl("1", sobjlists$nodal_involvement))] <- "N1"
sobjlists$N[which(sobjlists$N %in% c("N1", "N2", "N3"))] <- "N1-3"

p <- ggplot(unique(sobjlists[which(!is.na(sobjlists$N)),c(1,16,15)]), aes(x = N, y = percent_low, fill = N)) + 
  geom_boxplot(width = 0.5) +
  ylab("% TACSTD2+ cells") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("N0", "N1-3")), 
                     method="wilcox.test", label="p.format", color="black")+ 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     limits = c(0, 110))

ggsave("test.pdf", plot = p, width = 3, height = 2.5)
#ggsave("TROP2_N_corr_allBC_71422.pdf", plot = p, width = 3, height = 2.5)


# TROP2 clinical correlations (PATIENT LEVEL) --------------

sobjlists <- FetchData(object = cancer.epi, 
                       vars = c("samples",
                                "Patient", 
                                "BC.Subtype", 
                                "TROP2_ident", 
                                "Age", 
                                "Stage", 
                                "Grade",
                                "Tumor.Size", 
                                "nodal_involvement", 
                                "TNM.Classification",
                                "Ki67"))

sobjlists <- sobjlists %>% dplyr::group_by(Patient, 
                                           Age,
                                           Stage,
                                           Grade,
                                           Tumor.Size,
                                           nodal_involvement,
                                           Ki67,
                                           BC.Subtype, 
                                           TNM.Classification,
                                           TROP2_ident) %>% 
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 

sobjlists$percent_low <- NA
for (i in unique(sobjlists$Patient)) {
  sobjlists$percent_low[which(sobjlists$Patient == i)] <- sobjlists$percent[which((sobjlists$Patient == i) & (sobjlists$TROP2_ident == "low"))]
}

sobjlists <- sobjlists[order(match(sobjlists$TROP2_ident, c("high", "med", "low"))),]
sobjlists <- sobjlists[order(sobjlists$percent_low, decreasing = F),] # order by TROP2
sobjlists <- sobjlists[order(sobjlists$BC.Subtype),] # order by BC Subtype
sobjlists$Patient <- factor(sobjlists$Patient, levels = unique(sobjlists$Patient))
sobjlists$percent_low <- 100 - sobjlists$percent_low

# Clinical variables relationship

#sobjlists <- sobjlists[which(sobjlists$BC.Subtype == "TNBC"),]
colnames(sobjlists)
age_linreg <- unique(sobjlists[,c(1,2,14)])
age_linreg$Age <- as.numeric(age_linreg$Age)
age_linreg <- age_linreg[which(!is.na(age_linreg$Age)),]

p <- ggscatter(age_linreg, x = "Age", y = "percent_low",
               add = "reg.line",
               conf.int = TRUE,
               add.params = list(color = "blue", fill = "lightgray"),
               xlab = "Age",
               ylab = "% TACSTD2+ cells") + theme_classic() +
  stat_cor(method = "pearson") + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     limits = c(0, 120))
ggsave("reintegrated_TROP2_age_scatter_allBC_PATIENTLEVEL_111822.pdf", plot = p, width = 3, height = 2.5)
ggsave("TROP2_age_scatter_allBC_PATIENTLEVEL_111622.pdf", plot = p, width = 3, height = 2.5)
#ggsave("TROP2_age_scatter_allBC_PATIENTLEVEL_101022.pdf", plot = p, width = 3, height = 2.5)

sobjlists$Age[which(sobjlists$Age < 45)] <- "<45yo"
sobjlists$Age[which((sobjlists$Age >= 45) & (sobjlists$Age <= 65))] <- "45-65yo"
sobjlists$Age[which(sobjlists$Age > 65)] <- ">65yo"
sobjlists$Age <- factor(sobjlists$Age, levels = c("<45yo", "45-65yo", ">65yo"))

p <- ggplot(unique(sobjlists[,c(1,2,14)]), aes(x = Age, y = percent_low, fill = Age)) + 
  geom_boxplot(width = 0.5) +
  ylab("% TACSTD2+ cells") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("45-65yo", ">65yo"),
                                        c("<45yo", "45-65yo"),
                                        c("<45yo", ">65yo")), 
                     method="wilcox.test", label="p.format", color="black")

ggsave("TROP2_age_corr_allBC_PATIENTLEVEL_111622.pdf", plot = p, width = 4, height = 3.5)
#ggsave("TROP2_age_corr_allBC_PATIENTLEVEL_101022.pdf", plot = p, width = 4, height = 3.5)

sobjlists$N <- NA
sobjlists$N[which(grepl("N0", sobjlists$TNM.Classification))] <- "N0"
sobjlists$N[which(grepl("N1", sobjlists$TNM.Classification))] <- "N1"
sobjlists$N[which(grepl("N2", sobjlists$TNM.Classification))] <- "N2"
sobjlists$N[which(grepl("N3", sobjlists$TNM.Classification))] <- "N3"
sobjlists$N[which(grepl("no", sobjlists$nodal_involvement))] <- "N0"
sobjlists$N[which(grepl("1", sobjlists$nodal_involvement))] <- "N1"
sobjlists$N[which(sobjlists$N %in% c("N1", "N2", "N3"))] <- "N1-3"

colnames(sobjlists)
p <- ggplot(unique(sobjlists[which(!is.na(sobjlists$N)),c(1,15,14)]), aes(x = N, y = percent_low, fill = N)) + 
  geom_boxplot(width = 0.5) +
  ylab("% TACSTD2+ cells") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("N0", "N1-3")), 
                     method="wilcox.test", label="p.format", color="black")+ 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     limits = c(0, 110))

ggsave("reintegrated_TROP2_N_corr_allBC_PATIENTLEVEL_111822.pdf", plot = p, width = 3, height = 2.5)
ggsave("TROP2_N_corr_allBC_PATIENTLEVEL_111622.pdf", plot = p, width = 3, height = 2.5)
#ggsave("TROP2_N_corr_allBC_PATIENTLEVEL_101022.pdf", plot = p, width = 3, height = 2.5)



# ----------
# ----------
# HER2 analysis -------------

cancer.epi$HER2_ident <- NA
cancer.epi$HER2_ident[(cancer.epi$ERBB2 <= quantile(cancer.epi$ERBB2, 0.975))] <- "med"
cancer.epi$HER2_ident[(cancer.epi$ERBB2 >= quantile(cancer.epi$ERBB2, 0.975))] <- "high"
cancer.epi$HER2_ident[(cancer.epi$ERBB2 <= 0)] <- "low"
# NAsubset <- cancer.epi@meta.data[which(is.na(cancer.epi$HER2_ident) == T),]
# summary(NAsubset$HER2)
#cancer.epi$HER2_ident[is.na(cancer.epi$HER2_ident)] <- "high"

Idents(cancer.epi) <- cancer.epi$HER2_ident

j.markers_DGE <- FindAllMarkers(cancer.epi,
                                slot = "data",
                                min.cells.group = 5,
                                min.pct = 0.2,
                                logfc.threshold = 0,
                                test.use = "MAST")

j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)

write.csv(j.markers_DGE, "reintegrated_HER2_DGEs_0.2minpct0logfcmincells5_111822_YESmito.csv")
#write.csv(j.markers_DGE, "HER2_DGEs_0.2minpct0logfcmincells5_111522_YESmito.csv")

j.markers_DGE_filtered <- DEG_Remove_mito(j.markers_DGE)
write.csv(j.markers_DGE_filtered, "reintegrated_HER2_DGEs_0.2minpct0logfcmincells5_111822_NOmito.csv")
#write.csv(j.markers_DGE_filtered, "HER2_DGEs_0.2minpct0logfcmincells5_111522_NOmito.csv")


j.markers_DGE_nothresh <- FindAllMarkers(cancer.epi,
                                slot = "data",
                                test.use = "MAST")

j.markers_DGE_nothresh$min.pct.diff <- abs(j.markers_DGE_nothresh$pct.1 - j.markers_DGE_nothresh$pct.2)

write.csv(j.markers_DGE_nothresh, "reintegrated_HER2_DGEs_NOTHRESH_111822.csv")
#write.csv(j.markers_DGE_nothresh, "HER2_DGEs_NOTHRESH_111522.csv")


# setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/061822/Fig1")
# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/CancerEpiFinal")
# write.csv(j.markers_DGE_nothresh, "HER2_DGEs_NOTHRESH_71422.csv")
# write.csv(j.markers_DGE, "HER2_DGEs_0.2minpct0logfcmincells5_71422.csv")
# 
# j.markers_DGE_filtered <- read.csv("HER2_DGEs_NOTHRESH_71422.csv", header = T, row.names = 1)
# # j.markers_DGE_filtered <- j.markers_DGE[which(j.markers_DGE$min.pct.diff > 0.1),]
# # j.markers_DGE_filtered <- j.markers_DGE_filtered[which(j.markers_DGE_filtered$p_val_adj < 0.05),]
# j.markers_DGE_filtered <- DEG_Remove_mito(j.markers_DGE_filtered)
# 
# setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only")
# j.markers_DGE_filtered <- read.csv("ERBB2_ABSNOLUTELYNOTHRESH.csv")
# j.markers_DGE_filtered <- DEG_Remove_mito(j.markers_DGE_filtered)
# write.csv(j.markers_DGE_filtered, "ERBB2_NOTHRESHnomito.csv")



# HER2 MA plot ----------

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig1")
j.markers_DGE <- read.csv("HER2_DGEs_NOTHRESH_111522.csv")
head(j.markers_DGE)
rownames(j.markers_DGE) <- j.markers_DGE$X
j.markers_DGE_filtered <- DEG_Remove_mito(j.markers_DGE)

Idents(cancer.epi) <- cancer.epi$HER2_ident
avgexp <- AverageExpression(object = cancer.epi, features = unique(j.markers_DGE_filtered$gene))
avgexp <- as.data.frame(avgexp[['RNA']])

#j.markers_DGE_filtered <- j.markers_DGE_filtered[which(j.markers_DGE_filtered$cluster == "high"),]
#j.markers_DGE_filtered <- j.markers_DGE_filtered[which(j.markers_DGE_filtered$cluster == "med"),]
j.markers_DGE_filtered <- j.markers_DGE_filtered[which(j.markers_DGE_filtered$cluster == "low"),]

avgexp <- avgexp[which(rownames(avgexp) %in% unique(j.markers_DGE_filtered$gene)),]

data <- cbind(j.markers_DGE_filtered, avgexp)
colnames(data)
#data_input <- data[,c(12,3,6,8)] #high
#data_input <- data[,c(11,3,6,8)] #med
data_input <- data[,c(10,3,6,8)] #low
colnames(data_input) <- c("baseMean", "log2FoldChange", "padj", "gene")

options(ggrepel.max.overlaps = 13)
p <- ggmaplot(data_input, 
              fdr = 0.05, fc = 1.5, 
              size = 1,
              palette = c("#B31B21", "#1465AC", "darkgray"),
              genenames = as.vector(data_input$gene),
              top = 100,
              font.label = c("plain", 11),
              select.top.method = "fc",
              ggtheme = ggplot2::theme_classic()) +
  xlim(0,7) + ylim(-2.5,4) + NoLegend()

#ggsave("HER2_MAplot_111822.pdf", plot = p, width = 4.3, height = 4)
#ggsave("HER2_MAplotMEDvall_111822.pdf", plot = p, width = 4.3, height = 4)
ggsave("HER2_MAplotLOWvall_111822.pdf", plot = p, width = 4.3, height = 4)


# ggsave("HER2_MAplot_71422.pdf", plot = p, width = 4.3, height = 4)
# ggsave("HER2_MAplotLOWvall_71822.pdf", plot = p, width = 4.3, height = 4)
# ggsave("HER2_MAplotNEGvall_71822.pdf", plot = p, width = 4.3, height = 4)

# HER2 GO analysis -----------------------------

library(org.Hs.eg.db)
# library(EnsDb.Hsapiens.v79)
# orgdb = "EnsDb.Hsapiens.v79"
orgdb = "org.Hs.eg.db"


dfsample <- read.csv("HER2_DGEs_0.2minpct0logfcmincells5_111522_NOmito.csv", header = T, row.names = 1)
dfsample <- DEG_Remove_mito(dfsample)
#dfsample <- dfsample[which(dfsample$avg_log2FC >= 0.1),]
dfsample <- dfsample[which(abs(dfsample$min.pct.diff) >= 0.1),]
#dfsample <- dfsample[which(dfsample$avg_log2FC > 0),]
dfsample <- split(dfsample$gene, dfsample$cluster)

dfsample$high <- bitr(dfsample$high, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$med <- bitr(dfsample$med, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
dfsample$low <- bitr(dfsample$low, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)

genelist <- list("high" = dfsample$high$ENTREZID,
                 "med" = dfsample$med$ENTREZID,
                 "low" = dfsample$low$ENTREZID)

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


p <- dotplot(GOclusterplot, includeAll = F)
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

ggsave("HER2_hallmark_orghsdb_minpctdiff0.1_111822.pdf", plot = p, width = 15, height = 10, units = "cm")


# ggsave("HER2_hallmark_71022.pdf", plot = p, width = 15, height = 10, units = "cm")
# ggsave("HER2_hallmark_logfc0.1_71122.pdf", plot = p, width = 15, height = 10, units = "cm")
ggsave("HER2_hallmark_orghsdb_minpctdiff0.1_71122.pdf", plot = p, width = 15, height = 10, units = "cm")
#ggsave("HER2_hallmark_minpctdiff0.1_71122.pdf", plot = p, width = 15, height = 10, units = "cm")
#ggsave("test.pdf", plot = p, width = 15, height = 10, units = "cm")


# Barplot of HER2 high cells by sample ------------

sobjlists <- FetchData(object = cancer.epi, 
                       vars = c("samples",
                                "Patient", 
                                "BC.Subtype", 
                                "HER2_ident"))

sobjlists <- sobjlists %>% dplyr::group_by(samples, 
                                           Patient, 
                                           BC.Subtype, 
                                           HER2_ident) %>% 
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 

sobjlists$percent_low <- NA
for (i in unique(sobjlists$samples)) {
  sobjlists$percent_low[which(sobjlists$samples == i)] <- sobjlists$percent[which((sobjlists$samples == i) & (sobjlists$HER2_ident == "low"))]
}

sobjlists <- sobjlists[order(match(sobjlists$HER2_ident, c("high", "med", "low"))),]
sobjlists <- sobjlists[order(sobjlists$percent_low, decreasing = F),] # order by HER2
sobjlists <- sobjlists[order(sobjlists$BC.Subtype),] # order by BC Subtype
sobjlists$samples <- factor(sobjlists$samples, levels = unique(sobjlists$samples))

bp <- ggplot(sobjlists, 
             aes(x = samples, 
                 y = percent, 
                 group = as.factor(BC.Subtype),
                 fill = HER2_ident)) +
  scale_fill_manual(name = "ERBB2 expression", 
                    values = rev(c("light grey", viridis(2))), 
                    unique(sobjlists$HER2_ident)) +
  geom_bar(stat = "identity", width = 0.93) +
  geom_text(aes(label = paste0(round(100 - sobjlists$percent_low), "%"), 
                y = 120, 
                fill = NULL), 
            angle = 90,
            colour = "#666666",
            size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 12)) + 
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, 135), 
                     breaks = seq(0,100,25)) + 
  facet_nested( ~ BC.Subtype, 
                scales = "free", 
                space = "free", 
                switch = "x" 
  ) +
  ylab("% Cells") + xlab("Samples (IHC Subtype)") +
  theme(strip.text.x = element_text(size = 14, face = "italic"), 
        strip.background = element_rect(colour = "#FFFFFF", 
                                        size = 1.5, 
                                        fill = "#EEEEEE"), 
        panel.spacing.x = unit(-0.1, "lines"))

ggsave("HER2_sample_barplot_111622.pdf", plot = bp, width = 14, height = 2)

ggsave("HER2_sample_barplot_111622.pdf", plot = bp, width = 14, height = 2) 
#ggsave("HER2_sample_barplot_71422.pdf", plot = bp, width = 14, height = 2)
#ggsave("HER2_sample_barplot_71022.pdf", plot = bp, width = 14, height = 2)

# HER2 correlation analysis ---------

Idents(cancer.epi) <- cancer.epi$HER2_ident

colnames(cancer.epi@meta.data)
data_tot <- cancer.epi@meta.data[,c(4,14,150, 137:149)]
data_tot <- cancer.epi@meta.data[,c(4,14,138, 102:114)]
#data_tot <- cancer.epi@meta.data[,c(4,14,93:105,107)] ##preprint

data <- data_tot %>% dplyr::group_by(HER2_ident) %>% 
  dplyr::summarise(ESR1 = mean(ESR1),
                   ERBB2 = mean(ERBB2),
                   ERBB3 = mean(ERBB3),
                   PIK3CA = mean(PIK3CA),
                   NTRK = mean(NTRK),
                   CD274 = mean(CD274),
                   EGFR = mean(EGFR),
                   FGFR = mean(FGFR),
                   TACSTD2 = mean(TACSTD2),
                   CDK = mean(CDK),
                   AR = mean(AR),
                   NECTIN2 = mean(NECTIN2),
                   LAG3 = mean(LAG3)) 

data <- as.data.frame(data)
rownames(data) <- data$HER2_ident
data <- data[,-1]

data <- t(apply(data, 2, function(x) (x-mean(x))/sd(x)))
#data <- data[,c(1,3,2,4)]

data_pval <- vector() 
data_corval <- vector() 
colnames(data_tot)
for (i in c(3:15)) {
  j <- cor.test(data_tot$ERBB2, data_tot[,i], 
                alternative = "two.sided",
                method = "pearson",
                exact = NULL, conf.level = 0.95)
  j.p <- j$p.value
  j.c <- j$estimate #correlation coefficient
  names(j.p) <- colnames(data_tot)[i]
  names(j.c) <- colnames(data_tot)[i]
  data_pval <- append(data_pval, j.p)
  data_corval <- append(data_corval, j.c)
  
}

data_pval
data_corval

data <- t(data) 
data <- data[,-which(colnames(data) == "ERBB2")]

p <- Heatmap(data, 
             cluster_rows = F, 
             column_split = 3, 
             row_split = c("high", "med", "low"),
             heatmap_legend_param = list(
               title = "Z-score \nAvg Expression",
               at = c(-2, 0, 2)),
             column_title_gp = gpar(col = c("red", "blue"), fontsize = 0),
             row_title_gp = gpar(fontsize = 0)
)
ggsave("reintegrated_HER2_corr_heatmap_81822.pdf", plot = as.ggplot(p), width = 7.5, height = 3)  
ggsave("HER2_corr_heatmap_81822.pdf", plot = as.ggplot(p), width = 7.5, height = 3)  

ggsave("HER2_corr_heatmap_71422.pdf", plot = as.ggplot(p), width = 7.5, height = 3)  

# HER2 clinical correlations (SAMPLE LEVEL) --------------

sobjlists <- FetchData(object = cancer.epi, 
                       vars = c("samples",
                                "Patient", 
                                "BC.Subtype", 
                                "HER2_ident", 
                                "Age", 
                                "Stage", 
                                "Grade",
                                "Tumor.Size", 
                                "nodal_involvement", 
                                "TNM.Classification",
                                "Ki67"))

sobjlists <- sobjlists %>% dplyr::group_by(samples, 
                                           Patient, 
                                           Age,
                                           Stage,
                                           Grade,
                                           Tumor.Size,
                                           nodal_involvement,
                                           Ki67,
                                           BC.Subtype, 
                                           TNM.Classification,
                                           HER2_ident) %>% 
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 

sobjlists$percent_low <- NA
for (i in unique(sobjlists$samples)) {
  sobjlists$percent_low[which(sobjlists$samples == i)] <- sobjlists$percent[which((sobjlists$samples == i) & (sobjlists$HER2_ident == "low"))]
}

sobjlists <- sobjlists[order(match(sobjlists$HER2_ident, c("high", "med", "low"))),]
sobjlists <- sobjlists[order(sobjlists$percent_low, decreasing = F),] # order by HER2
sobjlists <- sobjlists[order(sobjlists$BC.Subtype),] # order by BC Subtype
sobjlists$samples <- factor(sobjlists$samples, levels = unique(sobjlists$samples))
sobjlists$percent_low <- 100 - sobjlists$percent_low

# Clinical variables relationship

sobjlists <- sobjlists[which(sobjlists$BC.Subtype != "HER2+"),]

colnames(sobjlists)
age_linreg <- unique(sobjlists[,c(1,3,15)])
age_linreg$Age <- as.numeric(age_linreg$Age)
age_linreg <- age_linreg[which(!is.na(age_linreg$Age)),]

p <- ggscatter(age_linreg, x = "Age", y = "percent_low",
               add = "reg.line",
               conf.int = TRUE,
               add.params = list(color = "blue", fill = "lightgray"),
               xlab = "Age",
               ylab = "% ERBB2+ cells") +
  stat_cor(method = "pearson") #+ 
# scale_y_continuous(breaks = c(0, 25, 50),
#                    limits = c(0, 50))
ggsave("test.pdf", plot = p, width = 3, height = 2.5)
#ggsave("HER2_age_scatter_nonHER2_71422.pdf", plot = p, width = 3, height = 2.5)


sobjlists$Age[which(sobjlists$Age < 45)] <- "<45yo"
sobjlists$Age[which((sobjlists$Age >= 45) & (sobjlists$Age <= 65))] <- "45-65yo"
sobjlists$Age[which(sobjlists$Age > 65)] <- ">65yo"
sobjlists$Age <- factor(sobjlists$Age, levels = c("<45yo", "45-65yo", ">65yo"))

colnames(sobjlists)
p <- ggplot(unique(sobjlists[,c(1,3,15)]), aes(x = Age, y = percent_low, fill = Age)) + 
  geom_boxplot(width = 0.5) +
  ylab("% ERBB2+ cells") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("45-65yo", ">65yo"),
                                        c("<45yo", "45-65yo"),
                                        c("<45yo", ">65yo")), 
                     method="wilcox.test", label="p.format", color="black")

ggsave("HER2_age_corr_nonHER2_71422.pdf", plot = p, width = 4, height = 3.5)
#ggsave("HER2_age_corr_nonHER2_71022.pdf", plot = p, width = 4, height = 3.5)

# sobjlists$Stage[which(sobjlists$Stage == "NA")] <- NA
# sobjlists$Grade[which(sobjlists$Grade == "NA")] <- NA

sobjlists$N <- NA
sobjlists$N[which(grepl("N0", sobjlists$TNM.Classification))] <- "N0"
sobjlists$N[which(grepl("N1", sobjlists$TNM.Classification))] <- "N1"
sobjlists$N[which(grepl("N2", sobjlists$TNM.Classification))] <- "N2"
sobjlists$N[which(grepl("N3", sobjlists$TNM.Classification))] <- "N3"
sobjlists$N[which(grepl("no", sobjlists$nodal_involvement))] <- "N0"
sobjlists$N[which(grepl("1", sobjlists$nodal_involvement))] <- "N1"
sobjlists$N[which(sobjlists$N %in% c("N1", "N2", "N3"))] <- "N1-3"

colnames(sobjlists)
p <- ggplot(unique(sobjlists[which(!is.na(sobjlists$N)),c(1,16,15)]), aes(x = N, y = percent_low, fill = N)) + 
  geom_boxplot(width = 0.5) +
  ylab("% ERBB2+ cells") +
  ylim(0,33) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("N0", "N1-3")), 
                     method="wilcox.test", label="p.format", color="black")

ggsave("test.pdf", plot = p, width = 3, height = 2.5)
#ggsave("HER2_N_corr_nonHER2_71422.pdf", plot = p, width = 3, height = 2.5)



# HER2 clinical correlations (PATIENT LEVEL) --------------

sobjlists <- FetchData(object = cancer.epi, 
                       vars = c("samples",
                                "Patient", 
                                "BC.Subtype", 
                                "HER2_ident", 
                                "Age", 
                                "Stage", 
                                "Grade",
                                "Tumor.Size", 
                                "nodal_involvement", 
                                "TNM.Classification",
                                "Ki67"))

sobjlists <- sobjlists %>% dplyr::group_by(Patient, 
                                           Age,
                                           Stage,
                                           Grade,
                                           Tumor.Size,
                                           nodal_involvement,
                                           Ki67,
                                           BC.Subtype, 
                                           TNM.Classification,
                                           HER2_ident) %>% 
  dplyr::summarise(Nb = n()) %>% 
  dplyr::mutate(C = sum(Nb)) %>% 
  dplyr::mutate(percent = Nb/C*100) 

sobjlists$percent_low <- NA
for (i in unique(sobjlists$Patient)) {
  sobjlists$percent_low[which(sobjlists$Patient == i)] <- sobjlists$percent[which((sobjlists$Patient == i) & (sobjlists$HER2_ident == "low"))]
}

sobjlists <- sobjlists[order(match(sobjlists$HER2_ident, c("high", "med", "low"))),]
sobjlists <- sobjlists[order(sobjlists$percent_low, decreasing = F),] # order by HER2
sobjlists <- sobjlists[order(sobjlists$BC.Subtype),] # order by BC Subtype
sobjlists$Patient <- factor(sobjlists$Patient, levels = unique(sobjlists$Patient))
sobjlists$percent_low <- 100 - sobjlists$percent_low

# Clinical variables relationship

sobjlists <- sobjlists[which(sobjlists$BC.Subtype != "HER2+"),]

colnames(sobjlists)
age_linreg <- unique(sobjlists[,c(1,2,14)])
age_linreg$Age <- as.numeric(age_linreg$Age)
age_linreg <- age_linreg[which(!is.na(age_linreg$Age)),]

p <- ggscatter(age_linreg, x = "Age", y = "percent_low",
               add = "reg.line",
               conf.int = TRUE,
               add.params = list(color = "blue", fill = "lightgray"),
               xlab = "Age",
               ylab = "% ERBB2+ cells") +
  stat_cor(method = "pearson") #+ 
# scale_y_continuous(breaks = c(0, 25, 50),
#                    limits = c(0, 50))
ggsave("reintegrated_HER2_age_scatter_nonHER2_PATIENTLEVEL_111822.pdf", plot = p, width = 3, height = 2.5)
ggsave("HER2_age_scatter_nonHER2_PATIENTLEVEL_111822.pdf", plot = p, width = 3, height = 2.5)
#ggsave("HER2_age_scatter_nonHER2_PATIENTLEVEL_101022.pdf", plot = p, width = 3, height = 2.5) ##preprint
#ggsave("HER2_age_scatter_nonHER2_71022.pdf", plot = p, width = 3, height = 2.5)


sobjlists$Age[which(sobjlists$Age < 45)] <- "<45yo"
sobjlists$Age[which((sobjlists$Age >= 45) & (sobjlists$Age <= 65))] <- "45-65yo"
sobjlists$Age[which(sobjlists$Age > 65)] <- ">65yo"
sobjlists$Age <- factor(sobjlists$Age, levels = c("<45yo", "45-65yo", ">65yo"))

colnames(sobjlists)
p <- ggplot(unique(sobjlists[,c(1,2,14)]), aes(x = Age, y = percent_low, fill = Age)) + 
  geom_boxplot(width = 0.5) +
  ylab("% ERBB2+ cells") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("45-65yo", ">65yo"),
                                        c("<45yo", "45-65yo"),
                                        c("<45yo", ">65yo")), 
                     method="wilcox.test", label="p.format", color="black")

ggsave("HER2_age_corr_nonHER2_PATIENTLEVEL_101022.pdf", plot = p, width = 4, height = 3.5)
#ggsave("HER2_age_corr_nonHER2_71022.pdf", plot = p, width = 4, height = 3.5)

# sobjlists$Stage[which(sobjlists$Stage == "NA")] <- NA
# sobjlists$Grade[which(sobjlists$Grade == "NA")] <- NA

sobjlists$N <- NA
sobjlists$N[which(grepl("N0", sobjlists$TNM.Classification))] <- "N0"
sobjlists$N[which(grepl("N1", sobjlists$TNM.Classification))] <- "N1"
sobjlists$N[which(grepl("N2", sobjlists$TNM.Classification))] <- "N2"
sobjlists$N[which(grepl("N3", sobjlists$TNM.Classification))] <- "N3"
sobjlists$N[which(grepl("no", sobjlists$nodal_involvement))] <- "N0"
sobjlists$N[which(grepl("1", sobjlists$nodal_involvement))] <- "N1"
sobjlists$N[which(sobjlists$N %in% c("N1", "N2", "N3"))] <- "N1-3"

colnames(sobjlists)
p <- ggplot(unique(sobjlists[which(!is.na(sobjlists$N)),c(1,15,14)]), aes(x = N, y = percent_low, fill = N)) + 
  geom_boxplot(width = 0.5) +
  ylab("% ERBB2+ cells") +
  ylim(0,40) + #33) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  NoLegend() +
  stat_compare_means(comparisons = list(c("N0", "N1-3")), 
                     method="wilcox.test", label="p.format", color="black")

ggsave("reintegrated_HER2_N_corr_nonHER2_PATIENTLEVEL_111822.pdf", plot = p, width = 3, height = 2.5)
ggsave("HER2_N_corr_nonHER2_PATIENTLEVEL_111822.pdf", plot = p, width = 3, height = 2.5)
#ggsave("HER2_N_corr_nonHER2_PATIENTLEVEL_101022.pdf", plot = p, width = 3, height = 2.5)
#ggsave("HER2_N_corr_nonHER2_71022.pdf", plot = p, width = 3, height = 2.5)




#downstream TROP2 enrichment ======================

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/CancerEpiFinal")
TROP2_downstream <- read.csv(file = "BIOGRID-GENE-110248-4.4.211.tab3_TROP2.csv")
head(TROP2_downstream)

TROP2_downstream.A <- TROP2_downstream$Official.Symbol.Interactor.A
TROP2_downstream.A <- TROP2_downstream.A[TROP2_downstream.A != "TACSTD2"]

TROP2_downstream.B <- TROP2_downstream$Official.Symbol.Interactor.B
TROP2_downstream.B <- TROP2_downstream.B[TROP2_downstream.B != "TACSTD2"]

TROP2_downstream_list <- c(TROP2_downstream.B,TROP2_downstream.B)
TROP2_downstream_list_unique <- unique(TROP2_downstream_list)
length(TROP2_downstream_list_unique)


setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/CancerEpiFinal")
TROP2_downstream <- read.csv(file = "TACSTD2_NCBI_interactions.csv")
head(TROP2_downstream)
TROP2_downstream$Gene <- sub("\xca", "", TROP2_downstream$Gene)
TROP2_downstream_interactants <- TROP2_downstream$Gene
TROP2_downstream_interactants <- TROP2_downstream_interactants[TROP2_downstream_interactants != ""]


TROP2.DEGS <- read.csv("TROP2_DGEs_71422_logfc0mincell5minpct0.2.csv", header = T, row.names = 1)
head(TROP2.DEGS)
TROP2.DEGS <- DEG_Remove_mito(TROP2.DEGS)
TROP2.DEGS.up <- TROP2.DEGS[which(TROP2.DEGS$avg_log2FC >= 0),]

TROP2.high <- TROP2.DEGS.up[TROP2.DEGS.up$cluster == "high",]$gene
TROP2.high.downstream <- TROP2.high[which(TROP2.high %in% TROP2_downstream_list_unique)]
TROP2.high.downstream
(length(unique(TROP2.high.downstream)) / length(unique(TROP2.high))) * 100


TROP2.med <- TROP2.DEGS.up[TROP2.DEGS.up$cluster == "med",]$gene
TROP2.med.downstream <- TROP2.med[which(TROP2.med %in% TROP2_downstream_list_unique)]
TROP2.med.downstream
(length(unique(TROP2.med.downstream)) / length(unique(TROP2.med))) * 100

TROP2.low <- TROP2.DEGS.up[TROP2.DEGS.up$cluster == "low",]$gene
TROP2.low.downstream <- TROP2.low[which(TROP2.low %in% TROP2_downstream_list_unique)]
TROP2.low.downstream
(length(unique(TROP2.low.downstream)) / length(unique(TROP2.low))) * 100


#downstream HER2 enrichment ======================

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/CancerEpiFinal")
HER2_downstream <- read.csv(file = "BIOGRID-GENE-108376-4.4.211.tab3_HER2.csv")
head(HER2_downstream)

HER2_downstream.A <- HER2_downstream$Official.Symbol.Interactor.A
HER2_downstream.A <- HER2_downstream.A[HER2_downstream.A != "ERBB2"]

HER2_downstream.B <- HER2_downstream$Official.Symbol.Interactor.B
HER2_downstream.B <- HER2_downstream.B[HER2_downstream.B != "ERBB2"]

HER2_downstream_list <- c(HER2_downstream.B,HER2_downstream.B)
HER2_downstream_list_unique <- unique(HER2_downstream_list)
length(HER2_downstream_list_unique)


setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/CancerEpiFinal")
HER2_downstream <- read.csv(file = "ERBB2_NCBI_interactions.csv")
head(HER2_downstream)
HER2_downstream$Gene <- sub("\xca", "", HER2_downstream$Gene)
HER2_downstream_interactants <- HER2_downstream$Gene
HER2_downstream_interactants <- HER2_downstream_interactants[HER2_downstream_interactants != ""]



HER2.DEGS <- read.csv("HER2_DGEs_71022.csv", header = T, row.names = 1)
head(HER2.DEGS)
HER2.DEGS.up <- HER2.DEGS[which(HER2.DEGS$avg_log2FC >= 0),]

HER2.high <- HER2.DEGS.up[HER2.DEGS.up$cluster == "high",]$gene
HER2.high.downstream <- HER2.high[which(HER2.high %in% HER2_downstream_list_unique)]
HER2.high.downstream
(length(unique(HER2.high.downstream)) / length(unique(HER2.high))) * 100

HER2.med <- HER2.DEGS.up[HER2.DEGS.up$cluster == "med",]$gene
HER2.med.downstream <- HER2.med[which(HER2.med %in% HER2_downstream_list_unique)]
HER2.med.downstream
(length(unique(HER2.med.downstream)) / length(unique(HER2.med))) * 100

HER2.low <- HER2.DEGS.up[HER2.DEGS.up$cluster == "low",]$gene
HER2.low.downstream <- HER2.low[which(HER2.low %in% HER2_downstream_list_unique)]
HER2.low.downstream
(length(unique(HER2.low.downstream)) / length(unique(HER2.low))) * 100


