setwd("/endosome/work/InternalMedicine/s437775/BreastCancer_Integrated/Analysis/Fig2")
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# This will generate all of the sub-figures for Figure 2, -------
# and is the exploratory analysis done on cancer epithelial cells
# based on expression of 2 clinical targets, TACSTD2/TROP2 and ERBB2/HER2. --------

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

# Functions ------------------------------------------------------

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

setwd("/endosome/work/InternalMedicine/s437775/BreastCancer_Integrated/PrimaryBreastAtlas")
cancer.epi <- readRDS("primarybreastatlas_allcells_v1.1_022723.rds")
cancer.epi <- subset(cancer.epi, subset = Cell_Type_Annotation == "Cancer Epithelial Cells")

DefaultAssay(cancer.epi) <- "RNA"
cancer.epi <- NormalizeData(cancer.epi, assay = "RNA")
head(cancer.epi@meta.data)

# Stratify cells by TACSTD2 expression level ------------

cancer.epi$TROP2_ident <- NA
cancer.epi$TROP2_ident[(cancer.epi$TACSTD2 >= quantile(cancer.epi$TACSTD2,0.9))] <- "high"
cancer.epi$TROP2_ident[(cancer.epi$TACSTD2 <= quantile(cancer.epi$TACSTD2,0.9))] <- "med"
cancer.epi$TROP2_ident[(cancer.epi$TACSTD2 <= 0)] <- "low"

# Stratify cells by ERBB2 expression level --------------
cancer.epi$HER2_ident <- NA
cancer.epi$HER2_ident[(cancer.epi$ERBB2 <= quantile(cancer.epi$ERBB2, 0.975))] <- "med"
cancer.epi$HER2_ident[(cancer.epi$ERBB2 >= quantile(cancer.epi$ERBB2, 0.975))] <- "high"
cancer.epi$HER2_ident[(cancer.epi$ERBB2 <= 0)] <- "low"

# -------------
# -------------
# Fig. 2A: Barplot of TACSTD2 high cells by sample --------------------------

sobjlists <- FetchData(object = cancer.epi, 
                       vars = c("samples", "Patient", "BC.Subtype", "TROP2_ident"))

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
                switch = "x") +
  ylab("% Cells") + xlab("Samples (IHC Subtype)") +
  theme(strip.text.x = element_text(size = 14, face = "italic"), 
        strip.background = element_rect(colour = "#FFFFFF", 
                                        size = 1.5, 
                                        fill = "#EEEEEE"), 
        panel.spacing.x = unit(-0.1, "lines"))

ggsave("TROP2_sample_barplot_112322.pdf", plot = bp, width = 14, height = 2)

# -------------
# -------------
# Fig. 2B: Barplot of ERBB2 high cells by sample --------------------------

sobjlists <- FetchData(object = cancer.epi, 
                       vars = c("samples", "Patient", "BC.Subtype", "TROP2_ident"))

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
                switch = "x") +
  ylab("% Cells") + xlab("Samples (IHC Subtype)") +
  theme(strip.text.x = element_text(size = 14, face = "italic"), 
        strip.background = element_rect(colour = "#FFFFFF", 
                                        size = 1.5, 
                                        fill = "#EEEEEE"), 
        panel.spacing.x = unit(-0.1, "lines"))

ggsave("TROP2_sample_barplot_112322.pdf", plot = bp, width = 14, height = 2)



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

# -------------
# -------------
# Fig. 2C: TACSTD2 and clinical target correlation analysis ---------

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

# -------------
# -------------

# Fig. 2D: ERBB2 and clinical target correlation analysis ---------

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

# -------------
# -------------
# Fig. 2E: DEG analysis for TACSTD2 --------------------------------

Idents(cancer.epi) <- cancer.epi$TROP2_ident

setwd("/endosome/work/InternalMedicine/s437775/BreastCancer_Integrated/Analysis/Fig1")

j.markers_DGE <- FindAllMarkers(cancer.epi,
                                slot = "data",
                                min.cells.group = 5,
                                min.pct = 0.2,
                                logfc.threshold = 0,
                                test.use = "MAST")

j.markers_DGE$min.pct.diff <- abs(j.markers_DGE$pct.1 - j.markers_DGE$pct.2)
j.markers_DGE <- DEG_Remove_mito(j.markers_DGE)
write.csv(j.markers_DGE, "TACSTD2_DGEs_nomitogenes.csv")

# MA plot for DEGs ------------

setwd("/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/Analysis/Fig1")
j.markers_DGE <- read.csv("TROP2_DGEs_ABSOLUTELYNOTHRESH_111522.csv")
head(j.markers_DGE)
rownames(j.markers_DGE) <- j.markers_DGE$X
j.markers_DGE_filtered <- DEG_Remove_mito(j.markers_DGE)

Idents(cancer.epi) <- cancer.epi$TROP2_ident
avgexp <- AverageExpression(object = cancer.epi, features = unique(j.markers_DGE_filtered$gene))
avgexp <- as.data.frame(avgexp[['RNA']])

j.markers_DGE_filtered <- j.markers_DGE_filtered[which(j.markers_DGE_filtered$cluster == "low"),]

avgexp <- avgexp[which(rownames(avgexp) %in% unique(j.markers_DGE_filtered$gene)),]

data <- cbind(j.markers_DGE_filtered, avgexp)
colnames(data)
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

# GO analysis for DEGs -----------------------------

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


# -------------
# -------------
# Fig. 2F: DEG analysis for ERBB2 ----------------------------

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



# MA plot for DEGs ----------

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

# GO analysis for DEGs -----------------------------

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



# -------------
# Fig. 2G: Nodal status correlation for TACSTD2 ------------------

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
for (i in unique(sobjlists$samples)) {
  sobjlists$percent_low[which(sobjlists$samples == i)] <- sobjlists$percent[which((sobjlists$samples == i) & (sobjlists$TROP2_ident == "low"))]
}

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

ggsave("TROP2_N_corr_allBC_71422.pdf", plot = p, width = 3, height = 2.5)

# ----------
# ----------
# Fig. 2H: Nodal status correlation for ERBB2 -----------------

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

# Exclude HER2 clinical subtype patients from analysis 
sobjlists <- sobjlists[which(sobjlists$BC.Subtype != "HER2+"),]

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

# -------------
# -------------
# Fig. 2I/J: Cancer epithelial cell heterogeneity exploration ----------
# Prepare cancer epithelial object -------------------------------------------------------

DefaultAssay(cancer.epi) <- "RNA"
cancer.epi <- NormalizeData(cancer.epi, assay = "RNA")
cancer.epi <- FindVariableFeatures(cancer.epi, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)
all.genes <- rownames(cancer.epi)
cancer.epi <- ScaleData(cancer.epi, features = all.genes)

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

# Quantify ITTH (similarity scores) --------------------------------------------------------------

# Housekeeping genes to exclude from similarity analysis 

# read in list of human housekeeping genes from https://housekeeping.unicamp.br/
housekeeping_genes_cancer.epi <- read.csv("Housekeeping_GenesHuman.csv", sep =";")[,2]
housekeeping_genes_cancer.epi <- unique(housekeeping_genes_cancer.epi[which(housekeeping_genes_cancer.epi %in% rownames(cancer.epi))])
housekeeping_genes_cancer.epi <- which(rownames(cancer.epi) %in% housekeeping_genes_cancer.epi)

# Calculate tumor similarity matrices and save to .csv ------------------------
tumor_score <- data.frame()
tumor_score_sc50 <- data.frame() 
tumor_score_GM <- data.frame()

# generate similarity matrix for all tumor samples 
for (i in samples) {
  j <- as.data.frame(GetAssayData(object = subset(cancer.epi, subset = samples == i), 
                                  assay = "RNA", 
                                  slot = "data"))
  
  # calculate similarity matrix and save to file (for all non-housekeeping genes and SC50 genes only)
  k <- simil(j,
             housekeeping_genes_cancer.epi,
             paste0(gsub(" ", "",i), "_", gsub(" ", "",i),"_global.rds"),
             "corr")
  
  k2 <- simil_sc50(j,
                   paste0(gsub(" ", "",i), "_", gsub(" ", "",i),"_sc50.rds"),
                   "corr")

    # append similarity matrix quantile and mean to tumor_score dataframe
  tumor_score <- rbind(tumor_score,
                       data.frame(name = i, value = as.list(k)))
  
  tumor_score_sc50 <- rbind(tumor_score_sc50,
                            data.frame(name = i, value = as.list(k2)))
}

# save tumor_score to .csv file 
write.csv(tumor_score, "cancerepi_globalITH_all.csv")
write.csv(tumor_score_sc50, "cancerepi_sc50ITH_all.csv")

# ROGUE all cancer epi --------

as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

expr <- as_matrix(GetAssayData(cancer.epi, slot = "data", assay = "RNA"))
meta <- cancer.epi@meta.data
expr <- matr.filter(expr, min.cells = 10, min.genes = 10)
ent.res <- SE_fun(expr)

pdf("cancerepi_all_SEplot.pdf", width = 5, height = 5)
SEplot(ent.res)
dev.off()

rogue.value <- CalculateRogue(ent.res, 
                              platform = "UMI")
rogue.value

rogue.res <- rogue(expr, 
                   labels = meta$BC.Subtype,
                   samples = meta$samples, 
                   platform = "UMI",
                   span = 0.6)
rogue.res

pdf("cancerepi_all_ROGUEboxplot.pdf", width = 5, height = 5)
p
dev.off()

mydata <- rogue.res
for (i in c(1:dim(rogue.res)[1])) {
  mydata[i,1] <- rownames(rogue.res)[i]
  if (length(which(!is.na(rogue.res[i,]))) > 0) {
    mydata[i,2] <- rogue.res[i, which(!is.na(rogue.res[i,]))]
    mydata[i,3] <- colnames(rogue.res)[which(!is.na(rogue.res[i,]))]
  }
}
colnames(mydata) <- c("Patient", "ROGUE", "BC.Subtype")

p <- ggplot(mydata, aes(x = BC.Subtype, y = ROGUE, fill = BC.Subtype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(comparisons = list(c("HER2+", "HR+"), 
                                        c("HR+", "TNBC"), 
                                        c("HER2+", "TNBC")),
                     method="wilcox.test", color="black",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")))

# Create dataframe of tumor ITH scores to plot --------------------------------

scores <- data.frame(unique(sobjlists[,c(1,3)]))
rownames(scores) <- scores$samples

# read in global and sc50 ITH scores
tumor_score <- read.csv("cancerepi_globalITH_all.csv", stringsAsFactors = FALSE, na.strings = "unknown")[,-1]
rownames(tumor_score) <- tumor_score$name

tumor_score_sc50 <- read.csv("cancerepi_sc50ITH_all.csv", stringsAsFactors = FALSE, na.strings = "unknown")[,-1]
rownames(tumor_score_sc50) <- tumor_score_sc50$name

tumor_score_rogue <- as.data.frame(rogue.res)
colnames(tumor_score_rogue) <- c("ROGUE", "name")

scores <- dplyr::left_join(scores, tumor_score_rogue, by = c("samples" = "name"))
scores <- dplyr::left_join(scores, tumor_score, by = c("samples" = "name"))
scores <- dplyr::left_join(scores, tumor_score_sc50, by = c("samples" = "name"))

# create dataframe of scores for plotting
scores_graph <- scores[,c(1, 2, 3, 9, 15)]
scores_graph[,c(4,5)] <- (-1)*scores_graph[,c(4,5)]
scores_graph <- as.data.frame(t(scores_graph)[-1,])
rownames(scores_graph) <- c("BC.Subtype", "ROGUE", "global     ", "sc50")
colnames(scores_graph) <- scores$samples

# normalize ITH scores
scores_graph_rescale <- scores_graph[1,]
temp <- mutate_all(scores_graph[-1,], function(x) as.numeric(as.character(x)))
temp <- t(apply(temp, 1, function(x)(1-((x-min(x))/(max(x)-min(x))))))
scores_graph_rescale <- rbind(scores_graph_rescale, temp)

# order results by ITH score
scores_graph <- scores_graph[-1,]
scores_graph <- mutate_all(scores_graph, as.numeric)

scores_graph_rescale <- scores_graph_rescale[,order(as.numeric(scores_graph_rescale[2,]),
                                                    decreasing = T)]
scores_graph_rescale <- scores_graph_rescale[,order(as.character(scores_graph_rescale[1,]),
                                                    decreasing = F)]
scores_graph_rescale <- scores_graph_rescale[-1,]
scores_graph_rescale <- 1 - mutate_all(scores_graph_rescale, as.numeric)
scores_graph_rescale[scores_graph_rescale > 1] <- 1
scores_graph_rescale[scores_graph_rescale < 0] <- 0

for (i in 1:dim(scores_graph_rescale)[2]) {
  score <- max(sobjlists[which(sobjlists$samples == colnames(scores_graph_rescale)[i]),]$percent)
  scores_graph_rescale[3,i] <- (100 - score)/100
}

scores_graph_rescale[4,] <- abs(scores_graph_rescale[1,] - scores_graph_rescale[3,])
scores_graph_rescale[4,which(scores_graph_rescale[4,] > 0.5)] <- 1
scores_graph_rescale[4,which(scores_graph_rescale[4,] <= 0.5)] <- 0

# Create similarity heatmap ---------------------------------------------------

# lots of custom formatting in the parameters
mydata <- scores_graph_rescale[c(1:4),]

hm <- as.ggplot(ComplexHeatmap::pheatmap(mydata, 
                                         #col = rev(brewer.pal(n = 11, name = "RdBu")), 
                                         col = rocket(500),
                                         breaks = seq(min(mydata[1,]),
                                                      max(mydata[1,]), 
                                                      length.out = 11),
                                         legend = TRUE, 
                                         cellheight = 30,
                                         fontsize = 12,
                                         show_rownames = TRUE, show_colnames = FALSE, 
                                         cluster_rows = FALSE, cluster_cols = FALSE, 
                                         treeheight_col = 0, treeheight_row = 0, 
                                         border_color = "white", 
                                         border_gp = gpar(col = "black", lty = 2),
                                         gaps_col = c(12, 41), 
                                         heatmap_legend_param = list(legend_direction = "horizontal", 
                                                                     legend_width = unit(6.5, "cm"),
                                                                     at = c(min(mydata[1,]), 
                                                                            max(mydata[1,])),
                                                                     labels = c("                    Least\n                    Heterogeneous", 
                                                                                "Most                    \nHeterogeneous                    "), 
                                                                     labels_gp = gpar(fontsize = 12),
                                                                     title_gp = gpar(fontsize = 14, fontface = "plain"),
                                                                     title_position = "topcenter", 
                                                                     title ="Heterogeneity Score"))) + 
  theme(plot.margin = unit(c(0,1,-0.5,2), "cm"), 
        legend.position = "bottom")

# Create stacked % bar plot by subtype label ----------------------------------

# reorder sobjlists 
sobjlists <- sobjlists[order(match(sobjlists$samples, colnames(mydata))), ]
sobjlists$samples <- factor(sobjlists$samples, levels = unique(sobjlists$samples))

sc50_order <- c("Basal_SC", 
                "Her2E_SC",
                "LumA_SC", 
                "LumB_SC")
sobjlists <- sobjlists[order(match(sobjlists$samples, sc50_order)), ]
sobjlists$sc50.Pred <- factor(sobjlists$sc50.Pred, levels = c("Basal_SC", 
                                                              "Her2E_SC",
                                                              "LumA_SC", 
                                                              "LumB_SC"))

# create barplot (again lots of custom formatting)
bp <- ggplot(sobjlists, 
             aes(x = samples, 
                 y = percent, 
                 group = as.factor(BC.Subtype),
                 fill = sc50.Pred)) +
  geom_bar(stat = "identity", width = 0.93) +
  ylab("% Cells (SC50)") + xlab("Samples (IHC Subtype)") +
  theme_minimal() + 
  theme(plot.margin = unit(c(1.2,6.03,-2,0.65), "cm"), 
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) + 
  scale_fill_manual(values = c("#76b5c5", "#ad76c5", "#c58676", "#8dc576"), 
                    sc50_order) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.1)) + 
  theme(legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14)) +
  facet_nested( ~ BC.Subtype, 
                scales = "free", 
                space = "free", 
                switch = "x" 
  ) +
  theme(strip.text.x = element_text(size = 14, face = "italic"), 
        strip.background = element_rect(colour = "#FFFFFF", 
                                        size = 1.5, 
                                        fill = "#EEEEEE"), 
        panel.spacing.x = unit(-0.1, "lines"))

# Combine bar graph and heatmaps ----------------------------------------------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun/simil")

ggsave("cancerepi_barplot_globalITH_ROGUE.pdf", 
       plot = grid.arrange(bp, hm,
                           nrow = 2, ncol = 1, 
                           widths = c(200), 
                           heights = c(10, 15)), 
       width = 60, height = 9.5, units = "cm") # combine bar graph and heatmap with specified graph dimensions


# -------------
# -------------
