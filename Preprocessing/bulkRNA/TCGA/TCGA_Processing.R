


library(TFEA.ChIP)
library(org.Hs.eg.db)
library(TCGAbiolinks)
library(survminer)
library(survival)
library(openxlsx)
library(edgeR)
library(RTCGA)
library(RTCGA.clinical)
library(pheatmap)
library(cetcolor)
library(readxl)

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
library(glmGamPoi)
library(sctransform)#, lib.loc="/opt/R/4.0.2/lib64/R/library")
library(matrixStats)#, lib.loc="/opt/R/4.0.2/lib64/R/library")
library(sparseMatrixStats)
library(DESeq2)
library(genefu)




#downloading TCGA data ==================


#https://www.youtube.com/watch?v=oyAn5B-vLus
GDCprojects = getGDCprojects()
GDCprojects[c("project_id", "name")]



#https://github.com/Sijin-ZhangLab/PanMyeloid/blob/master/Survival_analysis_function.R
clinical <- read.xlsx("TCGA-CDR-SupplementalTableS1.xlsx") #https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018
clinical <- clinical[,c("type","bcr_patient_barcode","gender","race","age_at_initial_pathologic_diagnosis","ajcc_pathologic_tumor_stage","OS","OS.time","DSS","DSS.time","DFI","DFI.time","PFI","PFI.time")]
colnames(clinical) <- c("type","sample","gender","race","age","stage","OS","OS.time","DSS","DSS.time","DFI","DFI.time","PFI","PFI.time")
rownames(clinical) <- clinical$sample


### generate survival dataset
cancers <- c("LUAD","LUSC","KICH","KIRC","KIRP","BRCA","STAD","SKCM","LIHC","COAD","READ","ESCA","OV","PAAD","THCA","UCEC","DLBC")
for(cancerType in cancers){
  ###cancerType <- "COAD"
  #https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/query.html <- for more info on options
  query <- GDCquery(project = paste0("TCGA-",cancerType), 
                    experimental.strategy = "RNA-Seq",
                    data.category = "Transcriptome Profiling", 
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts", 
                    legacy = FALSE) #to get one aligned using hg38
  GDCdownload(query)
  Rnaseq <- GDCprepare(query)
  query_results <- getResults(query)

  saveRDS(Rnaseq,paste0(cancerType,"_SummExpObj.rds"))
  saveRDS(query_results,paste0(cancerType,"_results.rds"))
  
  count_matrix <- assay(Rnaseq)
  samplesNT <- TCGAquery_SampleTypes(barcode = colnames(count_matrix),typesample = c("NT"))
  samplesTP <- TCGAquery_SampleTypes(barcode = colnames(count_matrix),typesample = c("TP"))
  gexp <- count_matrix[,samplesTP]
  
  clinical_subset <- clinical[clinical$type == cancerType,]
  colnames(gexp) <-  sapply(strsplit(colnames(gexp),'-'),function(x) paste0(x[1:3],collapse="-"))
  intersect_sample <- intersect(clinical_subset$sample, colnames(gexp))
  gexp_subset <- gexp[,intersect_sample]
  clinical_subset <- clinical_subset[clinical_subset$sample %in% intersect_sample,]
  
  tmp <- list(expression=gexp_subset,clinical=clinical_subset)
  
  saveRDS(tmp,paste0(cancerType,".rds"))
}

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/TCGA_Analysis/GDCdata")
BRCA <- readRDS(file = "BRCA_SummExpObj.rds")
setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/TCGA_Analysis")

#https://rdrr.io/bioc/TCGAbiolinks/man/TCGAanalyze_SurvivalKM.html and
#https://costalab.ukaachen.de/open_data/Bioinformatics_Analysis_in_R_2019/BIAR_D3/handout.html
# clinical_patient_Cancer <- GDCquery_clinic("TCGA-BRCA","clinical")
# #also pick with high immune infiltrate
# clinical_BRCAsub <- clinical_patient_Cancer[clinical_patient_Cancer$gender=="female", #& clinical_patient_Cancer$prior_treatment=="No", 
#                                             colnames(clinical_patient_Cancer)]
# head(clinical_BRCAsub)
# table(clinical_BRCAsub$gender)
# table(clinical_BRCAsub$prior_treatment)

colnames(colData(BRCA)) #metadata

#BRCA.subset <- BRCA[,((BRCA$gender == "female") & is.na(BRCA$gender == "female") == F)]
BRCA.subset <- BRCA
dim(BRCA.subset)
dim(BRCA)
table(BRCA$gender)
table(BRCA.subset$gender)

head(colData(BRCA)$treatments)
table(BRCA.subset$gender)


#https://costalab.ukaachen.de/open_data/Bioinformatics_Analysis_in_R_2019/BIAR_D3/handout.html#6_gene_expression_and_survival
#dataBRCAcomplete <- log2(assay(BRCA))
dataBRCAcomplete <- assay(BRCA.subset)


geneinfo_df <- rowData(BRCA.subset)

dataBRCAcomplete <- as.data.frame(dataBRCAcomplete)
dataBRCAcomplete$ensembl_gene_id <- rownames(dataBRCAcomplete)
dataBRCAcomplete <- dataBRCAcomplete[,c(1223, 1:1222)]
dataBRCAcomplete[1:5,1:5]

head(geneinfo_df) 
geneinfo_df$ConvertedGenes <- NA
geneinfo_df<- geneinfo_df[,c(2, 4, 1, 3)]
head(geneinfo_df)


dataBRCAcomplete <- merge(geneinfo_df, dataBRCAcomplete, by = "ensembl_gene_id")
dataBRCAcomplete <- as.data.frame(dataBRCAcomplete)
dataBRCAcomplete[1:5,1:5]


#converting genes =====================================

library(limma)
library(org.Hs.eg.db)
oldgenes <- dataBRCAcomplete[,"external_gene_name"]
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

dataBRCAcomplete$ConvertedGenes <- usegenes
new_mat <- dataBRCAcomplete

new_mat <- new_mat[which(!(new_mat$ConvertedGenes %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- dataBRCAcomplete[which(dataBRCAcomplete$ConvertedGenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  k$ConvertedGenes <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}

rownames(new_mat) <- new_mat[,"ConvertedGenes"]
dataBRCAcomplete <- new_mat


dataBRCAcomplete_genedf <- dataBRCAcomplete[,1:4]
head(dataBRCAcomplete_genedf)

rownames(dataBRCAcomplete) <- dataBRCAcomplete$ConvertedGenes
dataBRCAcomplete <- dataBRCAcomplete[,-c(1:4)]

dataBRCAcomplete <- as.data.frame(sapply(dataBRCAcomplete, as.numeric))
rownames(dataBRCAcomplete) <- dataBRCAcomplete_genedf$external_gene_name
colnames(dataBRCAcomplete) <- gsub(".", "-", colnames(dataBRCAcomplete), fixed=TRUE)

saveRDS(dataBRCAcomplete, "TCGA_genesconverted_72722.rds")
saveRDS(dataBRCAcomplete_genedf, "TCGA_genesinfo_withconverted_72722.rds")


