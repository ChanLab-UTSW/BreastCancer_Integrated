

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This is the first step in processing the scRNA-seq  
# primary breast datasets.
# The original datasets can be found at the following GEO numbers/links*:
# *Also found in Table 1.


# Azizi (immune; BC01-BC08): GSE114725; Paper: https://www.sciencedirect.com/science/article/pii/S0092867418307232?via%3Dihub
# Azizi (T; BC09-BC11): GSE114724; Paper: https://www.sciencedirect.com/science/article/pii/S0092867418307232?via%3Dihub
# Karaayvaz: GSE118389; Paper: https://www.nature.com/articles/s41467-018-06052-0
# Pal: GSE161529; Paper: https://www.embopress.org/doi/full/10.15252/embj.2020107333
# Qian: https://lambrechtslab.sites.vib.be/en/pan-cancer-blueprint-tumour-microenvironment-0 
#      Qian's Paper: https://www.nature.com/articles/s41422-020-0355-0 
# Savas: GSE110686; Paper: https://www.nature.com/articles/s41591-018-0078-7
# Wu (2020): https://singlecell.broadinstitute.org/single_cell/study/SCP1106/stromal-cell-diversity-associated-with-immune-evasion-in-human-triple-negative-breast-cancer
#      Wu's Paper: https://www.embopress.org/doi/full/10.15252/embj.2019104063
# Wu (2021): GSE176078; Paper: https://www.nature.com/articles/s41588-021-00911-1#Fig1
# Xu: GSE180286; Paper: https://www.nature.com/articles/s41389-021-00355-6


# The output of this R script will be the input 
# of "DoubletFinder_Primary.R"

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





# ================================================================== ======
#Azizi  ============================================================ =========================================================================
# ================================================================== ======
# Loading Raw Data into RStudio ---------------------------------- 

filePaths = getGEOSuppFiles("GSE114727") 
tarF <- list.files(path = "./GSE114727/", pattern = "*.tar", full.names = TRUE) 
tarF
untar(tarF, exdir = "./GSE114727/") 
gzipF <- list.files(path = "./GSE114727/", pattern = "*.gz", full.names = TRUE) 
ldply(.data = gzipF, .fun = gunzip) 

list.files(path = "./GSE114727/", full.names = TRUE)
list.files(path = "./GSE114727/", pattern = "\\.csv$",full.names = TRUE)
list.files(path = "./GSE114727/", pattern = "\\.mtx$",full.names = TRUE)
list.files(path = "./GSE114727/", pattern = "*.genes.tsv$", full.names = TRUE) 
list.files(path = "./GSE114727/", pattern = "*.barcodes.tsv$", full.names = TRUE)

#3' specific (BC01-BC08)
filePaths = getGEOSuppFiles("GSE114725") 
tarF <- list.files(path = "./GSE114725/", pattern = "*.tar", full.names = TRUE) 
tarF
untar(tarF, exdir = "./GSE114725/") 
gzipF <- list.files(path = "./GSE114725/", pattern = "*.gz", full.names = TRUE) 
ldply(.data = gzipF, .fun = gunzip) 
list.files(path = "./GSE114725/", full.names = TRUE)

#5' specific (BC09-BC11)
filePaths = getGEOSuppFiles("GSE114724") 
tarF <- list.files(path = "./GSE114724/", pattern = "*.tar", full.names = TRUE) 
tarF
untar(tarF, exdir = "./GSE114724/") 
gzipF <- list.files(path = "./GSE114724/", pattern = "*.gz", full.names = TRUE) 
ldply(.data = gzipF, .fun = gunzip) 
list.files(path = "./GSE114724/", full.names = TRUE)
list.files(path = "./GSE114724/", pattern = "\\.mtx$",full.names = TRUE)
list.files(path = "./GSE114724/", pattern = "*.genes.tsv$", full.names = TRUE) 
list.files(path = "./GSE114724/", pattern = "*.barcodes.tsv$", full.names = TRUE) 



# BC01-BC08 raw matrix -----------------------------------------------------------

BC01BC08fullmatold <- read.csv(file = './GSE114725//GSE114725_rna_raw.csv')
BC01BC08fullmatold[1:5,1:5]

#flips rows and columns
BC01BC08fullmatflip <- data.frame(t(BC01BC08fullmatold[-1]))
colnames(BC01BC08fullmatflip) <- BC01BC08fullmatold[,1]
BC01BC08fullmatflip[1:5,1:5]

#adds cell name to patient name and removes pdata from main matrix
colnames(BC01BC08fullmatflip) <-paste(colnames(BC01BC08fullmatflip), BC01BC08fullmatflip[4,], sep = "_")
BC01BC08fullmatflip <- BC01BC08fullmatflip[-c(1,2,3,4),]
BC01BC08fullmatflip[1:5,1:5]

#takes pdata information from original matrix
BC01BC08fullmat_pdat <- data.frame("cells" = colnames(BC01BC08fullmatflip), "index" = BC01BC08fullmatold[,2], "replicate" = BC01BC08fullmatold[,3])
BC01BC08fullmat_pdat[1:3,1:3]

BC01BC08fullmat_pdat$samples <- BC01BC08fullmat_pdat$cells
BC01BC08fullmat_pdat$samples <- sub("\\_.*","",BC01BC08fullmat_pdat$samples)
BC01BC08fullmat_pdat$samples <- paste(BC01BC08fullmat_pdat$samples, BC01BC08fullmat_pdat$index, sep = "_")
BC01BC08fullmat_pdat$samples <- paste(BC01BC08fullmat_pdat$samples, BC01BC08fullmat_pdat$replicate, sep = "_")

#changes flipped name to simpler name 
BC01BC08fullmat <- BC01BC08fullmatflip
BC01BC08fullmat[1:5,1:5]

# BC1 matrix -----------------------------------------------------------

BC1mat <- grep(pattern =c("^BC1") , x = colnames(BC01BC08fullmat), value = TRUE)
BC1mat = BC01BC08fullmat[,grepl(c("^BC1"),colnames(BC01BC08fullmat))]
BC1mat[1:5,1:5]

BC1pdat <- grep(pattern =c("^BC1") , x = BC01BC08fullmat_pdat[,1], value = TRUE)
BC1pdat = BC01BC08fullmat_pdat[grepl(c("^BC1"),BC01BC08fullmat_pdat[,1]),]

#as.data.frame(BC1pdat)
BC1pdat[1:3,1:5]

BC1pdat$Patient <- "BC1"
BC1pdat$Grade <- "1"
BC1$pdat$Stage <- "NA"
BC1pdat$Molecular.Subtype <- "ER+PR+"
BC1pdat$Gene.Coverage <- "3'"
BC1pdat$RNA.Type <-"mRNA"
BC1pdat$Sequencer <- "Illumina HiSeq 2500"
BC1pdat$Library.Preparation <- "inDrop v2"
BC1pdat$Capture.Method <- "inDrop v2"

BC1pdat$BC.Subtype <- "HR+"
BC1pdat$BRCA.Status <- "Negative"
BC1pdat$Grade <- "1"
BC1pdat$Stage <- "NA"
BC1pdat$Age <- "38"
BC1pdat$Gender <- "Female"
BC1pdat$Ethnicity <- "NA"
BC1pdat$Menopause <- "Pre"
BC1pdat$Tumor.Size <- "1 cm"
BC1pdat$Histology <- "Ductal"
BC1pdat$Treatment.Status <- "Naive"

colnames(BC1pdat)[2] <- "Tissue.Source"
BC1pdat$Tissue.Source[which(BC1pdat$Tissue.Source == "BLOOD")] <- "Blood"
BC1pdat$Tissue.Source[which(BC1pdat$Tissue.Source == "NORMAL")] <- "Adjacent Normal Breast"
BC1pdat$Tissue.Source[which(BC1pdat$Tissue.Source == "TUMOR")] <- "Primary Tumor"
BC1pdat[1:3,1:3]

BC1pdat$Status <- BC1pdat$Tissue.Source
table(BC1pdat$Tissue.Source)
BC1pdat$Status[which(BC1pdat$Status == "Blood")] <- "Peripheral"
BC1pdat$Status[which(BC1pdat$Status == "Adjacent Normal Breast")] <- "Adjacent Normal"
BC1pdat$Status[which(BC1pdat$Status == "Primary Tumor")] <- "Primary"



# BC2 matrix -----------------------------------------------------------

BC2mat <- grep(pattern =c("^BC2") , x = colnames(BC01BC08fullmat), value = TRUE)
BC2mat = BC01BC08fullmat[,grepl(c("^BC2"),colnames(BC01BC08fullmat))]
BC2mat[1:5,1:5]

BC2pdat <- grep(pattern =c("^BC2") , x = BC01BC08fullmat_pdat[,1], value = TRUE)
BC2pdat = BC01BC08fullmat_pdat[grepl(c("^BC2"),BC01BC08fullmat_pdat[,1]),]
BC2pdat[1:3,1:3]

BC2pdat$Patient <- "BC2"

BC2pdat$Molecular.Subtype <- "ER+"
BC2pdat$Gene.Coverage <- "3'"
BC2pdat$RNA.Type <-"mRNA"
BC2pdat$Sequencer <- "Illumina HiSeq 2500"
BC2pdat$Library.Preparation <- "inDrop v2"
BC2pdat$Capture.Method <- "inDrop v2"

BC2pdat$BC.Subtype <- "HR+"
BC2pdat$BRCA.Status <- "NA"
BC2pdat$Grade <- "2"
BC2pdat$Stage <- "NA"
BC2pdat$Age <- "60"
BC2pdat$Gender <- "Female"
BC2pdat$Ethnicity <- "NA"
BC2pdat$Menopause <- "Post"
BC2pdat$Tumor.Size <- "3 cm"
BC2pdat$Histology <- "Ductal"
BC2pdat$Treatment.Status <- "Naive"

colnames(BC2pdat)[2] <- "Tissue.Source"
BC2pdat$Tissue.Source[which(BC2pdat$Tissue.Source == "LYMPHNODE")] <- "Lymph Node"
BC2pdat$Tissue.Source[which(BC2pdat$Tissue.Source == "NORMAL")] <- "Adjacent Normal Breast"
BC2pdat$Tissue.Source[which(BC2pdat$Tissue.Source == "TUMOR")] <- "Primary Tumor"
BC2pdat[1:3,1:3]

BC2pdat$Status <- BC2pdat$Tissue.Source
table(BC2pdat$Tissue.Source)
BC2pdat$Status[which(BC2pdat$Status == "Lymph Node")] <- "Metastatic"
BC2pdat$Status[which(BC2pdat$Status == "Adjacent Normal Breast")] <- "Adjacent Normal"
BC2pdat$Status[which(BC2pdat$Status == "Primary Tumor")] <- "Primary"



# BC3 matrix -----------------------------------------------------------

BC3mat <- grep(pattern =c("^BC3") , x = colnames(BC01BC08fullmat), value = TRUE)
BC3mat = BC01BC08fullmat[,grepl(c("^BC3"),colnames(BC01BC08fullmat))]
BC3mat[1:5,1:5]

BC3pdat <- grep(pattern =c("^BC3") , x = BC01BC08fullmat_pdat[,1], value = TRUE)
BC3pdat = BC01BC08fullmat_pdat[grepl(c("^BC3"),BC01BC08fullmat_pdat[,1]),]
BC3pdat[1:3,1:3]

BC3pdat$Patient <- "BC3"

BC3pdat$Molecular.Subtype <- "TNBC"
BC3pdat$Gene.Coverage <- "3'"
BC3pdat$RNA.Type <-"mRNA"
BC3pdat$Sequencer <- "Illumina HiSeq 2500"
BC3pdat$Library.Preparation <- "inDrop v2"
BC3pdat$Capture.Method <- "inDrop v2"

BC3pdat$BC.Subtype <- "TNBC"
BC3pdat$BRCA.Status <- "Negative"
BC3pdat$Grade <- "3"
BC3pdat$Stage <- "NA"
BC3pdat$Age <- "43"
BC3pdat$Gender <- "Female"
BC3pdat$Ethnicity <- "NA"
BC3pdat$Menopause <- "Pre"
BC3pdat$Tumor.Size <- "1.5 cm"
BC3pdat$Histology <- "Ductal"
BC3pdat$Treatment.Status <- "Naive"
colnames(BC3pdat)[2] <- "Tissue.Source"
BC3pdat$Tissue.Source[which(BC3pdat$Tissue.Source == "NORMAL")] <- "Adjacent Normal Breast"
BC3pdat$Tissue.Source[which(BC3pdat$Tissue.Source == "TUMOR")] <- "Primary Tumor"
BC3pdat[1:3,1:3]

BC3pdat$Status <- BC3pdat$Tissue.Source
table(BC3pdat$Tissue.Source)
BC3pdat$Status[which(BC3pdat$Status == "Normal")] <- "Adjacent Normal"
BC3pdat$Status[which(BC3pdat$Status == "Primary Tumor")] <- "Primary"



# BC4 matrix -----------------------------------------------------------

BC4mat <- grep(pattern =c("^BC4") , x = colnames(BC01BC08fullmat), value = TRUE)
BC4mat = BC01BC08fullmat[,grepl(c("^BC4"),colnames(BC01BC08fullmat))]
BC4mat[1:5,1:5]

BC4pdat <- grep(pattern =c("^BC4") , x = BC01BC08fullmat_pdat[,1], value = TRUE)
BC4pdat = BC01BC08fullmat_pdat[grepl(c("^BC4"),BC01BC08fullmat_pdat[,1]),]
BC4pdat[1:3,1:3]

BC4pdat$Patient <- "BC4"
BC4pdat$Molecular.Subtype <- "ER+PR+"
BC4pdat$Gene.Coverage <- "3'"
BC4pdat$RNA.Type <-"mRNA"
BC4pdat$Sequencer <- "Illumina HiSeq 2500"
BC4pdat$Library.Preparation <- "inDrop v2"
BC4pdat$Capture.Method <- "inDrop v2"

BC4pdat$BC.Subtype <- "HR+"
BC4pdat$BRCA.Status <- "Negative"
BC4pdat$Grade <- "1"
BC4pdat$Stage <- "NA"
BC4pdat$Age <- "52"
BC4pdat$Gender <- "Female"
BC4pdat$Ethnicity <- "NA"
BC4pdat$Menopause <- "Pre"
BC4pdat$Tumor.Size <- "2.1 cm"
BC4pdat$Histology <- "Ductal"
BC4pdat$Treatment.Status <- "Naive"
colnames(BC4pdat)[2] <- "Tissue.Source"
BC4pdat$Tissue.Source[which(BC4pdat$Tissue.Source == "BLOOD")] <- "Blood"
BC4pdat$Tissue.Source[which(BC4pdat$Tissue.Source == "TUMOR")] <- "Primary Tumor"
BC4pdat[1:3,1:3]

BC4pdat$Status <- BC4pdat$Tissue.Source
table(BC4pdat$Tissue.Source)
BC4pdat$Status[which(BC4pdat$Status == "Blood")] <- "Peripheral"
BC4pdat$Status[which(BC4pdat$Status == "Primary Tumor")] <- "Primary"



# BC5 matrix -----------------------------------------------------------

BC5mat <- grep(pattern =c("^BC5") , x = colnames(BC01BC08fullmat), value = TRUE)
BC5mat = BC01BC08fullmat[,grepl(c("^BC5"),colnames(BC01BC08fullmat))]
BC5mat[1:5,1:5]

BC5pdat <- grep(pattern =c("^BC5") , x = BC01BC08fullmat_pdat[,1], value = TRUE)
BC5pdat = BC01BC08fullmat_pdat[grepl(c("^BC5"),BC01BC08fullmat_pdat[,1]),]
BC5pdat[1:3,1:3]

BC5pdat$Patient <- "BC5"

BC5pdat$Molecular.Subtype <- "TNBC"
BC5pdat$Gene.Coverage <- "3'"
BC5pdat$RNA.Type <-"mRNA"
BC5pdat$Sequencer <- "Illumina HiSeq 2500"
BC5pdat$Library.Preparation <- "inDrop v2"
BC5pdat$Capture.Method <- "inDrop v2"

BC5pdat$BC.Subtype <- "TNBC"
BC5pdat$BRCA.Status <- "NA"
BC5pdat$Grade <- "3"
BC5pdat$Stage <- "NA"
BC5pdat$Age <- "78"
BC5pdat$Gender <- "Female"
BC5pdat$Ethnicity <- "NA"
BC5pdat$Menopause <- "Post"
BC5pdat$Tumor.Size <- "2 cm"
BC5pdat$Histology <- "Ductal"
BC5pdat$Treatment.Status <- "Naive"
colnames(BC5pdat)[2] <- "Tissue.Source"
BC5pdat$Tissue.Source[which(BC5pdat$Tissue.Source == "TUMOR")] <- "Primary Tumor"
BC5pdat[1:3,1:3]

BC5pdat$Status <- BC5pdat$Tissue.Source
table(BC5pdat$Tissue.Source)
BC5pdat$Status[which(BC5pdat$Status == "Primary Tumor")] <- "Primary"


# BC6 matrix -----------------------------------------------------------

BC6mat <- grep(pattern =c("^BC6") , x = colnames(BC01BC08fullmat), value = TRUE)
BC6mat = BC01BC08fullmat[,grepl(c("^BC6"),colnames(BC01BC08fullmat))]
BC6mat[1:5,1:5]

BC6pdat <- grep(pattern =c("^BC6") , x = BC01BC08fullmat_pdat[,1], value = TRUE)
BC6pdat = BC01BC08fullmat_pdat[grepl(c("^BC6"),BC01BC08fullmat_pdat[,1]),]
BC6pdat[1:3,1:3]

BC6pdat$Patient <- "BC6"

BC6pdat$Molecular.Subtype <- "ER+"
BC6pdat$Gene.Coverage <- "3'"
BC6pdat$RNA.Type <-"mRNA"
BC6pdat$Sequencer <- "Illumina HiSeq 2500"
BC6pdat$Library.Preparation <- "inDrop v2"
BC6pdat$Capture.Method <- "inDrop v2"

BC6pdat$BC.Subtype <- "HR+"
BC6pdat$BRCA.Status <- "NA"
BC6pdat$Grade <- "2"
BC6pdat$Stage <- "NA"
BC6pdat$Age <- "58"
BC6pdat$Gender <- "Female"
BC6pdat$Ethnicity <- "NA"
BC6pdat$Menopause <- "Post"
BC6pdat$Tumor.Size <- "1.3 cm"
BC6pdat$Histology <- "Ductal"
BC6pdat$Treatment.Status <- "Naive"

colnames(BC6pdat)[2] <- "Tissue.Source"
BC6pdat$Tissue.Source[which(BC6pdat$Tissue.Source == "TUMOR")] <- "Primary Tumor"
BC6pdat[1:3,1:3]

BC6pdat$Status <- BC6pdat$Tissue.Source
table(BC6pdat$Tissue.Source)
BC6pdat$Status[which(BC6pdat$Status == "Primary Tumor")] <- "Primary"


# BC7 matrix -----------------------------------------------------------

BC7mat <- grep(pattern =c("^BC7") , x = colnames(BC01BC08fullmat), value = TRUE)
BC7mat = BC01BC08fullmat[,grepl(c("^BC7"),colnames(BC01BC08fullmat))]
BC7mat[1:5,1:5]

BC7pdat <- grep(pattern =c("^BC7") , x = BC01BC08fullmat_pdat[,1], value = TRUE)
BC7pdat = BC01BC08fullmat_pdat[grepl(c("^BC7"),BC01BC08fullmat_pdat[,1]),]
BC7pdat[1:3,1:3]

BC7pdat$Patient <- "BC7"

BC7pdat$Molecular.Subtype <- "HER2+"
BC7pdat$Gene.Coverage <- "3'"
BC7pdat$RNA.Type <-"mRNA"
BC7pdat$Sequencer <- "Illumina HiSeq 2500"
BC7pdat$Library.Preparation <- "inDrop v2"
BC7pdat$Capture.Method <- "inDrop v2"

BC7pdat$BC.Subtype <- "HER2+"
BC7pdat$BRCA.Status <- "NA"
BC7pdat$Grade <- "3"
BC7pdat$Stage <- "NA"
BC7pdat$Age <- "65"
BC7pdat$Gender <- "Female"
BC7pdat$Ethnicity <- "NA"
BC7pdat$Menopause <- "Post"
BC7pdat$Tumor.Size <- "1.2 cm"
BC7pdat$Histology <- "Ductal"
BC7pdat$Treatment.Status <- "Naive"

colnames(BC7pdat)[2] <- "Tissue.Source"
BC7pdat$Tissue.Source[which(BC7pdat$Tissue.Source == "NORMAL")] <- "Adjacent Normal Breast"
BC7pdat$Tissue.Source[which(BC7pdat$Tissue.Source == "TUMOR")] <- "Primary Tumor"
BC7pdat[1:3,1:3]

BC7pdat$Status <- BC7pdat$Tissue.Source
table(BC7pdat$Tissue.Source)
BC7pdat$Status[which(BC7pdat$Status == "Adjacent Normal Breast")] <- "Adjacent Normal"
BC7pdat$Status[which(BC7pdat$Status == "Primary Tumor")] <- "Primary"


# BC8 matrix -----------------------------------------------------------

BC8mat <- grep(pattern =c("^BC8") , x = colnames(BC01BC08fullmat), value = TRUE)
BC8mat = BC01BC08fullmat[,grepl(c("^BC8"),colnames(BC01BC08fullmat))]
BC8mat[1:5,1:5]

BC8pdat <- grep(pattern =c("^BC8") , x = BC01BC08fullmat_pdat[,1], value = TRUE)
BC8pdat = BC01BC08fullmat_pdat[grepl(c("^BC8"),BC01BC08fullmat_pdat[,1]),]
BC8pdat[1:3,1:3]

BC8pdat$Patient <- "BC8"

BC8pdat$Molecular.Subtype <- "TNBC"
BC8pdat$Gene.Coverage <- "3'"
BC8pdat$RNA.Type <-"mRNA"
BC8pdat$Sequencer <- "Illumina HiSeq 2500"
BC8pdat$Library.Preparation <- "inDrop v2"
BC8pdat$Capture.Method <- "inDrop v2"

BC8pdat$BC.Subtype <- "TNBC"
BC8pdat$BRCA.Status <- "NA"
BC8pdat$Grade <- "2"
BC8pdat$Stage <- "NA"
BC8pdat$Age <- "72"
BC8pdat$Gender <- "Female"
BC8pdat$Ethnicity <- "NA"
BC8pdat$Menopause <- "Post"
BC8pdat$Tumor.Size <- "1.3 cm"
BC8pdat$Histology <- "Ductal"
BC8pdat$Treatment.Status <- "Naive"

colnames(BC8pdat)[2] <- "Tissue.Source"
BC8pdat$Tissue.Source[which(BC8pdat$Tissue.Source == "TUMOR")] <- "Primary Tumor"
BC8pdat[1:3,1:3]

BC8pdat$Status <- BC8pdat$Tissue.Source
table(BC8pdat$Tissue.Source)
BC8pdat$Status[which(BC8pdat$Status == "Primary Tumor")] <- "Primary"


#BC9_1 - TCELL ---------------------------------------------------------

BC9mat_T1 <- readMM(file = './GSE114724//GSM3148575_BC09_TUMOR1_matrix.mtx')
BC9mat_T1 <- as.matrix(BC9mat_T1)

BC9genes_T1 <- read.table(file = './GSE114724//GSM3148575_BC09_TUMOR1_genes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
rownames(BC9mat_T1) <- BC9genes_T1[,2]

BC9bar_T1 <- read.table(file = './GSE114724//GSM3148575_BC09_TUMOR1_barcodes.tsv', 
                        sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(BC9mat_T1) <- BC9bar_T1[,1]

colnames(BC9mat_T1) <- paste(colnames(BC9mat_T1), "BC9_T1", sep = "_")
BC9pdat_T1 <- data.frame("cells" = colnames(BC9mat_T1), "replicate" = "1")
BC9pdat_T1$samples <- BC9pdat_T1$cells
BC9pdat_T1$samples <- sub(".*?_", "", BC9pdat_T1$samples)
BC9pdat_T1$samples <- paste(BC9pdat_T1$samples, BC9mat_T1$replicate, sep = "_")

BC9pdat_T1$Patient <- "BC9_1"

BC9pdat_T1$Molecular.Subtype <- "ER+PR+"
BC9pdat_T1$Gene.Coverage <- "5'"
BC9pdat_T1$RNA.Type <-"mRNA"
BC9pdat_T1$Sequencer <- "Illumina NextSeq 500"
BC9pdat_T1$Library.Preparation <- "10x Genomics Chromium Single Cell Immune Profiling Solution"
BC9pdat_T1$Capture.Method <- "10X Genomics Chromium Single Cell 5' Library and Bead Kit"
BC9pdat_T1$Tissue.Source <- "Primary Tumor"
BC9pdat_T1$BC.Subtype <- "HR+"
BC9pdat_T1$Status <- "Primary"
BC9pdat_T1$BRCA.Status <- "NA"
BC9pdat_T1$Grade <- "2"
BC9pdat_T1$Stage <- "NA"
BC9pdat_T1$Age <- "65"
BC9pdat_T1$Gender <- "Female"
BC9pdat_T1$Ethnicity <- "NA"
BC9pdat_T1$Menopause <- "Post"
BC9pdat_T1$Tumor.Size <- "1.7 cm"
BC9pdat_T1$Histology <- "Lobular"
BC9pdat_T1$Treatment.Status <- "Naive"


#BC9_2 - TCELL ---------------------------------------------------------

BC9mat_T2 <- readMM(file = './GSE114724//GSM3148576_BC09_TUMOR2_matrix.mtx')
BC9mat_T2 <- as.matrix(BC9mat_T2)

BC9genes_T2 <- read.table(file = './GSE114724//GSM3148576_BC09_TUMOR2_genes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
rownames(BC9mat_T2) <- BC9genes_T2[,2]

BC9bar_T2 <- read.table(file = './GSE114724//GSM3148576_BC09_TUMOR2_barcodes.tsv', 
                        sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(BC9mat_T2) <- BC9bar_T2[,1]

colnames(BC9mat_T2) <- paste(colnames(BC9mat_T2), "BC9_T2", sep = "_")
BC9pdat_T2 <- data.frame("cells" = colnames(BC9mat_T2), "replicate" = "2")
BC9pdat_T2$samples <- BC9pdat_T2$cells
BC9pdat_T2$samples <- sub(".*?_", "", BC9pdat_T2$samples)
BC9pdat_T2$samples <- paste(BC9pdat_T2$samples, BC9pdat_T2$replicate, sep = "_")

BC9pdat_T2$Patient <- "BC9_2"
BC9pdat_T2$Molecular.Subtype <- "ER+PR+"
BC9pdat_T2$Gene.Coverage <- "5'"
BC9pdat_T2$RNA.Type <-"mRNA"
BC9pdat_T2$Sequencer <- "Illumina NextSeq 500"
BC9pdat_T2$Library.Preparation <- "10X Genomics Chromium Single Cell 5' Library and Bead Kit"
BC9pdat_T2$Capture.Method <- "10X Genomics Chromium Single Cell 5' Library and Bead Kit"
BC9pdat_T2$Tissue.Source <- "Primary Tumor"
BC9pdat_T2$BC.Subtype <- "HR+"
BC9pdat_T2$Status <- "Primary"
BC9pdat_T2$BRCA.Status <- "NA"
BC9pdat_T2$Grade <- "2"
BC9pdat_T2$Stage <- "NA"
BC9pdat_T2$Age <- "65"
BC9pdat_T2$Gender <- "Female"
BC9pdat_T2$Ethnicity <- "NA"
BC9pdat_T2$Menopause <- "Post"
BC9pdat_T2$Tumor.Size <- "1.7 cm"
BC9pdat_T2$Histology <- "Lobular"
BC9pdat_T2$Treatment.Status <- "NA"

#BC10 - TCELL ---------------------------------------------------------

BC10mat_T <- readMM(file = './GSE114724//GSM3148577_BC10_TUMOR1_matrix.mtx')
BC10mat_T <- as.matrix(BC10mat_T)

BC10genes_T <- read.table(file = './GSE114724//GSM3148577_BC10_TUMOR1_genes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
rownames(BC10mat_T) <- BC10genes_T[,2]

BC10bar_T <- read.table(file = './GSE114724//GSM3148577_BC10_TUMOR1_barcodes.tsv', 
                        sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(BC10mat_T) <- BC10bar_T[,1]

colnames(BC10mat_T) <- paste(colnames(BC10mat_T), "BC10_T1", sep = "_")
BC10pdat_T <- data.frame("cells" = colnames(BC10mat_T), "replicate" = "1")
BC10pdat_T$samples <- BC10pdat_T$cells
BC10pdat_T$samples <- sub(".*?_", "", BC10pdat_T$samples)
BC10pdat_T$samples <- paste(BC10pdat_T$samples, BC10pdat_T$replicate, sep = "_")

BC10pdat_T$Patient <- "BC10"
BC10pdat_T$Molecular.Subtype <- "TNBC"
BC10pdat_T$Gene.Coverage <- "5'"
BC10pdat_T$RNA.Type <-"mRNA"
BC10pdat_T$Sequencer <- "Illumina NextSeq 500"
BC10pdat_T$Library.Preparation <- "10X Genomics Chromium Single Cell 5' Library and Bead Kit"
BC10pdat_T$Capture.Method <- "10X Genomics Chromium Single Cell 5' Library and Bead Kit"
BC10pdat_T$Tissue.Source <- "Primary Tumor"
BC10pdat_T$BC.Subtype <- "TNBC"
BC10pdat_T$Status <- "Primary"
BC10pdat_T$BRCA.Status <- "Negative"
BC10pdat_T$Grade <- "3"
BC10pdat_T$Stage <- "NA"
BC10pdat_T$Age <- "73"
BC10pdat_T$Gender <- "Female"
BC10pdat_T$Ethnicity <- "NA"
BC10pdat_T$Menopause <- "Post"
BC10pdat_T$Tumor.Size <- "1.5 cm"
BC10pdat_T$Histology <- "Ductal"
BC10pdat_T$Treatment.Status <- "NA"

#BC11_1 - TCELL ---------------------------------------------------------

BC11mat_T1 <- readMM(file = './GSE114724//GSM3148578_BC11_TUMOR1_matrix.mtx')
BC11mat_T1 <- as.matrix(BC11mat_T1)

BC11genes_T1 <- read.table(file = './GSE114724//GSM3148578_BC11_TUMOR1_genes.tsv', 
                           sep = '\t', header = FALSE, stringsAsFactors = FALSE)
rownames(BC11mat_T1) <- BC11genes_T1[,2]

BC11bar_T1 <- read.table(file = './GSE114724//GSM3148578_BC11_TUMOR1_barcodes.tsv', 
                         sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(BC11mat_T1) <- BC11bar_T1[,1]

colnames(BC11mat_T1) <- paste(colnames(BC11mat_T1), "BC11_T1", sep = "_")
BC11pdat_T1 <- data.frame("cells" = colnames(BC11mat_T1), "replicate" = "1")
BC11pdat_T1$samples <- BC11pdat_T1$cells
BC11pdat_T1$samples <- sub(".*?_", "", BC11pdat_T1$samples)
BC11pdat_T1$samples <- paste(BC11pdat_T1$samples, BC11pdat_T1$replicate, sep = "_")

BC11pdat_T1$Patient <- "BC11_1"

BC11pdat_T1$Molecular.Subtype <- "HER2+"
BC11pdat_T1$Gene.Coverage <- "5'"
BC11pdat_T1$RNA.Type <-"mRNA"
BC11pdat_T1$Sequencer <- "Illumina NextSeq 500"
BC11pdat_T1$Library.Preparation <- "10X Genomics Chromium Single Cell 5' Library and Bead Kit"
BC11pdat_T1$Capture.Method <- "10X Genomics Chromium Single Cell 5' Library and Bead Kit"
BC11pdat_T1$Tissue.Source <- "Primary Tumor"
BC11pdat_T1$BC.Subtype <- "HER2+"
BC11pdat_T1$Status <- "Primary"
BC11pdat_T1$BRCA.Status <- "NA"
BC11pdat_T1$Grade <- "3"
BC11pdat_T1$Stage <- "NA"
BC11pdat_T1$Age <- "50"
BC11pdat_T1$Gender <- "Female"
BC11pdat_T1$Ethnicity <- "NA"
BC11pdat_T1$Menopause <- "Post"
BC11pdat_T1$Tumor.Size <- "2.2 cm"
BC11pdat_T1$Histology <- "Ductal"
BC11pdat_T1$Treatment.Status <- "NA"

#BC11_2 - TCELL ---------------------------------------------------------

BC11mat_T2 <- readMM(file = './GSE114724//GSM3148579_BC11_TUMOR2_matrix.mtx')
BC11mat_T2 <- as.matrix(BC11mat_T2)

BC11genes_T2 <- read.table(file = './GSE114724//GSM3148579_BC11_TUMOR2_genes.tsv', 
                           sep = '\t', header = FALSE, stringsAsFactors = FALSE)
rownames(BC11mat_T2) <- BC11genes_T2[,2]

BC11bar_T2 <- read.table(file = './GSE114724//GSM3148579_BC11_TUMOR2_barcodes.tsv', 
                         sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(BC11mat_T2) <- BC11bar_T2[,1]

colnames(BC11mat_T2) <- paste(colnames(BC11mat_T2), "BC11_T2", sep = "_")
BC11pdat_T2 <- data.frame("cells" = colnames(BC11mat_T2), "replicate" = "2")
BC11pdat_T2$samples <- BC11pdat_T2$cells
BC11pdat_T2$samples <- sub(".*?_", "", BC11pdat_T2$samples)
BC11pdat_T2$samples <- paste(BC11pdat_T2$samples, BC11pdat_T1$replicate, sep = "_")

BC11pdat_T2$Patient <- "BC11_2"

BC11pdat_T2$Molecular.Subtype <- "HER2+"
BC11pdat_T2$Gene.Coverage <- "5'"
BC11pdat_T2$RNA.Type <-"mRNA"
BC11pdat_T2$Sequencer <- "Illumina NextSeq 500"
BC11pdat_T2$Library.Preparation <- "10X Genomics Chromium Single Cell 5' Library and Bead Kit"
BC11pdat_T2$Capture.Method <- "10X Genomics Chromium Single Cell 5' Library and Bead Kit"
BC11pdat_T2$Tissue.Source <- "Primary Tumor"
BC11pdat_T2$BC.Subtype <- "HER2+"
BC11pdat_T2$Status <- "Primary"
BC11pdat_T2$BRCA.Status <- "NA"
BC11pdat_T2$Grade <- "3"
BC11pdat_T2$Stage <- "NA"
BC11pdat_T2$Age <- "50"
BC11pdat_T2$Gender <- "Female"
BC11pdat_T2$Ethnicity <- "NA"
BC11pdat_T2$Menopause <- "Post"
BC11pdat_T2$Tumor.Size <- "2.2 cm"
BC11pdat_T2$Histology <- "Ductal"
BC11pdat_T2$Treatment.Status <- "NA"

# joined (immune) =================================================================

joined_immune <- cbind(BC1mat, BC2mat, BC3mat, BC4mat, BC5mat, 
                       BC6mat, BC7mat, BC8mat)
joined_immune <- as.data.frame(joined_immune)

library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(joined_immune)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- joined_immune[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- joined_immune[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}

joined_immune <- new_mat
joined_immune[1:5,1:5]

# joined (T) =================================================================


joined_T <- cbind(BC9mat_T1, BC9mat_T2, BC10mat_T, BC11mat_T1, BC11mat_T2)
joined_T <- as.data.frame(joined_T)

library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(joined_T)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- joined_T[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- joined_T[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}

joined_T <- new_mat

#pdata ==============================================

pdat_immune <- rbind(BC1pdat, BC2pdat, BC3pdat, BC4pdat, BC5pdat,
                     BC6pdat, BC7pdat, BC8pdat)

pdat_T <- rbind(BC9pdat_T1, BC9pdat_T2, BC10pdat_T, BC11pdat_T1, BC11pdat_T2)



#sobj object creation =====================================================================

#immune _________________________________________________________________
rownames(pdat_immune) <- pdat_immune$cells
fdat_immune <- toupper(as.matrix(rownames(joined_immune)))
#_________________________________________________________________
rownames(fdat_immune) <- fdat_immune[,1]
fdat_immune <- data.frame(fdat_immune)
common_colnames <- "gene_short_name"
colnames(fdat_immune) <- common_colnames
rownames(joined_immune) <- rownames(fdat_immune)
#_________________________________________________________________
sobj_pre1 <- CreateSeuratObject(counts = joined_immune)
pdat_immune <- as.data.frame(pdat_immune)
sobj_pre1 <-AddMetaData(sobj_pre1,metadata=pdat_immune)
head(sobj_pre1@meta.data)
#_________________________________________________________________
sobj_pre1[["RNA"]]@meta.features<-fdat_immune
head(sobj_pre1[["RNA"]]@meta.features)
slotNames(sobj_pre1[["RNA"]])


#T cells _________________________________________________________________
rownames(pdat_T) <- pdat_T$cells
fdat_T <- toupper(as.matrix(rownames(joined_T)))
#_________________________________________________________________
rownames(fdat_T) <- fdat_T[,1]
fdat_T <- data.frame(fdat_T)
common_colnames <- "gene_short_name"
colnames(fdat_T) <- common_colnames
rownames(joined_T) <- rownames(fdat_T)
#_________________________________________________________________
sobj_pre2 <- CreateSeuratObject(counts = joined_T)
sobj_pre2 <-AddMetaData(sobj_pre2,metadata=pdat_T)
head(sobj_pre2@meta.data)
#_________________________________________________________________
sobj_pre2[["RNA"]]@meta.features<-fdat_T
head(sobj_pre2[["RNA"]]@meta.features)
slotNames(sobj_pre2[["RNA"]])


sobj1 <- sobj_pre1
sobj2 <- sobj_pre2


saveRDS(sobj2, file = "AziziTcells_unfiltered_62622.rds")


sobj1 <- subset(x = sobj1, subset = Tissue.Source == "Blood", invert = TRUE)
sobj1 <- subset(x = sobj1, subset = Tissue.Source == "Adjacent Normal Breast", invert = TRUE)
sobj1_prim <- subset(x = sobj1, subset = Tissue.Source == "Lymph Node", invert = TRUE)



saveRDS(sobj1_prim, file = "AziziPRIMARYimmune_unfiltered_62622.rds")


sobj1 <- sobj1_prim



#cell cycle ========================================================================
#from: https://satijalab.org/seurat/articles/cell_cycle_vignette.html

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sobj1 <- CellCycleScoring(sobj1, s.features = s.genes, g2m.features = g2m.genes)
sobj1$CC.Difference <- sobj1$S.Score - sobj1$G2M.Score
# view cell cycle scores and phase assignments
head(sobj1[[]])

sobj2 <- CellCycleScoring(sobj2, s.features = s.genes, g2m.features = g2m.genes)
sobj2$CC.Difference <- sobj2$S.Score - sobj2$G2M.Score
# view cell cycle scores and phase assignments
head(sobj2[[]])


# Preprocessing and Filtering (SCTransform) --------------------------------------------------------------------

mito.genes <- grep(pattern = "^MT\\.", x = rownames(sobj1@assays[["RNA"]]), value = TRUE)
print(mito.genes)
sobj1 <- PercentageFeatureSet(sobj1, pattern = "^MT\\.", col.name = "percent.mt")

#blood cell contamination
hemo.genes <- grep(pattern = "^HB[^(P)]", x = rownames(sobj1@assays[["RNA"]]), value = TRUE)
print(hemo.genes)
sobj1 <- PercentageFeatureSet(sobj1, "^HB[^(P)]", col.name = "percent.hb")
platelet.genes <- grep(pattern = "^PECAM1|^PF4\\b", x = rownames(sobj1@assays[["RNA"]]), value = TRUE)
print(platelet.genes)

sobj1 <- PercentageFeatureSet(sobj1, "^PECAM1|^PF4\\b", col.name = "percent.platelet")

nCount95 <- quantile(sobj1@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
nCount5 <- quantile(sobj1@meta.data$nCount_RNA, 0.05)#05)  
nFeat95 <- quantile(sobj1@meta.data$nFeature_RNA, 0.95) #0.9) 
nFeat5 <- quantile(sobj1@meta.data$nFeature_RNA, 0.05) #0.1) 
ptMt5 <- quantile(sobj1@meta.data$percent.mt, 0.05)
ptMt90 <- quantile(sobj1@meta.data$percent.mt, 0.95)
ptHb95 <- quantile(sobj1@meta.data$percent.hb, 0.95)



sobj1 <- subset(x = sobj1, subset = nFeature_RNA > nFeat5 & nFeature_RNA < nFeat95)
sobj1 <- subset(x = sobj1, subset = percent.mt > ptMt5 & percent.mt < ptMt90)
sobj1 <- subset(x = sobj1, subset = nCount_RNA > nCount5 & nCount_RNA < nCount95)
sobj1 <- subset(x = sobj1, subset = percent.hb < ptHb95)

summary(sobj1@meta.data$percent.mt)
summary(sobj1@meta.data$percent.hb)
dim(sobj1)
dim(sobj1_prim)
dim(sobj2)


colnames(sobj1@meta.data)


saveRDS(sobj1, file = "AziziPRIMimmune_62622.rds")




mito.genes <- grep(pattern = "^MT\\.", x = rownames(sobj2@assays[["RNA"]]), value = TRUE)
print(mito.genes)
sobj2 <- PercentageFeatureSet(sobj2, pattern = "^MT\\.", col.name = "percent.mt")

hemo.genes <- grep(pattern = "^HB[^(P)]", x = rownames(sobj2@assays[["RNA"]]), value = TRUE)
print(hemo.genes)
sobj2 <- PercentageFeatureSet(sobj2, "^HB[^(P)]", col.name = "percent.hb")
platelet.genes <- grep(pattern = "^PECAM1|^PF4\\b", x = rownames(sobj2@assays[["RNA"]]), value = TRUE)
print(platelet.genes)

sobj2 <- PercentageFeatureSet(sobj2, "^PECAM1|^PF4\\b", col.name = "percent.platelet")

nCount95 <- quantile(sobj2@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
nCount5 <- quantile(sobj2@meta.data$nCount_RNA, 0.05)#05)  
nFeat95 <- quantile(sobj2@meta.data$nFeature_RNA, 0.95) #0.9) 
nFeat5 <- quantile(sobj2@meta.data$nFeature_RNA, 0.05) #0.1) 
ptMt5 <- quantile(sobj2@meta.data$percent.mt, 0.05)
ptMt90 <- quantile(sobj2@meta.data$percent.mt, 0.95)
ptHb95 <- quantile(sobj2@meta.data$percent.hb, 0.95)



sobj2 <- subset(x = sobj2, subset = nFeature_RNA > nFeat5 & nFeature_RNA < nFeat95)
sobj2 <- subset(x = sobj2, subset = percent.mt > ptMt5 & percent.mt < ptMt90)
sobj2 <- subset(x = sobj2, subset = nCount_RNA > nCount5 & nCount_RNA < nCount95)
sobj2 <- subset(x = sobj2, subset = percent.hb < ptHb95)

summary(sobj2@meta.data$percent.mt)
summary(sobj2@meta.data$percent.hb)
dim(sobj2)
dim(sobj2)



saveRDS(sobj2, file = "AziziT_filtered_62622.rds")



# ================================================================== =======
# Karaayvaz ======================================================== ====================================
# ================================================================== ======

# Loading Raw Data into RStudio ---------------------------------- 

filePaths = getGEOSuppFiles("GSE118389") 
tarF <- list.files(path = "./GSE118389/", pattern = "*.tar", full.names = TRUE) 
untar(tarF, exdir = "./GSE118389/") 
gzipF <- list.files(path = "./GSE118389/", pattern = "*.gz", full.names = TRUE) 
ldply(.data = gzipF, .fun = gunzip) 

list.files(path = "./GSE118389/", full.names = TRUE)


# full matrix ----------------------------------------------------------


fullmat <- read.table(file = './GSE118389//GSE118389_counts_rsem.txt', 
                      sep = '\t', header = TRUE, stringsAsFactors = FALSE)
fullmat[1:5,1:5]
dim(fullmat)
class(fullmat)


# PT039 matrix -----------------------------------------------------------
PT039mat <- grep(pattern =c("^PT039") , x = colnames(fullmat), value = TRUE)
PT039mat = fullmat[,grepl(c("^PT039"),colnames(fullmat))]
dim(PT039mat)
PT039mat[1:5,1:5]


PT039pdat <- data.frame("cells" = colnames(PT039mat), 
                        "lymphovascular_invasion" = "no", 
                        "nodal_involvement" = "no",
                        "BRCA_status" = "Negative", "Tissue Source" = "Primary Tumor", 
                        "BC Subtype" = "TNBC", "Status" = "Primary",
                        "Molecular.Subtype" = "TNBC", "Stage" = "NA", "Grade" = "3",
                        "Tumor.Size" = "9.5", "Age" = "64", "Gender" = "Female", 
                        "Parity" = "NA", "Gravidity" = "NA", "Menopause" = "NA",
                        "Treatment.Status" = "Naive")
#PT039_batch <- list(c("B1", "B3", "B8", "B9"))
#PT039pdat$batch <- PT039_batch
PT039pdat$samples <- PT039pdat$cells
PT039pdat$samples <- gsub("(.+?_.+?)_.*" ,"\\1", PT039pdat$samples)

PT039pdat$Patient <- "PT039"
PT039pdat$Gene.Coverage <- "Full-Length"
PT039pdat$RNA.Type <-"mRNA"
PT039pdat$Sequencer <- "Illumina Genome Analyzer"
PT039pdat$Library.Preparation <- "Smart-Seq2"
PT039pdat$Capture.Method <- "Smart-Seq2"


# PT058 matrix -----------------------------------------------------------
PT058mat <- grep(pattern =c("^PT058") , x = colnames(fullmat), value = TRUE)
PT058mat = fullmat[,grepl(c("^PT058"),colnames(fullmat))]
dim(PT058mat)
PT058mat[1:5,1:5]


PT058pdat <- data.frame("cells" = colnames(PT058mat), 
                        "lymphovascular_invasion" = "yes", 
                        "nodal_involvement" = "no",
                        "BRCA_status" = "NA", "Tissue Source" = "Primary Tumor", 
                        "BC Subtype" = "TNBC", "Status" = "Primary",
                        "Molecular.Subtype" = "TNBC", "Stage" = "NA", "Grade" = "3",
                        "Tumor.Size" = "1.5", "Age" = "46", "Gender" = "Female",
                        "Parity" = "NA", "Gravidity" = "NA", "Menopause" = "NA",
                        "Treatment.Status" = "Naive")
PT058pdat$samples <- PT058pdat$cells
PT058pdat$samples <- gsub("(.+?_.+?)_.*" ,"\\1", PT058pdat$samples)

#PT039pdat$batch <- "PT039_batch"B2
PT058pdat$Patient <- "PT058"
PT058pdat$Gene.Coverage <- "Full-Length"
PT058pdat$RNA.Type <-"mRNA"
PT058pdat$Sequencer <- "Illumina Genome Analyzer"
PT058pdat$Library.Preparation <- "Smart-Seq2"
PT058pdat$Capture.Method <- "Smart-Seq2"

# PT081 matrix -----------------------------------------------------------
PT081mat <- grep(pattern =c("^PT081") , x = colnames(fullmat), value = TRUE)
PT081mat = fullmat[,grepl(c("^PT081"),colnames(fullmat))]
dim(PT081mat)
PT081mat[1:5,1:5]


PT081pdat <- data.frame("cells" = colnames(PT081mat), 
                        "lymphovascular_invasion" = "no", 
                        "nodal_involvement" = "no",
                        "BRCA_status" = "Negative", "Tissue Source" = "Primary Tumor", 
                        "BC Subtype" = "TNBC", "Status" = "Primary",
                        "Molecular.Subtype" = "TNBC", "Stage" = "NA", "Grade" = "3",
                        "Tumor.Size" = "2.3", "Age" = "45", "Gender" = "Female",
                        "Parity" = "NA", "Gravidity" = "NA", "Menopause" = "NA",
                        "Treatment.Status" = "Naive")
PT081pdat$samples <- PT081pdat$cells
PT081pdat$samples <- gsub("(.+?_.+?)_.*" ,"\\1", PT081pdat$samples)

#PT081_batch <- list(c("B3", "B4"))
#PT081pdat$batch <- PT081_batch

PT081pdat$Patient <- "PT081"
PT081pdat$Gene.Coverage <- "Full-Length"
PT081pdat$RNA.Type <-"mRNA"
PT081pdat$Sequencer <- "Illumina Genome Analyzer"
PT081pdat$Library.Preparation <- "Smart-Seq2"
PT081pdat$Capture.Method <- "Smart-Seq2"


# PT084 matrix -----------------------------------------------------------
PT084mat <- grep(pattern =c("^PT084") , x = colnames(fullmat), value = TRUE)
PT084mat = fullmat[,grepl(c("^PT084"),colnames(fullmat))]
dim(PT084mat)
PT084mat[1:5,1:5]


PT084pdat <- data.frame("cells" = colnames(PT084mat), 
                        "lymphovascular_invasion" = "yes", 
                        "nodal_involvement" = "yes(1)",
                        "BRCA_status" = "BRCA2+", "Tissue Source" = "Primary Tumor", 
                        "BC Subtype" = "TNBC", "Status" = "Primary",
                        "Molecular.Subtype" = "TNBC", "Stage" = "NA", "Grade" = "3",
                        "Tumor.Size" = "2.3", "Age" = "52", "Gender" = "Female",
                        "Parity" = "NA", "Gravidity" = "NA", "Menopause" = "NA",
                        "Treatment.Status" = "Naive")
#PT084_batch <- list(c("B5", "B6", "B8"))
#PT084pdat$batch <- PT084_batch
PT084pdat$samples <- PT084pdat$cells
PT084pdat$samples <- gsub("(.+?_.+?)_.*" ,"\\1", PT084pdat$samples)

PT084pdat$Patient <- "PT084"
PT084pdat$Gene.Coverage <- "Full-Length"
PT084pdat$RNA.Type <-"mRNA"
PT084pdat$Sequencer <- "Illumina Genome Analyzer"
PT084pdat$Library.Preparation <- "Smart-Seq2"
PT084pdat$Capture.Method <- "Smart-Seq2"


# PT089 matrix -----------------------------------------------------------
PT089mat <- grep(pattern =c("^PT089") , x = colnames(fullmat), value = TRUE)
PT089mat = fullmat[,grepl(c("^PT089"),colnames(fullmat))]
dim(PT089mat)
PT089mat[1:5,1:5]


PT089pdat <- data.frame("cells" = colnames(PT089mat), 
                        "lymphovascular_invasion" = "no", 
                        "nodal_involvement" = "no",
                        "BRCA_status" = "Negative", "Tissue Source" = "Primary Tumor", 
                        "BC Subtype" = "TNBC", "Status" = "Primary",
                        "Molecular.Subtype" = "TNBC", "Stage" = "NA", "Grade" = "2",
                        "Tumor.Size" = "1.5", "Age" = "44", "Gender" = "Female", 
                        "Parity" = "NA", "Gravidity" = "NA", "Menopause" = "NA",
                        "Treatment.Status" = "Naive")
#PT089_batch <- list(c("B3", "B5", "B6"))
#PT089pdat$batch <- PT089_batch
PT089pdat$samples <- PT089pdat$cells
PT089pdat$samples <- gsub("(.+?_.+?)_.*" ,"\\1", PT089pdat$samples)

PT089pdat$Patient <- "PT089"

PT089pdat$Gene.Coverage <- "Full-Length"
PT089pdat$RNA.Type <-"mRNA"
PT089pdat$Sequencer <- "Illumina Genome Analyzer"
PT089pdat$Library.Preparation <- "Smart-Seq2"
PT089pdat$Capture.Method <- "Smart-Seq2"

# PT126 matrix -----------------------------------------------------------
PT126mat <- grep(pattern =c("^PT126") , x = colnames(fullmat), value = TRUE)
PT126mat = fullmat[,grepl(c("^PT126"),colnames(fullmat))]
dim(PT126mat)
PT126mat[1:5,1:5]


PT126pdat <- data.frame("cells" = colnames(PT126mat), 
                        "lymphovascular_invasion" = "yes", 
                        "nodal_involvement" = "yes(1)",
                        "BRCA_status" = "BRCA2_VUS", "Tissue Source" = "Primary Tumor", 
                        "BC Subtype" = "TNBC", "Status" = "Primary",
                        "Molecular.Subtype" = "TNBC", "Stage" = "NA", "Grade" = "3",
                        "Tumor.Size" = "1.7", "Age" = "42", "Gender" = "Female",
                        "Parity" = "NA", "Gravidity" = "NA", "Menopause" = "NA",
                        "Treatment.Status" = "Naive")
#PT039pdat$batch <- "B7"
PT126pdat$samples <- PT126pdat$cells
PT126pdat$samples <- gsub("(.+?_.+?)_.*" ,"\\1", PT126pdat$samples)

PT126pdat$Patient <- "PT126"

PT126pdat$Gene.Coverage <- "Full-Length"
PT126pdat$RNA.Type <-"mRNA"
PT126pdat$Sequencer <- "Illumina Genome Analyzer"
PT126pdat$Library.Preparation <- "Smart-Seq2"
PT126pdat$Capture.Method <- "Smart-Seq2"



#Karaayvaz Joined ===============================

joined <- cbind(PT039mat, PT058mat, PT081mat, PT084mat, PT089mat, PT126mat)
joined <- as.data.frame(joined)


library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(joined)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- joined[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- joined[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}

joined <- new_mat

#pdata ====================================

pdat_joined <- rbind(PT039pdat, PT058pdat, PT081pdat, PT084pdat, PT089pdat, PT126pdat)
head(pdat_joined, 5)


rownames(pdat_joined) <- pdat_joined$cells 
fdat <- as.data.frame(rownames(joined))
colnames(fdat) <- "gene_short_name"
rownames(fdat) <- fdat$gene_short_name
#rownames(fdat) = make.names(fdat$V1, unique=TRUE)
#rownames(joined) <- rownames(fdat)



joined <- as.data.frame(joined)

#cell type label ======================

celltype <- read.csv(file = "Karaayvaz cell type labels.csv")
celltype[1:5,1:2]
rownames(celltype) <- celltype[,1]
dim(celltype)
celltype <- celltype[1:1112,"V2", drop = FALSE]
celltype[1:5,1:1]


#sobj object =====================================
sobj_pre <- CreateSeuratObject(counts = joined, project = "Karaayvaz")
sobj_pre <-AddMetaData(sobj_pre,metadata=pdat_joined)
sobj_pre <- AddMetaData(sobj_pre, metadata = celltype, col.name = "celltype_minor")
head(sobj_pre@meta.data)
#gene name input
sobj_pre[["RNA"]]@meta.features<-fdat
head(sobj_pre[["RNA"]]@meta.features)

sobj_pre$orig.ident <- "Karaayvaz"

sobj <- sobj_pre

saveRDS(sobj, file = "Karaayvaz_unfiltered_62622.rds")


#cell cycle ========================================================================
#from: https://satijalab.org/seurat/articles/cell_cycle_vignette.html

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes)
sobj$CC.Difference <- sobj$S.Score - sobj$G2M.Score
# view cell cycle scores and phase assignments
head(sobj[[]])

# Preprocessing and Filtering (SCTransform) --------------------------------------------------------------------

mito.genes <- grep(pattern = "^MT-", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(mito.genes)
sobj <- PercentageFeatureSet(sobj, pattern = "^MT", col.name = "percent.mt")

#ribosomal genes
ribo.genes <- grep(pattern = "^RP[SL]", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(ribo.genes)
sobj <- PercentageFeatureSet(sobj, "^RP[SL]", col.name = "percent.ribo")

#blood cell contamination
hemo.genes <- grep(pattern = "^HB[^(P)]", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(hemo.genes)
platelet.genes <- grep(pattern = "^PECAM1|^PF4\\b", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(platelet.genes)
sobj <- PercentageFeatureSet(sobj, "^HB[^(P)]", col.name = "percent.hb")
sobj <- PercentageFeatureSet(sobj, "^PECAM1|^PF4\\b", col.name = "percent.platelet")

nCount95 <- quantile(sobj@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
nCount5 <- quantile(sobj@meta.data$nCount_RNA, 0.05)#05)  
nFeat95 <- quantile(sobj@meta.data$nFeature_RNA, 0.95) #0.9) 
nFeat5 <- quantile(sobj@meta.data$nFeature_RNA, 0.05) #0.1) 
ptMt5 <- quantile(sobj@meta.data$percent.mt, 0.05)
ptMt90 <- quantile(sobj@meta.data$percent.mt, 0.95)
ptHb95 <- quantile(sobj@meta.data$percent.hb, 0.95)

sobj <- subset(x = sobj, subset = nFeature_RNA > nFeat5 & nFeature_RNA < nFeat95)
sobj <- subset(x = sobj, subset = percent.mt > ptMt5 & percent.mt < ptMt90)
sobj <- subset(x = sobj, subset = nCount_RNA > nCount5 & nCount_RNA < nCount95)
sobj <- subset(x = sobj, subset = percent.hb < ptHb95)

summary(sobj@meta.data$percent.mt)
summary(sobj@meta.data$percent.hb)
dim(sobj)
dim(sobj_pre)


colnames(sobj@meta.data)
sobj$BMI <- "NA"
sobj$Histology <- "NA"


sobj$Ki67 <- sobj$Patient
sobj$Ki67[which(sobj$Ki67 == "PT039")] <- "35% (IHC)"
sobj$Ki67[which(sobj$Ki67 == "PT058")] <- "<10% (IHC)"
sobj$Ki67[which(sobj$Ki67 == "PT081")] <- ">50% (IHC)"
sobj$Ki67[which(sobj$Ki67 == "PT084")] <- "50% (IHC)"
sobj$Ki67[which(sobj$Ki67 == "PT089")] <- "15% (IHC)"
sobj$Ki67[which(sobj$Ki67 == "PT126")] <- "65% (IHC)"

saveRDS(sobj, file = "Karaayvaz_filtered_62622.rds")


# ================================================================== ======

# Pal ============================================================== =======================
# ================================================================== ======

# Loading Raw Data into RStudio ---------------------------------- 

filePaths = getGEOSuppFiles("GSE161529") 
tarF <- list.files(path = "./GSE161529/", pattern = "*.tar", full.names = TRUE) 
ldply(.data = tarF, .fun = untar, exdir = "./GSE161529/")
gzipF <- list.files(path = "./GSE161529/", pattern = "*.gz", full.names = TRUE) 
ldply(.data = gzipF, .fun = gunzip) 

list.files(path = "./GSE161529/", full.names = TRUE)
tarF <- list.files(path = "./GSE161529/", pattern = "*.tar", full.names = TRUE) 
ldply(.data = tarF, .fun = untar, exdir = "./GSE161529/")


list.files(path = "./GSE161529/", pattern = "\\.mtx$",full.names = TRUE)
list.files(path = "./GSE161529/", pattern = "*.features.tsv$", full.names = TRUE) 
list.files(path = "./GSE161529/", pattern = "*.barcodes.tsv$", full.names = TRUE) 



#TN-0126-Tot (TNBC total, NA, NA) ----------------------------------------

N0126_TNBCTot.mat <- readMM(file = './GSE161529//GSM4909281_TN-MH0126-matrix.mtx')
N0126_TNBCTot.mat <- as.matrix(N0126_TNBCTot.mat)
N0126_TNBCTot.mat[1:5,1:5]
dim(N0126_TNBCTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0126_TNBCTot.mat) <- genes[,2]
N0126_TNBCTot.mat[1:5,1:5]

N0126_TNBCTot.bar <- read.table(file = './GSE161529//GSM4909281_TN-MH0126-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0126_TNBCTot.mat) <- N0126_TNBCTot.bar[,1]
N0126_TNBCTot.mat[1:5,1:5]

colnames(N0126_TNBCTot.mat) <- paste(colnames(N0126_TNBCTot.mat), "0126_TNBCTot", sep = "_")
N0126_TNBCTot.pdat <- data.frame("cells" = colnames(N0126_TNBCTot.mat),"samples" = "0126_TNBC", "Patient" = "0126","Molecular Subtype" = "TNBC",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "TNBC", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                                 "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "64",
                                 "Histology" = "NA", "Tumor.Size" = "6.4", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "NA", "Age.At.Menopause" = "NA")




#TN-0135-Tot (TNBC total, NA, NA) ----------------------------------------

N0135_TNBCTot.mat <- readMM(file = './GSE161529//GSM4909282_TN-MH0135-matrix.mtx')
N0135_TNBCTot.mat <- as.matrix(N0135_TNBCTot.mat)
N0135_TNBCTot.mat[1:5,1:5]
dim(N0135_TNBCTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0135_TNBCTot.mat) <- genes[,2]
N0135_TNBCTot.mat[1:5,1:5]

N0135_TNBCTot.bar <- read.table(file = './GSE161529//GSM4909282_TN-MH0135-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0135_TNBCTot.mat) <- N0135_TNBCTot.bar[,1]
N0135_TNBCTot.mat[1:5,1:5]

colnames(N0135_TNBCTot.mat) <- paste(colnames(N0135_TNBCTot.mat), "0135_TNBCTot", sep = "_")
N0135_TNBCTot.pdat <- data.frame("cells" = colnames(N0135_TNBCTot.mat),"samples" = "0135_TNBC", "Patient" = "0135","Molecular Subtype" = "TNBC",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "TNBC", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                                 "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "61",
                                 "Histology" = "NA", "Tumor.Size" = "2.2", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "NA", "Age.At.Menopause" = "NA")





#TN-0106-Tot (TNBC total, NA, NA) ----------------------------------------

N0106_TNBCTot.mat <- readMM(file = './GSE161529//GSM4909283_TN-SH0106-matrix.mtx')
N0106_TNBCTot.mat <- as.matrix(N0106_TNBCTot.mat)
N0106_TNBCTot.mat[1:5,1:5]
dim(N0106_TNBCTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0106_TNBCTot.mat) <- genes[,2]
N0106_TNBCTot.mat[1:5,1:5]

N0106_TNBCTot.bar <- read.table(file = './GSE161529//GSM4909283_TN-SH0106-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0106_TNBCTot.mat) <- N0106_TNBCTot.bar[,1]
N0106_TNBCTot.mat[1:5,1:5]

colnames(N0106_TNBCTot.mat) <- paste(colnames(N0106_TNBCTot.mat), "0106_TNBCTot", sep = "_")
N0106_TNBCTot.pdat <- data.frame("cells" = colnames(N0106_TNBCTot.mat),"samples" = "0106_TNBC", "Patient" = "0106","Molecular Subtype" = "TNBC",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "TNBC", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA",  
                                 "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "65",
                                 "Histology" = "NA", "Tumor.Size" = "2.5", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "NA", "Age.At.Menopause" = "NA")






#TN-0114-Tot (TNBC total, NA, NA) ----------------------------------------

N0114_TNBCTot.mat <- readMM(file = './GSE161529//GSM4909284_TN-MH0114-T2-matrix.mtx')
N0114_TNBCTot.mat <- as.matrix(N0114_TNBCTot.mat)
N0114_TNBCTot.mat[1:5,1:5]
dim(N0114_TNBCTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0114_TNBCTot.mat) <- genes[,2]
N0114_TNBCTot.mat[1:5,1:5]

N0114_TNBCTot.bar <- read.table(file = './GSE161529//GSM4909284_TN-MH0114-T2-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0114_TNBCTot.mat) <- N0114_TNBCTot.bar[,1]
N0114_TNBCTot.mat[1:5,1:5]

colnames(N0114_TNBCTot.mat) <- paste(colnames(N0114_TNBCTot.mat), "0114_TNBCTot", sep = "_")
N0114_TNBCTot.pdat <- data.frame("cells" = colnames(N0114_TNBCTot.mat),"samples" = "0114_TNBC", "Patient" = "0114","Molecular Subtype" = "TNBC",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "TNBC", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                                 "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "84",
                                 "Histology" = "NA", "Tumor.Size" = "1.7", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "NA", "Age.At.Menopause" = "NA")




#HER2-0308-Tot (HER2+ total, NA, NA) ----------------------------------------

N0308_HER2Tot.mat <- readMM(file = './GSE161529//GSM4909289_HER2-AH0308-matrix.mtx')
N0308_HER2Tot.mat <- as.matrix(N0308_HER2Tot.mat)
N0308_HER2Tot.mat[1:5,1:5]
dim(N0308_HER2Tot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0308_HER2Tot.mat) <- genes[,2]
N0308_HER2Tot.mat[1:5,1:5]

N0308_HER2Tot.bar <- read.table(file = './GSE161529//GSM4909289_HER2-AH0308-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0308_HER2Tot.mat) <- N0308_HER2Tot.bar[,1]
N0308_HER2Tot.mat[1:5,1:5]

colnames(N0308_HER2Tot.mat) <- paste(colnames(N0308_HER2Tot.mat), "0308_Her2Tot", sep = "_")
N0308_HER2Tot.pdat <- data.frame("cells" = colnames(N0308_HER2Tot.mat),"samples" = "0308_Her2", "Patient" = "0308","Molecular Subtype" = "HER2+",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HER2+", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                                 "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "32",
                                 "Histology" = "NA", "Tumor.Size" = "2.0", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "NA", "Age.At.Menopause" = "NA")



#HER2-0337-Tot (HER2+ total, NA, NA) ----------------------------------------

N0337_HER2Tot.mat <- readMM(file = './GSE161529//GSM4909290_HER2-PM0337-matrix.mtx')
N0337_HER2Tot.mat <- as.matrix(N0337_HER2Tot.mat)
N0337_HER2Tot.mat[1:5,1:5]
dim(N0337_HER2Tot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0337_HER2Tot.mat) <- genes[,2]
N0337_HER2Tot.mat[1:5,1:5]

N0337_HER2Tot.bar <- read.table(file = './GSE161529//GSM4909290_HER2-PM0337-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0337_HER2Tot.mat) <- N0337_HER2Tot.bar[,1]
N0337_HER2Tot.mat[1:5,1:5]

colnames(N0337_HER2Tot.mat) <- paste(colnames(N0337_HER2Tot.mat), "0337_Her2Tot", sep = "_")
N0337_HER2Tot.pdat <- data.frame("cells" = colnames(N0337_HER2Tot.mat),"samples" = "0337_Her2", "Patient" = "0337","Molecular Subtype" = "HER2+",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HER2+", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                                 "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "66",
                                 "Histology" = "NA", "Tumor.Size" = "6.7", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "NA", "Age.At.Menopause" = "NA")




#HER2-0031-Tot (HER2+ total, NA, NA) ----------------------------------------

N0031_HER2Tot.mat <- readMM(file = './GSE161529//GSM4909291_HER2-MH0031-matrix.mtx')
N0031_HER2Tot.mat <- as.matrix(N0031_HER2Tot.mat)
N0031_HER2Tot.mat[1:5,1:5]
dim(N0031_HER2Tot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0031_HER2Tot.mat) <- genes[,2]
N0031_HER2Tot.mat[1:5,1:5]

N0031_HER2Tot.bar <- read.table(file = './GSE161529//GSM4909291_HER2-MH0031-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0031_HER2Tot.mat) <- N0031_HER2Tot.bar[,1]
N0031_HER2Tot.mat[1:5,1:5]

colnames(N0031_HER2Tot.mat) <- paste(colnames(N0031_HER2Tot.mat), "0031_Her2Tot", sep = "_")
N0031_HER2Tot.pdat <- data.frame("cells" = colnames(N0031_HER2Tot.mat),"samples" = "0031_Her2", "Patient" = "0031","Molecular Subtype" = "HER2+",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HER2+", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA",  
                                 "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "47",
                                 "Histology" = "NA", "Tumor.Size" = "1.8", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "NA", "Age.At.Menopause" = "NA")





#HER2-0069-Tot (HER2+ total, NA, NA) ----------------------------------------

N0069_HER2Tot.mat <- readMM(file = './GSE161529//GSM4909292_HER2-MH0069-matrix.mtx')
N0069_HER2Tot.mat <- as.matrix(N0069_HER2Tot.mat)
N0069_HER2Tot.mat[1:5,1:5]
dim(N0069_HER2Tot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0069_HER2Tot.mat) <- genes[,2]
N0069_HER2Tot.mat[1:5,1:5]

N0069_HER2Tot.bar <- read.table(file = './GSE161529//GSM4909292_HER2-MH0069-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0069_HER2Tot.mat) <- N0069_HER2Tot.bar[,1]
N0069_HER2Tot.mat[1:5,1:5]

colnames(N0069_HER2Tot.mat) <- paste(colnames(N0069_HER2Tot.mat), "0069_Her2Tot", sep = "_")
N0069_HER2Tot.pdat <- data.frame("cells" = colnames(N0069_HER2Tot.mat),"samples" = "0069_Her2", "Patient" = "0069","Molecular Subtype" = "HER2+",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HER2+", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                                 "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "71",
                                 "Histology" = "NA", "Tumor.Size" = "2.7", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "NA", "Age.At.Menopause" = "NA")






#HER2-0161-Tot (HER2+ total, NA, NA) ----------------------------------------

N0161_HER2Tot.mat <- readMM(file = './GSE161529//GSM4909293_HER2-MH0161-matrix.mtx')
N0161_HER2Tot.mat <- as.matrix(N0161_HER2Tot.mat)
N0161_HER2Tot.mat[1:5,1:5]
dim(N0161_HER2Tot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0161_HER2Tot.mat) <- genes[,2]
N0161_HER2Tot.mat[1:5,1:5]

N0161_HER2Tot.bar <- read.table(file = './GSE161529//GSM4909293_HER2-MH0161-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0161_HER2Tot.mat) <- N0161_HER2Tot.bar[,1]
N0161_HER2Tot.mat[1:5,1:5]

colnames(N0161_HER2Tot.mat) <- paste(colnames(N0161_HER2Tot.mat), "0161_Her2Tot", sep = "_")
N0161_HER2Tot.pdat <- data.frame("cells" = colnames(N0161_HER2Tot.mat),"samples" = "0161_Her2", "Patient" = "0161","Molecular Subtype" = "HER2+",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HER2+", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                                 "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "80",
                                 "Histology" = "NA", "Tumor.Size" = "4.5", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "NA", "Age.At.Menopause" = "NA")




#HER2-0176-Tot (HER2+ total, NA, NA) ----------------------------------------

N0176_HER2Tot.mat <- readMM(file = './GSE161529//GSM4909293_HER2-MH0161-matrix.mtx')
N0176_HER2Tot.mat <- as.matrix(N0176_HER2Tot.mat)
N0176_HER2Tot.mat[1:5,1:5]
dim(N0176_HER2Tot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0176_HER2Tot.mat) <- genes[,2]
N0176_HER2Tot.mat[1:5,1:5]

N0176_HER2Tot.bar <- read.table(file = './GSE161529//GSM4909293_HER2-MH0161-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0176_HER2Tot.mat) <- N0176_HER2Tot.bar[,1]
N0176_HER2Tot.mat[1:5,1:5]

colnames(N0176_HER2Tot.mat) <- paste(colnames(N0176_HER2Tot.mat), "0176_Her2Tot", sep = "_")
N0176_HER2Tot.pdat <- data.frame("cells" = colnames(N0176_HER2Tot.mat),"samples" = "0176_Her2", "Patient" = "0176","Molecular Subtype" = "HER2+",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HER2+", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                                 "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "60",
                                 "Histology" = "NA", "Tumor.Size" = "2.0", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "NA", "Age.At.Menopause" = "NA")





#PR-0319-Tot (HR+ total, NA, NA) ----------------------------------------

N0319_HRTot.mat <- readMM(file = './GSE161529//GSM4909295_ER-AH0319-matrix.mtx')
N0319_HRTot.mat <- as.matrix(N0319_HRTot.mat)
N0319_HRTot.mat[1:5,1:5]
dim(N0319_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0319_HRTot.mat) <- genes[,2]
N0319_HRTot.mat[1:5,1:5]

N0319_HRTot.bar <- read.table(file = './GSE161529//GSM4909295_ER-AH0319-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0319_HRTot.mat) <- N0319_HRTot.bar[,1]
N0319_HRTot.mat[1:5,1:5]

colnames(N0319_HRTot.mat) <- paste(colnames(N0319_HRTot.mat), "0319_HRTot", sep = "_")
N0319_HRTot.pdat <- data.frame("cells" = colnames(N0319_HRTot.mat),"samples" = "0319_HR", "Patient" = "0319","Molecular Subtype" = "PR+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "58",
                               "Histology" = "NA", "Tumor.Size" = "2.7", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")






#ER-0001-Tot (HR+ total, NA, NA) ----------------------------------------

N0001_HRTot.mat <- readMM(file = './GSE161529//GSM4909296_ER-MH0001-matrix.mtx')
N0001_HRTot.mat <- as.matrix(N0001_HRTot.mat)
N0001_HRTot.mat[1:5,1:5]
dim(N0001_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0001_HRTot.mat) <- genes[,2]
N0001_HRTot.mat[1:5,1:5]

N0001_HRTot.bar <- read.table(file = './GSE161529//GSM4909296_ER-MH0001-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0001_HRTot.mat) <- N0001_HRTot.bar[,1]
N0001_HRTot.mat[1:5,1:5]

colnames(N0001_HRTot.mat) <- paste(colnames(N0001_HRTot.mat), "0001_HRTot", sep = "_")
N0001_HRTot.pdat <- data.frame("cells" = colnames(N0001_HRTot.mat),"samples" = "0001_HR", "Patient" = "0001","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "58",
                               "Histology" = "NA", "Tumor.Size" = "3.2", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")




#ER-0125-Tot (HR+ total, NA, NA) ----------------------------------------

N0125_HRTot.mat <- readMM(file = './GSE161529//GSM4909297_ER-MH0125-matrix.mtx')
N0125_HRTot.mat <- as.matrix(N0125_HRTot.mat)
N0125_HRTot.mat[1:5,1:5]
dim(N0125_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0125_HRTot.mat) <- genes[,2]
N0125_HRTot.mat[1:5,1:5]

N0125_HRTot.bar <- read.table(file = './GSE161529//GSM4909297_ER-MH0125-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0125_HRTot.mat) <- N0125_HRTot.bar[,1]
N0125_HRTot.mat[1:5,1:5]

colnames(N0125_HRTot.mat) <- paste(colnames(N0125_HRTot.mat), "0125_HRTot", sep = "_")
N0125_HRTot.pdat <- data.frame("cells" = colnames(N0125_HRTot.mat),"samples" = "0125_HR", "Patient" = "0125","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "2","BRCA.Status" = "Negative", "Age" = "45",
                               "Histology" = "NA", "Tumor.Size" = "4.8", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")





#ER-0360-Tot (HR+ total, NA, NA) ----------------------------------------

N0360_HRTot.mat <- readMM(file = './GSE161529//GSM4909298_ER-PM0360-matrix.mtx')
N0360_HRTot.mat <- as.matrix(N0360_HRTot.mat)
N0360_HRTot.mat[1:5,1:5]
dim(N0360_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0360_HRTot.mat) <- genes[,2]
N0360_HRTot.mat[1:5,1:5]

N0360_HRTot.bar <- read.table(file = './GSE161529//GSM4909298_ER-PM0360-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0360_HRTot.mat) <- N0360_HRTot.bar[,1]
N0360_HRTot.mat[1:5,1:5]

colnames(N0360_HRTot.mat) <- paste(colnames(N0360_HRTot.mat), "0360_HRTot", sep = "_")
N0360_HRTot.pdat <- data.frame("cells" = colnames(N0360_HRTot.mat),"samples" = "0360_HR", "Patient" = "0360","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "2","BRCA.Status" = "Negative", "Age" = "70",
                               "Histology" = "NA", "Tumor.Size" = "5.0", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Lobular Carcinoma", "Age.At.Menopause" = "NA")






#ER-0114_T3-Tot (HR+ total, NA, NA) ----------------------------------------

N0114T3_HRTot.mat <- readMM(file = './GSE161529//GSM4909299_ER-MH0114-T3-matrix.mtx')
N0114T3_HRTot.mat <- as.matrix(N0114T3_HRTot.mat)
N0114T3_HRTot.mat[1:5,1:5]
dim(N0114T3_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0114T3_HRTot.mat) <- genes[,2]
N0114T3_HRTot.mat[1:5,1:5]

N0114T3_HRTot.bar <- read.table(file = './GSE161529//GSM4909299_ER-MH0114-T3-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0114T3_HRTot.mat) <- N0114T3_HRTot.bar[,1]
N0114T3_HRTot.mat[1:5,1:5]

colnames(N0114T3_HRTot.mat) <- paste(colnames(N0114T3_HRTot.mat), "0114_T3_HRTot", sep = "_")
N0114T3_HRTot.pdat <- data.frame("cells" = colnames(N0114T3_HRTot.mat),"samples" = "0114_T3_HR", "Patient" = "0114_T3","Molecular Subtype" = "ER+",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                                 "Stage" = "NA", "Grade" = "1","BRCA.Status" = "Negative", "Age" = "84",
                                 "Histology" = "NA", "Tumor.Size" = "2.0", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")







#ER-0032-Tot (HR+ total, NA, NA) ----------------------------------------

N0032_HRTot.mat <- readMM(file = './GSE161529//GSM4909300_ER-MH0032-matrix.mtx')
N0032_HRTot.mat <- as.matrix(N0032_HRTot.mat)
N0032_HRTot.mat[1:5,1:5]
dim(N0032_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0032_HRTot.mat) <- genes[,2]
N0032_HRTot.mat[1:5,1:5]

N0032_HRTot.bar <- read.table(file = './GSE161529//GSM4909300_ER-MH0032-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0032_HRTot.mat) <- N0032_HRTot.bar[,1]
N0032_HRTot.mat[1:5,1:5]

colnames(N0032_HRTot.mat) <- paste(colnames(N0032_HRTot.mat), "0032_HRTot", sep = "_")
N0032_HRTot.pdat <- data.frame("cells" = colnames(N0032_HRTot.mat),"samples" = "0032_HR", "Patient" = "0032","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "55",
                               "Histology" = "NA", "Tumor.Size" = "9.0", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")





#ER-0042-Tot (HR+ total, NA, NA) ----------------------------------------

N0042_HRTot.mat <- readMM(file = './GSE161529//GSM4909301_ER-MH0042-matrix.mtx')
N0042_HRTot.mat <- as.matrix(N0042_HRTot.mat)
N0042_HRTot.mat[1:5,1:5]
dim(N0042_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0042_HRTot.mat) <- genes[,2]
N0042_HRTot.mat[1:5,1:5]

N0042_HRTot.bar <- read.table(file = './GSE161529//GSM4909301_ER-MH0042-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0042_HRTot.mat) <- N0042_HRTot.bar[,1]
N0042_HRTot.mat[1:5,1:5]

colnames(N0042_HRTot.mat) <- paste(colnames(N0042_HRTot.mat), "0042_HRTot", sep = "_")
N0042_HRTot.pdat <- data.frame("cells" = colnames(N0042_HRTot.mat),"samples" = "0042_HR", "Patient" = "0042","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "2","BRCA.Status" = "Negative", "Age" = "58",
                               "Histology" = "NA", "Tumor.Size" = "1.8", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")




#ER-0025-Tot (HR+ total, NA, NA) ----------------------------------------

N0025_HRTot.mat <- readMM(file = './GSE161529//GSM4909302_ER-MH0025-matrix.mtx')
N0025_HRTot.mat <- as.matrix(N0025_HRTot.mat)
N0025_HRTot.mat[1:5,1:5]
dim(N0025_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0025_HRTot.mat) <- genes[,2]
N0025_HRTot.mat[1:5,1:5]

N0025_HRTot.bar <- read.table(file = './GSE161529//GSM4909302_ER-MH0025-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0025_HRTot.mat) <- N0025_HRTot.bar[,1]
N0025_HRTot.mat[1:5,1:5]

colnames(N0025_HRTot.mat) <- paste(colnames(N0025_HRTot.mat), "0025_HRTot", sep = "_")
N0025_HRTot.pdat <- data.frame("cells" = colnames(N0025_HRTot.mat),"samples" = "0025_HR", "Patient" = "0025","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "2","BRCA.Status" = "Negative", "Age" = "52",
                               "Histology" = "NA", "Tumor.Size" = "2.3", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")



#ER-0151-Tot (HR+ total, NA, NA) ----------------------------------------

N0151_HRTot.mat <- readMM(file = './GSE161529//GSM4909303_ER-MH0151-matrix.mtx')
N0151_HRTot.mat <- as.matrix(N0151_HRTot.mat)
N0151_HRTot.mat[1:5,1:5]
dim(N0151_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0151_HRTot.mat) <- genes[,2]
N0151_HRTot.mat[1:5,1:5]

N0151_HRTot.bar <- read.table(file = './GSE161529//GSM4909303_ER-MH0151-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0151_HRTot.mat) <- N0151_HRTot.bar[,1]
N0151_HRTot.mat[1:5,1:5]

colnames(N0151_HRTot.mat) <- paste(colnames(N0151_HRTot.mat), "0151_HRTot", sep = "_")
N0151_HRTot.pdat <- data.frame("cells" = colnames(N0151_HRTot.mat),"samples" = "0151_HR", "Patient" = "0151","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "2","BRCA.Status" = "Negative", "Age" = "49",
                               "Histology" = "NA", "Tumor.Size" = "3.5", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")


#ER-0163-Tot (HR+ total, NA, NA) ----------------------------------------

N0163_HRTot.mat <- readMM(file = './GSE161529//GSM4909304_ER-MH0163-matrix.mtx')
N0163_HRTot.mat <- as.matrix(N0163_HRTot.mat)
N0163_HRTot.mat[1:5,1:5]
dim(N0163_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0163_HRTot.mat) <- genes[,2]
N0163_HRTot.mat[1:5,1:5]

N0163_HRTot.bar <- read.table(file = './GSE161529//GSM4909304_ER-MH0163-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0163_HRTot.mat) <- N0163_HRTot.bar[,1]
N0163_HRTot.mat[1:5,1:5]

colnames(N0163_HRTot.mat) <- paste(colnames(N0163_HRTot.mat), "0163_HRTot", sep = "_")
N0163_HRTot.pdat <- data.frame("cells" = colnames(N0163_HRTot.mat),"samples" = "0163_HR", "Patient" = "0163","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "45",
                               "Histology" = "NA", "Tumor.Size" = "4.5", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")



#ER-00297C-Tot (HR+ total, NA, NA) ----------------------------------------

N00297C_HRTot.mat <- readMM(file = './GSE161529//GSM4909305_ER-MH0029-7C-matrix.mtx')
N00297C_HRTot.mat <- as.matrix(N00297C_HRTot.mat)
N00297C_HRTot.mat[1:5,1:5]
dim(N00297C_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N00297C_HRTot.mat) <- genes[,2]
N00297C_HRTot.mat[1:5,1:5]

N00297C_HRTot.bar <- read.table(file = './GSE161529//GSM4909305_ER-MH0029-7C-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N00297C_HRTot.mat) <- N00297C_HRTot.bar[,1]
N00297C_HRTot.mat[1:5,1:5]

colnames(N00297C_HRTot.mat) <- paste(colnames(N00297C_HRTot.mat), "0029_7C_HRTot", sep = "_")
N00297C_HRTot.pdat <- data.frame("cells" = colnames(N00297C_HRTot.mat),"samples" = "0029_7C_HR", "Patient" = "0029","Molecular Subtype" = "ER+",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                                 "Stage" = "NA", "Grade" = "1","BRCA.Status" = "Negative", "Age" = "59",
                                 "Histology" = "NA", "Tumor.Size" = "0.8", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")




#ER-00299C-Tot (HR+ total, NA, NA) ----------------------------------------

N00299C_HRTot.mat <- readMM(file = './GSE161529//GSM4909306_ER-MH0029-9C-matrix.mtx')
N00299C_HRTot.mat <- as.matrix(N00299C_HRTot.mat)
N00299C_HRTot.mat[1:5,1:5]
dim(N00299C_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N00299C_HRTot.mat) <- genes[,2]
N00299C_HRTot.mat[1:5,1:5]

N00299C_HRTot.bar <- read.table(file = './GSE161529//GSM4909306_ER-MH0029-9C-barcodes.tsv', 
                                sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N00299C_HRTot.mat) <- N00299C_HRTot.bar[,1]
N00299C_HRTot.mat[1:5,1:5]

colnames(N00299C_HRTot.mat) <- paste(colnames(N00299C_HRTot.mat), "0029_9C_HRTot", sep = "_")
N00299C_HRTot.pdat <- data.frame("cells" = colnames(N00299C_HRTot.mat),"samples" = "0029_9C_HR", "Patient" = "0029","Molecular Subtype" = "ER+",
                                 "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                                 "Capture Method" = "10X Genomics Chromium",
                                 "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                                 "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                                 "Stage" = "NA", "Grade" = "1","BRCA.Status" = "Negative", "Age" = "59",
                                 "Histology" = "NA", "Tumor.Size" = "3.8", "Treatment.Status" = "Naive",
                                 "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")





#ER-0040-Tot (HR+ total, NA, NA) ----------------------------------------

N0040_HRTot.mat <- readMM(file = './GSE161529//GSM4909307_ER-MH0040-matrix.mtx')
N0040_HRTot.mat <- as.matrix(N0040_HRTot.mat)
N0040_HRTot.mat[1:5,1:5]
dim(N0040_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0040_HRTot.mat) <- genes[,2]
N0040_HRTot.mat[1:5,1:5]

N0040_HRTot.bar <- read.table(file = './GSE161529//GSM4909307_ER-MH0040-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0040_HRTot.mat) <- N0040_HRTot.bar[,1]
N0040_HRTot.mat[1:5,1:5]

colnames(N0040_HRTot.mat) <- paste(colnames(N0040_HRTot.mat), "0040_HRTot", sep = "_")
N0040_HRTot.pdat <- data.frame("cells" = colnames(N0040_HRTot.mat),"samples" = "0040_HR", "Patient" = "0040","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "69",
                               "Histology" = "NA", "Tumor.Size" = "5.0", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")



#ER-0043-Tot (HR+ total, NA, NA) ----------------------------------------

N0043_HRTot.mat <- readMM(file = './GSE161529//GSM4909309_ER-MH0043-T-matrix.mtx')
N0043_HRTot.mat <- as.matrix(N0043_HRTot.mat)
N0043_HRTot.mat[1:5,1:5]
dim(N0043_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0043_HRTot.mat) <- genes[,2]
N0043_HRTot.mat[1:5,1:5]

N0043_HRTot.bar <- read.table(file = './GSE161529//GSM4909309_ER-MH0043-T-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0043_HRTot.mat) <- N0043_HRTot.bar[,1]
N0043_HRTot.mat[1:5,1:5]

colnames(N0043_HRTot.mat) <- paste(colnames(N0043_HRTot.mat), "0043_HRTot", sep = "_")
N0043_HRTot.pdat <- data.frame("cells" = colnames(N0043_HRTot.mat),"samples" = "0043_HR", "Patient" = "0043","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "55",
                               "Histology" = "NA", "Tumor.Size" = "5.5", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")




#ER-0056-Tot (HR+ total, NA, NA) ----------------------------------------

N0056_HRTot.mat <- readMM(file = './GSE161529//GSM4909311_ER-MH0056-T-matrix.mtx')
N0056_HRTot.mat <- as.matrix(N0056_HRTot.mat)
N0056_HRTot.mat[1:5,1:5]
dim(N0056_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0056_HRTot.mat) <- genes[,2]
N0056_HRTot.mat[1:5,1:5]

N0056_HRTot.bar <- read.table(file = './GSE161529//GSM4909311_ER-MH0056-T-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0056_HRTot.mat) <- N0056_HRTot.bar[,1]
N0056_HRTot.mat[1:5,1:5]

colnames(N0056_HRTot.mat) <- paste(colnames(N0056_HRTot.mat), "0056_HRTot", sep = "_")
N0056_HRTot.pdat <- data.frame("cells" = colnames(N0056_HRTot.mat),"samples" = "0056_HR", "Patient" = "0056","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "2","BRCA.Status" = "Negative", "Age" = "66",
                               "Histology" = "NA", "Tumor.Size" = "4.5", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")



#ER-0064-Tot (HR+ total, NA, NA) ----------------------------------------

N0064_HRTot.mat <- readMM(file = './GSE161529//GSM4909313_ER-MH0064-T-matrix.mtx')
N0064_HRTot.mat <- as.matrix(N0064_HRTot.mat)
N0064_HRTot.mat[1:5,1:5]
dim(N0064_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0064_HRTot.mat) <- genes[,2]
N0064_HRTot.mat[1:5,1:5]

N0064_HRTot.bar <- read.table(file = './GSE161529//GSM4909313_ER-MH0064-T-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0064_HRTot.mat) <- N0064_HRTot.bar[,1]
N0064_HRTot.mat[1:5,1:5]

colnames(N0064_HRTot.mat) <- paste(colnames(N0064_HRTot.mat), "0064_HRTot", sep = "_")
N0064_HRTot.pdat <- data.frame("cells" = colnames(N0064_HRTot.mat),"samples" = "0064_HR", "Patient" = "0064_HR","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "2","BRCA.Status" = "Negative", "Age" = "65",
                               "Histology" = "NA", "Tumor.Size" = "7.0", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")




#ER-0167-Tot (HR+ total, NA, NA) ----------------------------------------

N0167_HRTot.mat <- readMM(file = './GSE161529//GSM4909315_ER-MH0167-T-matrix.mtx')
N0167_HRTot.mat <- as.matrix(N0167_HRTot.mat)
N0167_HRTot.mat[1:5,1:5]
dim(N0167_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0167_HRTot.mat) <- genes[,2]
N0167_HRTot.mat[1:5,1:5]

N0167_HRTot.bar <- read.table(file = './GSE161529//GSM4909315_ER-MH0167-T-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0167_HRTot.mat) <- N0167_HRTot.bar[,1]
N0167_HRTot.mat[1:5,1:5]

colnames(N0167_HRTot.mat) <- paste(colnames(N0167_HRTot.mat), "0167_HRTot", sep = "_")
N0167_HRTot.pdat <- data.frame("cells" = colnames(N0167_HRTot.mat),"samples" = "0167_HR", "Patient" = "0167","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "83",
                               "Histology" = "NA", "Tumor.Size" = "3.4", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")





#ER-0173-Tot (HR+ total, NA, NA) ----------------------------------------

N0173_HRTot.mat <- readMM(file = './GSE161529//GSM4909317_ER-MH0173-T-matrix.mtx')
N0173_HRTot.mat <- as.matrix(N0173_HRTot.mat)
N0173_HRTot.mat[1:5,1:5]
dim(N0173_HRTot.mat)

genes <- read.table(file = './GSE161529//GSE161529_features.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(genes)
rownames(N0173_HRTot.mat) <- genes[,2]
N0173_HRTot.mat[1:5,1:5]

N0173_HRTot.bar <- read.table(file = './GSE161529//GSM4909317_ER-MH0173-T-barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(N0173_HRTot.mat) <- N0173_HRTot.bar[,1]
N0173_HRTot.mat[1:5,1:5]

colnames(N0173_HRTot.mat) <- paste(colnames(N0173_HRTot.mat), "0173_HRTot", sep = "_")
N0173_HRTot.pdat <- data.frame("cells" = colnames(N0173_HRTot.mat),"samples" = "0173_HR", "Patient" = "0173","Molecular Subtype" = "ER+",
                               "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium Single Cell 3'", 
                               "Capture Method" = "10X Genomics Chromium",
                               "Sequencer" = "Illumina NextSeq 500", "Tissue Source" = "Primary Tumor", "Population" = "Total","BC Subtype" = "HR+", "Status" = "Primary",
                               "Gender" = "Female", "Menopause" = "NA", "Parity" = "NA", 
                               "Stage" = "NA", "Grade" = "3","BRCA.Status" = "Negative", "Age" = "83",
                               "Histology" = "NA", "Tumor.Size" = "3.8", "Treatment.Status" = "Naive",
                               "Other.Clinical" = "Invasive Ductal Carcinoma NST", "Age.At.Menopause" = "NA")




#joined =============================

joined_prim <- cbind(N0126_TNBCTot.mat, N0135_TNBCTot.mat,
                     N0106_TNBCTot.mat, N0114_TNBCTot.mat,
                     N0308_HER2Tot.mat, N0337_HER2Tot.mat,
                     N0031_HER2Tot.mat, N0069_HER2Tot.mat,
                     N0161_HER2Tot.mat, N0176_HER2Tot.mat,
                     N0319_HRTot.mat, N0001_HRTot.mat, N0125_HRTot.mat,
                     N0360_HRTot.mat, N0114T3_HRTot.mat, N0032_HRTot.mat,
                     N0042_HRTot.mat, N0025_HRTot.mat, N0151_HRTot.mat,
                     N0163_HRTot.mat, N00297C_HRTot.mat, N00299C_HRTot.mat,
                     N0040_HRTot.mat, N0043_HRTot.mat, N0056_HRTot.mat,
                     N0064_HRTot.mat, N0167_HRTot.mat, N0173_HRTot.mat)

joined_prim <- as.data.frame(joined_prim)


library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(joined_prim)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- joined_prim[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- joined_prim[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}

joined_prim <- new_mat


pdat_prim <- rbind(N0126_TNBCTot.pdat, N0135_TNBCTot.pdat,
                   N0106_TNBCTot.pdat, N0114_TNBCTot.pdat,
                   N0308_HER2Tot.pdat, N0337_HER2Tot.pdat,
                   N0031_HER2Tot.pdat, N0069_HER2Tot.pdat,
                   N0161_HER2Tot.pdat, N0176_HER2Tot.pdat,
                   N0319_HRTot.pdat, N0001_HRTot.pdat, N0125_HRTot.pdat,
                   N0360_HRTot.pdat, N0114T3_HRTot.pdat, N0032_HRTot.pdat,
                   N0042_HRTot.pdat, N0025_HRTot.pdat, N0151_HRTot.pdat,
                   N0163_HRTot.pdat, N00297C_HRTot.pdat, N00299C_HRTot.pdat,
                   N0040_HRTot.pdat, N0043_HRTot.pdat, N0056_HRTot.pdat,
                   N0064_HRTot.pdat, N0167_HRTot.pdat, N0173_HRTot.pdat)

#seuratobj -------------------------------


#_________________________________________________________________
rownames(pdat_prim) <- pdat_prim$cells
fdat_prim <- toupper(as.matrix(rownames(joined_prim)))
#_________________________________________________________________
rownames(fdat_prim) <- fdat_prim[,1]
fdat_prim <- data.frame(fdat_prim)
common_colnames <- "gene_short_name"
colnames(fdat_prim) <- common_colnames
rownames(joined_prim) <- rownames(fdat_prim)
#_________________________________________________________________
sobj_prim <- CreateSeuratObject(counts = joined_prim)
rownames(pdat_prim) <- colnames(joined_prim)
sobj_prim <-AddMetaData(sobj_prim,metadata=pdat_prim)
head(sobj_prim@meta.data)
#_________________________________________________________________
sobj_prim[["RNA"]]@meta.features<-fdat_prim
head(sobj_prim[["RNA"]]@meta.features)
slotNames(sobj_prim[["RNA"]])

sobj_prim$Gravidity <- "NA"


saveRDS(sobj_prim, file = "Pal_unfiltered_62622.rds")

#cell cycle ========================================================================
#from: https://satijalab.org/seurat/articles/cell_cycle_vignette.html

sobj <- sobj_prim
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes)
sobj$CC.Difference <- sobj$S.Score - sobj$G2M.Score
# view cell cycle scores and phase assignments
head(sobj[[]])

# Preprocessing and Filtering (SCTransform) --------------------------------------------------------------------


mito.genes <- grep(pattern = "^MT\\.", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(mito.genes)
sobj <- PercentageFeatureSet(sobj, pattern = "^MT\\.", col.name = "percent.mt")

#blood cell contamination
hemo.genes <- grep(pattern = "^HB[^(P)]", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(hemo.genes)
sobj <- PercentageFeatureSet(sobj, "^HB[^(P)]", col.name = "percent.hb")
platelet.genes <- grep(pattern = "^PECAM1|^PF4\\b", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(platelet.genes)
sobj <- PercentageFeatureSet(sobj, "^PECAM1|^PF4\\b", col.name = "percent.platelet")


nCount95 <- quantile(sobj@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
nCount5 <- quantile(sobj@meta.data$nCount_RNA, 0.05)#05)  
nFeat95 <- quantile(sobj@meta.data$nFeature_RNA, 0.95) #0.9) 
nFeat5 <- quantile(sobj@meta.data$nFeature_RNA, 0.05) #0.1) 
ptMt5 <- quantile(sobj@meta.data$percent.mt, 0.05)
ptMt90 <- quantile(sobj@meta.data$percent.mt, 0.65) #85 for prim and norm
ptHb95 <- quantile(sobj@meta.data$percent.hb, 0.95)

sobj <- subset(x = sobj, subset = nFeature_RNA > nFeat5 & nFeature_RNA < nFeat95)
sobj <- subset(x = sobj, subset = percent.mt > ptMt5 & percent.mt < ptMt90)
sobj <- subset(x = sobj, subset = nCount_RNA > nCount5 & nCount_RNA < nCount95)
sobj <- subset(x = sobj, subset = percent.hb < ptHb95)

summary(sobj@meta.data$percent.mt)
summary(sobj@meta.data$percent.hb)
dim(sobj)
dim(sobj_prim)
dim(sobj_norm)
dim(sobj_met)

saveRDS(sobj, file = "Pal_PRIM_filtered_62622.rds")

# =================================================================== ======

# Qian ============================================================== =======================

# =================================================================== ======

#data loading ========================================


mat <- readMM(file = './export/BC_counts/matrix.mtx')
mat <- as.matrix(mat)
mat[1:5,1:5]
dim(mat)

genes <- read.table(file = './export/BC_counts/genes.tsv', 
                    sep = '\t', header = FALSE, stringsAsFactors = FALSE)
head(genes)
rownames(mat) <- genes[,1]
mat[1:5,1:5]

bar <- read.table(file = './export/BC_counts/barcodes.tsv', 
                  sep = '\t', header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- bar[,1]
mat[1:5,1:5]

mat.dat <- as.data.frame(mat)


library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(mat.dat)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- mat.dat[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- mat.dat[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}


joined <- new_mat

#metadata ----------------------------------------------------


meta <- read.csv(file = './2103-Breastcancer_metadata.csv')
head(meta)

colnames(meta) <- c("cells","nGene", "nUMI", "CellFromTumor","Patient", "TumorType", "TumorSite", "celltype_minor")
head(meta)
meta$RNA.Type <- "Total RNA"
meta$Sequencer <- "Illumina NextSeq, Illumina HiSeq4000, Illumina NovaSeq6000"
meta$Tissue.Source <- "Primary Tumor"
meta$Status <- "Primary"
meta$Gender <- "Female"
meta$Gene.Coverage <- "5'"
meta$Library.Preparation <- "10X Genomics Chromium v2 5'"
meta$Capture.Method <- "10X Genomics Chromium v2 5'"
meta$Treatment.Status <- "Naive"
meta$Parity <- "NA"
meta$Gravidity <- "NA"
meta$Menopause <- "NA"
meta$Ethnicity <- "NA"
meta$BRCA.Status <- "Negative"
meta$Tumor.Size <- "NA"
meta$samples <- meta$Patient



#put patient into each column, then do the "which" thingy for rest of pdata
meta$Patient[which(meta$Patient == "41")] <- "BC_1"
meta$Patient[which(meta$Patient == "42")] <- "BC_2"
meta$Patient[which(meta$Patient == "43")] <- "BC_3"
meta$Patient[which(meta$Patient == "44")] <- "BC_4"
meta$Patient[which(meta$Patient == "45")] <- "BC_5"
meta$Patient[which(meta$Patient == "46")] <- "BC_6"
meta$Patient[which(meta$Patient == "47")] <- "BC_7"
meta$Patient[which(meta$Patient == "48")] <- "BC_8"
meta$Patient[which(meta$Patient == "49")] <- "BC_9"
meta$Patient[which(meta$Patient == "50")] <- "BC_10"
meta$Patient[which(meta$Patient == "51")] <- "BC_11"
meta$Patient[which(meta$Patient == "52")] <- "BC_12"
meta$Patient[which(meta$Patient == "53")] <- "BC_13"
meta$Patient[which(meta$Patient == "54")] <- "BC_14"

meta$Stage <- meta$Patient
meta$Stage[which(meta$Stage == "BC_1")] <- "3"
meta$Stage[which(meta$Stage == "BC_2")] <- "3"
meta$Stage[which(meta$Stage == "BC_3")] <- "3"
meta$Stage[which(meta$Stage == "BC_4")] <- "3"
meta$Stage[which(meta$Stage == "BC_5")] <- "3"
meta$Stage[which(meta$Stage == "BC_6")] <- "3"
meta$Stage[which(meta$Stage == "BC_7")] <- "3"
meta$Stage[which(meta$Stage == "BC_8")] <- "3"
meta$Stage[which(meta$Stage == "BC_9")] <- "3"
meta$Stage[which(meta$Stage == "BC_10")] <- "3"
meta$Stage[which(meta$Stage == "BC_11")] <- "3"
meta$Stage[which(meta$Stage == "BC_12")] <- "3"
meta$Stage[which(meta$Stage == "BC_13")] <- "3"
meta$Stage[which(meta$Stage == "BC_14")] <- "2"


meta$Histology <- meta$Patient
meta$Histology[which(meta$Histology == "BC_1")] <- "Invasive Ductal Carcinoma"
meta$Histology[which(meta$Histology == "BC_2")] <- "Invasive Ductal Carcinoma"
meta$Histology[which(meta$Histology == "BC_3")] <- "Invasive Ductal Carcinoma"
meta$Histology[which(meta$Histology == "BC_4")] <- "Invasive Ductal Carcinoma"
meta$Histology[which(meta$Histology == "BC_5")] <- "Invasive Apocrine Carcinoma"
meta$Histology[which(meta$Histology == "BC_6")] <- "Invasive Ductal Carcinoma"
meta$Histology[which(meta$Histology == "BC_7")] <- "Metaplastic Carcinoma"
meta$Histology[which(meta$Histology == "BC_8")] <- "Invasive Ductal Carcinoma"
meta$Histology[which(meta$Histology == "BC_9")] <- "Invasive Ductal Carcinoma"
meta$Histology[which(meta$Histology == "BC_10")] <- "Invasive Ductal Carcinoma"
meta$Histology[which(meta$Histology == "BC_11")] <- "Invasive Ductal Carcinoma"
meta$Histology[which(meta$Histology == "BC_12")] <- "Invasive Ductal Carcinoma"
meta$Histology[which(meta$Histology == "BC_13")] <- "Invasive Ductal Carcinoma"
meta$Histology[which(meta$Histology == "BC_14")] <- "Invasive Ductal Carcinoma"

meta$Molecular.Subtype <- meta$Patient
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_1")] <- "HER2+"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_2")] <- "TNBC"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_3")] <- "TNBC"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_4")] <- "TNBC"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_5")] <- "TNBC"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_6")] <- "HER2+"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_7")] <- "TNBC"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_8")] <- "HER2+"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_9")] <- "Luminal-HER2+"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_10")] <- "TNBC"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_11")] <- "TNBC"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_12")] <- "TNBC"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_13")] <- "Luminal B"
meta$Molecular.Subtype[which(meta$Molecular.Subtype == "BC_14")] <- "Luminal A"

meta$BC.Subtype <- meta$Patient
meta$BC.Subtype[which(meta$BC.Subtype == "BC_1")] <- "HER2+"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_2")] <- "TNBC"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_3")] <- "TNBC"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_4")] <- "TNBC"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_5")] <- "TNBC"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_6")] <- "HER2+"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_7")] <- "TNBC"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_8")] <- "HER2+"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_9")] <- "HER2+"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_10")] <- "TNBC"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_11")] <- "TNBC"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_12")] <- "TNBC"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_13")] <- "HR+"
meta$BC.Subtype[which(meta$BC.Subtype == "BC_14")] <- "HR+"

meta$TNM.Classification <- meta$Patient
meta$TNM.Classification[which(meta$TNM.Classification == "BC_1")] <- "pT1cN0M0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_2")] <- "pT2N0M0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_3")] <- "pT3N0M0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_4")] <- "pT1N0M0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_5")] <- "pT2N1aM0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_6")] <- "pT2N0M0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_7")] <- "pT3N0M0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_8")] <- "pT3N0M0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_9")] <- "pT2N1M0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_10")] <- "pT2N0M0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_11")] <- "pT1cN0M0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_12")] <- "pT2NxM0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_13")] <- "pT2N1M0"
meta$TNM.Classification[which(meta$TNM.Classification == "BC_14")] <- "pT2N1M0"

meta$Age <- meta$Patient
meta$Age[which(meta$Age == "BC_1")] <- "40-45"
meta$Age[which(meta$Age == "BC_2")] <- "50-55"
meta$Age[which(meta$Age == "BC_3")] <- "50-55"
meta$Age[which(meta$Age == "BC_4")] <- "86-90"
meta$Age[which(meta$Age == "BC_5")] <- "60-65"
meta$Age[which(meta$Age == "BC_6")] <- "70-75"
meta$Age[which(meta$Age == "BC_7")] <- "80-85"
meta$Age[which(meta$Age == "BC_8")] <- "60-65"
meta$Age[which(meta$Age == "BC_9")] <- "40-45"
meta$Age[which(meta$Age == "BC_10")] <- "50-55"
meta$Age[which(meta$Age == "BC_11")] <- "46-50"
meta$Age[which(meta$Age == "BC_12")] <- "76-80"
meta$Age[which(meta$Age == "BC_13")] <- "60-65"
meta$Age[which(meta$Age == "BC_14")] <- "46-50"

#Seurat object =============================================

#_________________________________________________________________
rownames(meta) <- meta$cells
fdat <- toupper(as.matrix(rownames(mat.dat)))
fdat <- as.data.frame(fdat)
#_________________________________________________________________
rownames(fdat) <- fdat[,1]
fdat <- data.frame(fdat)
common_colnames <- c("gene_short_name")
colnames(fdat) <- common_colnames
#_________________________________________________________________
sobj_pre <- CreateSeuratObject(counts = mat.dat)
rownames(meta) <- colnames(mat.dat)
sobj_pre <-AddMetaData(sobj_pre,metadata=meta)
head(sobj_pre@meta.data)
#_________________________________________________________________
sobj_pre[["RNA"]]@meta.features<-fdat
head(sobj_pre[["RNA"]]@meta.features)
slotNames(sobj_pre[["RNA"]])

sobj_pre$Grade <- "NA"

saveRDS(sobj_pre, file = "pancancer_unfiltered_62622.rds")

sobj <- sobj_pre

#cell cycle ========================================================================
#from: https://satijalab.org/seurat/articles/cell_cycle_vignette.html

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes)
sobj$CC.Difference <- sobj$S.Score - sobj$G2M.Score
# view cell cycle scores and phase assignments
head(sobj[[]])

# Preprocessing and Filtering (SCTransform) --------------------------------------------------------------------


mito.genes <- grep(pattern = "^MT-", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(mito.genes)
sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")

#blood cell contamination
hemo.genes <- grep(pattern = "^HB[^(P)]", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(hemo.genes)
sobj <- PercentageFeatureSet(sobj, "^HB[^(P)]", col.name = "percent.hb")
platelet.genes <- grep(pattern = "^PECAM1|^PF4\\b", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(platelet.genes)
sobj <- PercentageFeatureSet(sobj, "^PECAM1|^PF4\\b", col.name = "percent.platelet")

nCount95 <- quantile(sobj@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
nCount5 <- quantile(sobj@meta.data$nCount_RNA, 0.05)#05)  
nFeat95 <- quantile(sobj@meta.data$nFeature_RNA, 0.95) #0.9) 
nFeat5 <- quantile(sobj@meta.data$nFeature_RNA, 0.05) #0.1) 
ptMt5 <- quantile(sobj@meta.data$percent.mt, 0.05)
ptMt90 <- quantile(sobj@meta.data$percent.mt, 0.85)
ptHb95 <- quantile(sobj@meta.data$percent.hb, 0.95)

sobj <- subset(x = sobj, subset = nFeature_RNA > nFeat5 & nFeature_RNA < nFeat95)
sobj <- subset(x = sobj, subset = percent.mt > ptMt5 & percent.mt < ptMt90)
sobj <- subset(x = sobj, subset = nCount_RNA > nCount5 & nCount_RNA < nCount95)
sobj <- subset(x = sobj, subset = percent.hb < ptHb95)

summary(sobj@meta.data$percent.mt)
summary(sobj@meta.data$percent.hb)
dim(sobj)
dim(sobj_pre)

sobj$samples <- sobj$Patient

saveRDS(sobj, file = "Qian_filtered_62622.rds")


# =================================================================== ======

# Savas ============================================================= ============
# =================================================================== ======

# Loading Raw Data into RStudio ---------------------------------- 

filePaths = getGEOSuppFiles("GSE110686") 
tarF <- list.files(path = "./GSE110686/", pattern = "*.tar", full.names = TRUE) 
untar(tarF, exdir = "./GSE110686/") 
gzipF <- list.files(path = "./GSE110686/", pattern = "*.gz", full.names = TRUE) 
ldply(.data = gzipF, .fun = gunzip) 

list.files(path = "./GSE110686/", pattern = "\\.mtx$",full.names = TRUE)
list.files(path = "./GSE110686/", pattern = "*.genes.tsv$", full.names = TRUE) 
list.files(path = "./GSE110686/", pattern = "*.barcodes.tsv$", full.names = TRUE) 


# Tumor Infiltrating Lymphocytes 20 Channel 1 ----------------------------------------------------------
matrix20_1 <- readMM(file = './GSE110686//GSM3011853_tils20_channel1_matrix.mtx')
matrix20_1 <- as.matrix(matrix20_1)
#matrix20_1[1:5,1:5]
#dim(matrix20_1)

genes20_1 <- read.table(file = './GSE110686//GSM3011853_tils20_channel1_genes.tsv', 
                        sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(genes20_1)
rownames(matrix20_1) <- genes20_1[,2]
#matrix20_1[1:5,1:5]

barcodes20_1 <- read.table(file = './GSE110686//GSM3011853_tils20_channel1_barcodes.tsv', 
                           sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes20_1)

colnames(matrix20_1) <- barcodes20_1[,1]
#matrix20_1[1:5,1:5]

colnames(matrix20_1) <- paste(colnames(matrix20_1), "TIL20_1", sep = "_")
pdat20_1 <- data.frame("samples" = colnames(matrix20_1), "Patient" = "TIL20", "Batch" = "20_1", "Molecular Subtype" = "TNBC",
                       "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium v2", 
                       "Capture Method" = "10X Genomics Chromium v2",
                       "Sequencer" = "Illumina HiSeq 2500", "Tissue Source" = "Primary Tumor", 
                       "BC Subtype" = "TNBC", "Status" = "Primary",
                       "Stage" = "NA", "Grade" = "NA",
                       "Tumor.Size" = "1.2-1.1", "Age" = "NA", "Gender" = "Female", 
                       "Parity" = "NA", "Gravidity" = "NA", "Menopause" = "NA",
                       "Treatment.Status" = "Naive", "BMI" = "NA", "Histology" = "NA")




# Tumor Infiltrating Lymphocytes 20 Channel 2 ----------------------------------------------------------

matrix20_2 <- readMM(file = './GSE110686//GSM3011853_tils20_channel2_matrix.mtx')
matrix20_2 <- as.matrix(matrix20_2)
#matrix20_2[1:5,1:5]
#dim(matrix20_2)

genes20_2 <- read.table(file = './GSE110686//GSM3011853_tils20_channel2_genes.tsv', 
                        sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(genes20_2)
rownames(matrix20_2) <- genes20_2[,2]
#matrix20_2[1:5,1:5]

barcodes20_2 <- read.table(file = './GSE110686//GSM3011853_tils20_channel2_barcodes.tsv', 
                           sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes20_2)
colnames(matrix20_2) <- barcodes20_2[,1]
#matrix20_2[1:5,1:5]

colnames(matrix20_2) <- paste(colnames(matrix20_2), "TIL20_2", sep = "_")
pdat20_2 <- data.frame("samples" = colnames(matrix20_2), "Patient" = "TIL20", "Batch" = "20_2", "Molecular Subtype" = "TNBC",
                       "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium v2", 
                       "Capture Method" = "10X Genomics Chromium v2",
                       "Sequencer" = "Illumina HiSeq 2500", "Tissue Source" = "Primary Tumor", 
                       "BC Subtype" = "TNBC", "Status" = "Primary",
                       "Stage" = "NA", "Grade" = "NA",
                       "Tumor.Size" = "1.2-1.1", "Age" = "NA", "Gender" = "Female", 
                       "Parity" = "NA", "Gravidity" = "NA", "Menopause" = "NA",
                       "Treatment.Status" = "Naive", "BMI" = "NA", "Histology" = "NA")



# Tumor Infiltrating Lymphocytes 32 ----------------------------------------------------------

matrix32 <- readMM(file = './GSE110686//GSM3011854_tils32_matrix.mtx')
matrix32 <- as.matrix(matrix32)
#matrix32[1:5,1:5]
#dim(matrix32)

genes32 <- read.table(file = './GSE110686//GSM3011854_tils32_genes.tsv', 
                      sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(genes32)
rownames(matrix32) <- genes32[,2]
#matrix32[1:5,1:5]

barcodes32 <- read.table(file = './GSE110686//GSM3011854_tils32_barcodes.tsv', 
                         sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes32)
colnames(matrix32) <- barcodes32[,1]
#matrix32[1:5,1:5]

colnames(matrix32) <- paste(colnames(matrix32), "TIL32", sep = "_")
pdat32 <- data.frame("samples" = colnames(matrix32), "Patient" = "TIL32", "Batch" = "32_1", "Molecular Subtype" = "TNBC",
                     "RNA Type" = "mRNA", "Gene Coverage" = "3'", "Library Preparation" = "10X Genomics Chromium v2", 
                     "Capture Method" = "10X Genomics Chromium v2",
                     "Sequencer" = "Illumina HiSeq 2500", "Tissue Source" = "Primary Tumor", 
                     "BC Subtype" = "TNBC", "Status" = "Primary",
                     "Stage" = "NA", "Grade" = "NA",
                     "Tumor.Size" = "1.2-1.1", "Age" = "NA", "Gender" = "Female", 
                     "Parity" = "NA", "Gravidity" = "NA", "Menopause" = "NA",
                     "Treatment.Status" = "Naive", "BMI" = "NA", "Histology" = "NA")


# Joined Matrices ----------------------------------------------------------

#case 1 = TIL20
#case 2 = TIL32
joined <- cbind(matrix20_1, matrix20_2, matrix32)
joined <- as.data.frame(joined)
#joined <- cbind(matrix20, matrix32)
dim(joined)


library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(joined)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- joined[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- joined[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}

joined <- new_mat

#pdat ===========================================

pdat_joined <- rbind(pdat20_1, pdat20_2, pdat32)


#_________________________________________________________________
rownames(pdat_joined) <- pdat_joined$samples
fdat_joined <- toupper(as.matrix(rownames(joined)))
#_________________________________________________________________
rownames(fdat_joined) <- fdat_joined[,1]
fdat_joined <- data.frame(fdat_joined)
common_colnames <- "gene_short_name"
colnames(fdat_joined) <- common_colnames
rownames(joined) <- rownames(fdat_joined)
#_________________________________________________________________
sobj_pre <- CreateSeuratObject(counts = joined)
sobj_pre <-AddMetaData(sobj_pre,metadata=pdat_joined)
head(sobj_pre@meta.data)
#_________________________________________________________________
sobj_pre[["RNA"]]@meta.features<-fdat_joined
head(sobj_pre[["RNA"]]@meta.features)
slotNames(sobj_pre[["RNA"]])

sobj_pre$celltype_minor <- "T Cells"
sobj_pre$Library.Preparation <- "Illumina TruSeq RNA Sample Prep Kit v2"

saveRDS(sobj_pre, file = "Savas_unfiltered_62622.rds")


sobj <- sobj_pre

#cell cycle ========================================================================
#from: https://satijalab.org/seurat/articles/cell_cycle_vignette.html

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes)
sobj$CC.Difference <- sobj$S.Score - sobj$G2M.Score
# view cell cycle scores and phase assignments
head(sobj[[]])

# Preprocessing and Filtering (SCTransform) --------------------------------------------------------------------

mito.genes <- grep(pattern = "^MT\\.", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(mito.genes)
sobj <- PercentageFeatureSet(sobj, pattern = "^MT\\.", col.name = "percent.mt")

#ribosomal genes
ribo.genes <- grep(pattern = "^RP[SL]", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(ribo.genes)
sobj <- PercentageFeatureSet(sobj, "^RP[SL]", col.name = "percent.ribo")

#blood cell contamination
hemo.genes <- grep(pattern = "^HB[^(P)]", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(hemo.genes)
platelet.genes <- grep(pattern = "^PECAM1|^PF4\\b", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(platelet.genes)
sobj <- PercentageFeatureSet(sobj, "^HB[^(P)]", col.name = "percent.hb")
sobj <- PercentageFeatureSet(sobj, "^PECAM1|^PF4\\b", col.name = "percent.platelet")

nCount95 <- quantile(sobj@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
nCount5 <- quantile(sobj@meta.data$nCount_RNA, 0.05)#05)  
nFeat95 <- quantile(sobj@meta.data$nFeature_RNA, 0.95) #0.9) 
nFeat5 <- quantile(sobj@meta.data$nFeature_RNA, 0.05) #0.1) 
ptMt5 <- quantile(sobj@meta.data$percent.mt, 0.05)
ptMt90 <- quantile(sobj@meta.data$percent.mt, 0.95)
ptHb95 <- quantile(sobj@meta.data$percent.hb, 0.95)

sobj <- subset(x = sobj, subset = nFeature_RNA > nFeat5 & nFeature_RNA < nFeat95)
sobj <- subset(x = sobj, subset = percent.mt > ptMt5 & percent.mt < ptMt90)
sobj <- subset(x = sobj, subset = nCount_RNA > nCount5 & nCount_RNA < nCount95)
sobj <- subset(x = sobj, subset = percent.hb < ptHb95)

summary(sobj@meta.data$percent.mt)
summary(sobj@meta.data$percent.hb)
dim(sobj)
dim(sobj_pre)

saveRDS(sobj, file = "Savas_filtered_62622.rds")

# =================================================================== ======

# Wu (2020) ========================================================= ============
# =================================================================== ======




#fullmat ------------------------------------------------

fullmat <- read.csv("/Wu Files/Wu_EMBO_countr_matrix.csv")
fullmat[1:5,1:5]

rownames(fullmat) <- fullmat[,1]
fullmat[1:5,1:5]
#actual like full matrix 
fullmetadata <- read.csv("/Wu Files/Wu_EMBO_metadata.csv")
fullmetadata[1:5,1:5]
fullmetadata <- fullmetadata[-1,]
fullmetadata <- fullmetadata[,c(1,2,4,7)]
head(fullmetadata)
fullmetadata$patientID[which(fullmetadata$patientID == "P1")] <- "Patient 1"
fullmetadata$patientID[which(fullmetadata$patientID == "P2")] <- "Patient 2"
fullmetadata$patientID[which(fullmetadata$patientID == "P3")] <- "Patient 3"
fullmetadata$patientID[which(fullmetadata$patientID == "P4")] <- "Patient 4"
fullmetadata$patientID[which(fullmetadata$patientID == "P5")] <- "Patient 5"

colnames(fullmetadata) <- c("cells","patientID", "Patient", "celltype_minor")
fullmetadata$samples <- fullmetadata$Patient
fullmetadata$Molecular.Subtype <- "TNBC"
fullmetadata$Pathologic.Stage <- "3"
fullmetadata$Gene.Coverage <- "3'"
fullmetadata$RNA.Type <- "mRNA"
fullmetadata$Library.Preparation <- "Human Tumor Dissociation Kit (Miltenyi Biotech)"
fullmetadata$Capture.Method <- "10X Genomics Single Cell 3' v2"
fullmetadata$Sequencer <- "Illumina NextSeq 500"
fullmetadata$Tissue.Source <- "Primary Tumor"
fullmetadata$BC.Subtype <- "TNBC"
fullmetadata$Status <- "Primary"



#P1 (CID44041)-------------------------------------------------------------------------

P1mat <- grep(pattern =c("^CID44041") , x = colnames(fullmat), value = TRUE)
P1mat = fullmat[,grepl(c("^CID44041"),colnames(fullmat))]
P1mat[1:5,1:5]

#P2 (CID44971)-------------------------------------------------------------------------

P2mat <- grep(pattern =c("^CID44971") , x = colnames(fullmat), value = TRUE)
P2mat = fullmat[,grepl(c("^CID44971"),colnames(fullmat))]
P2mat[1:5,1:5]

#P3 (CID44991)-------------------------------------------------------------------------

P3mat <- grep(pattern =c("^CID44991") , x = colnames(fullmat), value = TRUE)
P3mat = fullmat[,grepl(c("^CID44991"),colnames(fullmat))]
P3mat[1:5,1:5]

#P4 (CID4513)-------------------------------------------------------------------------

P4mat <- grep(pattern =c("^CID4513") , x = colnames(fullmat), value = TRUE)
P4mat = fullmat[,grepl(c("^CID4513"),colnames(fullmat))]
P4mat[1:5,1:5]


#P5 (CID4515)-------------------------------------------------------------------------

P5mat <- grep(pattern =c("^CID4515") , x = colnames(fullmat), value = TRUE)
P5mat = fullmat[,grepl(c("^CID4515"),colnames(fullmat))]
P5mat[1:5,1:5]


# joined -------------------------------------

joined <- cbind(P1mat, P2mat, P3mat, P4mat, P5mat)
joined <- as.data.frame(joined)

library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(joined)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- joined[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- joined[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}


joined <- new_mat

#create seurat object ===================================


#_________________________________________________________________
rownames(fullmetadata) <- fullmetadata$cells
fdat_joined <- toupper(as.matrix(rownames(joined)))
fullmetadata <- fullmetadata[,-1]
#_________________________________________________________________
rownames(fdat_joined) <- fdat_joined[,1]
fdat_joined <- data.frame(fdat_joined)
common_colnames <- "gene_short_name"
colnames(fdat_joined) <- common_colnames
rownames(joined) <- rownames(fdat_joined)
#_________________________________________________________________
sobj_pre <- CreateSeuratObject(counts = joined)
sobj_pre <-AddMetaData(sobj_pre,metadata=fullmetadata)
head(sobj_pre@meta.data)
#_________________________________________________________________
sobj_pre[["RNA"]]@meta.features<-fdat_joined
head(sobj_pre[["RNA"]]@meta.features)
slotNames(sobj_pre[["RNA"]])

colnames(sobj_pre@meta.data) <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "patientID", "Patient",
                                  "celltype_minor", "samples", "Molecular.Subtype", "Grade", "Gene.Coverage",
                                  "RNA.Type", "Library.Preparation", "Capture.Method", "Sequencer", "Tissue.Source",
                                  "BC.Subtype", "Status")

sobj_pre$Gender <- "Female"
sobj_pre$Stage <- "NA"
sobj_pre$Ethnicity <- "NA"
sobj_pre$Parity <- "NA"
sobj_pre$Gravidity <- "NA"
sobj_pre$Menopause <- "NA"
sobj_pre$BMI <- "NA"
sobj_pre$Tumor.Size <- "NA"
sobj_pre$Histology <- "NA"
sobj_pre$BRCA.Status <- "NA"

sobj_pre$Age <- sobj_pre$Patient
sobj_pre$Age[which(sobj_pre$Age == "Patient 1")] <- "35"
sobj_pre$Age[which(sobj_pre$Age == "Patient 2")] <- "49"
sobj_pre$Age[which(sobj_pre$Age == "Patient 3")] <- "47"
sobj_pre$Age[which(sobj_pre$Age == "Patient 4")] <- "73"
sobj_pre$Age[which(sobj_pre$Age == "Patient 5")] <- "67"


sobj_pre$Ki67 <- sobj_pre$Patient
sobj_pre$Ki67[which(sobj_pre$Ki67 == "Patient 1")] <- "70%"
sobj_pre$Ki67[which(sobj_pre$Ki67 == "Patient 2")] <- "40%"
sobj_pre$Ki67[which(sobj_pre$Ki67 == "Patient 3")] <- "NA"
sobj_pre$Ki67[which(sobj_pre$Ki67 == "Patient 4")] <- "75%"
sobj_pre$Ki67[which(sobj_pre$Ki67 == "Patient 5")] <- "60%"

sobj_pre$Treatment.Status <- sobj_pre$Patient
sobj_pre$Treatment.Status[which(sobj_pre$Treatment.Status == "Patient 1")] <- "Naive"
sobj_pre$Treatment.Status[which(sobj_pre$Treatment.Status == "Patient 2")] <- "Naive"
sobj_pre$Treatment.Status[which(sobj_pre$Treatment.Status == "Patient 3")] <- "Niave"
sobj_pre$Treatment.Status[which(sobj_pre$Treatment.Status == "Patient 4")] <- "Treated (AC)"
sobj_pre$Treatment.Status[which(sobj_pre$Treatment.Status == "Patient 5")] <- "Naive"

saveRDS(sobj_pre, file = "Wu_unfiltered_62622.rds")

sobj <- sobj_pre

#cell cycle ========================================================================
#from: https://satijalab.org/seurat/articles/cell_cycle_vignette.html

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes)
sobj$CC.Difference <- sobj$S.Score - sobj$G2M.Score
# view cell cycle scores and phase assignments
head(sobj[[]])

# Preprocessing and Filtering (SCTransform) --------------------------------------------------------------------


mito.genes <- grep(pattern = "^MT-", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(mito.genes)
sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")

#blood cell contamination
hemo.genes <- grep(pattern = "^HB[^(P)]", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(hemo.genes)
sobj <- PercentageFeatureSet(sobj, "^HB[^(P)]", col.name = "percent.hb")
platelet.genes <- grep(pattern = "^PECAM1|^PF4\\b", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(platelet.genes)
sobj <- PercentageFeatureSet(sobj, "^PECAM1|^PF4\\b", col.name = "percent.platelet")

nCount95 <- quantile(sobj@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
nCount5 <- quantile(sobj@meta.data$nCount_RNA, 0.05)#05)  
nFeat95 <- quantile(sobj@meta.data$nFeature_RNA, 0.95) #0.9) 
nFeat5 <- quantile(sobj@meta.data$nFeature_RNA, 0.05) #0.1) 
ptMt5 <- quantile(sobj@meta.data$percent.mt, 0.05)
ptMt90 <- quantile(sobj@meta.data$percent.mt, 0.95)
ptHb95 <- quantile(sobj@meta.data$percent.hb, 0.95)

sobj <- subset(x = sobj, subset = nFeature_RNA > nFeat5 & nFeature_RNA < nFeat95)
sobj <- subset(x = sobj, subset = percent.mt > ptMt5 & percent.mt < ptMt90)
sobj <- subset(x = sobj, subset = nCount_RNA > nCount5 & nCount_RNA < nCount95)
sobj <- subset(x = sobj, subset = percent.hb < ptHb95)

summary(sobj@meta.data$percent.mt)
summary(sobj@meta.data$percent.hb)
dim(sobj)
dim(sobj_pre)

saveRDS(sobj, file = "OldWu_filtered_62622.rds")


# =================================================================== ======

# Wu (2021) ========================================================= ============
# =================================================================== ======





# Loading Raw Data into RStudio ---------------------------------- 

filePaths = getGEOSuppFiles("GSE176078") 
tarF <- list.files(path = "./GSE176078/", pattern = "*.tar", full.names = TRUE) 
ldply(.data = tarF, .fun = untar, exdir = "./GSE176078/")
gzipF <- list.files(path = "./GSE176078/", pattern = "*.gz", full.names = TRUE) 
ldply(.data = gzipF, .fun = gunzip) 

list.files(path = "./GSE176078/", full.names = TRUE)
tarF <- list.files(path = "./GSE176078/", pattern = "*.tar", full.names = TRUE) 
ldply(.data = tarF, .fun = untar, exdir = "./GSE176078/")


list.files(path = "./GSE176078//CID3586", pattern = "\\.mtx$",full.names = TRUE)
list.files(path = "./GSE176078/", pattern = "*.genes.tsv$", full.names = TRUE) 
list.files(path = "./GSE176078/", pattern = "*.barcodes.tsv$", full.names = TRUE) 



# CID3941 - ER+ ----------------------------------------------------------

CID3941.mat <- readMM(file = './GSE176078//CID3941/count_matrix_sparse.mtx')
CID3941.mat <- as.matrix(CID3941.mat)
CID3941.mat[1:5,1:5]
dim(CID3941.mat)

CID3941.genes <- read.table(file = './GSE176078//CID3941/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID3941.genes)
rownames(CID3941.mat) <- CID3941.genes[,1]
CID3941.mat[1:5,1:5]

CID3941.bar <- read.table(file = './GSE176078//CID3941/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID3941.mat) <- CID3941.bar[,1]
CID3941.mat[1:5,1:5]


CID3941.meta <- read.csv(file = './GSE176078//CID3941/metadata.csv')
head(CID3941.meta)
CID3941.pdat <- CID3941.meta[,c(5,6,7,8,9)]

CID3941.pdat$Patient <- CID3941.meta$orig.ident
CID3941.pdat$Molecular.Subtype <- "ER+"
CID3941.pdat$RNA.Type <- "Total RNA"
CID3941.pdat$Gene.Coverage <- "3' and 5'"
CID3941.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3941.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3941.pdat$Sequencer <- "Illumina NextSeq 500"
CID3941.pdat$Tissue.Source <- "Primary Tumor"
CID3941.pdat$BC.Subtype <- "HR+"
CID3941.pdat$Status <- "Primary"
CID3941.pdat$Treatment.Status <- "Naive"

CID3941.pdat$Grade <- "2"
CID3941.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID3941.pdat$Other.Clinical <- "Multifocal tumour with associated high grade DCIS"
CID3941.pdat$TNM.Classification <- "pT1c, N1a, Mx"
CID3941.pdat$Ki67 <- "0.1"
CID3941.pdat$BRCA.Status <- "NA"
CID3941.pdat$Tumor.Size <- "NA"
CID3941.pdat$Age <- "50"
CID3941.pdat$Gender <- "Female"
CID3941.pdat$Parity <- "NA"
CID3941.pdat$Gravidity <- "NA"
CID3941.pdat$BMI <- "NA"
CID3941.pdat$Menopause <- "NA"
CID3941.pdat$Ethnicity <- "NA"

# CID3948 - ER+ ----------------------------------------------------------


CID3948.mat <- readMM(file = './GSE176078//CID3948/count_matrix_sparse.mtx')
CID3948.mat <- as.matrix(CID3948.mat)
CID3948.mat[1:5,1:5]
dim(CID3948.mat)

CID3948.genes <- read.table(file = './GSE176078//CID3948/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID3948.genes)
rownames(CID3948.mat) <- CID3948.genes[,1]
CID3948.mat[1:5,1:5]

CID3948.bar <- read.table(file = './GSE176078//CID3948/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID3948.mat) <- CID3948.bar[,1]
CID3948.mat[1:5,1:5]


CID3948.meta <- read.csv(file = './GSE176078//CID3948/metadata.csv')
head(CID3948.meta)
CID3948.pdat <- CID3948.meta[,c(5,6,7,8,9)]

CID3948.pdat$Patient <- CID3948.meta$orig.ident
CID3948.pdat$Molecular.Subtype <- "ER+"
CID3948.pdat$RNA.Type <- "Total RNA"
CID3948.pdat$Gene.Coverage <- "3' and 5'"
CID3948.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3948.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3948.pdat$Sequencer <- "Illumina NextSeq 500"
CID3948.pdat$Tissue.Source <- "Primary Tumor"
CID3948.pdat$BC.Subtype <- "HR+"
CID3948.pdat$Status <- "Primary"
CID3948.pdat$Treatment.Status <- "Naive"


CID3948.pdat$Grade <- "3"
CID3948.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID3948.pdat$Other.Clinical <- "Associated LCIS, with LVI and perineural invasion"
CID3948.pdat$TNM.Classification <- "pT2, N2a"
CID3948.pdat$Ki67 <- "~10%"
CID3948.pdat$BRCA.Status <- "NA"
CID3948.pdat$Tumor.Size <- "NA"
CID3948.pdat$Age <- "82"
CID3948.pdat$Gender <- "Female"
CID3948.pdat$Parity <- "NA"
CID3948.pdat$Gravidity <- "NA"
CID3948.pdat$BMI <- "NA"
CID3948.pdat$Menopause <- "NA"
CID3948.pdat$Ethnicity <- "NA"
# CID3963 - ER+ ----------------------------------------------------------

CID3963.mat <- readMM(file = './GSE176078//CID3963/count_matrix_sparse.mtx')
CID3963.mat <- as.matrix(CID3963.mat)
CID3963.mat[1:5,1:5]
dim(CID3963.mat)

CID3963.genes <- read.table(file = './GSE176078//CID3963/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID3963.genes)
rownames(CID3963.mat) <- CID3963.genes[,1]
CID3963.mat[1:5,1:5]

CID3963.bar <- read.table(file = './GSE176078//CID3963/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID3963.mat) <- CID3963.bar[,1]
CID3963.mat[1:5,1:5]


CID3963.meta <- read.csv(file = './GSE176078//CID3963/metadata.csv')
head(CID3963.meta)
CID3963.pdat <- CID3963.meta[,c(5,6,7,8,9)]

CID3963.pdat$Patient <- CID3963.meta$orig.ident
CID3963.pdat$Molecular.Subtype <- "ER+"
CID3963.pdat$RNA.Type <- "Total RNA"
CID3963.pdat$Gene.Coverage <- "3' and 5'"
CID3963.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3963.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3963.pdat$Sequencer <- "Illumina NextSeq 500"
CID3963.pdat$Tissue.Source <- "Primary Tumor"
CID3963.pdat$BC.Subtype <- "HR+"
CID3963.pdat$Status <- "Primary"
CID3963.pdat$Treatment.Status <- "Treated (AC, Paclitaxel, Herceptin (administered for Dx 3 years prior)"


CID3963.pdat$Grade <- "3"
CID3963.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID3963.pdat$Other.Clinical <- "Probable recurrence from 3 years prior"
CID3963.pdat$TNM.Classification <- "pT2, pN0, Mx, Stage IIA"
CID3963.pdat$Ki67 <- "0.43"
CID3963.pdat$BRCA.Status <- "NA"
CID3963.pdat$Tumor.Size <- "NA"
CID3963.pdat$Age <- "61"
CID3963.pdat$Gender <- "Female"
CID3963.pdat$Parity <- "NA"
CID3963.pdat$Gravidity <- "NA"
CID3963.pdat$BMI <- "NA"
CID3963.pdat$Menopause <- "NA"
CID3963.pdat$Ethnicity <- "NA"
# CID4040 - ER+ ----------------------------------------------------------

CID4040.mat <- readMM(file = './GSE176078//CID4040/count_matrix_sparse.mtx')
CID4040.mat <- as.matrix(CID4040.mat)
CID4040.mat[1:5,1:5]
dim(CID4040.mat)

CID4040.genes <- read.table(file = './GSE176078//CID4040/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4040.genes)
rownames(CID4040.mat) <- CID4040.genes[,1]
CID4040.mat[1:5,1:5]

CID4040.bar <- read.table(file = './GSE176078//CID4040/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4040.mat) <- CID4040.bar[,1]
CID4040.mat[1:5,1:5]


CID4040.meta <- read.csv(file = './GSE176078//CID4040/metadata.csv')
head(CID4040.meta)
CID4040.pdat <- CID4040.meta[,c(5,6,7,8,9)]

CID4040.pdat$Patient <- CID4040.meta$orig.ident
CID4040.pdat$Molecular.Subtype <- "ER+"
CID4040.pdat$RNA.Type <- "Total RNA"
CID4040.pdat$Gene.Coverage <- "3' and 5'"
CID4040.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4040.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4040.pdat$Sequencer <- "Illumina NextSeq 500"
CID4040.pdat$Tissue.Source <- "Primary Tumor"
CID4040.pdat$BC.Subtype <- "HR+"
CID4040.pdat$Status <- "Primary"
CID4040.pdat$Treatment.Status <- "Naive"


CID4040.pdat$Grade <- "3"
CID4040.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID4040.pdat$Other.Clinical <- "Associated high grade DCIS"
CID4040.pdat$TNM.Classification <- "pT2, N0"
CID4040.pdat$Ki67 <- ">50%"
CID4040.pdat$BRCA.Status <- "NA"
CID4040.pdat$Tumor.Size <- "NA"
CID4040.pdat$Age <- "57"
CID4040.pdat$Gender <- "Female"
CID4040.pdat$Parity <- "NA"
CID4040.pdat$Gravidity <- "NA"
CID4040.pdat$BMI <- "NA"
CID4040.pdat$Menopause <- "NA"
CID4040.pdat$Ethnicity <- "NA"

# CID4471 - ER+ - has normal epithelia ----------------------------------------------------------

CID4471.mat <- readMM(file = './GSE176078//CID4471/count_matrix_sparse.mtx')
CID4471.mat <- as.matrix(CID4471.mat)
CID4471.mat[1:5,1:5]
dim(CID4471.mat)

CID4471.genes <- read.table(file = './GSE176078//CID4471/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4471.genes)
rownames(CID4471.mat) <- CID4471.genes[,1]
CID4471.mat[1:5,1:5]

CID4471.bar <- read.table(file = './GSE176078//CID4471/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4471.mat) <- CID4471.bar[,1]
CID4471.mat[1:5,1:5]


CID4471.meta <- read.csv(file = './GSE176078//CID4471/metadata.csv')
head(CID4471.meta)
CID4471.pdat <- CID4471.meta[,c(5,6,7,8,9)]

CID4471.pdat$Patient <- CID4471.meta$orig.ident
CID4471.pdat$Molecular.Subtype <- "ER+"
CID4471.pdat$RNA.Type <- "Total RNA"
CID4471.pdat$Gene.Coverage <- "3' and 5'"
CID4471.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4471.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4471.pdat$Sequencer <- "Illumina NextSeq 500"
CID4471.pdat$Tissue.Source <- "Primary Tumor"
CID4471.pdat$BC.Subtype <- "HR+"
CID4471.pdat$Status <- "Primary"
CID4471.pdat$Treatment.Status <- "Naive"


CID4471.pdat$Grade <- "2"
CID4471.pdat$Histology <- "Invasive Lobular Carcinoma (ILC)"
CID4471.pdat$Other.Clinical <- "NA"
CID4471.pdat$TNM.Classification <- "pT3, pN0 (i+)"
CID4471.pdat$Ki67 <- "0.2"
CID4471.pdat$BRCA.Status <- "NA"
CID4471.pdat$Tumor.Size <- "NA"
CID4471.pdat$Age <- "55"
CID4471.pdat$Gender <- "Female"
CID4471.pdat$Parity <- "NA"
CID4471.pdat$Gravidity <- "NA"
CID4471.pdat$BMI <- "NA"
CID4471.pdat$Menopause <- "NA"
CID4471.pdat$Ethnicity <- "NA"

# CID4067 - ER+ ----------------------------------------------------------

CID4067.mat <- readMM(file = './GSE176078//CID4067/count_matrix_sparse.mtx')
CID4067.mat <- as.matrix(CID4067.mat)
CID4067.mat[1:5,1:5]
dim(CID4067.mat)

CID4067.genes <- read.table(file = './GSE176078//CID4067/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4067.genes)
rownames(CID4067.mat) <- CID4067.genes[,1]
CID4067.mat[1:5,1:5]

CID4067.bar <- read.table(file = './GSE176078//CID4067/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4067.mat) <- CID4067.bar[,1]
CID4067.mat[1:5,1:5]


CID4067.meta <- read.csv(file = './GSE176078//CID4067/metadata.csv')
head(CID4067.meta)
CID4067.pdat <- CID4067.meta[,c(5,6,7,8,9)]

CID4067.pdat$Patient <- CID4067.meta$orig.ident
CID4067.pdat$Molecular.Subtype <- "ER+"
CID4067.pdat$RNA.Type <- "Total RNA"
CID4067.pdat$Gene.Coverage <- "3' and 5'"
CID4067.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4067.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4067.pdat$Sequencer <- "Illumina NextSeq 500"
CID4067.pdat$Tissue.Source <- "Primary Tumor"
CID4067.pdat$BC.Subtype <- "HR+"
CID4067.pdat$Status <- "Primary"
CID4067.pdat$Treatment.Status <- "Naive"


CID4067.pdat$Grade <- "2"
CID4067.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID4067.pdat$Other.Clinical <- "Associated low grade DCIS and focal perineural invasion"
CID4067.pdat$TNM.Classification <- "pT2, N1(sn), Mx"
CID4067.pdat$Ki67 <- "3-4%"
CID4067.pdat$BRCA.Status <- "NA"
CID4067.pdat$Tumor.Size <- "NA"
CID4067.pdat$Age <- "85"
CID4067.pdat$Gender <- "Female"
CID4067.pdat$Parity <- "NA"
CID4067.pdat$Gravidity <- "NA"
CID4067.pdat$BMI <- "NA"
CID4067.pdat$Menopause <- "NA"
CID4067.pdat$Ethnicity <- "NA"

# CID4290A - ER+ - has normal epithelia ----------------------------------------------------------

CID4290A.mat <- readMM(file = './GSE176078//CID4290A/count_matrix_sparse.mtx')
CID4290A.mat <- as.matrix(CID4290A.mat)
CID4290A.mat[1:5,1:5]
dim(CID4290A.mat)

CID4290A.genes <- read.table(file = './GSE176078//CID4290A/count_matrix_genes.tsv', 
                             sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4290A.genes)
rownames(CID4290A.mat) <- CID4290A.genes[,1]
CID4290A.mat[1:5,1:5]

CID4290A.bar <- read.table(file = './GSE176078//CID4290A/count_matrix_barcodes.tsv', 
                           sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4290A.mat) <- CID4290A.bar[,1]
CID4290A.mat[1:5,1:5]


CID4290A.meta <- read.csv(file = './GSE176078//CID4290A/metadata.csv')
head(CID4290A.meta)
CID4290A.pdat <- CID4290A.meta[,c(5,6,7,8,9)]

CID4290A.pdat$Patient <- CID4290A.meta$orig.ident
CID4290A.pdat$Molecular.Subtype <- "ER+"
CID4290A.pdat$RNA.Type <- "Total RNA"
CID4290A.pdat$Gene.Coverage <- "3' and 5'"
CID4290A.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4290A.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4290A.pdat$Sequencer <- "Illumina NextSeq 500"
CID4290A.pdat$Tissue.Source <- "Primary Tumor"
CID4290A.pdat$BC.Subtype <- "HR+"
CID4290A.pdat$Status <- "Primary"
CID4290A.pdat$Treatment.Status <- "Naive"

CID4290A.pdat$Grade <- "2"
CID4290A.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID4290A.pdat$Other.Clinical <- "Locally advanced, skin and chest wall muscle involvement"
CID4290A.pdat$TNM.Classification <- "pT4b, Nx"
CID4290A.pdat$Ki67 <- "0.1"
CID4290A.pdat$BRCA.Status <- "NA"
CID4290A.pdat$Tumor.Size <- "NA"
CID4290A.pdat$Age <- "88"
CID4290A.pdat$Gender <- "Female"
CID4290A.pdat$Parity <- "NA"
CID4290A.pdat$Gravidity <- "NA"
CID4290A.pdat$BMI <- "NA"
CID4290A.pdat$Menopause <- "NA"
CID4290A.pdat$Ethnicity <- "NA"

# CID4398 - ER+ ----------------------------------------------------------

CID4398.mat <- readMM(file = './GSE176078//CID4398/count_matrix_sparse.mtx')
CID4398.mat <- as.matrix(CID4398.mat)
CID4398.mat[1:5,1:5]
dim(CID4398.mat)

CID4398.genes <- read.table(file = './GSE176078//CID4398/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4398.genes)
rownames(CID4398.mat) <- CID4398.genes[,1]
CID4398.mat[1:5,1:5]

CID4398.bar <- read.table(file = './GSE176078//CID4398/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4398.mat) <- CID4398.bar[,1]
CID4398.mat[1:5,1:5]


CID4398.meta <- read.csv(file = './GSE176078//CID4398/metadata.csv')
head(CID4398.meta)
CID4398.pdat <- CID4398.meta[,c(5,6,7,8,9)]

CID4398.pdat$Patient <- CID4398.meta$orig.ident
CID4398.pdat$Molecular.Subtype <- "ER+"
CID4398.pdat$RNA.Type <- "Total RNA"
CID4398.pdat$Gene.Coverage <- "3' and 5'"
CID4398.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4398.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4398.pdat$Sequencer <- "Illumina NextSeq 500"
CID4398.pdat$Tissue.Source <- "Primary Tumor"
CID4398.pdat$BC.Subtype <- "HR+"
CID4398.pdat$Status <- "Primary"
CID4398.pdat$Treatment.Status <- "Treated (Neoadjuvant FEC-D)"


CID4398.pdat$Grade <- "3"
CID4398.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID4398.pdat$Other.Clinical <- "Mixed morphology with associated high grade DCIS, extensive LVI and perineural invasion. RCB-III, minimal or no-response to chemotherapy"
CID4398.pdat$TNM.Classification <- "pT3, pN2a, pMx, Stage IIIA"
CID4398.pdat$Ki67 <- "0.75"
CID4398.pdat$BRCA.Status <- "NA"
CID4398.pdat$Tumor.Size <- "NA"
CID4398.pdat$Age <- "52"
CID4398.pdat$Gender <- "Female"
CID4398.pdat$Parity <- "NA"
CID4398.pdat$Gravidity <- "NA"
CID4398.pdat$BMI <- "NA"
CID4398.pdat$Menopause <- "NA"
CID4398.pdat$Ethnicity <- "NA"
# CID4530N - ER+ - has normal epithelia ----------------------------------------------------------

CID4530N.mat <- readMM(file = './GSE176078//CID4530N/count_matrix_sparse.mtx')
CID4530N.mat <- as.matrix(CID4530N.mat)
CID4530N.mat[1:5,1:5]
dim(CID4530N.mat)

CID4530N.genes <- read.table(file = './GSE176078//CID4530N/count_matrix_genes.tsv', 
                             sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4530N.genes)
rownames(CID4530N.mat) <- CID4530N.genes[,1]
CID4530N.mat[1:5,1:5]

CID4530N.bar <- read.table(file = './GSE176078//CID4530N/count_matrix_barcodes.tsv', 
                           sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4530N.mat) <- CID4530N.bar[,1]
CID4530N.mat[1:5,1:5]


CID4530N.meta <- read.csv(file = './GSE176078//CID4530N/metadata.csv')
head(CID4530N.meta)
CID4530N.pdat <- CID4530N.meta[,c(5,6,7,8,9)]

CID4530N.pdat$Patient <- CID4530N.meta$orig.ident
CID4530N.pdat$Molecular.Subtype <- "ER+"
CID4530N.pdat$RNA.Type <- "Total RNA"
CID4530N.pdat$Gene.Coverage <- "3' and 5'"
CID4530N.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4530N.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4530N.pdat$Sequencer <- "Illumina NextSeq 500"
CID4530N.pdat$Tissue.Source <- "Primary Tumor"
CID4530N.pdat$BC.Subtype <- "ER+"
CID4530N.pdat$Status <- "Primary"
CID4530N.pdat$Treatment.Status <- "Naive"

CID4530N.pdat$Grade <- "2"
CID4530N.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID4530N.pdat$Other.Clinical <- "Multifocal tumour with associated high grade DCIS and LVI"
CID4530N.pdat$TNM.Classification <- "pT3, pN3, pMx, Stage IIIA"
CID4530N.pdat$Ki67 <- "0.05"
CID4530N.pdat$BRCA.Status <- "NA"
CID4530N.pdat$Tumor.Size <- "NA"
CID4530N.pdat$Age <- "42"
CID4530N.pdat$Gender <- "Female"
CID4530N.pdat$Parity <- "NA"
CID4530N.pdat$Gravidity <- "NA"
CID4530N.pdat$BMI <- "NA"
CID4530N.pdat$Menopause <- "NA"
CID4530N.pdat$Ethnicity <- "NA"

# CID4535 - ER+ - has normal epithelia ----------------------------------------------------------

CID4535.mat <- readMM(file = './GSE176078//CID4535/count_matrix_sparse.mtx')
CID4535.mat <- as.matrix(CID4535.mat)
CID4535.mat[1:5,1:5]
dim(CID4535.mat)

CID4535.genes <- read.table(file = './GSE176078//CID4535/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4535.genes)
rownames(CID4535.mat) <- CID4535.genes[,1]
CID4535.mat[1:5,1:5]

CID4535.bar <- read.table(file = './GSE176078//CID4535/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4535.mat) <- CID4535.bar[,1]
CID4535.mat[1:5,1:5]


CID4535.meta <- read.csv(file = './GSE176078//CID4535/metadata.csv')
head(CID4535.meta)
CID4535.pdat <- CID4535.meta[,c(5,6,7,8, 9)]

CID4535.pdat$Patient <- CID4535.meta$orig.ident
CID4535.pdat$Molecular.Subtype <- "ER+"
CID4535.pdat$RNA.Type <- "Total RNA"
CID4535.pdat$Gene.Coverage <- "3' and 5'"
CID4535.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4535.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4535.pdat$Sequencer <- "Illumina NextSeq 500"
CID4535.pdat$Tissue.Source <- "Primary Tumor"
CID4535.pdat$BC.Subtype <- "ER+"
CID4535.pdat$Status <- "Primary"
CID4535.pdat$Treatment.Status <- "Naive"

CID4535.pdat$Grade <- "2"
CID4535.pdat$Histology <- "Invasive Lobular Carcinoma (ILC)"
CID4535.pdat$Other.Clinical <- "NA"
CID4535.pdat$TNM.Classification <- "pT2, pN0 (i+),Stage IIB"
CID4535.pdat$Ki67 <- "0.1"
CID4535.pdat$BRCA.Status <- "NA"
CID4535.pdat$Tumor.Size <- "NA"
CID4535.pdat$Age <- "47"
CID4535.pdat$Gender <- "Female"
CID4535.pdat$Parity <- "NA"
CID4535.pdat$Gravidity <- "NA"
CID4535.pdat$BMI <- "NA"
CID4535.pdat$Menopause <- "NA"
CID4535.pdat$Ethnicity <- "NA"

# CID4461 - ER+ ----------------------------------------------------------

CID4461.mat <- readMM(file = './GSE176078//CID4461/count_matrix_sparse.mtx')
CID4461.mat <- as.matrix(CID4461.mat)
CID4461.mat[1:5,1:5]
dim(CID4461.mat)

CID4461.genes <- read.table(file = './GSE176078//CID4461/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4461.genes)
rownames(CID4461.mat) <- CID4461.genes[,1]
CID4461.mat[1:5,1:5]

CID4461.bar <- read.table(file = './GSE176078//CID4461/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4461.mat) <- CID4461.bar[,1]
CID4461.mat[1:5,1:5]


CID4461.meta <- read.csv(file = './GSE176078//CID4461/metadata.csv')
head(CID4461.meta)
CID4461.pdat <- CID4461.meta[,c(5,6,7,8, 9)]

CID4461.pdat$Patient <- CID4461.meta$orig.ident
CID4461.pdat$Molecular.Subtype <- "ER+"
CID4461.pdat$RNA.Type <- "Total RNA"
CID4461.pdat$Gene.Coverage <- "3' and 5'"
CID4461.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4461.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4461.pdat$Sequencer <- "Illumina NextSeq 500"
CID4461.pdat$Tissue.Source <- "Primary Tumor"
CID4461.pdat$BC.Subtype <- "HR+"
CID4461.pdat$Status <- "Primary"
CID4461.pdat$Treatment.Status <- "Naive"

CID4461.pdat$Grade <- "2"
CID4461.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID4461.pdat$Other.Clinical <- "Associated intermediate to high grade DCIS, LVI and perineural invasion"
CID4461.pdat$TNM.Classification <- "pT3, N1a, Mx"
CID4461.pdat$Ki67 <- "0.15"
CID4461.pdat$BRCA.Status <- "NA"
CID4461.pdat$Tumor.Size <- "NA"
CID4461.pdat$Age <- "54"
CID4461.pdat$Gender <- "Female"
CID4461.pdat$Parity <- "NA"
CID4461.pdat$Gravidity <- "NA"
CID4461.pdat$BMI <- "NA"
CID4461.pdat$Menopause <- "NA"
CID4461.pdat$Ethnicity <- "NA"

# CID4463 - ER+ - has normal epithelia ----------------------------------------------------------

CID4463.mat <- readMM(file = './GSE176078//CID4463/count_matrix_sparse.mtx')
CID4463.mat <- as.matrix(CID4463.mat)
CID4463.mat[1:5,1:5]
dim(CID4463.mat)

CID4463.genes <- read.table(file = './GSE176078//CID4463/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4463.genes)
rownames(CID4463.mat) <- CID4463.genes[,1]
CID4463.mat[1:5,1:5]

CID4463.bar <- read.table(file = './GSE176078//CID4463/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4463.mat) <- CID4463.bar[,1]
CID4463.mat[1:5,1:5]


CID4463.meta <- read.csv(file = './GSE176078//CID4463/metadata.csv')
head(CID4463.meta)
CID4463.pdat <- CID4463.meta[,c(5,6,7,8,9)]

CID4463.pdat$Patient <- CID4463.meta$orig.ident
CID4463.pdat$Molecular.Subtype <- "ER+"
CID4463.pdat$RNA.Type <- "Total RNA"
CID4463.pdat$Gene.Coverage <- "3' and 5'"
CID4463.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4463.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4463.pdat$Sequencer <- "Illumina NextSeq 500"
CID4463.pdat$Tissue.Source <- "Primary Tumor"
CID4463.pdat$BC.Subtype <- "HR+"
CID4463.pdat$Status <- "Primary"
CID4463.pdat$Treatment.Status <- "Naive"


CID4463.pdat$Grade <- "2"
CID4463.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID4463.pdat$Other.Clinical <- "IDC with areas of lobular -like growth pattern, but is E-cadherin positive. Associated low through high grade DCIS and LVI"
CID4463.pdat$TNM.Classification <- "pT3, N1, Mx"
CID4463.pdat$Ki67 <- "0.5"
CID4463.pdat$BRCA.Status <- "NA"
CID4463.pdat$Tumor.Size <- "NA"
CID4463.pdat$Age <- "58"
CID4463.pdat$Gender <- "Female"
CID4463.pdat$Parity <- "NA"
CID4463.pdat$Gravidity <- "NA"
CID4463.pdat$BMI <- "NA"
CID4463.pdat$Menopause <- "NA"
CID4463.pdat$Ethnicity <- "NA"




# CID3586 - HER2+ - has normal epithelia ----------------------------------------------------------

CID3586.mat <- readMM(file = './GSE176078//CID3586/count_matrix_sparse.mtx')
CID3586.mat <- as.matrix(CID3586.mat)
CID3586.mat[1:5,1:5]
dim(CID3586.mat)

CID3586.genes <- read.table(file = './GSE176078//CID3586/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID3586.genes)
rownames(CID3586.mat) <- CID3586.genes[,1]
CID3586.mat[1:5,1:5]

CID3586.bar <- read.table(file = './GSE176078//CID3586/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID3586.mat) <- CID3586.bar[,1]
CID3586.mat[1:5,1:5]


CID3586.meta <- read.csv(file = './GSE176078//CID3586/metadata.csv')
head(CID3586.meta)
CID3586.pdat <- CID3586.meta[,c(5,6,7,8, 9)]

CID3586.pdat$Patient <- CID3586.meta$orig.ident
CID3586.pdat$Molecular.Subtype <- "HER2+/ER+"
CID3586.pdat$RNA.Type <- "Total RNA"
CID3586.pdat$Gene.Coverage <- "3' and 5'"
CID3586.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3586.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3586.pdat$Sequencer <- "Illumina NextSeq 500"
CID3586.pdat$Tissue.Source <- "Primary Tumor"
CID3586.pdat$BC.Subtype <- "HER2+"
CID3586.pdat$Status <- "Primary"
CID3586.pdat$Treatment.Status <- "Naive"

CID3586.pdat$Grade <- "3"
CID3586.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID3586.pdat$Other.Clinical <- "Multifocal tumour with associatied high grade DCIS and extensive LVI"
CID3586.pdat$TNM.Classification <- "pT(m)2, N2a"
CID3586.pdat$Ki67 <- "30-50%"
CID3586.pdat$BRCA.Status <- "NA"
CID3586.pdat$Tumor.Size <- "NA"
CID3586.pdat$Age <- "43"
CID3586.pdat$Gender <- "Female"
CID3586.pdat$Parity <- "NA"
CID3586.pdat$Gravidity <- "NA"
CID3586.pdat$BMI <- "NA"
CID3586.pdat$Menopause <- "NA"
CID3586.pdat$Ethnicity <- "NA"

# CID3838 - HER2+ ----------------------------------------------------------

CID3838.mat <- readMM(file = './GSE176078//CID3838/count_matrix_sparse.mtx')
CID3838.mat <- as.matrix(CID3838.mat)
CID3838.mat[1:5,1:5]
dim(CID3838.mat)

CID3838.genes <- read.table(file = './GSE176078//CID3838/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID3838.genes)
rownames(CID3838.mat) <- CID3838.genes[,1]
CID3838.mat[1:5,1:5]

CID3838.bar <- read.table(file = './GSE176078//CID3838/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID3838.mat) <- CID3838.bar[,1]
CID3838.mat[1:5,1:5]

CID3838.meta <- read.csv(file = './GSE176078//CID3838/metadata.csv')
head(CID3838.meta)
CID3838.pdat <- CID3838.meta[,c(5,6,7,8,9)]

CID3838.pdat$Patient <- CID3838.meta$orig.ident
CID3838.pdat$Molecular.Subtype <- "HER2+"
CID3838.pdat$RNA.Type <- "Total RNA"
CID3838.pdat$Gene.Coverage <- "3' and 5'"
CID3838.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3838.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3838.pdat$Sequencer <- "Illumina NextSeq 500"
CID3838.pdat$Tissue.Source <- "Primary Tumor"
CID3838.pdat$BC.Subtype <- "HER2+"
CID3838.pdat$Status <- "Primary"
CID3838.pdat$Treatment.Status <- "Naive"

CID3838.pdat$Grade <- "3"
CID3838.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID3838.pdat$Other.Clinical <- "Associated high grade DCIS"
CID3838.pdat$TNM.Classification <- "pT2, N1a"
CID3838.pdat$Ki67 <- "0.6"
CID3838.pdat$BRCA.Status <- "NA"
CID3838.pdat$Tumor.Size <- "NA"
CID3838.pdat$Age <- "49"
CID3838.pdat$Gender <- "Female"
CID3838.pdat$Parity <- "NA"
CID3838.pdat$Gravidity <- "NA"
CID3838.pdat$BMI <- "NA"
CID3838.pdat$Menopause <- "NA"
CID3838.pdat$Ethnicity <- "NA"

# CID3921 - HER2+ ----------------------------------------------------------

CID3921.mat <- readMM(file = './GSE176078//CID3921/count_matrix_sparse.mtx')
CID3921.mat <- as.matrix(CID3921.mat)
CID3921.mat[1:5,1:5]
dim(CID3921.mat)

CID3921.genes <- read.table(file = './GSE176078//CID3921/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID3921.genes)
rownames(CID3921.mat) <- CID3921.genes[,1]
CID3921.mat[1:5,1:5]

CID3921.bar <- read.table(file = './GSE176078//CID3921/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID3921.mat) <- CID3921.bar[,1]
CID3921.mat[1:5,1:5]


CID3921.meta <- read.csv(file = './GSE176078//CID3921/metadata.csv')
head(CID3921.meta)
CID3921.pdat <- CID3921.meta[,c(5,6,7,8, 9)]

CID3921.pdat$Patient <- CID3921.meta$orig.ident
CID3921.pdat$Molecular.Subtype <- "HER2+"
CID3921.pdat$RNA.Type <- "Total RNA"
CID3921.pdat$Gene.Coverage <- "3' and 5'"
CID3921.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3921.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3921.pdat$Sequencer <- "Illumina NextSeq 500"
CID3921.pdat$Tissue.Source <- "Primary Tumor"
CID3921.pdat$BC.Subtype <- "HER2+"
CID3921.pdat$Status <- "Primary"
CID3921.pdat$Treatment.Status <- "Naive"


CID3921.pdat$Grade <- "3"
CID3921.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID3921.pdat$Other.Clinical <- "Associated high grade DCIS and focal LVI"
CID3921.pdat$TNM.Classification <- "pT2, N2a (Stage IIIA)"
CID3921.pdat$Ki67 <- ">50%"
CID3921.pdat$BRCA.Status <- "NA"
CID3921.pdat$Tumor.Size <- "NA"
CID3921.pdat$Age <- "60"
CID3921.pdat$Gender <- "Female"
CID3921.pdat$Parity <- "NA"
CID3921.pdat$Gravidity <- "NA"
CID3921.pdat$BMI <- "NA"
CID3921.pdat$Menopause <- "NA"
CID3921.pdat$Ethnicity <- "NA"

# CID4066 - HER2+ - has normal epithelia ----------------------------------------------------------

CID4066.mat <- readMM(file = './GSE176078//CID4066/count_matrix_sparse.mtx')
CID4066.mat <- as.matrix(CID4066.mat)
CID4066.mat[1:5,1:5]
dim(CID4066.mat)

CID4066.genes <- read.table(file = './GSE176078//CID4066/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4066.genes)
rownames(CID4066.mat) <- CID4066.genes[,1]
CID4066.mat[1:5,1:5]

CID4066.bar <- read.table(file = './GSE176078//CID4066/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4066.mat) <- CID4066.bar[,1]
CID4066.mat[1:5,1:5]


CID4066.meta <- read.csv(file = './GSE176078//CID4066/metadata.csv')
head(CID4066.meta)
CID4066.pdat <- CID4066.meta[,c(5,6,7,8,9)]

CID4066.pdat$Patient <- CID4066.meta$orig.ident
CID4066.pdat$Molecular.Subtype <- "HER2+/ER+"
CID4066.pdat$RNA.Type <- "Total RNA"
CID4066.pdat$Gene.Coverage <- "3' and 5'"
CID4066.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4066.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4066.pdat$Sequencer <- "Illumina NextSeq 500"
CID4066.pdat$Tissue.Source <- "Primary Tumor"
CID4066.pdat$BC.Subtype <- "HER2+"
CID4066.pdat$Status <- "Primary"
CID4066.pdat$Treatment.Status <- "Treated (Neoadjuvant AC)"


CID4066.pdat$Grade <- "2"
CID4066.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID4066.pdat$Other.Clinical <- "Associated high grade DCIS and extensive LVI. RCB-III, minimal or no-response to chemotherapy"
CID4066.pdat$TNM.Classification <- "pT2 N2a Mx"
CID4066.pdat$Ki67 <- "0.3"
CID4066.pdat$BRCA.Status <- "NA"
CID4066.pdat$Tumor.Size <- "NA"
CID4066.pdat$Age <- "41"
CID4066.pdat$Gender <- "Female"
CID4066.pdat$Parity <- "NA"
CID4066.pdat$Gravidity <- "NA"
CID4066.pdat$BMI <- "NA"
CID4066.pdat$Menopause <- "NA"
CID4066.pdat$Ethnicity <- "NA"

# CID45171 - HER2+ ----------------------------------------------------------

CID45171.mat <- readMM(file = './GSE176078//CID45171/count_matrix_sparse.mtx')
CID45171.mat <- as.matrix(CID45171.mat)
CID45171.mat[1:5,1:5]
dim(CID45171.mat)

CID45171.genes <- read.table(file = './GSE176078//CID45171/count_matrix_genes.tsv', 
                             sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID45171.genes)
rownames(CID45171.mat) <- CID45171.genes[,1]
CID45171.mat[1:5,1:5]

CID45171.bar <- read.table(file = './GSE176078//CID45171/count_matrix_barcodes.tsv', 
                           sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID45171.mat) <- CID45171.bar[,1]
CID45171.mat[1:5,1:5]


CID45171.meta <- read.csv(file = './GSE176078//CID45171/metadata.csv')
head(CID45171.meta)
CID45171.pdat <- CID45171.meta[,c(5,6,7,8,9)]

CID45171.pdat$Patient <- CID45171.meta$orig.ident
CID45171.pdat$Molecular.Subtype <- "HER2+"
CID45171.pdat$RNA.Type <- "Total RNA"
CID45171.pdat$Gene.Coverage <- "3' and 5'"
CID45171.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID45171.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID45171.pdat$Sequencer <- "Illumina NextSeq 500"
CID45171.pdat$Tissue.Source <- "Primary Tumor"
CID45171.pdat$BC.Subtype <- "HER2+"
CID45171.pdat$Status <- "Primary"
CID45171.pdat$Treatment.Status <- "Naive"



CID45171.pdat$Grade <- "3"
CID45171.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID45171.pdat$Other.Clinical <- "NA"
CID45171.pdat$TNM.Classification <- "NA"
CID45171.pdat$Ki67 <- "0.8"
CID45171.pdat$BRCA.Status <- "NA"
CID45171.pdat$Tumor.Size <- "NA"
CID45171.pdat$Age <- "58"
CID45171.pdat$Gender <- "Female"
CID45171.pdat$Parity <- "NA"
CID45171.pdat$Gravidity <- "NA"
CID45171.pdat$BMI <- "NA"
CID45171.pdat$Menopause <- "NA"
CID45171.pdat$Ethnicity <- "NA"


# CID44041 - TNBC - has normal epithelia ----------------------------------------------------------

CID44041.mat <- readMM(file = './GSE176078//CID44041/count_matrix_sparse.mtx')
CID44041.mat <- as.matrix(CID44041.mat)
CID44041.mat[1:5,1:5]
dim(CID44041.mat)

CID44041.genes <- read.table(file = './GSE176078//CID44041/count_matrix_genes.tsv', 
                             sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID44041.genes)
rownames(CID44041.mat) <- CID44041.genes[,1]
CID44041.mat[1:5,1:5]

CID44041.bar <- read.table(file = './GSE176078//CID44041/count_matrix_barcodes.tsv', 
                           sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID44041.mat) <- CID44041.bar[,1]
CID44041.mat[1:5,1:5]


CID44041.meta <- read.csv(file = './GSE176078//CID44041/metadata.csv')
head(CID44041.meta)
CID44041.pdat <- CID44041.meta[,c(5,6,7,8,9)]

CID44041.pdat$Patient <- CID44041.meta$orig.ident
CID44041.pdat$Molecular.Subtype <- "TNBC"
CID44041.pdat$RNA.Type <- "Total RNA"
CID44041.pdat$Gene.Coverage <- "3' and 5'"
CID44041.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID44041.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID44041.pdat$Sequencer <- "Illumina NextSeq 500"
CID44041.pdat$Tissue.Source <- "Primary Tumor"
CID44041.pdat$BC.Subtype <- "TNBC"
CID44041.pdat$Status <- "Primary"
CID44041.pdat$Treatment.Status <- "Naive"

CID44041.pdat$Grade <- "3"
CID44041.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID44041.pdat$Other.Clinical <- "Associated high grade DCIS and focal LVI"
CID44041.pdat$TNM.Classification <- "pT2, N1a, Mx"
CID44041.pdat$Ki67 <- "0.7"
CID44041.pdat$BRCA.Status <- "NA"
CID44041.pdat$Tumor.Size <- "NA"
CID44041.pdat$Age <- "35"
CID44041.pdat$Gender <- "Female"
CID44041.pdat$Parity <- "NA"
CID44041.pdat$Gravidity <- "NA"
CID44041.pdat$BMI <- "NA"
CID44041.pdat$Menopause <- "NA"
CID44041.pdat$Ethnicity <- "NA"

# CID4465 - TNBC - has normal epithelia ----------------------------------------------------------

CID4465.mat <- readMM(file = './GSE176078//CID4465/count_matrix_sparse.mtx')
CID4465.mat <- as.matrix(CID4465.mat)
CID4465.mat[1:5,1:5]
dim(CID4465.mat)

CID4465.genes <- read.table(file = './GSE176078//CID4465/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4465.genes)
rownames(CID4465.mat) <- CID4465.genes[,1]
CID4465.mat[1:5,1:5]

CID4465.bar <- read.table(file = './GSE176078//CID4465/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4465.mat) <- CID4465.bar[,1]
CID4465.mat[1:5,1:5]


CID4465.meta <- read.csv(file = './GSE176078//CID4465/metadata.csv')
head(CID4465.meta)
CID4465.pdat <- CID4465.meta[,c(5,6,7,8,9)]

CID4465.pdat$Patient <- CID4465.meta$orig.ident
CID4465.pdat$Molecular.Subtype <- "TNBC"
CID4465.pdat$RNA.Type <- "Total RNA"
CID4465.pdat$Gene.Coverage <- "3' and 5'"
CID4465.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4465.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4465.pdat$Sequencer <- "Illumina NextSeq 500"
CID4465.pdat$Tissue.Source <- "Primary Tumor"
CID4465.pdat$BC.Subtype <- "TNBC"
CID4465.pdat$Status <- "Primary"
CID4465.pdat$Treatment.Status <- "Naive"

CID4465.pdat$Grade <- "3"
CID4465.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID4465.pdat$Other.Clinical <- "Basal phenotype - patchy CK5/6 and p63 positivity. Associated high grade DCIS at periphery of tumour mass"
CID4465.pdat$TNM.Classification <- "PT2, N0(sn) Mx"
CID4465.pdat$Ki67 <- "0.7"
CID4465.pdat$BRCA.Status <- "NA"
CID4465.pdat$Tumor.Size <- "NA"
CID4465.pdat$Age <- "54"
CID4465.pdat$Gender <- "Female"
CID4465.pdat$Parity <- "NA"
CID4465.pdat$Gravidity <- "NA"
CID4465.pdat$BMI <- "NA"
CID4465.pdat$Menopause <- "NA"
CID4465.pdat$Ethnicity <- "NA"


# CID3946 - TNBC ----------------------------------------------------------

CID3946.mat <- readMM(file = './GSE176078//CID3946/count_matrix_sparse.mtx')
CID3946.mat <- as.matrix(CID3946.mat)
CID3946.mat[1:5,1:5]
dim(CID3946.mat)

CID3946.genes <- read.table(file = './GSE176078//CID3946/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID3946.genes)
rownames(CID3946.mat) <- CID3946.genes[,1]
CID3946.mat[1:5,1:5]

CID3946.bar <- read.table(file = './GSE176078//CID3946/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID3946.mat) <- CID3946.bar[,1]
CID3946.mat[1:5,1:5]


CID3946.meta <- read.csv(file = './GSE176078//CID3946/metadata.csv')
head(CID3946.meta)
CID3946.pdat <- CID3946.meta[,c(5,6,7,8,9)]

CID3946.pdat$Patient <- CID3946.meta$orig.ident
CID3946.pdat$Molecular.Subtype <- "TNBC"
CID3946.pdat$RNA.Type <- "Total RNA"
CID3946.pdat$Gene.Coverage <- "3' and 5'"
CID3946.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3946.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID3946.pdat$Sequencer <- "Illumina NextSeq 500"
CID3946.pdat$Tissue.Source <- "Primary Tumor"
CID3946.pdat$BC.Subtype <- "TNBC"
CID3946.pdat$Status <- "Primary"
CID3946.pdat$Treatment.Status <- "Naive"

CID3946.pdat$Grade <- "3"
CID3946.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID3946.pdat$Other.Clinical <- "Basal phenotype. Reactive lymphoid inflitrate with germinal centres"
CID3946.pdat$TNM.Classification <- "pT2, N0, Mx"
CID3946.pdat$Ki67 <- "0.6"
CID3946.pdat$BRCA.Status <- "NA"
CID3946.pdat$Tumor.Size <- "NA"
CID3946.pdat$Age <- "52"
CID3946.pdat$Gender <- "Female"
CID3946.pdat$Parity <- "NA"
CID3946.pdat$Gravidity <- "NA"
CID3946.pdat$BMI <- "NA"
CID3946.pdat$Menopause <- "NA"
CID3946.pdat$Ethnicity <- "NA"

# CID4495 - TNBC ----------------------------------------------------------

CID4495.mat <- readMM(file = './GSE176078//CID4495/count_matrix_sparse.mtx')
CID4495.mat <- as.matrix(CID4495.mat)
CID4495.mat[1:5,1:5]
dim(CID4495.mat)

CID4495.genes <- read.table(file = './GSE176078//CID4495/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4495.genes)
rownames(CID4495.mat) <- CID4495.genes[,1]
CID4495.mat[1:5,1:5]

CID4495.bar <- read.table(file = './GSE176078//CID4495/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4495.mat) <- CID4495.bar[,1]
CID4495.mat[1:5,1:5]


CID4495.meta <- read.csv(file = './GSE176078//CID4495/metadata.csv')
head(CID4495.meta)
CID4495.pdat <- CID4495.meta[,c(5,6,7,8,9)]

CID4495.pdat$Patient <- CID4495.meta$orig.ident
CID4495.pdat$Molecular.Subtype <- "TNBC"
CID4495.pdat$RNA.Type <- "Total RNA"
CID4495.pdat$Gene.Coverage <- "3' and 5'"
CID4495.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4495.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4495.pdat$Sequencer <- "Illumina NextSeq 500"
CID4495.pdat$Tissue.Source <- "Primary Tumor"
CID4495.pdat$BC.Subtype <- "TNBC"
CID4495.pdat$Status <- "Primary"
CID4495.pdat$Treatment.Status <- "Naive"

CID4495.pdat$Grade <- "3"
CID4495.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID4495.pdat$Other.Clinical <- "Medullary features"
CID4495.pdat$TNM.Classification <- "pT1c, pN0"
CID4495.pdat$Ki67 <- "0.8"
CID4495.pdat$BRCA.Status <- "NA"
CID4495.pdat$Tumor.Size <- "NA"
CID4495.pdat$Age <- "63"
CID4495.pdat$Gender <- "Female"
CID4495.pdat$Parity <- "NA"
CID4495.pdat$Gravidity <- "NA"
CID4495.pdat$BMI <- "NA"
CID4495.pdat$Menopause <- "NA"
CID4495.pdat$Ethnicity <- "NA"

# CID44971 - TNBC - has normal epithelia ----------------------------------------------------------

CID44971.mat <- readMM(file = './GSE176078//CID44971/count_matrix_sparse.mtx')
CID44971.mat <- as.matrix(CID44971.mat)
CID44971.mat[1:5,1:5]
dim(CID44971.mat)

CID44971.genes <- read.table(file = './GSE176078//CID44971/count_matrix_genes.tsv', 
                             sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID44971.genes)
rownames(CID44971.mat) <- CID44971.genes[,1]
CID44971.mat[1:5,1:5]

CID44971.bar <- read.table(file = './GSE176078//CID44971/count_matrix_barcodes.tsv', 
                           sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID44971.mat) <- CID44971.bar[,1]
CID44971.mat[1:5,1:5]


CID44971.meta <- read.csv(file = './GSE176078//CID44971/metadata.csv')
head(CID44971.meta)
CID44971.pdat <- CID44971.meta[,c(5,6,7,8,9)]

CID44971.pdat$Patient <- CID44971.meta$orig.ident
CID44971.pdat$Molecular.Subtype <- "TNBC"
CID44971.pdat$RNA.Type <- "Total RNA"
CID44971.pdat$Gene.Coverage <- "3' and 5'"
CID44971.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID44971.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID44971.pdat$Sequencer <- "Illumina NextSeq 500"
CID44971.pdat$Tissue.Source <- "Primary Tumor"
CID44971.pdat$BC.Subtype <- "TNBC"
CID44971.pdat$Status <- "Primary"
CID44971.pdat$Treatment.Status <- "Naive"

CID44971.pdat$Grade <- "3"
CID44971.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID44971.pdat$Other.Clinical <- "Highly atypical cells with circumscribed periphery, associated high grade DCIS and LVI. Accompanying lymphoid stroma"
CID44971.pdat$TNM.Classification <- "pT2, N1a, Mx"
CID44971.pdat$Ki67 <- "0.4"
CID44971.pdat$BRCA.Status <- "NA"
CID44971.pdat$Tumor.Size <- "NA"
CID44971.pdat$Age <- "49"
CID44971.pdat$Gender <- "Female"
CID44971.pdat$Parity <- "NA"
CID44971.pdat$Gravidity <- "NA"
CID44971.pdat$BMI <- "NA"
CID44971.pdat$Menopause <- "NA"
CID44971.pdat$Ethnicity <- "NA"

# CID44991 - TNBC - has normal epithelia ----------------------------------------------------------

CID44991.mat <- readMM(file = './GSE176078//CID44991/count_matrix_sparse.mtx')
CID44991.mat <- as.matrix(CID44991.mat)
CID44991.mat[1:5,1:5]
dim(CID44991.mat)

CID44991.genes <- read.table(file = './GSE176078//CID44991/count_matrix_genes.tsv', 
                             sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID44991.genes)
rownames(CID44991.mat) <- CID44991.genes[,1]
CID44991.mat[1:5,1:5]

CID44991.bar <- read.table(file = './GSE176078//CID44991/count_matrix_barcodes.tsv', 
                           sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID44991.mat) <- CID44991.bar[,1]
CID44991.mat[1:5,1:5]


CID44991.meta <- read.csv(file = './GSE176078//CID44991/metadata.csv')
head(CID44991.meta)
CID44991.pdat <- CID44991.meta[,c(5,6,7,8,9)]

CID44991.pdat$Patient <- CID44991.meta$orig.ident
CID44991.pdat$Molecular.Subtype <- "TNBC"
CID44991.pdat$RNA.Type <- "Total RNA"
CID44991.pdat$Gene.Coverage <- "3' and 5'"
CID44991.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID44991.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID44991.pdat$Sequencer <- "Illumina NextSeq 500"
CID44991.pdat$Tissue.Source <- "Primary Tumor"
CID44991.pdat$BC.Subtype <- "TNBC"
CID44991.pdat$Status <- "Primary"
CID44991.pdat$Treatment.Status <- "Naive"

CID44991.pdat$Grade <- "3"
CID44991.pdat$Histology <- "Invasive Ductal Carcinoma (IDC)"
CID44991.pdat$Other.Clinical <- "BRCA2 mutation"
CID44991.pdat$TNM.Classification <- "NA"
CID44991.pdat$Ki67 <- "60-70%"
CID44991.pdat$BRCA.Status <- "BRCA2"
CID44991.pdat$Tumor.Size <- "NA"
CID44991.pdat$Age <- "47"
CID44991.pdat$Gender <- "Female"
CID44991.pdat$Parity <- "NA"
CID44991.pdat$Gravidity <- "NA"
CID44991.pdat$BMI <- "NA"
CID44991.pdat$Menopause <- "NA"
CID44991.pdat$Ethnicity <- "NA"

# CID4513 - TNBC ----------------------------------------------------------

CID4513.mat <- readMM(file = './GSE176078//CID4513/count_matrix_sparse.mtx')
CID4513.mat <- as.matrix(CID4513.mat)
CID4513.mat[1:5,1:5]
dim(CID4513.mat)

CID4513.genes <- read.table(file = './GSE176078//CID4513/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4513.genes)
rownames(CID4513.mat) <- CID4513.genes[,1]
CID4513.mat[1:5,1:5]

CID4513.bar <- read.table(file = './GSE176078//CID4513/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4513.mat) <- CID4513.bar[,1]
CID4513.mat[1:5,1:5]


CID4513.meta <- read.csv(file = './GSE176078//CID4513/metadata.csv')
head(CID4513.meta)
CID4513.pdat <- CID4513.meta[,c(5,6,7,8,9)]

CID4513.pdat$Patient <- CID4513.meta$orig.ident
CID4513.pdat$Molecular.Subtype <- "TNBC (metaplastic)"
CID4513.pdat$RNA.Type <- "Total RNA"
CID4513.pdat$Gene.Coverage <- "3' and 5'"
CID4513.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4513.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4513.pdat$Sequencer <- "Illumina NextSeq 500"
CID4513.pdat$Tissue.Source <- "Primary Tumor"
CID4513.pdat$BC.Subtype <- "TNBC"
CID4513.pdat$Status <- "Primary"
CID4513.pdat$Treatment.Status <- "Treated (Neoadjuvant AC (4x), Paclitaxel (3x))"

CID4513.pdat$Grade <- "3"
CID4513.pdat$Histology <- "Medullary Breast Carcinoma (MDC)"
CID4513.pdat$Other.Clinical <- "Metaplastic, spindle cell carcinoma with areas of sarcomatous appearance and inflammatory infiltrate. LVI present.  RCB-II, partial pathological response to chemotherapy"
CID4513.pdat$TNM.Classification <- "pT3, pN0, Mx, Stage IIB"
CID4513.pdat$Ki67 <- "0.75"
CID4513.pdat$BRCA.Status <- "NA"
CID4513.pdat$Tumor.Size <- "NA"
CID4513.pdat$Age <- "73"
CID4513.pdat$Gender <- "Female"
CID4513.pdat$Parity <- "NA"
CID4513.pdat$Gravidity <- "NA"
CID4513.pdat$BMI <- "NA"
CID4513.pdat$Menopause <- "NA"
CID4513.pdat$Ethnicity <- "NA"

# CID4515 - TNBC - has normal epithelia ----------------------------------------------------------

CID4515.mat <- readMM(file = './GSE176078//CID4515/count_matrix_sparse.mtx')
CID4515.mat <- as.matrix(CID4515.mat)
CID4515.mat[1:5,1:5]
dim(CID4515.mat)

CID4515.genes <- read.table(file = './GSE176078//CID4515/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4515.genes)
rownames(CID4515.mat) <- CID4515.genes[,1]
CID4515.mat[1:5,1:5]

CID4515.bar <- read.table(file = './GSE176078//CID4515/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4515.mat) <- CID4515.bar[,1]
CID4515.mat[1:5,1:5]


CID4515.meta <- read.csv(file = './GSE176078//CID4515/metadata.csv')
head(CID4515.meta)
CID4515.pdat <- CID4515.meta[,c(5,6,7,8,9)]

CID4515.pdat$Patient <- CID4515.meta$orig.ident
CID4515.pdat$Molecular.Subtype <- "TNBC"
CID4515.pdat$RNA.Type <- "Total RNA"
CID4515.pdat$Gene.Coverage <- "3' and 5'"
CID4515.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4515.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4515.pdat$Sequencer <- "Illumina NextSeq 500"
CID4515.pdat$Tissue.Source <- "Primary Tumor"
CID4515.pdat$BC.Subtype <- "TNBC"
CID4515.pdat$Status <- "Primary"
CID4515.pdat$Treatment.Status <- "Naive"

CID4515.pdat$Grade <- "3"
CID4515.pdat$Histology <- "Ivasive Ductal Carcinoma (IDC)"
CID4515.pdat$Other.Clinical <- "Basal phenotype: CK5/6+ focal 40%, CK14+ focal 30%. Associated high grade DCIS and patchy lymphoid infiltrate"
CID4515.pdat$TNM.Classification <- "PpT1c, pN1, Mi, Stage IIA"
CID4515.pdat$Ki67 <- "0.6"
CID4515.pdat$BRCA.Status <- "NA"
CID4515.pdat$Tumor.Size <- "NA"
CID4515.pdat$Age <- "67"
CID4515.pdat$Gender <- "Female"
CID4515.pdat$Parity <- "NA"
CID4515.pdat$Gravidity <- "NA"
CID4515.pdat$BMI <- "NA"
CID4515.pdat$Menopause <- "NA"
CID4515.pdat$Ethnicity <- "NA"

# CID4523 - TNBC ----------------------------------------------------------

CID4523.mat <- readMM(file = './GSE176078//CID4523/count_matrix_sparse.mtx')
CID4523.mat <- as.matrix(CID4523.mat)
CID4523.mat[1:5,1:5]
dim(CID4523.mat)

CID4523.genes <- read.table(file = './GSE176078//CID4523/count_matrix_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dim(CID4523.genes)
rownames(CID4523.mat) <- CID4523.genes[,1]
CID4523.mat[1:5,1:5]

CID4523.bar <- read.table(file = './GSE176078//CID4523/count_matrix_barcodes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
#dim(barcodes2032)

colnames(CID4523.mat) <- CID4523.bar[,1]
CID4523.mat[1:5,1:5]


CID4523.meta <- read.csv(file = './GSE176078//CID4523/metadata.csv')
head(CID4523.meta)
CID4523.pdat <- CID4523.meta[,c(5,6,7,8,9)]

CID4523.pdat$Patient <- CID4523.meta$orig.ident
CID4523.pdat$Molecular.Subtype <- "TNBC (metaplastic)"
CID4523.pdat$RNA.Type <- "Total RNA"
CID4523.pdat$Gene.Coverage <- "3' and 5'"
CID4523.pdat$Library.Preparation <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4523.pdat$Capture.Method <- "10X Genomics Chromium Single-Cell v2 3' and 5' Chemistry Library"
CID4523.pdat$Sequencer <- "Illumina NextSeq 500"
CID4523.pdat$Tissue.Source <- "Primary Tumor"
CID4523.pdat$BC.Subtype <- "TNBC"
CID4523.pdat$Status <- "Primary"
CID4523.pdat$Treatment.Status <- "Treated ( Neoadjuvant AC (4x), Paclitaxel (1x))"


CID4523.pdat$Grade <- "3"
CID4523.pdat$Histology <- "Medullary Breast Carcinoma (MDC)"
CID4523.pdat$Other.Clinical <- "Metaplastic carcinoma with sebaceous differentiation. LVI present. RCB-II, partial pathological response to chemotherapy"
CID4523.pdat$TNM.Classification <- "pT2, pN0 (i+), pM0, Stage IIA"
CID4523.pdat$Ki67 <- "0.9"
CID4523.pdat$BRCA.Status <- "NA"
CID4523.pdat$Tumor.Size <- "NA"
CID4523.pdat$Age <- "52"
CID4523.pdat$Gender <- "Female"
CID4523.pdat$Parity <- "NA"
CID4523.pdat$Gravidity <- "NA"
CID4523.pdat$BMI <- "NA"
CID4523.pdat$Menopause <- "NA"
CID4523.pdat$Ethnicity <- "NA"



# Joined Matrices ----------------------------------------------------------

joined <- cbind(CID3586.mat, CID3838.mat, CID3921.mat, CID3941.mat, CID3946.mat,
                CID3948.mat, CID3963.mat, CID4040.mat, CID4066.mat, CID4067.mat,
                CID4290A.mat, CID4398.mat, CID44041.mat, CID4461.mat, CID4463.mat,
                CID4465.mat, CID4471.mat, CID4495.mat, CID44971.mat, CID44991.mat,
                CID4513.mat, CID4515.mat, CID45171.mat, CID4523.mat, CID4530N.mat,
                CID4535.mat)
dim(joined)

joined <- as.data.frame(joined)

library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(joined)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- joined[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- joined[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}


joined <- new_mat

#pdata =================================

pdat_joined <- rbind(CID3586.pdat, CID3838.pdat, CID3921.pdat, CID3941.pdat, CID3946.pdat,
                     CID3948.pdat, CID3963.pdat, CID4040.pdat, CID4066.pdat, CID4067.pdat,
                     CID4290A.pdat, CID4398.pdat, CID44041.pdat, CID4461.pdat, CID4463.pdat,
                     CID4465.pdat, CID4471.pdat, CID4495.pdat, CID44971.pdat, CID44991.pdat,
                     CID4513.pdat, CID4515.pdat, CID45171.pdat, CID4523.pdat, CID4530N.pdat,
                     CID4535.pdat)
length(unique(joined[,1]))==length(joined[,1])
length(unique(pdat_joined[,1]))==length(pdat_joined[,1])



#_________________________________________________________________
rownames(pdat_joined) <- pdat_joined$samples
fdat_joined <- toupper(as.matrix(rownames(joined)))
#_________________________________________________________________
rownames(fdat_joined) <- fdat_joined[,1]
fdat_joined <- data.frame(fdat_joined)
common_colnames <- c("gene_short_name")
colnames(fdat_joined) <- common_colnames
rownames(joined) <- rownames(fdat_joined)
#_________________________________________________________________
sobj_pre <- CreateSeuratObject(counts = joined)
rownames(pdat_joined) <- colnames(joined)
sobj_pre <-AddMetaData(sobj_pre,metadata=pdat_joined)
head(sobj_pre@meta.data)
#_________________________________________________________________
sobj_pre[["RNA"]]@meta.features<-fdat_joined
head(sobj_pre[["RNA"]]@meta.features)
slotNames(sobj_pre[["RNA"]])

sobj <- sobj_pre

sobj_prim <- subset(x = sobj, subset = celltype_major == "Normal Epithelial", invert = TRUE)

sobj_prim$samples <- sobj_prim$Patient

saveRDS(sobj_prim, file = "Wu2021_primonly_unfiltered_62622.rds")

sobj <- sobj_prim
#cell cycle ========================================================================
#from: https://satijalab.org/seurat/articles/cell_cycle_vignette.html

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes)
sobj$CC.Difference <- sobj$S.Score - sobj$G2M.Score
# view cell cycle scores and phase assignments
head(sobj[[]])

# Preprocessing and Filtering (SCTransform) --------------------------------------------------------------------


mito.genes <- grep(pattern = "^MT-", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(mito.genes)
sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")

#blood cell contamination
hemo.genes <- grep(pattern = "^HB[^(P)]", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(hemo.genes)
sobj <- PercentageFeatureSet(sobj, "^HB[^(P)]", col.name = "percent.hb")
platelet.genes <- grep(pattern = "^PECAM1|^PF4\\b", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(platelet.genes)
sobj <- PercentageFeatureSet(sobj, "^PECAM1|^PF4\\b", col.name = "percent.platelet")

nCount95 <- quantile(sobj@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
nCount5 <- quantile(sobj@meta.data$nCount_RNA, 0.05)#05)  
nFeat95 <- quantile(sobj@meta.data$nFeature_RNA, 0.95) #0.9) 
nFeat5 <- quantile(sobj@meta.data$nFeature_RNA, 0.05) #0.1) 
ptMt5 <- quantile(sobj@meta.data$percent.mt, 0.05)
ptMt95 <- quantile(sobj@meta.data$percent.mt, 0.95)
ptHb95 <- quantile(sobj@meta.data$percent.hb, 0.95)

sobj <- subset(x = sobj, subset = nFeature_RNA > nFeat5 & nFeature_RNA < nFeat95)
sobj <- subset(x = sobj, subset = percent.mt > ptMt5 & percent.mt < ptMt95)
sobj <- subset(x = sobj, subset = nCount_RNA > nCount5 & nCount_RNA < nCount95)
sobj <- subset(x = sobj, subset = percent.hb < ptHb95)

summary(sobj@meta.data$percent.mt)
summary(sobj@meta.data$percent.hb)
dim(sobj)
dim(sobj_prim)

DefaultAssay(sobj) <- "RNA"
saveRDS(sobj, file = "NewWu_PRIM_filtered_62622.rds")





# =================================================================== ======

# Xu ================================================================ ============
# =================================================================== ======








# Loading Raw Data into RStudio ---------------------------------- 

filePaths = getGEOSuppFiles("GSE180286") 
tarF <- list.files(path = "./GSE180286/", pattern = "*.tar", full.names = TRUE) 
ldply(.data = tarF, .fun = untar, exdir = "./GSE180286/")
gzipF <- list.files(path = "./GSE180286/", pattern = "*.gz", full.names = TRUE) 
ldply(.data = gzipF, .fun = gunzip) 

list.files(path = "./GSE180286/", full.names = TRUE)

list.files(path = "./GSE180286/", pattern = "\\.txt$",full.names = TRUE)


# P1_prim ---------------

P1prim.mat <- read.table(file = './GSE180286//GSM5457199_A2019-1.expression_matrix.txt', 
                         sep = '\t', header = TRUE, stringsAsFactors = FALSE)
dim(P1prim.mat)
P1prim.mat[1:5,1:5]

colnames(P1prim.mat) <- paste(colnames(P1prim.mat), "P1_Prim", sep = "_")
P1prim.pdat <- data.frame("cells" = colnames(P1prim.mat), "samples" = "P1prim","Patient" = "P1", "BC.Subtype" = "TNBC",
                          "Molecular Subtype" = "TNBC", "RNA Type" = "mRNA", "Gene Coverage" = "NA", 
                          "Library Preparation" = "Singleron GEXSCOPE Single Cell RNAseq Library Kit", 
                          "Capture Method" = "Singleron GEXSCOPE Single Cell RNAseq Library Kit", "Sequencer" = "Illumina NovaSeq", 
                          "Treatment.Status" = "Naive", "Gender" = "Female", "Age" = "NA", "Ki67" = "NA",
                          "TNM.Classification" = "T2N1M0", "Status" = "Primary")

head(P1prim.pdat)
rownames(P1prim.pdat) <- P1prim.pdat$cells

P1prim.mat <- as.data.frame(P1prim.mat)

library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(P1prim.mat)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- P1prim.mat[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- P1prim.mat[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}


P1prim.mat <- new_mat

P1prim <- CreateSeuratObject(P1prim.mat, meta.data = P1prim.pdat)


P1prim <- CreateSeuratObject(P1prim.mat, meta.data = P1prim.pdat)
head(P1prim@meta.data)
# P2_prim ---------------

P2prim.mat <- read.table(file = './GSE180286//GSM5457202_B2019-1.expression_matrix.txt', 
                         sep = '\t', header = TRUE, stringsAsFactors = FALSE)
dim(P2prim.mat)
P2prim.mat[1:5,1:5]

P2prim.mat <- as.data.frame(P2prim.mat)

library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(P2prim.mat)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- P2prim.mat[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- P2prim.mat[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}


P2prim.mat <- new_mat
colnames(P2prim.mat) <- paste(colnames(P2prim.mat), "P2_Prim", sep = "_")
P2prim.pdat <- data.frame("cells" = colnames(P2prim.mat), "samples" = "P2prim","Patient" = "P2", "BC.Subtype" = "HR+",
                          "Molecular Subtype" = "LumB", "RNA Type" = "mRNA", "Gene Coverage" = "NA", 
                          "Library Preparation" = "Singleron GEXSCOPE Single Cell RNAseq Library Kit", 
                          "Capture Method" = "Singleron GEXSCOPE Single Cell RNAseq Library Kit", "Sequencer" = "Illumina NovaSeq", 
                          "Treatment.Status" = "Naive", "Gender" = "Female", "Age" = "54", "Ki67" = "40%",
                          "TNM.Classification" = "T3N2M0", "Status" = "Primary")
rownames(P2prim.pdat) <- P2prim.pdat$cells



P2prim <- CreateSeuratObject(P2prim.mat, meta.data = P2prim.pdat)
# P3_prim ---------------

P3prim.mat <- read.table(file = './GSE180286//GSM5457205_C2020-1.expression_matrix.txt', 
                         sep = '\t', header = TRUE, stringsAsFactors = FALSE)
dim(P3prim.mat)
P3prim.mat[1:5,1:5]


P3prim.mat <- as.data.frame(P3prim.mat)


library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(P3prim.mat)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- P3prim.mat[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- P3prim.mat[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}


P3prim.mat <- new_mat

colnames(P3prim.mat) <- paste(colnames(P3prim.mat), "P3_Prim", sep = "_")
P3prim.pdat <- data.frame("cells" = colnames(P3prim.mat), "samples" = "P3prim","Patient" = "P3", "BC.Subtype" = "HER2+",
                          "Molecular Subtype" = "Her2", "RNA Type" = "mRNA", "Gene Coverage" = "NA", 
                          "Library Preparation" = "Singleron GEXSCOPE Single Cell RNAseq Library Kit", 
                          "Capture Method" = "Singleron GEXSCOPE Single Cell RNAseq Library Kit", "Sequencer" = "Illumina NovaSeq", 
                          "Treatment.Status" = "Naive", "Gender" = "Female", "Age" = "NA", "Ki67" = "NA",
                          "TNM.Classification" = "T1N1M0", "Status" = "Primary")

rownames(P3prim.pdat) <- P3prim.pdat$cells

P3prim <- CreateSeuratObject(P3prim.mat, meta.data = P3prim.pdat)

# P4_prim ---------------

P4prim.mat <- read.table(file = './GSE180286//GSM5457208_D2020-1.expression_matrix.txt', 
                         sep = '\t', header = TRUE, stringsAsFactors = FALSE)
dim(P4prim.mat)
P4prim.mat[1:5,1:5]

P4prim.mat <- as.data.frame(P4prim.mat)

library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(P4prim.mat)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- P4prim.mat[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- P4prim.mat[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}


P4prim.mat <- new_mat

colnames(P4prim.mat) <- paste(colnames(P4prim.mat), "P4_Prim", sep = "_")
P4prim.pdat <- data.frame("cells" = colnames(P4prim.mat), "samples" = "P4prim","Patient" = "P4", "BC.Subtype" = "TNBC",
                          "Molecular Subtype" = "TNBC", "RNA Type" = "mRNA", "Gene Coverage" = "NA", 
                          "Library Preparation" = "Singleron GEXSCOPE Single Cell RNAseq Library Kit", 
                          "Capture Method" = "Singleron GEXSCOPE Single Cell RNAseq Library Kit", "Sequencer" = "Illumina NovaSeq", 
                          "Treatment.Status" = "Naive", "Gender" = "Female", "Age" = "NA", "Ki67" = "NA",
                          "TNM.Classification" = "T2N2M0", "Status" = "Primary")

rownames(P4prim.pdat) <- P4prim.pdat$cells


P4prim <- CreateSeuratObject(P4prim.mat, meta.data = P4prim.pdat)

# P5_prim ---------------

P5prim.mat <- read.table(file = './GSE180286//GSM5457211_E2020-1.expression_matrix.txt', 
                         sep = '\t', header = TRUE, stringsAsFactors = FALSE)
dim(P5prim.mat)
P5prim.mat[1:5,1:5]

P5prim.mat <- as.data.frame(P5prim.mat)

library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(P5prim.mat)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)

new_mat <- P5prim.mat[!duplicated(usegenes),]
rownames(new_mat) <- unique(usegenes)
new_mat <- new_mat[which(!(rownames(new_mat) %in%
                             usegenes[duplicated(usegenes)])),]

for (i in unique(usegenes[duplicated(usegenes)])) {
  j <- P5prim.mat[which(usegenes == i),]
  j <- apply(j, 2, FUN = max)
  k <- data.frame(t(j))
  rownames(k) <- i
  colnames(k) <- colnames(new_mat)
  new_mat <- rbind(new_mat,k)
}


P5prim.mat <- new_mat

colnames(P5prim.mat) <- paste(colnames(P5prim.mat), "P5_Prim", sep = "_")
P5prim.pdat <- data.frame("cells" = colnames(P5prim.mat), "samples" = "P5prim","Patient" = "P5", "BC.Subtype" = "HER2+",
                          "Molecular Subtype" = "Her2", "RNA Type" = "mRNA", "Gene Coverage" = "NA", 
                          "Library Preparation" = "Singleron GEXSCOPE Single Cell RNAseq Library Kit", 
                          "Capture Method" = "Singleron GEXSCOPE Single Cell RNAseq Library Kit", "Sequencer" = "Illumina NovaSeq", 
                          "Treatment.Status" = "Naive", "Gender" = "Female", "Age" = "NA", "Ki67" = "NA",
                          "TNM.Classification" = "T2N1M0", "Status" = "Primary")

rownames(P5prim.pdat) <- P5prim.pdat$cells

P5prim <- CreateSeuratObject(P5prim.mat, meta.data = P5prim.pdat)


# Joined Matrices ----------------------------------------------------------

sobj_prim <- merge(x = P1prim, y = c(P2prim, P3prim, P4prim, P5prim),  
                   add.cell.ids = c("P1", "P2", "P3", "P4", "P5"))

saveRDS(sobj_prim, file = "Xu_PRIM_unfiltered_62622.rds")


sobj <- sobj_prim


#cell cycle ========================================================================
#from: https://satijalab.org/seurat/articles/cell_cycle_vignette.html

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes)
sobj$CC.Difference <- sobj$S.Score - sobj$G2M.Score
# view cell cycle scores and phase assignments
head(sobj[[]])

# Preprocessing and Filtering (SCTransform) --------------------------------------------------------------------


mito.genes <- grep(pattern = "^MT-", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(mito.genes)
sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")

#blood cell contamination
hemo.genes <- grep(pattern = "^HB[^(P)]", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(hemo.genes)
sobj <- PercentageFeatureSet(sobj, "^HB[^(P)]", col.name = "percent.hb")
platelet.genes <- grep(pattern = "^PECAM1|^PF4\\b", x = rownames(sobj@assays[["RNA"]]), value = TRUE)
print(platelet.genes)
sobj <- PercentageFeatureSet(sobj, "^PECAM1|^PF4\\b", col.name = "percent.platelet")

nCount95 <- quantile(sobj@meta.data$nCount_RNA, 0.95)#95) # calculate value in the 95th percentile)
nCount5 <- quantile(sobj@meta.data$nCount_RNA, 0.05)#05)  
nFeat95 <- quantile(sobj@meta.data$nFeature_RNA, 0.95) #0.9) 
nFeat5 <- quantile(sobj@meta.data$nFeature_RNA, 0.05) #0.1) 
ptMt5 <- quantile(sobj@meta.data$percent.mt, 0.05)
ptMt95 <- quantile(sobj@meta.data$percent.mt, 0.80) #80 for prim
ptHb95 <- quantile(sobj@meta.data$percent.hb, 0.95)

sobj <- subset(x = sobj, subset = nFeature_RNA > nFeat5 & nFeature_RNA < nFeat95)
sobj <- subset(x = sobj, subset = percent.mt > ptMt5 & percent.mt < ptMt95)
sobj <- subset(x = sobj, subset = nCount_RNA > nCount5 & nCount_RNA < nCount95)
sobj <- subset(x = sobj, subset = percent.hb < ptHb95)

summary(sobj@meta.data$percent.mt)
summary(sobj@meta.data$percent.hb)
dim(sobj)
dim(sobj_prim)

saveRDS(sobj, file = "Xu_PRIM_filtered_62622.rds")



