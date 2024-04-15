

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






# Loading Raw Data into RStudio (mRNA) ---------------------------------- 

#tar -xvf GSE194040_RAW.tar
#mv *.gz GSE194040

list.files(path = "./GSE194040/", full.names = TRUE)
mafiles <- list.files(path = "./GSE194040/", pattern = "\\+.txt$",full.names = TRUE)


#read in microarray data -------------------------------------------------

test.sample <- read.table(file = './GSE194040//GSM5826887_ISPY2_999733_GPL20078_SingleChannel_FullGenome+.txt')

mafiles <- list.files(path = "./GSE194040/", pattern = "\\+.txt$",full.names = TRUE)

library(limma)
x <- read.maimages(mafiles, source = "agilent", green.only = T, 
                   other.columns = c("ProbeName",	"gMeanSignal",	
                                     "gBGMeanSignal",	"gMedianSignal",	
                                     "gBGMedianSignal",	"gNormalizedSignal (log2)",
                                     "gIsFeatNonUnifOL"))

#structure of class object: https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/EList.html
str(x)
head(x$genes) #probe names - need to convert to gene names
colnames(x$E)<-gsub("^.*\\ISPY2_","",colnames(x$E))
colnames(x$E) <- gsub("_.*", "", colnames(x$E))
colnames(x$E) <- paste("X", colnames(x$E), sep = "")
head(x$E) # is values from gMedianSignal column

colnames(x$Eb)<-gsub("^.*\\ISPY2_","",colnames(x$Eb))
colnames(x$Eb) <- gsub("_.*", "", colnames(x$Eb))
colnames(x$Eb) <- paste("X", colnames(x$Eb), sep = "")
head(x$Eb) #is values from gBGMedianSignal

head(x$targets) #patient/sample names
names(x$other) #all values from original sample files
class(x$other$gMedianSignal) #all values from original sample files
class(x)

corrected.values <- read.table(file = './GSE194040/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt')
head(corrected.values)

x$other$corrected <- corrected.values
head(x$other$corrected)
head(x$E)

# Add gene annotation ==================

probe.gene.annot <- read.table(file = './GSE194040/GSE194040_ProbeAnnotation_ISPY2Edit_GPL30493_wasGPL16233.txt', 
                               sep = '\t', header = T, stringsAsFactors = FALSE)
dim(probe.gene.annot)
probe.gene.annot[1:5,1:4]
probes <- x$genes
head(probes)
probegenedat <- merge(probes, probe.gene.annot, by ="ProbeName")
head(probegenedat)

x$genes <- probegenedat
head(x$genes) 


# Add metadata ==================

#download tableS2
meta <- read.csv("mmc3.csv", header = T)
meta$Patient.Identifier <- paste("X", meta$Patient.Identifier, sep = "")
head(meta)

patientID <- x$targets
patientID$Patient.Identifier <- patientID$FileName
rownames(patientID) <- 1:nrow(patientID)
patientID$Patient.Identifier<-gsub("^.*\\ISPY2_","",patientID$Patient.Identifier)
patientID$Patient.Identifier <- gsub("_.*", "", patientID$Patient.Identifier)
patientID$Patient.Identifier <- paste("X", patientID$Patient.Identifier, sep = "")
head(patientID)

patientID.withmeta <- merge(patientID, meta, by = "Patient.Identifier")
head(patientID.withmeta)

x$targets <- patientID.withmeta
head(x$targets)

# processed version ======================

#just to get the Elist object version 
#already normalized it so I'll just replace it
y <- backgroundCorrect(x, method = "normexp")
y <- normalizeBetweenArrays(y, method = "quantile")
class(y)

head(y$E)
y$other$bgCorrectNormBtn <- y$E 
y$E <- y$other$corrected #just put their correction in the E slot

head(y$E) #their corrected values
head(y$targets) #patient metadata
head(y$genes) #gene metadata
head(y$source) ##from agilent
names(y$other) #includes RAW SAMPLE VALUES in columns "ProbeName" throught "gIsFeatNonUnifOL"

rownames(y$E)
test <- y$E["ATP5C1",]
test[1:7,1:5]

which(y$genes == "HIST1H2BI", arr.ind=TRUE)



saveRDS(y, file = "iSPY_microarrayobject.rds")


##Extra =========================================== ================
## ================================================ ========

##convert gene names =======================================================

y <- readRDS("iSPY_microarrayobject.rds")


library(limma)
library(org.Hs.eg.db)
oldgenes <- rownames(y$E)
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes<-ifelse(is.na(newgenes),oldgenes,newgenes)
head(oldgenes)
head(newgenes)
head(usegenes)

which(usegenes == "ATP5F1C")

joined <- y$E

investigate <- usegenes[duplicated(usegenes)]
investigate

test <- joined[usegenes == "H2BC5",]
test[1:8,1:5]
test_trans <- t(test)
test_trans <- as.data.frame(test_trans)
test_trans[1:5,1:2]
HIST1H2BI <- summary(test_trans$C2orf27A)
HIST1H2BI <- summary(test_trans$C2orf27B)


#fit =====================================

fit <- read.table(file = './GSE194040/GSE194040_ISPY2_LinearFit_forCombiningExp_GPL30493_wasGPL16233.txt', 
                  sep = '\t', header = T, stringsAsFactors = FALSE)
head(fit)

topTable(fit, coef = 4, n =20)

fit2 <- read.table(file = './GSE194040/GSE194040_ISPY2_LinearFit_forCombiningExp_20078.txt', 
                   sep = '\t', header = T, stringsAsFactors = FALSE)

head(fit2)
