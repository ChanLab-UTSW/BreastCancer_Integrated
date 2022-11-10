



#dataset: https://www.cancerrxgene.org/downloads/bulk_download 

#NOTE!!!! 
## https://github.com/CancerRxGene/gdscIC50/tree/master/vignettes
## ^ shows description of raw data files, as well as code for processing it

##python package for analysis: https://gdsctools.readthedocs.io/en/master/index.html

#libraries =====
library(RCurl)


#download files (NOTE, NOT SURE THIS WORKS, BUT GIVES YOU IDEA OF WHAT EACH FILE IS) USE INDIVIDUAL LOADING CODE TO GET ACTUAL FILE =====

#https://felixfan.github.io/download-files/ <- actually worked
#https://www.youtube.com/watch?v=EBfx1L16qlM

save_path = "/project/InternalMedicine/Chan_lab/shared/IntegratedBreast_s204665/ReviewerAdditions/Drug_Celllines/drug_orig_dataset_files/drug_screen_files/"

dataset_dir = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.4/"
setwd(save_path)

filenames = getURL(dataset_dir, ftp.use.epsv = FALSE, dirlistonly = TRUE)
filenames <- strsplit(filenames, "\n")
filenames = unlist(filenames)
filenames

#filenames.working <- filenames[c(7,9:10,12:13)]

for (ind_file in filenames) #filenames.working) {
{
  
  # csv <- getURL(paste0(dataset_dir, ind_file), ssl.verifypeer = FALSE)
  # test <- read.csv(textConnection(csv))
  # head(test)
  # 
  # write.csv(test, ind_file)
  
  download.file(paste(dataset_dir, ind_file, sep = ""), paste(save_path, "/", ind_file,
                                                      sep = ""))
}


#ind file (USE THIS CODE TO GET ACTUAL .csv FILES) ----------

filenames #<- gets you individual file name list

url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.4/GDSC1_fitted_dose_response_24Jul22.csv"
url = "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.4/GDSC2_fitted_dose_response_24Jul22.csv"

csv <- getURL(url, ssl.verifypeer = FALSE)
test <- read.csv(textConnection(csv))
head(test)

#write.csv(test, "USETHIS_actual_GDSC1_fitted.csv")
write.csv(test, "USETHIS_actual_GDSC2_fitted.csv")

test.2 <- read.csv("USETHIS_actual_GDSC1_fitted.csv")
#test.2 <- read.csv("./GDSCtools_mobems/GF_BRCA_mobem.csv")

unique(sort(test.2$TCGA_DESC))

#preview files =====


test.2 <- read.csv("GDSC1_fitted_dose_response_24Jul22.csv")
test.2 <- read.csv("./GDSCtools_mobems/GF_BRCA_mobem.csv")

