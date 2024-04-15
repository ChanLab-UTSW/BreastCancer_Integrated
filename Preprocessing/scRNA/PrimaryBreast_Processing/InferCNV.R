# SETUP -----------------------------------------------------------------------
# Libraries -------------------------------------------------------------------

library(BiocManager)
library(plyr)
library(dplyr) 
library(Matrix)
library(Seurat)
library(ggplot2)
library(Rcpp)
library(RcppZiggurat)
library(Rfast)

#BiocManager::install("infercnv")
library(infercnv)

# ---------
# ---------
# Load in Seurat objects -------------------------------------------------------

setwd("/project/InternalMedicine/Chan_lab/shared/FinalObjects/Primary_Only/")
combo.reference <- readRDS("PrimObject_FINALnoZgenedem_71222.rds") 
DefaultAssay(combo.reference) <- "RNA"

cancer.epi <- readRDS("cancerepi_withGEmetadata_111422.rds")
DefaultAssay(cancer.epi) <- "RNA"

# ---------
# ---------
# Prepare functions ---------------------------------

prepare_infercnv_metadata <- function(seurat_object, subset_data = FALSE, count_df, for_infercnv=TRUE) {
  temp_metadata <- data.frame(cell_ids = names(Idents(seurat_object)), 
                              cell_type = Idents(seurat_object), 
                              stringsAsFactors = F)
  temp_metadata <- temp_metadata[order(temp_metadata$cell_type),]
  
  # remove unassigned cells and CAFs (as they may be malignant):
  temp_metadata <- temp_metadata[grep("[u,U]nknown|[u,U]nassigned|Fibroblasts", 
                                      temp_metadata$cell_type, invert=T),]
  
  # remove T cells
  temp_metadata <- temp_metadata[grep("[t,T][-_][c,C]ell", temp_metadata$cell_type, invert=T),]
    
  # only include cells present in count_df:
  print(paste0("No cells in metadata df before filtering for those in count df = ", 
               nrow(temp_metadata)))
  temp_metadata <- temp_metadata[temp_metadata$cell_ids %in% colnames(count_df),]
  
  print(paste0("No cells in metadata df after filtering for those in count df = ", 
               nrow(temp_metadata)))
    
  if (for_infercnv) {
    # label cells in clusters of < 2 cells as 'outliers' as these will break InferCNV:
    temp_metadata$cell_type <- as.character(temp_metadata$cell_type)
    cluster_list <- split(temp_metadata, temp_metadata$cell_type)
    metadata_outliers_labelled <- do.call(
      "rbind", lapply(cluster_list, function(x) {
        if (nrow(x) < 2) {
          x$cell_type <- gsub("_[0-9].*$", "_outlier", x$cell_type)
        }
        return(x)
      })
    )
    
    # remove outlier cell types with <2 cells:
    cluster_list2 <- split(metadata_outliers_labelled, 
                           metadata_outliers_labelled$cell_type)
    i=1
    metadata_final <- do.call(
      "rbind", lapply(cluster_list2, function(x) {
        if (nrow(x) < 2) {
          print(paste0(cluster_list[[i]]$cell_type[1], 
                       " removed as contained <2 cells"))
          i <<- i+1
          return(NULL)
        } else {
          i <<- i+1
          return(x)
        }
      })
    )
  } else {
    metadata_final <- temp_metadata
  }
  rownames(metadata_final) <- metadata_final$cell_ids 
  
  # record number per cell type:
  number_per_cell_type <- as.data.frame(table(metadata_final$cell_type))
  
  # create and label results list:
  result_list <- list(metadata_final, number_per_cell_type, seurat_object)
  names(result_list) <- c("metadata", "number_per_group", "seurat")
  
  return(result_list)
} 

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

# ---------
# ---------
# Create gene order file -------

# start from Trinity CTAT/CVN files
gene_order <- read.table("gencode_v21_gen_pos.complete.txt")
gene_order$V1 <- purrr::map_df(strsplit(gene_order$V1, "\\|"), ~as.data.frame(t(.)))[,1]

oldgenes <- gene_order$V1
newgenes <- alias2SymbolTable(oldgenes, species = "Hs")
usegenes <- ifelse(is.na(newgenes),oldgenes,newgenes)
gene_order$V1 <- usegenes

gene_order$V1 <- gsub("\\.", "-", gene_order$V1)
gene_order$V1 <- toupper(gene_order$V1)

length(rownames(seurat_10X)[which(!rownames(seurat_10X) %in% gene_order$V1)])

gene_order_nodup <- gene_order[-which(duplicated(gene_order$V1)),]

write.table(gene_order, "infercnv_gene_order.txt", sep="\t", col.names = F, row.names = F)
write.table(gene_order_nodup, "infercnv_gene_order_nodup.txt", sep="\t", col.names = F, row.names = F)

# ---------
# ---------
# Run inferCNV by dataset ------- 

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/inferCNV")

for (i in unique(combo.reference$orig.ident)[6]) {
  # generate input matrix and metadata files -------
  
  # load seurat object:
  seurat_10X <- subset(combo.reference, subset = orig.ident == i)
  i <- paste0(i, "_2")
  seurat_10X <- subset(seurat_10X, subset = samples %in% unique(seurat_10X$samples)[8:14])
  
  seurat_10X <- nochemo_pre
  Idents(seurat_10X) <- paste0(seurat_10X$celltype_final, "_", seurat_10X$samples)
  i <- "Bassez"
  
  # create raw matrix input file and subset if necessary:
  count_df <- GetAssayData(seurat_10X , slot = "counts")
  print(
    paste0(
      "Dimensions of count df = ", paste(as.character(dim(count_df)), collapse=",")
    )
  )
  
  # create metadata df:
  infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, 
                                                 subset_data=subset_data, 
                                                 count_df, for_infercnv=T)
  seurat_10X <- infercnv_metadata$seurat
  print(paste0("Cell types are: ", unique(infercnv_metadata$metadata$cell_type)))
  
  # only keep cells in metadata df:
  print(paste0("No cells in count df before filtering for those in metadata df = ", 
               ncol(count_df)))
  count_df <- count_df[,colnames(count_df) %in% infercnv_metadata$metadata$cell_ids]
  print(paste0("No cells in count df after filtering for those in metadata df = ", 
               ncol(count_df)))
  
  # generate cluster metric plots for epithelial cluster:
  epithelial_clusters <- grep("Epithelial", unique(infercnv_metadata$metadata$cell_type), value=T)
  print(paste0("Epithelial cluster = ", epithelial_clusters))
  
  png(paste0(i, "_metrics_by_epithelial_cluster.png"),
      width=14, height=8, res=300, units='in')
  temp_violinplot <- VlnPlot(
    object = seurat_10X,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    pt.size = 1.5,
    idents = epithelial_clusters)
  print(temp_violinplot)
  dev.off()
  
  # remove CAFs and T cells  from analysis
  infercnv_metadata$metadata <- infercnv_metadata$metadata[
    grep("Fibroblasts", infercnv_metadata$metadata$cell_type, invert=T),
    ]
  
  infercnv_metadata$metadata <- infercnv_metadata$metadata[
    grep("T Cells", infercnv_metadata$metadata$cell_type, invert=T),
    ]
  
  # collapse all stromal cells into 'stromal' cell type:
  infercnv_metadata$metadata$cell_type[
    grep("Epithelial", infercnv_metadata$metadata$cell_type, invert=T)
    ] <- paste0("Stromal_", infercnv_metadata$metadata$cell_type[
      grep("Epithelial", infercnv_metadata$metadata$cell_type, invert=T)])
  
  infercnv_metadata$metadata$cell_type[
    grep("Epithelial", infercnv_metadata$metadata$cell_type, invert=T)
    ] <- gsub("^([^_]*_[^_]*).*", "\\1", infercnv_metadata$metadata$cell_type[grep("Epithelial", infercnv_metadata$metadata$cell_type, invert=T)])

  # write all input files for inferCNV ------
  
  # if no epithelial clusters present, abort:
  if (length(epithelial_clusters) < 1) {
    print(paste0("No epithelial/myoepithelial clusters detected, aborting..."))
  } else {
    # write count, metadata files and new seurat object:
    print("Creating inferCNV raw counts file...")
    
    count_df <- as_matrix(count_df)
    write.table(count_df, paste0(i, "_input_matrix.txt"), 
                quote=F, sep="\t", col.names=T, row.names=T)
    
    write.table(infercnv_metadata$metadata, paste0(i, "_metadata.txt"), 
                quote=F, sep="\t", col.names=F, row.names=F)
    
    write.table(infercnv_metadata$number_per_group, paste0(i, "_number_per_group.txt"), 
                quote=F, col.names=F, row.names=F, sep="\t")
    
    #saveRDS(seurat_10X, "seurat_object_annotated.Rdata")
  }
    
  # define normals for inferCNV reference (stromal cells, no fibroblasts) ------
  
  normals <- grep(
    "Epithelial", 
    unique(infercnv_metadata$metadata$cell_type[
      infercnv_metadata$metadata$cell_ids %in% colnames(count_df)
      ]), value=T, 
    invert=T
  )
  
  # run inferCNV -----
  
  if (length(epithelial_clusters) >= 1) {
    print(paste0("Normal is: ", normals))
    
    print("Creating inferCNV object...")
    raw_path <- paste0(i, "_input_matrix.txt")
    annotation_path <- paste0(i, "_metadata.txt")
    gene_path <- "infercnv_gene_order_nodup.txt"
    
    initial_infercnv_object <- CreateInfercnvObject(
      raw_counts_matrix=raw_path,
      annotations_file=annotation_path,
      delim="\t",
      gene_order_file=gene_path,
      ref_group_names=normals
    )
      
    library(parallel)
    numcores <- detectCores()
    
    new_gene_order = data.frame()
    for (chr_name in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                       "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                       "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) {
      new_gene_order = rbind(new_gene_order, 
                             initial_infercnv_object@gene_order[which(initial_infercnv_object@gene_order[["chr"]]
                                                              == chr_name) , , drop=FALSE])
    }
    
    names(new_gene_order) <- c("chr", "start", "stop")
    initial_infercnv_object@gene_order = new_gene_order
    initial_infercnv_object@expr.data = initial_infercnv_object@expr.data[rownames(new_gene_order), , drop=FALSE]

    library(igraph)
    library(reticulate)
    library(leiden) 
    
    print("InferCNV object created, running inferCNV...")
    system.time(
      infercnv_output <- try(
        infercnv::run(
          initial_infercnv_object,
          num_threads=numcores-1,
          HMM = T,
          out_dir= paste0("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/inferCNV/", i),
          cutoff=0.1,
          window_length=101,
          max_centered_threshold=3,
          cluster_by_groups=T,
          plot_steps=F,
          denoise=T,
          sd_amplifier=2,
          analysis_mode = "subclusters"
        )
      )
    )
  }
  
  # remove temporary files ----- 
  
  out_dir <= "/work/InternalMedicine/s437775/simil/simil_cancerepi/091322"
  system(paste0("rm ", out_dir, "/*infercnv_obj"))
  system(paste0("rm ", out_dir, "/*dat"))
  system(paste0("rm ", out_dir, "/*preliminary*"))
  system(paste0("rm ", out_dir, "/run.final.infercnv_obj"))
  
  
}

# ---------
# ---------
# Add metadata to cancer.epi object ------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/inferCNV")
 
CNVmetadata <- data.frame()
for (i in c("Wu", "Xu", "Qian", "Karaayvaz", 
            "Pal_Prim_1", "Pal_Prim_2", "Pal_Prim_3", "Pal_Prim_4", 
            "Wu2021prim_1", "Wu2021prim_2", "Wu2021prim_3")) {
  
  infercnv_output_meta <- as.data.frame(t(read.table(paste0(i, "/infercnv.observations.txt"))))
  infercnv_output_meta$cell <- "Epithelial"
  
  infercnv_ref_meta <- as.data.frame(t(read.table(paste0(i, "/infercnv.references.txt"))))
  infercnv_ref_meta$cell <- "Reference"
  
  heatmap <- rbind(infercnv_output_meta, infercnv_ref_meta)
  heatmap[,1:(dim(heatmap)[2]-1)] <- as.data.frame(scales::rescale(as.matrix(heatmap[,1:(dim(heatmap)[2]-1)]), c(-1,1)))
  heatmap$totalCNV <- rowMeans(heatmap[,1:(dim(heatmap)[2]-1)]^2)
  
  pdf(paste0(i, "_CNA_epithelial_density_plot.pdf"))
  print(ggplot(heatmap, aes(x = totalCNV, 
                      color = factor(cell, levels = c("Reference", "Epithelial")), 
                      fill = factor(cell, levels = c("Reference", "Epithelial")))) + 
    geom_histogram(alpha=0.5, position="identity") + theme_classic() + theme_bw() + 
    scale_color_discrete(name = "Cell") + scale_fill_discrete(name = "Cell"))
  dev.off()
  
  #heatmap_add <- as.data.frame(heatmap$totalCNV)
  heatmap_add <- as.data.frame(t(heatmap))
  rownames(heatmap_add) <- colnames(heatmap)
  heatmap_add$join <- rownames(heatmap_add) 
  #colnames(heatmap_add) <- c("totalCNV")
  
  if (i == "Wu") {
    CNVmetadata <- heatmap_add
  } else {
    CNVmetadata <- full_join(CNVmetadata, heatmap_add, by = "join")
  }
}

#write.csv(CNVmetadata, "CNV_metadata.csv")
write.csv(CNVmetadata, "CNV_metadata_Bassez.csv")

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/inferCNV")
CNVmetadata <- read.csv("CNV_metadata_Bassez.csv", row.names = 1)

CNVmetadata <- heatmap
rownames(CNVmetadata) <- gsub("\\.", "-", rownames(CNVmetadata))
cancer.epi <- AddMetaData(cancer.epi, CNVmetadata)#, col.name = "totalCNV")

pdf("totalCNV_UMAP_Bassez.pdf")
FeaturePlot(cancer.epi, feature = "totalCNV", 
            min.cutoff = 0,
            max.cutoff = 0.05)
dev.off()

for (i in 1:10) { 
  p <- FeaturePlot(cancer.epi, feature = paste0("GM",i), min.cutoff = 0, max.cutoff = 2)
  ggsave(file = paste0("GE", i, "_UMAP.pdf"),
         plot = as.ggplot(p))
}
# ---------
# ---------
# CNV analysis -------

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/inferCNV")
cancer.epi <- readRDS("cancerepi_withCNVmetadata_102922.rds")

avg_CNV <- data.frame()
for (i in unique(cancer.epi$samples)) {
  subset <- subset(cancer.epi, subset = samples == i)
  #subset_CNV <- subset@meta.data[,c(135:4942,4944:15186)]
  subset_CNV <- subset@meta.data[,c(57:7367)]
  subset_CNV <- purrr::map_df(subset_CNV, ~ as.numeric(gsub(" ", "", .x)))
  subset_CNV <- colMeans(subset_CNV, na.rm = T)
  subset_CNV <- as.data.frame(t(subset_CNV))
  avg_CNV <- rbind(avg_CNV, subset_CNV)
}

rownames(avg_CNV) <- unique(cancer.epi$samples)
write.csv(avg_CNV, "cancerepi_avgCNVpersample.csv")

avg_CNV <- read.csv("cancerepi_avgCNVpersample.csv", header = T, row.names = 1)

cancerepi_cor_CNV <- data.frame()
for (i in unique(cancer.epi$samples)) {
  subset <- subset(cancer.epi, subset = samples == i)
  #subset_CNV <- subset@meta.data[,c(135:4942,4944:15186)]
  subset_CNV <- subset@meta.data[,c(57:7367)]
  subset_CNV <- purrr::map_df(subset_CNV, ~ as.numeric(gsub(" ", "", .x)))
  subset_CNV <- as.data.frame(subset_CNV)
  rownames(subset_CNV) <- rownames(subset@meta.data)
  subset_CNV <- t(subset_CNV)
  if(length(which(!is.na(avg_CNV[which(rownames(avg_CNV) == i),]))) > 0) {
    cor_CNV <- cor(t(avg_CNV[which(rownames(avg_CNV) == i),]), subset_CNV, use = "pairwise.complete.obs")
    rownames(cor_CNV) <- c("corCNV")
    cancerepi_cor_CNV <- rbind(cancerepi_cor_CNV, t(cor_CNV))
  }
}

cancer.epi <- AddMetaData(cancer.epi, cancerepi_cor_CNV, col.name = "corCNV")

cancer.epi <- readRDS("cancerepi_withCNVmetadata_102922.rds")

cancer.epi$totalCNV <- as.numeric(cancer.epi$totalCNV)
cancer.epi$corCNV <- as.numeric(cancer.epi$corCNV)
cancer.epi$cancer <- NA
cancer.epi$cancer[which((cancer.epi$totalCNV < 0.02) & (cancer.epi$corCNV < 0.4))] <- "N"
cancer.epi$cancer[which((cancer.epi$totalCNV >= 0.02) | (cancer.epi$corCNV >= 0.4))] <- "Y"

setwd("/endosome/work/InternalMedicine/s437775/simil/simil_cancerepi/091322/110322_rerun")
saveRDS(cancer.epi, "Bassez_withCNVmetadata_111022.rds")

pdf("cancerepi_CNVscatter_cancerepi.pdf", width = 6, height = 4)
ggscatter(cancer.epi@meta.data, x = "totalCNV", y = "corCNV", 
          color = "cancer",
          size = 0.1,
          xlab = "CNV Signal", ylab = "CNV Correlation", 
          xlim = c(0,0.4), ylim = c(-0.3, 1)) + 
  geom_hline(yintercept = 0.4, linetype = 'dotted', col = 'red') + 
  geom_vline(xintercept = 0.02, linetype = 'dotted', col = 'red') + 
  scale_color_discrete(name = 'Classification', 
                       labels = c("Malignant (Cancer)","Normal (Non-cancer)")) + 
  guides(shape = guide_none())
dev.off()

noncancerepi <- subset(cancer.epi, subset = cancer == "NA")

pdf("cancerepi_UMAP_cancerepi.pdf", width = 5, height = 5)
DimPlot(cancer.epi, group.by = "cancer")
dev.off()

# ---------
# ---------
