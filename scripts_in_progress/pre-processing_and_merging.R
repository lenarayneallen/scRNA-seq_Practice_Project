library(Seurat)
library(tidyverse)
library(GEOquery)
library(dplyr)
library(DoubletFinder)
library(Matrix)
library(fields)
library(KernSmooth)
library(ROCR)
library(devtools)
library(remotes)
library(ggplot2)
library(patchwork)

#read in counts ---------------------------------------------------------------
E1_YF <- ReadMtx(mtx = "GSM5910784_Case1-YF_matrix.mtx.gz",
                       features = "GSM5910784_Case1-YF_features.tsv.gz",
                       cells = "GSM5910784_Case1-YF_barcodes.tsv.gz")


E1_ZY <- ReadMtx(mtx = "GSM5910785_Case1-ZY_matrix.mtx.gz",
                       features = "GSM5910785_Case1-ZY_features.tsv.gz",
                       cells = "GSM5910785_Case1-ZY_barcodes.tsv.gz")

E2_YF <- ReadMtx(mtx = "GSM5910787_Case2-YF_matrix.mtx.gz",
                       features = "GSM5910787_Case2-YF_features.tsv.gz",
                       cells = "GSM5910787_Case2-YF_barcodes.tsv.gz")

E2_ZY <- ReadMtx(mtx = "GSM5910788_Case2-ZY_matrix.mtx.gz",
                       features = "GSM5910788_Case2-ZY_features.tsv.gz",
                       cells = "GSM5910788_Case2-ZY_barcodes.tsv.gz")

E2_ZC <- ReadMtx(mtx = "GSM5910786_Case2-ZC_matrix.mtx.gz",
                       features = "GSM5910786_Case2-ZC_features.tsv.gz",
                       cells = "GSM5910786_Case2-ZC_barcodes.tsv.gz")

E3_YF <- ReadMtx(mtx = "GSM5910789_Case3-YF_matrix.mtx.gz",
                       features = "GSM5910789_Case3-YF_features.tsv.gz",
                       cells = "GSM5910789_Case3-YF_barcodes.tsv.gz")

E3_ZY <- ReadMtx(mtx = "GSM5910790_Case3-ZY_matrix.mtx.gz",
                       features = "GSM5910790_Case3-ZY_features.tsv.gz",
                       cells = "GSM5910790_Case3-ZY_barcodes.tsv.gz")

E4_ZY <- ReadMtx(mtx = "GSM5910791_Case4-ZY_matrix.mtx.gz",
                       features = "GSM5910791_Case4-ZY_features.tsv.gz",
                       cells = "GSM5910791_Case4-ZY_barcodes.tsv.gz")

counts_list <- list(E1_YF, E1_ZY, E2_YF, E2_ZY, E2_ZC, E3_YF, E3_ZY, E4_ZY)
sample_names <-c('E1_YF', 'E1_ZY', 'E2_YF', 'E2_ZY', 'E2_ZC', 'E3_YF', 'E3_ZY', 'E4_ZY')

#set the Seurat assay version to v3 for compatibility with DoubletFinder later on
options(Seurat.object.assay.version = "v3")
seurat_objects <- list()

#create seurat objects ---------------------------------------------------------
for (count in counts_list){
  seurat_object <- CreateSeuratObject(count)
  seurat_objects <- append(seurat_objects, seurat_object) 
}
names(seurat_objects) <- c('E1_YF_ser', 'E1_ZY_ser', 'E2_YF_ser', 'E2_ZY_ser', 'E2_ZC_ser', 'E3_YF_ser', 'E3_ZY_ser', 'E4_ZY_ser')



#initial filtering--------------------------------------------------------------
filtered_objects <- list()

for (seurat_obj in seurat_objects){
  #create column of metadata indicating percent of mitochondrial genes per cell
  seurat_obj$PerMito <- PercentageFeatureSet(seurat_obj, pattern = '^MT-')
  
  #filter, keeping only cells expressing more than 200 genes, with less than 20% 
  #of mitochondrial reads
  seurat_obj_filtered <- subset(seurat_obj, subset = nFeature_RNA > 200 &
                                  PerMito < 20)
  #append to list of filtered objects
  filtered_objects <- append(filtered_objects, seurat_obj_filtered)
  
}

names(filtered_objects) <- c('E1_YF_fil', 'E1_ZY_fil', 'E2_YF_fil', 'E2_ZY_fil', 'E2_ZC_fil', 'E3_YF_fil', 'E3_ZY_fil', 'E4_ZY_fil')

E1_YF_fil <- filtered_objects[[1]]
E1_ZY_fil <- filtered_objects[[2]]
E2_YF_fil <- filtered_objects[[3]]
E2_ZY_fil <- filtered_objects[[4]]
E2_ZC_fil <- filtered_objects[[5]]
E3_YF_fil <- filtered_objects[[6]]
E3_ZY_fil <- filtered_objects[[7]]
E4_ZY_fil <- filtered_objects[[8]]

fils <- list(E1_YF_fil, E1_ZY_fil, E2_YF_fil, E2_ZY_fil, E2_ZC_fil, E3_YF_fil, E3_ZY_fil, E4_ZY_fil)

#doubletfinder------------------------------------------------------------------
#pre-process data
pp_objects <- list()
for (object in fils){
  seurat_pp <- NormalizeData(object)
  seurat_pp <- FindVariableFeatures(object = seurat_pp, selection.method = "vst", nfeatures = 2000)
  seurat_pp <- ScaleData(object = seurat_pp)
  seurat_pp <- RunPCA(object = seurat_pp)
  seurat_pp <- FindNeighbors(object = seurat_pp, dims = 1:30)
  seurat_pp <- FindClusters(object = seurat_pp)
  seurat_pp <- RunUMAP(object = seurat_pp, dims = 1:30)
  pp_objects <- append(pp_objects, seurat_pp)
}


names(pp_objects) <- c('E1_YF_pp', 'E1_ZY_pp', 'E2_YF_pp', 'E2_ZY_pp', 'E2_ZC_pp', 'E3_YF_pp', 'E3_ZY_pp', 'E4_ZY_pp')

E1_YF_pp <- pp_objects[[1]]
E1_ZY_pp <- pp_objects[[2]]
E2_YF_pp <- pp_objects[[3]]
E2_ZY_pp <- pp_objects[[4]]
E2_ZC_pp <- pp_objects[[5]]
E3_YF_pp <- pp_objects[[6]]
E3_ZY_pp <- pp_objects[[7]]
E4_ZY_pp <- pp_objects[[8]]

pps <- c(E1_YF_pp, E1_ZY_pp, E2_YF_pp, E2_ZY_pp, E2_ZC_pp, E3_YF_pp, E3_ZY_pp, E4_ZY_pp)

#identify pK
pK_list <- c()
for (object in pp_objects){
  #estimate pK value for each object/sample (no ground truth method)
  sweep_res <- paramSweep_v3(object, PCs = 1:30, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  
  ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    geom_point() +
    geom_line()
  
  #choose the pK value that corresponds with the maximum BCMetric
  pK <- bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  
  pK <- as.numeric(as.character(pK[[1]]))
  pK_list <- append(pK_list, pK)
}


#calculating nExp for each sample individually due to differing
#multiplet rates; taking homotypic doublets into account
nExp_adj_list <- c()

E1_YF_anno <- E1_YF_pp@meta.data$seurat_clusters
E1_YF_homotypic <- modelHomotypic(E1_YF_anno)
E1_YF_nExp <- round(0.064*nrow(E1_YF_pp@meta.data))
E1_YF_nExp_adj <- round(E1_YF_nExp*(1-E1_YF_homotypic))
nExp_adj_list <- append(nExp_adj_list, E1_YF_nExp_adj)

E1_ZY_anno <- E1_ZY_pp@meta.data$seurat_clusters
E1_ZY_homotypic <- modelHomotypic(E1_ZY_anno)
E1_ZY_nExp <- round(0.064*nrow(E1_ZY_pp@meta.data))
E1_ZY_nExp_adj <- round(E1_ZY_nExp*(1-E1_ZY_homotypic))
nExp_adj_list <- append(nExp_adj_list, E1_ZY_nExp_adj)

E2_YF_anno <- E2_YF_pp@meta.data$seurat_clusters
E2_YF_homotypic <- modelHomotypic(E2_YF_anno)
E2_YF_nExp <- round(0.08*nrow(E2_YF_pp@meta.data))
E2_YF_nExp_adj <- round(E2_YF_nExp*(1-E2_YF_homotypic))
nExp_adj_list <- append(nExp_adj_list, E2_YF_nExp_adj)

E2_ZY_anno <- E2_ZY_pp@meta.data$seurat_clusters
E2_ZY_homotypic <- modelHomotypic(E2_ZY_anno)
E2_ZY_nExp <- round(0.072*nrow(E2_ZY_pp@meta.data))
E2_ZY_nExp_adj <- round(E2_ZY_nExp*(1-E2_ZY_homotypic))
nExp_adj_list <- append(nExp_adj_list, E2_ZY_nExp_adj)

E2_ZC_anno <- E2_ZC_pp@meta.data$seurat_clusters
E2_ZC_homotypic <- modelHomotypic(E2_ZC_anno)
E2_ZC_nExp <- round(0.048*nrow(E2_ZC_pp@meta.data))
E2_ZC_nExp_adj <- round(E2_ZC_nExp*(1-E2_ZC_homotypic))
nExp_adj_list <- append(nExp_adj_list, E2_ZC_nExp_adj)

E3_YF_anno <- E3_YF_pp@meta.data$seurat_clusters
E3_YF_homotypic <- modelHomotypic(E3_YF_anno)
E3_YF_nExp <- round(0.064*nrow(E3_YF_pp@meta.data))
E3_YF_nExp_adj <- round(E3_YF_nExp*(1-E3_YF_homotypic))
nExp_adj_list <- append(nExp_adj_list, E3_YF_nExp_adj)

E3_ZY_anno <- E3_ZY_pp@meta.data$seurat_clusters
E3_ZY_homotypic <- modelHomotypic(E3_ZY_anno)
E3_ZY_nExp <- round(0.056*nrow(E3_ZY_pp@meta.data))
E3_ZY_nExp_adj <- round(E3_ZY_nExp*(1-E3_ZY_homotypic))
nExp_adj_list <- append(nExp_adj_list, E3_ZY_nExp_adj)

E4_ZY_anno <- E4_ZY_pp@meta.data$seurat_clusters
E4_ZY_homotypic <- modelHomotypic(E4_ZY_anno)
E4_ZY_nExp <- round(0.008*nrow(E4_ZY_pp@meta.data))
E4_ZY_nExp_adj <- round(E4_ZY_nExp*(1-E4_ZY_homotypic))
nExp_adj_list <- append(nExp_adj_list, E4_ZY_nExp_adj)

#run doubletfinder on each sample individually 
for (i in 1:8){
  dbf <- doubletFinder_v3(pps[[i]],
                          PCs = 1:30, 
                          pN = 0.25,
                          pK = pK_list[[i]], 
                          nExp = nExp_adj_list[[i]], 
                          reuse.pANN = FALSE, 
                          sct = FALSE)
  assign(paste0(sample_names[[i]], "_dbf"), dbf)
}

dbfs <- c(E1_YF_dbf, E1_ZY_dbf, E2_YF_dbf, E2_ZY_dbf, E2_ZC_dbf, E3_YF_dbf, E3_ZY_dbf, E4_ZY_dbf)

#renaming doubletfinder output column as "doublet_ind" for all samples
E1_YF_dbf$doublet_ind <- E1_YF_dbf$DF.classifications_0.25_0.28_471
E1_YF_dbf$DF.classifications_0.25_0.28_471 <- NULL

E1_ZY_dbf$doublet_ind <- E1_ZY_dbf$DF.classifications_0.25_0.3_502
E1_ZY_dbf$DF.classifications_0.25_0.3_502 <- NULL

E2_YF_dbf$doublet_ind <- E2_YF_dbf$DF.classifications_0.25_0.3_821
E2_YF_dbf$DF.classifications_0.25_0.3_821 <- NULL

E2_ZY_dbf$doublet_ind <- E2_ZY_dbf$DF.classifications_0.25_0.3_587
E2_ZY_dbf$DF.classifications_0.25_0.3_587 <- NULL

E2_ZC_dbf$doublet_ind <- E2_ZC_dbf$DF.classifications_0.25_0.21_281
E2_ZC_dbf$DF.classifications_0.25_0.21_281 <- NULL

E3_YF_dbf$doublet_ind <- E3_YF_dbf$DF.classifications_0.25_0.2_472
E3_YF_dbf$DF.classifications_0.25_0.2_472 <- NULL

E3_ZY_dbf$doublet_ind <- E3_ZY_dbf$DF.classifications_0.25_0.29_396
E3_ZY_dbf$DF.classifications_0.25_0.29_396 <- NULL

E4_ZY_dbf$doublet_ind <- E4_ZY_dbf$DF.classifications_0.25_0.3_10
E4_ZY_dbf$DF.classifications_0.25_0.3_10 <- NULL


dbfs <- list(E1_YF_dbf, E1_ZY_dbf, E2_YF_dbf, E2_ZY_dbf, E2_ZC_dbf, E3_YF_dbf, E3_ZY_dbf, E4_ZY_dbf)

#selecting singlets only from each object
for (i in 1:8){
  object <- dbfs[[i]]
  singlets <- subset(object, cells = WhichCells(object, expression = doublet_ind == "Singlet"))
  assign(paste0(sample_names[[i]], "_singlets"), singlets)
}

singlets <- list(E1_YF_singlets, E1_ZY_singlets, E2_YF_singlets, E2_ZY_singlets, E2_ZC_singlets, E3_YF_singlets, E3_ZY_singlets, E4_ZY_singlets)

#selecting only singlets from original filtered but NOT PRE PROCESSED Seurat object
for (i in 1:8){
  sing_obj <- singlets[[i]]
  fil_obj <- fils[[i]]
  fil_obj$cellname <- row.names(fil_obj@meta.data)
  singlet_cells <- row.names(sing_obj@meta.data)
  seurat_filt_singlets <- subset(fil_obj, subset = cellname %in% singlet_cells)
  assign(paste0(sample_names[[i]], "_filt_singlets"), seurat_filt_singlets)
}                 

filt_singlets <- list(E1_YF_filt_singlets, E1_ZY_filt_singlets, E2_YF_filt_singlets, E2_ZY_filt_singlets, E2_ZC_filt_singlets, E3_YF_filt_singlets, E3_ZY_filt_singlets, E4_ZY_filt_singlets)

#change assay type to V5 for each object
for (i in 1:8){
  object <- filt_singlets[[i]]
  object@assays$RNA5 <- as(object = object@assays$RNA, Class = "Assay5")
  object@active.assay <- "RNA5"
  object[["RNA"]] <- NULL
  assign(paste0(sample_names[[i]], "_filt_no_doublets_5"), object)
}
#removing unneccesary objects for memory purposes 
remove(fils, filtered_objects, pp_objects, pps)
remove(seurat_objects)
remove(E1_YF_filt_singlets, E1_ZY_filt_singlets, E2_YF_filt_singlets, E2_ZY_filt_singlets, E2_ZC_filt_singlets, E3_YF_filt_singlets, E3_ZY_filt_singlets, E4_ZY_filt_singlets)
remove(E1_YF_singlets, E1_ZY_singlets, E2_YF_singlets, E2_ZY_singlets, E2_ZC_singlets, E3_YF_singlets, E3_ZY_singlets, E4_ZY_singlets)
remove(E1_YF_dbf, E1_ZY_dbf, E2_YF_dbf, E2_ZY_dbf, E2_ZC_dbf, E3_YF_dbf, E3_ZY_dbf, E4_ZY_dbf)
remove(E1_YF_pp, E1_ZY_pp, E2_YF_pp, E2_ZY_pp, E2_ZC_pp, E3_YF_pp, E3_ZY_pp, E4_ZY_pp)
remove(E1_YF_fil, E1_ZY_fil, E2_YF_fil, E2_ZY_fil, E2_ZC_fil, E3_YF_fil, E3_ZY_fil, E4_ZY_fil)
remove(seurat_filt_singlets,fil_obj, seurat_pp, seurat_object, seruat_obj)

###MERGING AND INTEGRATION --------------------------------------------------                 
#merge Seurat objects
#edit cell id names to be more intuitive/match the paper
merged_seurat_obj <- merge(E1_YF_filt_no_doublets_5,
                           y = c(E1_ZY_filt_no_doublets_5,
                                 E2_YF_filt_no_doublets_5,
                                 E2_ZY_filt_no_doublets_5,
                                 E2_ZC_filt_no_doublets_5,
                                 E3_YF_filt_no_doublets_5,
                                 E3_ZY_filt_no_doublets_5,
                                 E4_ZY_filt_no_doublets_5),
                           add.cell.ids = c('P1PT_P1_PT', 'P1HM_P1_HM', 'P2PT_P2_PT', 'P2HM_P2_HM', 'P2NT_P2_NT', 'P3PT_P3_PT', 'P3HM_P3_HM', 'P4HM_P4_HM'),
                           project = 'pdac')

#save merged object
saveRDS(merged_seurat_obj, "merged.rds")

#adding sample column
merged_seurat_obj$sample <- rownames(merged_seurat_obj@meta.data)
merged_seurat_obj@meta.data <- separate(merged_seurat_obj@meta.data, col = 'sample', into = c('sample','patient', 'type', 'barcode'),
                                        sep = '_')


#log normalization
merged_seurat_obj <- NormalizeData(merged_seurat_obj,
                                   normalization.method = "LogNormalize")
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj)
merged_seurat_obj <- ScaleData(merged_seurat_obj, 
                                vars.to.regress = c("PerMito", "nCount_RNA"))

obj.list <- SplitObject(merged_seurat_obj, split.by = 'sample')


features <- SelectIntegrationFeatures(object.list = obj.list)
