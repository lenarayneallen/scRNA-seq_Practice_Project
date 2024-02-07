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
library(gridExtra)
library(patchwork)
library(cowplot)
library(AnnotationHub)
library(BiocManager)
library(BiocFileCache)
library(ensembldb)
library(BPCells)
library(rcartocolor)

### PROCESSING MERGED OBJECT FOR INTEGRATION -------------------------------------------
#read in merged object
merged <- readRDS("merged.rds")

#create a metadata column indicating sample and split by sample, patient, type, and barcode
merged$sample <- rownames(merged@meta.data)
merged@meta.data <- separate(merged@meta.data, col = 'sample', into = c('sample','patient', 'type', 'barcode'),
                                        sep = '_')

#standard preprocessing workflow
merged <- NormalizeData(object = merged)
merged <- FindVariableFeatures(object = merged)
merged <- ScaleData(object = merged)
merged <- RunPCA(object = merged)
ElbowPlot(merged)
merged <- FindNeighbors(object = merged, dims = 1:30)
merged <- FindClusters(object = merged)
merged <- RunUMAP(object = merged, dims = 1:30)

#printing UMAPS grouped by patient, type, and sample to check for batch effects
print(plot_by_patient <- DimPlot(merged, reduction = 'umap', group.by = 'patient'))
print(plot_by_type <- DimPlot(merged, reduction = 'umap', group.by = 'type',
                        cols = c('red', 'green', 'blue')))
print(plot_by_sample <- DimPlot(merged, reduction = 'umap', group.by = 'sample'))

#splitting object into layers by sample
objects_split <- SplitObject(merged, split.by = 'sample')

#log normalize, find variable features, and scale data
#regress out variation due to UMI counts and % mitochondrial reads
for(i in 1:length(objects_split)){
  objects_split[[i]] <- NormalizeData(objects_split[[i]],
                                 normalization.method = "LogNormalize")
  objects_split[[i]] <- FindVariableFeatures(objects_split[[i]])
  objects_split[[i]] <- ScaleData(objects_split[[i]], 
                                 vars.to.regress = c("PerMito", "nCount_RNA"))
}


###INTEGRATION & CLUSTERING-------------------------------------------------------------
#choose top 2000 variable features for integration
features <- SelectIntegrationFeatures(object.list = objects_split)

#find integration anchors
anchors <- FindIntegrationAnchors(object.list = objects_split,
                                 anchor.features = features)
#save anchors for future use/convenience
saveRDS(anchors, 'anchors.rds')
saveRDS(objects_split, 'split.rds')


#integrate and save integrated object
integrated <- IntegrateData(anchorset = anchors)
saveRDS(integrated, 'integrated.rds')


#re-run preprocessing steps
integrated <- ScaleData(object = integrated)
integrated <- RunPCA(object = integrated)
integrated <- RunUMAP(object = integrated, dims = 1:30)

#find neighbors and perform clustering at a resolution of 0.8
integrated <- FindNeighbors(object = integrated, dims = 1:30)
integrated <- FindClusters(object = integrated, 
                           resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

#save integrated and clustered object
saveRDS(integrated,'integrated_clustered.rds')
Idents(object = integrated) <- "integrated_snn_res.0.8"


###CLUSTRING QUALITY CONTROL & VISUALIZATION--------------------------------------------
#visualize clusters on integrated dataset
print(DimPlot(integrated, 
              reduction = "umap", 
              label = TRUE, 
              label.size = 6))

#UMAP visualization: cells per cluster of each sample
integrated_per_cluster_per_sample <- DimPlot(integrated, 
                          label = TRUE, 
                          split.by = "sample", 
                          ncol = 4) + NoLegend()

print(integrated_per_cluster_per_sample)


#UMAP visualization
int_umap_by_sample <- DimPlot(object = integrated, reduction = 'umap', group.by = 'sample') + labs(title = "integrated UMAP by sample")
int_umap_by_patient <- DimPlot(object = integrated, reduction = 'umap', group.by = 'patient') + labs(title = "integrated UMAP by patient")
print(int_umap_by_sample)
print(int_umap_by_patient)

###FINDING MARKERS & PRE-ANNOTATION------------------------------
#finding all markers
markers <- Find(object = integrated, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)
saveRDS(markers, 'all_markers.rds')

allmarkers <-readRDS('all_markers.rds')

#retrieve annotations from annotation hub and merge with all markers
ah <- AnnotationHub()
human_ens <- query(ah, 
             c("Homo sapiens", "EnsDb"), 
             ignore.case = TRUE)
id <- human_ens %>%
      mcols() %>%
      rownames() %>%
      tail(n = 1)

edb <- ah[[id]]

annotations <- genes(edb,
                     return.type = "data.frame")

annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

allmarkers.annotated <- left_join(x = allmarkers, y = unique(annotations[, c("gene_name", "description")]), 
                                  by = c("gene" = "gene_name"))

#create separate dataframes for each cluster containing top 20 genes
#compare top 20 genes to Panglao DB
for(i in 0:23){
  x <- allmarkers.annotated[allmarkers.annotated$cluster == i,]
  x <- top_n(x,
             n = 20, 
             wt = avg_log2FC)
  assign(paste0('cluster_', i, '_annotated'), x)
}

saveRDS(allmarkers.annotated, 'allmarkers.annotated.rds')
