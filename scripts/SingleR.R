library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)
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
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(SingleCellExperiment)
library(BPCells)
library(ggsci)

#PREPROCESSING ------------------------------------------
#load in reference genome
ref <- celldex::HumanPrimaryCellAtlasData(ensembl = TRUE)
View(as.data.frame(colData(ref)))

#read in integrated and clustered seurat object; select resolution
integrated_clustered <- readRDS('integrated_clustered.rds')
Idents(object = integrated_clustered) <- "integrated_snn_res.0.8"

#get NORMALIZED counts from the RNA layer (not the integrated layer)
#using BPCells to convert normalized counts for each sample to on-disk matrices
write_matrix_dir(mat = integrated_clustered[["RNA5"]]$data.1.1, 
                 dir = '/dir/zhangcounts/data.1.1')
write_matrix_dir(mat = integrated_clustered[["RNA5"]]$data.2.2, 
                 dir = '/dir/zhangcounts/data.2.2')
write_matrix_dir(mat = integrated_clustered[["RNA5"]]$data.3.3, 
                 dir = '/dir/zhangcounts/data.3.3')
write_matrix_dir(mat = integrated_clustered[["RNA5"]]$data.4.4, 
                 dir = '/dir/zhangcounts/data.4.4')
write_matrix_dir(mat = integrated_clustered[["RNA5"]]$data.5.5, 
                 dir = '/dir/zhangcounts/data.5.5')
write_matrix_dir(mat = integrated_clustered[["RNA5"]]$data.6.6, 
                 dir = '/dir/zhangcounts/data.6.6')
write_matrix_dir(mat = integrated_clustered[["RNA5"]]$data.7.7, 
                 dir = '/dir/zhangcounts/data.7.7')
write_matrix_dir(mat = integrated_clustered[["RNA5"]]$data.8.8, 
                 dir = '/dir/zhangcounts/data.8.8')

data.1.1.mat <- open_matrix_dir(dir = '/dir/zhangcounts/data.1.1')
data.2.2.mat <- open_matrix_dir(dir = '/dir/zhangcounts/data.2.2')
data.3.3.mat <- open_matrix_dir(dir = '/dir/zhangcounts/data.3.3')
data.4.4.mat <- open_matrix_dir(dir = '/dir/zhangcounts/data.4.4')
data.5.5.mat <- open_matrix_dir(dir = '/dir/zhangcounts/data.5.5')
data.6.6.mat <- open_matrix_dir(dir = '/dir/zhangcounts/data.6.6')
data.7.7.mat <- open_matrix_dir(dir = '/dir/zhangcounts/data.7.7')
data.8.8.mat <- open_matrix_dir(dir = '/dir/zhangcounts/data.8.8')

counts_list = c(data.1.1.mat, data.2.2.mat, data.3.3.mat, data.4.4.mat,
                data.5.5.mat, data.6.6.mat, data.7.7.mat,
                data.8.8.mat)

#Creating a new seurat object with only the on-disk matrices as counts in order to make joining layers less memory intensive
joined_object <- CreateSeuratObject(counts = counts_list, 
                                    meta.data = integrated_clustered@meta.data,
                                    project = "ZHANG")
#joining layers
joined_object <- JoinLayers(joined_object)

#converting the on-disk matrices back to a dgCMatrix
joined_object[["RNA"]]$counts <- as(object = joined_object[["RNA"]]$counts, Class = "dgCMatrix")

#getting counts
counts <- GetAssayData(joined_object, layer = 'counts')
gc()
#ensure counts and reference both use ensembl IDs so that SingleR will work properly

require(EnsDb.Hsapiens.v86)
ens <- mapIds(EnsDb.Hsapiens.v86,
              keys = rownames(counts),
              column = 'GENEID',
              keytype = 'SYMBOL')
all(rownames(counts) == names(ens))


keep <- !is.na(ens)
ens <- ens[keep]
counts <- counts[keep,]
rownames(counts) <- ens

#RUNNING SingleR -------------------------------------------------------------
#run SingleR at the single cell level
predictions <- SingleR(test = counts, 
                       ref = ref, 
                       labels = ref$label.main)

saveRDS(predictions, file = "predictions_norm.RDS")

#add single R cell identity predictions (single cell level) back to original integrated Seurat object
predictions_norm <- readRDS("predictions_norm.RDS")
integrated_clustered$singleR.labels <- predictions_norm$labels[match(rownames(integrated_clustered@meta.data), rownames(predictions_norm))]

saveRDS(integrated_clustered, 'int_clust_norm_singleR.RDS')



#run singleR at the CLUSTER level
predictions_by_cluster <- SingleR(test = counts, 
                                  ref = ref,
                                  labels = ref$label.main,
                                  clusters = joined_object@meta.data[["integrated_snn_res.0.8.x"]])

saveRDS(predictions_by_cluster, file = "predictions_by_cluster.RDS")

#add single R cell identity predictions (cluster level) back to joined Seurat object
joined_object[["SingleR.cluster.labels"]] <- predictions_by_cluster$labels[match(joined_object[[]][["integrated_snn_res.0.8.x"]], rownames(predictions_by_cluster))]


#add single R cell identity predictions (cluster level) back to integrated object
integrated_clustered$backup.rownames <- row.names(integrated_clustered@meta.data)
integrated_clustered@meta.data <- base::merge(integrated_clustered@meta.data, joined_object@meta.data, by = "row.names")
row.names(integrated_clustered@meta.data) <- integrated_clustered$Row.names
saveRDS(integrated_clustered, 'int_norm_singleR_bycluster.RDS')



#PLOTTING ----------------------------------
#by single cell
#plot UMAP grouping by singleR predictions at single cell resolution
DimPlot(integrated_clustered, reduction = 'umap', group.by = 'singleR.labels') +
  ggtitle("SingleR predicted identities by cell")

#make bar plot showing most common cell type per cluster at a resolution of 0.8
table(integrated_clustered$singleR.labels, integrated_clustered$integrated_snn_res.0.8.x)

ggplot(integrated_clustered@meta.data, aes(x = integrated_snn_res.0.8.x, fill = singleR.labels)) +
  geom_bar(position = "fill") +
  labs(color = "SingleR labels by cell", x = "cluster") +
  scale_fill_igv()

#by cluster
#plot by SingleR cluster at cluster resolution
DimPlot(integrated_clustered, reduction = 'umap', group.by = 'SingleR.cluster.labels') +
  ggtitle("SingleR predicted identities by cluster")
  

