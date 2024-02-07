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
library(RColorBrewer)
library(rcartocolor)

###read in file------------------------------------------------------------------------
integrated_clustered_final <- readRDS('integrated_clustered_final.RDS')

#set identities at 0.8 resolution
Idents(object = integrated_clustered_final) <- "integrated_snn_res.0.8.x"

###"MANUAL" ANNOTATION: visualization of marker genes-----------------------------------
#generate feature plots of markers for cell types identified in paper but not predicted by SingleR

#feature plots of known ductal cell markers
FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = "CFTR", 
            order = TRUE, 
            min.cutoff = 'q10',
            label = TRUE)

FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = "TFF1", 
            order = TRUE, 
            min.cutoff = 'q10',
            label = TRUE)

FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = "KRT19", 
            order = TRUE, 
            min.cutoff = 'q10',
            label = TRUE)

#feature plots of MKI67+ ductal cell markers
FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = "MKI67", 
            order = TRUE, 
            min.cutoff = 'q10',
            label = TRUE)

#feature plots of known mast cell markers
FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = "TPSAB1", 
            order = TRUE, 
            min.cutoff = 'q10',
            label = TRUE)

FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = "KIT", 
            order = TRUE, 
            min.cutoff = 'q10',
            label = TRUE)

#feature plots of known acinar cell markers
FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = "PRSS1", 
            order = TRUE, 
            min.cutoff = 'q10',
            label = TRUE)

FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = "KLK1", 
            order = TRUE, 
            min.cutoff = 'q10',
            label = TRUE)

#feature plots of known plasma cell markers
FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = "MZB1", 
            order = TRUE, 
            min.cutoff = 'q10',
            label = TRUE)

#feature plots of known fibroblast markers
FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = "COL1A1", 
            order = TRUE, 
            min.cutoff = 'q10',
            label = TRUE)


###"MANUAL" ANNOTATION: renaming identities --------------------------------------------
#Identites were determined by SingleR results and above feature plots of markers
integrated_clustered_final$orig.cluster <- Idents(integrated_clustered_final)
integrated_clustered_final <- RenameIdents(object = integrated_clustered_final, 
                                     "0" = "T cells",
                                     "1" = "T cells",
                                     "2" = "Ductal cells",
                                     "3" = "Mac/Mono",
                                     "4" = "Ductal cells",
                                     "5" = "Acinar cells",
                                     "6" = "Fibroblasts",
                                     "7" = "NK cells",
                                     "8" = "T cells",
                                     "9" = "Ductal cells",
                                     "10" = "Mac/Mono",
                                     "11" = "T cells",
                                     "12" = "Neutrophils",
                                     "13" = "T cells",
                                     "14" = "MKI67+ cells",
                                     "15" = "Mast cells",
                                     "16" = "Fibroblasts",
                                     "17" = "Mac/Mono", 
                                     "18" = "Mac/Mono", 
                                     "19" = "B cells", 
                                     "20" = "Plasma cells",
                                     "21" = "Endothelial cells",
                                     "22" = "MKI67+ cells",
                                     "23" = "Mac/Mono")


#save object with final identities
saveRDS(integrated_clustered_final, "integrated_clustered_final_idents.RDS")


#create UMAP plot with final identities
DimPlot(integrated_clustered_final, 
        reduction = "umap",
        label = TRUE, 
        label.size = 3)


###FURTHER VISUALIZATION ---------------------------------------------------------------
#set default assay
DefaultAssay(integrated_clustered_final) <- "integrated"
integrated_clustered_final$celltype <- Idents(integrated_clustered_final)


#show table of number of cells per identity among all conditions
table(Idents(integrated_clustered_final))


#plot bar plot of cell type abundance by condition
celltype_v_sampletype <- table(integrated_clustered_final$celltype, integrated_clustered_final$type.x)
cvs_df_perc <- as.data.frame(celltype_v_sampletype * 100 / rowSums(celltype_v_sampletype))

names(cvs_df_perc)[names(cvs_df_perc) == "Var1"] <- "celltype"
names(cvs_df_perc)[names(cvs_df_perc) == "Var2"] <- "condition"
names(cvs_df_perc)[names(cvs_df_perc) == "Freq"] <- "percent"

ggplot(cvs_df_perc, aes(x = celltype, y = percent, fill = condition)) +
  geom_bar(stat = 'identity', position = "dodge") +
  labs(colour = "condition", y = "percentage of cells", x = "cell type") +
  scale_fill_manual(values = c("HM" = "firebrick", "NT" = "royalblue3", "PT" = "tan2")) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("cell type abundance by condition")

#stacked bar plot of cell type % per sample condition
cvs_df <- as.data.frame(celltype_v_sampletype)

names(cvs_df)[names(cvs_df) == "Var1"] <- "celltype"
names(cvs_df)[names(cvs_df) == "Var2"] <- "condition"
names(cvs_df)[names(cvs_df) == "Freq"] <- "frequency"

ggplot(cvs_df, aes(x = condition, y = frequency, fill = celltype)) +
  geom_bar(stat = 'identity', position = "fill") +
  labs(colour = "cell type", y = "percentage of cells") +
  scale_fill_carto_d(name = "cell type", palette = "Bold") +
  theme_light() +
  ggtitle("Cell Type % per Sample Condition")

