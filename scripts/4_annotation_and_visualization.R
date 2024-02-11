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
library(Rmisc)
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

FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = c("CFTR", "TFF1", "MKI67"), 
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

FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = c("TPSAB1", "KIT",  "TPSB2", "MS4A2"), 
            order = TRUE, 
            min.cutoff = 'q10',
            label = TRUE) 

#feature plots of known acinar cell markers
FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = c("PRSS1", "KLK1"), 
            order = TRUE, 
            pt.size = 0.5,
            min.cutoff = 'q10',
            label = TRUE)

#feature plot fo known endocrine cell markers
FeaturePlot(integrated_clustered_final, 
            reduction = "umap",
            features = "INS", 
            order = TRUE, 
            min.cutoff = 'q10',
            label = FALSE)

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
                                     "5" = "Acinar/Endocrine",
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
        label.size = 4)


###FURTHER VISUALIZATION ---------------------------------------------------------------
#set default assay
DefaultAssay(integrated_clustered_final) <- "integrated"
integrated_clustered_final$celltype <- Idents(integrated_clustered_final)


#show table of number of cells per identity among all conditions
table(Idents(integrated_clustered_final))


#plot bar plot of cell type abundance by condition
celltype_v_sample <- table(integrated_clustered_final$celltype, integrated_clustered_final$sample.x)
celltype_v_sample_df <- as.data.frame(celltype_v_sample * 100 / rowSums(celltype_v_sample))

names(celltype_v_sample_df)[names(celltype_v_sample_df) == "Var1"] <- "celltype"
names(celltype_v_sample_df)[names(celltype_v_sample_df) == "Var2"] <- "sample"
names(celltype_v_sample_df)[names(celltype_v_sample_df) == "Freq"] <- "percent"

cst <- celltype_v_sample_df %>%
  mutate(type = case_when(
      grepl("HM", sample) ~ "HM",
      grepl("NT", sample) ~ "NT",
      grepl("PT", sample) ~ "PT"
    )
  )

summary <- summarySE(cst, measurevar = "percent", groupvars = c("celltype", "type"))

ggplot(summary, aes(x = celltype, y = percent, fill = type)) +
  geom_bar(stat = 'identity', position = "dodge") + 
  geom_errorbar(aes(ymin = percent-se, ymax = percent+se),
                width = .2,
                position = position_dodge(.9)) +
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
  

