# About this Analysis
### As a practice exercise, I used Seurat to analyze publicly available scRNA-seq data from Zhang et. al's 2023 study "Single cell transcriptomic analyses implicate an immunosuppressive tumor microenvironment in pancreatic cancer liver metastasis".

Zhang, S., Fang, W., Zhou, S. et al. Single cell transcriptomic analyses implicate an immunosuppressive tumor microenvironment in pancreatic cancer liver metastasis. Nat Commun 14, 5123 (2023). https://doi.org/10.1038/s41467-023-40727-7

**Link to the publication:** https://www.nature.com/articles/s41467-023-40727-7


## Introduction:
Here I will provide an overview of this scRNA-seq analysis, attempt to explain my methodology, and compare my results to those of the authors. As I am a beginner (and I have not been able to consult anyone with more knowledge on the subject for help),  some of the steps that I have taken may not be considered “best practices”; this analysis is exclusively a learning exercise and should not be taken as instructional materials. 

## Experimental Background: 
In this study, the authors sought to understand the cellular composition and microenvironment of primary Pancreatic Ductal Adenocarcinoma (PDAC) tumors and PDAC hepatic metastases through the use of scRNA-seq. Eight clinical samples from four patients were analyzed, and the samples derived from each patient are as follows:

![image](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/5cb378fb-3477-4305-88ee-7ff39aa5802c)



## Author's analysis methods:
For comparison purposes, here are the methods as described in the Zhang et al. publication: 

### "_scRNA-seq data processing, cluster annotation and data integration_"

"_The 10x Chromium single-cell RNA sequencing (scRNA-seq) data were processed using CellRanger (v3.1.0; 10x Genomics) for alignment, barcode assignment and unique molecular identifier (UMI) counting (using the genome reference set GRCh38-3.0.0). Filtered count matrices were converted to sparse matrices using the Seurat package (v3.2.3)55, and cells expressing less than 200 genes as well as with more than 20% mitochondrial reads, were excluded from the downstream analysis. The ‘doubletFinder_v3’ method from the DoubletFinder package (v2.0.3)56 was applied for additional cell filtering. Filtered data were then log normalized and scaled, with cell–cell variation due to UMI counts and percent mitochondrial reads regressed out._"

"_To avoid batch effects among samples and experiments, integration of single-cell data was performed using Seurat’s canonical correlation analysis (CCA) integration method. A total of 2000 features for anchoring (the ‘FindIntegrationAnchors’ function) and 30 dimensions for alignment (‘IntegrateData’) were used. Cell clustering was performed by ‘FindClusters’ function at a resolution of 0.8 and the top 20 genes were used to define cell identity. Dimensionality reduction was performed with ‘RunUMAP’ function and visualized by Uniform Manifold Approximation and Projection (UMAP). For subgroup cell clustering, cells of different types were extracted separately and clustered by their respective first 30 principal components (PCs) using different resolutions based on visual inspection._"

### "_Identification of signature genes_"
"_We applied the ‘FindAllMarkers’ function in Seurat to identify specific genes for each cell subset. For the selection of marker genes specific to each cell cluster/subset, we calculated the log2 fold change (log2FC) between two groups (a cell cluster/subset vs. other cells) using the ‘FindMarkers’ function with the Wilcoxon rank-sum test (default parameters)._"

## My analysis:
### Pre-processing and merging:
_**Creating the Seurat Object**_
- After reading in my counts matrices and creating the seurat objects for each sample, I oriented myself with the structure of the objects. It was immediately important to understand the differences between nFeature_RNA and nCount_RNA. Whereas nFeature_RNA provides the number of unique genes per individual cell, nCount_RNA provides the number of total molecules per individual cell. 
- Once I created the object, I filtered out cells expressing a percentage of mitochondrial genes greater than 20% and cells expressing less than 200 features (genes). This helps to remove cells that are low quality or may be dying.
- I then pre-processed each individual seurat object with the standard seurat pre-processing workflow. I created a new pre-processed object for each sample while retaining the filtered and unprocessed object for later use after running doubletfinder. 

_**Doubletfinder**_

- Like the authors, I utilized a tool called doubletfinder to identify and remove potential doublets from my seurat objects. Doublets are artifacts created when two cells are erroneously sequenced in a single reaction and thus counted as one cell.
- As it appears to be best practice to run doubletfinder with unmerged and/or unintegrated data, I ran doubletfinder on each sample individually. This is why I created a different seurat object for each sample in this initial preprocessing stage. 
- Doubletfinder requires the calculation of two parameters before running: pK and pExp. I utilized the no ground truth method to calculate the pK parameter (neighborhood size), and estimated the pExp for each seurat object/sample using a chart from 10x Genomics that defines approximate multiplet rates per number of cells loaded. 

- https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-


- Though the authors’ methodology indicates that they normalized their counts after running doubletfinder rather than before, doubletfinder required me to normalize my counts prior to running doubletfinder. To work around this, I created a new object after running the standard seurat pre-processing workflow (which includes normalization). I ran doubletfinder on this new object, identified the doublets, and then removed them from the original, un-normalized object.

_**Merging**_
- After running doubletfinder, I then merged the seurat objects for each individual sample into one seurat object.
  
### Integration and clustering:
_**Integration**_
- To prepare for integration, I re-ran the standard pre-processing workflow as defined above. I visualized the UMAP plots of the object, grouping by patient, sample, and sample type to assess for any bias, and I split the merged seurat object into layers by sample.
- I regressed out variation due to UMI counts and percent mitochondrial reads.
- As the authors noted, I selected the top 2000 variable features and used the default 30 dimensions for integration.

_**Clustering**_
- Though I performed clustering at a range of resolutions for curiosity's sake, I proceeded with analysis with the clustering results achieved at a resolution of 0.8 (as indicated by the authors).
- UMAP of integrated cells clustered at 0.8 resolution:
![post_integration_final](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/e4d05198-a84d-4c21-9239-8a4b49eae362)

- To assess for bias before proceeding with marker identification, I then visualized the integrated and clustered UMAP above in a variety of different ways:
  - _Split by sample:_
![int_per_cluster_per_sample](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/369988e5-4501-4cc9-9f61-444b7204cc1c)

  - _Grouped by sample:_

    ![int_umap_by_sample](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/cfcc2d56-b8ad-481a-9665-1dab8bffe73c)

  - _Grouped by patient:_
  
    ![int_umap_by_patient](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/b93ca81d-42c3-40f8-8705-fe17e0d84338)

_**Marker Identification**_
- I used the FindAllMarkers() function to identify the top 20 markers for each cluster. I annotated these markers using AnnotationHub and explored these markers using PanglaoDB.

### SingleR:
- Here my analysis branches from that of the authors; whereas the authors did not employ this tool, I decided to utilize SingleR for reference-based annotation.
- I utilized the HumanPrimaryCellAtlasData as a reference
- SingleR can provide cell-wise (default) or cluster-wise annotation. I generated SingleR predictions using both of these methods:
  - Cluster-wise annotation (faster!):
   ![singler_pred_ids_by_cluster](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/50b68000-78d4-4d14-bf90-be059495c4a0)

  - Single-cell-wise annotation:
![singleR_pred_by_cell](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/d47c15e4-80cc-4954-b3df-87b9248ea44e)

- Ultimately, I found the single-cell-wise resolution to be the most informative and nuanced; to understand the most prevalent cell types in each cluster, I generated a bar plot showing the proportion of SingleR-identified cell types in each cluster:
![singleRlabels_per_cluster_cellwise](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/88e75780-e0f6-494c-abf2-bad837fa0903)

### Annotation and visualization:
- I noticed that while mostly representative of the author's findings, SingleR did not identify a few cell types that were identified by the authors (mast cells, ductal cells, MKI67+ ductal cells, acinar cells, endocrine cells, fibroblasts, and plasma cells). In turn, SingleR identified some cell types that were not observed in the author's figure (Tissue Stem Cells, Neutrophils, Epithelial Cells, CMPs).
- To see if I could identify the groups noted by the authors but not by singleR, I generated feature plots of known markers for those groups:
    - _Ductal cell and MKI67 ductal cell markers_
    ![DUCTAL_MARKERS](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/c4e8a609-2b1f-41e3-a442-7bbb4e5cc0f0)
    - _Mast cell markers_
    ![MASTCELL_MARKERS](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/70bb922d-2416-4eb5-96a6-8c6eb904b557)
    - _Acinar and Endocrine Markers_
    ![ACINAR_MARKERS](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/9c248066-a71c-4606-9746-7aa87f08a779)
    ![ENDOCRINE](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/d281eb9d-f3af-400a-b9a6-94acddf26ea9)
    - _Plasma cell markers_
    ![plasma_MZB1_final](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/377c121a-08d4-4397-910c-afc35a4d7a01)
    - _Fibroblast markers_
    ![fibroblasts_COL1A1_final](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/ec3c055d-ffd1-491f-82fe-0d4ce034fdca)

- Starting with the Single-R annotations and making adjustments based off of the feature plots above, I assigned cluster identities as indicated in the below figure. It is important to note that this annotation is a guesstimate on my part. I understand that I am "working backwards" by looking specifically for markers of cell populations already identified by the authors, and that my annotation process was not particularly systematic, but I wanted to attempt it as an exercise nonetheless. 
![acinarendocrine](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/301b1dd6-7dc3-459e-97b3-a7d4dfd6c231)

- Here is my final annotated plot compared to that of the authors (Zhang et al. Figure 1C):
![ductal cell markers (2)](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/e83f4724-08db-4500-aa82-e52a78f0cbdd)

- From here, I was able to loosely recreate a bar plot from the paper (Zhang et al. Figure 1E) examining the abundance of each cell type by sample condition. Error bars for the pancreatic tumor (PT) and hepatic metastasis (HM) conditions represent the standard error of the mean; there are no error bars for the normal tissue (NT) condition as there is only one normal tissue sample. 

  <img src="https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/abb2dfa8-7531-4ad1-b665-37f0c4fbe27b" width=60% height=60%>
  <img src="https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/a89a669c-416e-470d-a3a3-b1d4725e2348" width=30% height=30%> 

- While I am unsure of their statistical significance, several trends observed in my plot are consistent with the authors' plot. These trends include:
  - increased abundance of acinar and endocrine cells (termed "secretory cells" in the paper) in the NT sample compared to the PT and HM samples
  - increased abundance of plasma cells in PT samples compared to HM samples
  - a lack of MKI67 cells in the NT sample
  - increased abundance of macrophage/monocyte cells in HM sample compared to the PT and NT samples
  - increased abundance of fibroblasts in PT samples compared to the HM and NT samples

- Despite these trends, some aspects of my results indicate that my analysis is in some way flawed. Notably, there is some presence of acinar/endocrine cells in the hepatic metastases according to my analysis. As acinar and endocrine cells are pancreatic tissue cells, they should theoretically not exist in the hepatic metastases. Furthermore, the authors did not identify acinar or endocrine cells in the hepatic metastases. This tells me that I likely misidentified the acinar/endocrine cluster as cluster 5; though this identification was made based on the feature plots of known acinar and endocrine marker genes (pictured above), these genes mapped onto cluster 5 relatively loosely.

- I was also able to create a stacked bar plot representing the cell type composition of each sample type:
  ![Cell Type PERC per Sample Condition FINAL](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/54db39cc-e2a0-4bd0-a49d-a37f980a9e15)




