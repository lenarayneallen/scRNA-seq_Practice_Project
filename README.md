# About this Analysis
### As a practice exercise, I used Seurat to analyze publicly available scRNA-seq data from Zhang et. al's 2023 study "Single cell transcriptomic analyses implicate an immunosuppressive tumor microenvironment in pancreatic cancer liver metastasis".

Zhang, S., Fang, W., Zhou, S. et al. Single cell transcriptomic analyses implicate an immunosuppressive tumor microenvironment in pancreatic cancer liver metastasis. Nat Commun 14, 5123 (2023). https://doi.org/10.1038/s41467-023-40727-7

**Link to the publication:** https://www.nature.com/articles/s41467-023-40727-7


## Introduction:
Here I will provide an overview of this scRNA-seq analysis, attempt to explain my methodology, and compare my results to those of the authors. As I am a beginner (and I have not been able to consult anyone with more knowledge on the subject for help), I understand that some of the steps that I have taken may not be considered “best practices”; this analysis is exclusively a learning exercise and should not be taken as instructional materials. 


## Author's analysis methods:
For comparison purposes, here are the methods as described in the Zhang et al. publication: 

### "_scRNA-seq data processing, cluster annotation and data integration_"

"_The 10x Chromium single-cell RNA sequencing (scRNA-seq) data were processed using CellRanger (v3.1.0; 10x Genomics) for alignment, barcode assignment and unique molecular identifier (UMI) counting (using the genome reference set GRCh38-3.0.0). Filtered count matrices were converted to sparse matrices using the Seurat package (v3.2.3)55, and cells expressing less than 200 genes as well as with more than 20% mitochondrial reads, were excluded from the downstream analysis. The ‘doubletFinder_v3’ method from the DoubletFinder package (v2.0.3)56 was applied for additional cell filtering. Filtered data were then log normalized and scaled, with cell–cell variation due to UMI counts and percent mitochondrial reads regressed out._"

"_To avoid batch effects among samples and experiments, integration of single-cell data was performed using Seurat’s canonical correlation analysis (CCA) integration method. A total of 2000 features for anchoring (the ‘FindIntegrationAnchors’ function) and 30 dimensions for alignment (‘IntegrateData’) were used. Cell clustering was performed by ‘FindClusters’ function at a resolution of 0.8 and the top 20 genes were used to define cell identity. Dimensionality reduction was performed with ‘RunUMAP’ function and visualized by Uniform Manifold Approximation and Projection (UMAP). For subgroup cell clustering, cells of different types were extracted separately and clustered by their respective first 30 principal components (PCs) using different resolutions based on visual inspection._"

### "_Identification of signature genes_"
"_We applied the ‘FindAllMarkers’ function in Seurat to identify specific genes for each cell subset. For the selection of marker genes specific to each cell cluster/subset, we calculated the log2 fold change (log2FC) between two groups (a cell cluster/subset vs. other cells) using the ‘FindMarkers’ function with the Wilcoxon rank-sum test (default parameters)._"


## Experimental Background: 



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

  _split by sample_
  ![int_per_cluster_per_sample](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/369988e5-4501-4cc9-9f61-444b7204cc1c)

  _grouped by sample_

  ![int_umap_by_sample](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/cfcc2d56-b8ad-481a-9665-1dab8bffe73c)

  _grouped by patient_
  
![int_umap_by_patient](https://github.com/lenarayneallen/Seurat_Practice_Project/assets/124638335/b93ca81d-42c3-40f8-8705-fe17e0d84338)


### SingleR:


### Annotation and visualization:
