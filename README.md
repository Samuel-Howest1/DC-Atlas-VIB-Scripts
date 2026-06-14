# DC-Atlas-VIB-Scripts

This is a repository for my internship in the VIB creating a DC cell Atlas. This repository contains all scripts related to this project.

This project use renv to keep package versions consistent between devices

# Content

The 3 main folders (Reports,Rshiny tool,Scripts) contain the following:

## Scripts

the Scripts folder contains 2 main sub-folders: "Processing of Data" and "Integration and Multimodel".

The "Processing of Data" folder contains all script used for the reports. It contains 3 scripts in total: a Pre-processing script, that filters out outliers cells, finds doublets and annotates cells. a Hashing script that identifies each hash-tag and find the doublets and negatives. and finally a Post-processing script that is used to manually correct unknown/suspicious cells, either correcting the annotation or removing doublets.

The "Integration and Multimodel" folder contains script to create the merged datasets and integrated them together for the cell atlas. the First script that should be used is the Merge_datasets.R to combine all Seurat objects form post-processing.

After which there is the option to use one of 3 integration methods marked as Integration\_"Method".R. Each one can be scored on the batch mixing and cell separation by the integration_Scoring.R script

Additionally there are some script that are either used to explore the data like Biomarkers.R that is used to find all biomarkes of the integrated Seuratobject for each cluster and Featue_plot_loop.R that is used the create feature plot for every orig.ident based on biomarker gene of Early or late Immature cDC1s and Early or Late Mature cDC1s.

Finally the Fix_ADT_Doublets script was used to correct an issue with the used dataset whereby the same ADT tag is used but with a different annotation creating doublet that cause for biases during the creation of the multimodal object.

The Last script that should be used is the Multi_model.R script that create a Seurat object with Weighted Nearest Neighbor technique to base cell placement in the UMAP on their Cite-seq profile.

## Reports

The Reports folder contains all the reports for the processed datasets. This means only what was done to processes the data. The datasets and the biological result are not present. Additionally not all are published and are kept private by the VIB until the paper is complete.

## Rshiny tool

Before the use of the Cell Atlas the script Create_Visualization_SeuratObject.R should be used. It remove all unnecessary data of the Seurat object that is not needed for the Cell Atlas thus creating a smaller Seurat Object for easier use of the Cell Atlas

The Rshiny tool contains 3 main pages used for analysis of the data. The first page is "Gene plots" Shows the whole umap graph gene expression based on the chosen gene and a violin plot for every cell type. Additionally it has a graph to look at the metadata column of the Seurat Object.

The second page is the "Cell Metadata" page which shows the gene expression, clusters, cell types, orig.idents, treatment, experiment metadata column for a chosen gene. You can filter based on all these metadata column how you want and the number of graphs dynamically change based on how you filter.

Finally there is a "Gene Comparison" page where you can look at the correlation between 2 genes.
