# Spatial-Transcriptomics-pipeline
Pipeline for processing high-resolution spatial transcriptomics data

In this github we present a pipeline to process high-resolution spatial transcriptomics data after obtaining the cell segmentation.
This pipeline consists on 3 steps: Clustering, Stitching, and Cell typing.

1. Clustering is done in R for this step you will need to install the following libraries:
dplyr, Seurat, patchwork, ggplot2, sctransform, SeuratData, SeuratDisk and SeuratObject.

2. Stitching is done in Python, for this step you will need to install the following libraries:
stitch2d, pandas, numpy, scanpy, anndata.

3. Cell typing is also done in Python, for this step you will need to install the following libraries:
pandas, numpy, scanpy, anndata, sklearn, seaborn, sklearn

Download the input data from the data folder. The gene expression file is called 'Lung5_Rep1_exprMat_file.csv' and the metadata file is called 'Lung5_Rep1_metadata_file.csv'. The expression matrix is zipped, you will have to unzip it before starting.

The images for the stitching section can be found in two folders:
1. CellComposite the morphology images stained with the morphology antibodies for membrane (CD298) and mitochondria, as well as DAPI stain for nuclei.
2. CellOverlay the cell segmented images.
*Both folders need to be stitched.

Download the results of the pipline from the following link https://drive.google.com/drive/folders/1jr-BvK8LVTIIf0HBPw30KyQLenentMVa?usp=share_link
For downloading the results from the clustering refer to the file 'Lung5-1-RNA-new-sct.h5ad'

The data was obtained from:
He, S., Bhatt, R., Brown, C. et al. High-plex imaging of RNA and proteins at subcellular resolution in fixed tissue by spatial molecular imaging. Nat Biotechnol 40, 1794–1806 (2022). https://doi-org.libweb.lib.utsa.edu/10.1038/s41587-022-01483-z

