#Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(SeuratData)
library(SeuratDisk)
library(SeuratObject)

#Load counts
counts<-read.csv(file='C:/Users/Wally/Documents/Research/Research-Lung/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/Lung5_Rep1_exprMat_file.csv',header=TRUE,row.names = 1)
#Load metadata
meta<-read.csv(file='C:/Users/Wally/Documents/Research/Research-Lung/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/Lung5_Rep1_metadata_file.csv',header=TRUE,row.names = 1)

#Create the Seurat object
so<-CreateSeuratObject(counts=counts,meta.data = meta,project='sct-lung5_1',min.cells=3,min.features = 0)


#Generate metrics and annotate mitochondrial counts
so[['percent.mt']]<-PercentageFeatureSet(so,pattern="^MT")

#Visualize QC metrics as a violin plot
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

#Filter to keep the cells with highest quality
so <- subset(so, subset = nCount_RNA > 250 & percent.mt < 5)

#Visualize QC metrics as a violin plot
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

#Normalization method
so <- SCTransform(so,assay = "RNA",  verbose = FALSE)
#Principal component analysis
so <- RunPCA(so, verbose = FALSE)
#UMAP
so <- RunUMAP(so, dims = 1:40, verbose = FALSE,n.neighbors=5, min.dist = 0.001)

#Find the nearest neighbors
so <- FindNeighbors(so, dims = 1:30, verbose = FALSE)
#Fin the clusters. Modify the resolution to get more or less clusters
so <- FindClusters(so, verbose = FALSE, resolution = 1.5)
#Plot the clusters
DimPlot(so, label = TRUE) + NoLegend()

#Save object
SaveH5Seurat(so, "C:/Users/Wally/Documents/Lung-Workshop/Lung5-1-RNA-new-sc1t.h5seurat")
Convert(source ='C:/Users/Wally/Documents/Lung-Workshop/Lung5-1-RNA-new-sc1t.h5seurat',dest= 'h5ad',assay = 'RNA') 
