###############################################################################################
# Pre-processing for data used to generate Figure 7, Figure S7 and S8 from Wang et al. Nature Communications 2022.

# This script is: STEP 1 of 4

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

#set working directory ====
setwd("")

#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
LG343_80.counts <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
LG343_80 <- CreateSeuratObject(counts = LG343_80.counts, project = "LG343_80", min.cells = 3, min.features = 200)
LG343_80[["Condition"]] = c('IKKbKO; Ctrl')

#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG343_80[["percent.mt"]] <- PercentageFeatureSet(object = LG343_80, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG343_80

VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)


#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
all
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
ElbowPlot(all)
all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

DimPlot(all, reduction = "umap", label = T)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG343_80_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
saveRDS(all,"LG343_80_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG343_80_PCA.rds")
all

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*9629) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_732", sct = FALSE)
#visualizing clusters and multiplet cells====
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
all <- RunTSNE(object = all, dims = 1:15, do.fast=TRUE, perplexity=100, max.iter=10000)
DimPlot(object = all, reduction = 'umap')

#calculate nearest neighbors
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
#DimPlot(object = all, reduction = 'tsne')
DimPlot(all, reduction = "umap", label = T)

Idents(object = all) <- "DF.classifications_0.25_0.005_732" #visualizing the singlet vs doublet cells
DimPlot(object = all, reduction = 'umap')

all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
singlets<-saveRDS(singlets,"LG343_80_singlets.rds")
singlets<-readRDS("LG343_80_singlets.rds")

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
ElbowPlot(singlets,ndims=50)
singlets <- RunTSNE(object = singlets, dims = 1:15, do.fast=TRUE, perplexity=100, max.iter=10000)
DimPlot(object = singlets, reduction = 'umap')

#calculate nearest neighbors
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
DimPlot(object = singlets, reduction = 'umap')

#finding markers of your clusters and annotate cells====
singlets.markers <- FindAllMarkers(object = singlets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(singlets.markers,"LG343_80_cluster_markers_singlets.csv")

singlets<-saveRDS(singlets,"LG343_80_singlets_PCA.rds")

###################################################
#for loading Cell Ranger counts:
LG343_78.counts <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
LG343_78 <- CreateSeuratObject(counts = LG343_78.counts, project = "IKKbCtrl_2", min.cells = 3, min.features = 200)
LG343_78[["Condition"]] = c('IKKbCtrl')
rm(LG343_78.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG343_78[["percent.mt"]] <- PercentageFeatureSet(object = LG343_78, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG343_78
pdf("LG343_78_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG343_78_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG343_78_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG343_78_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG343_78_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG343_78_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG343_78_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*10267) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_780", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG343_78_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG343_78_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_780" #visualizing the singlet vs doublet cells
pdf("LG343_78_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"LG343_78_after_doublet_detection.rds")

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG343_78_singlets.rds")

singlets<-readRDS("LG343_78_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG343_78_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG343_78_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG343_78_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG343_78_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG343_78_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG343_78_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG343_78_singlets_PCA.rds")
#######################################
#for loading Cell Ranger counts:
LG345_83.counts <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
LG345_83 <- CreateSeuratObject(counts = LG345_83.counts, project = "IKKb_WT_1", min.cells = 3, min.features = 200)
LG345_83[["Condition"]] = c('IKKb_WT')
rm(LG345_83.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG345_83[["percent.mt"]] <- PercentageFeatureSet(object = LG345_83, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG345_83
pdf("LG345_83_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG345_83_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG345_83_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG345_83_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG345_83_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG345_83_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG345_83_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9321) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_643", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG345_83_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG345_83_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_643" #visualizing the singlet vs doublet cells
pdf("LG345_83_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_552" #visualizing the singlet vs doublet cells
pdf("LG345_83_3_UMAP_singlets_doublets_3.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"LG345_83_after_doublet_detection.rds")

all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG345_83_singlets.rds")

singlets<-readRDS("LG345_83_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG345_83_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG345_83_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_83_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG345_83_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG345_83_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_83_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG345_83_singlets_PCA.rds")
###########################################
#for loading Cell Ranger counts:
LG345_131.counts <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
LG345_131 <- CreateSeuratObject(counts = LG345_131.counts, project = "IKKb_WT_2", min.cells = 3, min.features = 200)
LG345_131[["Condition"]] = c('IKKb_WT')
rm(LG345_131.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG345_131[["percent.mt"]] <- PercentageFeatureSet(object = LG345_131, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG345_131
pdf("LG345_131_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG345_131_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG345_131_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG345_131_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG345_131_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG345_131_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG345_131_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9839) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_679", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG345_131_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG345_131_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_679" #visualizing the singlet vs doublet cells
pdf("LG345_131_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"LG345_131_after_doublet_detection.rds")

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG345_131_singlets.rds")

singlets<-readRDS("LG345_131_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG345_131_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG345_131_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_131_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG345_131_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG345_131_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_131_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG345_131_singlets_PCA.rds")
########################################################
#for loading Cell Ranger counts:
LG343_107.counts <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
LG343_107 <- CreateSeuratObject(counts = LG343_107.counts, project = "LG343_107", min.cells = 3, min.features = 200)
LG343_107[["Condition"]] = c('IKKbKO; Ctrl; PS19/+')

#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG343_107[["percent.mt"]] <- PercentageFeatureSet(object = LG343_107, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG343_107

VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
all
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
ElbowPlot(all)
all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

DimPlot(all, reduction = "umap", label = T)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG343_107_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
saveRDS(all,"LG343_107_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG343_107_PCA.rds")
all

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*11255) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_855", sct = FALSE)
#visualizing clusters and multiplet cells====
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
all <- RunTSNE(object = all, dims = 1:15, do.fast=TRUE, perplexity=100, max.iter=10000)
DimPlot(object = all, reduction = 'umap')

#calculate nearest neighbors
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
#DimPlot(object = all, reduction = 'tsne')
DimPlot(all, reduction = "umap", label = T)

Idents(object = all) <- "DF.classifications_0.25_0.005_855" #visualizing the singlet vs doublet cells
DimPlot(object = all, reduction = 'umap')

all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
singlets<-saveRDS(singlets,"LG343_107_singlets.rds")
singlets<-readRDS("LG343_107_singlets.rds")

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
ElbowPlot(singlets,ndims=50)
singlets <- RunTSNE(object = singlets, dims = 1:15, do.fast=TRUE, perplexity=100, max.iter=10000)
DimPlot(object = singlets, reduction = 'umap')

#calculate nearest neighbors
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
DimPlot(object = singlets, reduction = 'umap')

#finding markers of your clusters and annotate cells====
singlets.markers <- FindAllMarkers(object = singlets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(singlets.markers,"LG343_107_cluster_markers_singlets.csv")

singlets<-saveRDS(singlets,"LG343_107_singlets_PCA.rds")
##############################################################
#for loading Cell Ranger counts:
LG343_79.counts <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
LG343_79 <- CreateSeuratObject(counts = LG343_79.counts, project = "IKKbCtrl;P301S+_2", min.cells = 3, min.features = 200)
LG343_79[["Condition"]] = c('IKKbCtrl;P301S+')
rm(LG343_79.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG343_79[["percent.mt"]] <- PercentageFeatureSet(object = LG343_79, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG343_79
pdf("LG343_79_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG343_79_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG343_79_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG343_79_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG343_79_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG343_79_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG343_79_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.091*12207) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_1111", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG343_79_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG343_79_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.01_1111" #visualizing the singlet vs doublet cells
pdf("LG343_79_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"LG343_79_after_doublet_detection.rds")

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG343_79_singlets.rds")

singlets<-readRDS("LG343_79_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG343_79_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG343_79_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG343_79_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG343_79_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG343_79_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG343_79_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG343_79_singlets_PCA.rds")
#########################################################
#for loading Cell Ranger counts:
LG345_98.counts <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
LG345_98 <- CreateSeuratObject(counts = LG345_98.counts, project = "IKKbWT;P301S+_1", min.cells = 3, min.features = 200)
LG345_98[["Condition"]] = c('IKKbWT;P301S+')
rm(LG345_98.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG345_98[["percent.mt"]] <- PercentageFeatureSet(object = LG345_98, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG345_98
pdf("LG345_98_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG345_98_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG345_98_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG345_98_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG345_98_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG345_98_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG345_98_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.084*11145) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_936", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG345_98_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG345_98_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_936" #visualizing the singlet vs doublet cells
pdf("LG345_98_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"LG345_98_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG345_98_singlets.rds")

singlets<-readRDS("LG345_98_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG345_98_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG345_98_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_98_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG345_98_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG345_98_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_98_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG345_98_singlets_PCA.rds")
#######################################################
#for loading Cell Ranger counts:
LG343_126.counts <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
LG343_126 <- CreateSeuratObject(counts = LG343_126.counts, project = "IKKbKO;P301S+_1", min.cells = 3, min.features = 200)
LG343_126[["Condition"]] = c('IKKbKO;P301S+')
rm(LG343_126.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG343_126[["percent.mt"]] <- PercentageFeatureSet(object = LG343_126, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG343_126
pdf("LG343_126_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG343_126_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG343_126_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG343_126_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG343_126_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG343_126_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG343_126_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*10578) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_804", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG343_126_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG343_126_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_804" #visualizing the singlet vs doublet cells
pdf("LG343_126_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"LG343_126_after_doublet_detection.rds")

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG343_126_singlets.rds")

singlets<-readRDS("LG343_126_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG343_126_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG343_126_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG343_126_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG343_126_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG343_126_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG343_126_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG343_126_singlets_PCA.rds")
#######################################################
#for loading Cell Ranger counts:
LG343_108.counts <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
LG343_108 <- CreateSeuratObject(counts = LG343_108.counts, project = "LG343_108", min.cells = 3, min.features = 200)
LG343_108[["Condition"]] = c('IKKbKO; Cre/+; PS19/+')

#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG343_108[["percent.mt"]] <- PercentageFeatureSet(object = LG343_108, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG343_108

VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)


#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
all
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
ElbowPlot(all)
all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

DimPlot(all, reduction = "umap", label = T)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG343_108_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
saveRDS(all,"LG343_108_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG343_108_PCA.rds")
all

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7500) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_405", sct = FALSE)
#visualizing clusters and multiplet cells====
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
all <- RunTSNE(object = all, dims = 1:15, do.fast=TRUE, perplexity=100, max.iter=10000)
DimPlot(object = all, reduction = 'umap')

#calculate nearest neighbors
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
#DimPlot(object = all, reduction = 'tsne')
DimPlot(all, reduction = "umap", label = T)

Idents(object = all) <- "DF.classifications_0.25_0.005_405" #visualizing the singlet vs doublet cells
DimPlot(object = all, reduction = 'umap')

all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG343_108_singlets.rds")

singlets<-readRDS("LG343_108_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG343_108_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG343_108_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG343_108_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG343_108_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG343_108_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG343_108_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG343_108_singlets_PCA.rds")

##############################################################
#for loading Cell Ranger counts:
LG345_84.counts <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
LG345_84 <- CreateSeuratObject(counts = LG345_84.counts, project = "IKKbCA;P301S+_1", min.cells = 3, min.features = 200)
LG345_84[["Condition"]] = c('IKKbCA;P301S+')
rm(LG345_84.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG345_84[["percent.mt"]] <- PercentageFeatureSet(object = LG345_84, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG345_84
pdf("LG345_84_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG345_84_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG345_84_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG345_84_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG345_84_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG345_84_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG345_84_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*10976) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_834", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG345_84_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG345_84_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_834" #visualizing the singlet vs doublet cells
pdf("LG345_84_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"LG345_84_after_doublet_detection.rds")

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG345_84_singlets.rds")

singlets<-readRDS("LG345_84_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG345_84_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG345_84_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_84_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG345_84_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG345_84_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_84_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG345_84_singlets_PCA.rds")
######################################################
#for loading Cell Ranger counts:
LG345_133.counts <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
LG345_133 <- CreateSeuratObject(counts = LG345_133.counts, project = "IKKbCA;P301S+_2", min.cells = 3, min.features = 200)
LG345_133[["Condition"]] = c('IKKbCA;P301S+')
rm(LG345_133.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG345_133[["percent.mt"]] <- PercentageFeatureSet(object = LG345_133, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG345_133
pdf("LG345_133_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG345_133_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 1000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG345_133_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG345_133_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG345_133_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG345_133_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG345_133_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.084*11293) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_949", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG345_133_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG345_133_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_949" #visualizing the singlet vs doublet cells
pdf("LG345_133_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"LG345_133_after_doublet_detection.rds")

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG345_133_singlets.rds")

singlets<-readRDS("LG345_133_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG345_133_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG345_133_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_133_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG345_133_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG345_133_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_133_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG345_133_singlets_PCA.rds")

