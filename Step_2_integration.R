###############################################################################################
# Pre-processing for data used to generate Figure 7, Figure S7 and S8 from Wang et al. Nature Communications 2022.

# This script is: STEP 2 of 4

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("")
LG343_80 <- readRDS(file = "LG343_80_singlets_PCA.rds")
LG343_78 <- readRDS(file = "LG343_78_singlets_PCA.rds")
LG343_107 <- readRDS(file = "LG343_107_singlets_PCA.rds")
LG343_126 <- readRDS(file = "LG343_126_singlets_PCA.rds")
LG343_79 <- readRDS(file = "LG343_79_singlets_PCA.rds")
LG343_108 <- readRDS(file = "LG343_108_singlets_PCA.rds")
LG345_83 <- readRDS(file = "LG345_83_singlets_PCA.rds")
LG345_131 <- readRDS(file = "LG345_131_singlets_PCA.rds")
LG345_98 <- readRDS(file = "LG345_98_singlets_PCA.rds")
LG345_84 <- readRDS(file = "LG345_84_singlets_PCA.rds")
LG345_133 <- readRDS(file = "LG345_133_singlets_PCA.rds")


setwd("")
Ctrl <- c(LG343_80, LG343_78, LG345_83, LG345_131)
Ctrl_P301S <- c(LG343_107, LG343_79, LG345_98)
IKKbKO_P301S <- c(LG343_126, LG343_108)
IKKbCA_P301S <- c(LG345_84, LG345_133)

anchors_Ctrl <- FindIntegrationAnchors(object.list = Ctrl, dims = 1:30)
Ctrl_integrated <- IntegrateData(anchorset = anchors_Ctrl, dims = 1:30)
anchors_Ctrl_P301S <- FindIntegrationAnchors(object.list = Ctrl_P301S, dims = 1:30)
Ctrl_P301S_integrated <- IntegrateData(anchorset = anchors_Ctrl_P301S, dims = 1:30)
anchors_IKKbKO_P301S <- FindIntegrationAnchors(object.list = IKKbKO_P301S, dims = 1:30)
IKKbKO_P301S_integrated <- IntegrateData(anchorset = anchors_IKKbKO_P301S, dims = 1:30)
anchors_IKKbCA_P301S <- FindIntegrationAnchors(object.list = IKKbCA_P301S, dims = 1:30)
IKKbCA_P301S_integrated <- IntegrateData(anchorset = anchors_IKKbCA_P301S, dims = 1:30)

LG343_all <- c(Ctrl_integrated,Ctrl_P301S_integrated,IKKbKO_P301S_integrated, IKKbCA_P301S_integrated)
anchors_all <- FindIntegrationAnchors(object.list = LG343_all, dims = 1:30)
Tau_integrated <- IntegrateData(anchorset = anchors_all, dims = 1:30)


pdf("Tau_QC.pdf", width=12, height=4)
Idents(Tau_integrated) <- "orig.ident"
VlnPlot(object = Tau_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
Idents(Tau_integrated) <- "Condition"
VlnPlot(object = Tau_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

DefaultAssay(Tau_integrated) <- 'integrated'

Tau_integrated <- ScaleData(Tau_integrated, verbose = FALSE)
Tau_integrated <- RunPCA(Tau_integrated, features = VariableFeatures(object = Tau_integrated), verbose = FALSE)

Tau_integrated <- FindNeighbors(Tau_integrated, dims = 1:20)
Tau_integrated <- FindClusters(Tau_integrated, resolution = 0.1)
Tau_integrated <- RunUMAP(Tau_integrated, dims = 1: 20)

str(Tau_integrated)

DefaultAssay(Tau_integrated) <- 'RNA'
Tau_integrated <- NormalizeData(Tau_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
Tau_integrated <- ScaleData(Tau_integrated, features = rownames(Tau_integrated))

pdf("Tau_integrated_umap.pdf", width=5, height=4)
DimPlot(Tau_integrated, reduction = 'umap', label = T)
dev.off()
pdf("Tau_integrated_umap_split_individual.pdf", width=10, height=9)
DimPlot(Tau_integrated, reduction = "umap", split.by = "orig.ident", label = T, ncol = 4)
dev.off()
pdf("Tau_integrated_umap_split_Condition.pdf", width=9, height=8)
DimPlot(Tau_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()

saveRDS(Tau_integrated, file = 'Tau_integrated_PCA_0.1.rds')

DefaultAssay(Tau_integrated) <- 'RNA'

Tau_markers <- FindAllMarkers(Tau_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(Tau_markers, "Tau_markers.csv")

DefaultAssay(Tau_integrated) <- 'RNA'
pdf("Tau_umap_test.pdf", width=8, height=6)
DimPlot(Tau_integrated, reduction = 'umap', label = T)
dev.off()

#Add marker genes

sig_EN<-c("Snap25","Slc17a7", "Nrgn","Gad1", "Gad2","Plp1", "Mbp", "Mobp", "Clu", "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r",
          "Pdgfra", "Vcan","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("Tau_annotation_combine.pdf", width=10, height=5)
DotPlot(object = Tau_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()
