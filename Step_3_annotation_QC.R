###############################################################################################
# Pre-processing for data used to generate Figure 7, Figure S7 and S8 from Wang et al. Nature Communications 2022.

# This script is: STEP 3 of 4

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################


library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)


setwd("")

Tau_integrated <- readRDS("Tau_integrated_PCA_0.1.rds")

Tau_integrated <- subset(Tau_integrated, idents = c("22","23"), invert = TRUE)

Tau_integrated <- RenameIdents(Tau_integrated,
                                 `0` = "oligodendrocytes", `1`="inhibitory neurons", `2`="excitatory neurons", `3`="excitatory neurons",
                                 `4`="excitatory neurons", `5`="astrocytes", `6`="unknown", `7`="excitatory neurons",
                                 `8`="excitatory neurons", `9`="inhibitory neurons", `10`="inhibitory neurons", `11`="microglia",
                                 `12`="OPCs", `13`="inhibitory neurons", `14`="inhibitory neurons", `15`="excitatory neurons",
                                 `16`="excitatory neurons", `17`="excitatory neurons", `18`="excitatory neurons", `19`="vascular cells",
                                 `20`="vascular cells", `21`="vascular cells"
)

pdf("Tau_umap_annotation_withLabel.pdf", width=6.5, height=4.5)
DimPlot(Tau_integrated, reduction = 'umap', label = T)
dev.off()

Tau_integrated$celltype.orig.ident <- paste(Idents(Tau_integrated), Tau_integrated$orig.ident, sep = "_")
Tau_integrated$celltype <- Idents(Tau_integrated)

saveRDS(Tau_integrated, file = "Tau_integrated_ready_4_DEGs.rds")

Idents(Tau_integrated)
Idents(Tau_integrated) <- "Condition"
Tau_integrated <- RenameIdents(Tau_integrated,
                               `IKKbCtrl` = "ikbkb+/+", `IKKbWT`="ikbkb+/+", `IKKbCtrl;P301S+`="ikbkb+/+;P301S+", `IKKbWT;P301S+`="ikbkb+/+;P301S+",
                               `IKKbKO;P301S+`="ikbkb-/-;P301S+", `IKKbCA;P301S+`="ikbkbCA;P301S+"
)

Tau_integrated$Genotype <- Idents(Tau_integrated)

# calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(Tau_integrated$Genotype,Tau_integrated$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
        geom_bar(stat="identity")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        xlab("Genotype")+
        ylab("Cell type ratio per genotype") + RotatedAxis()

ggsave("genotype_celltype_distribution.pdf",plot=last_plot(),path="~/Desktop/data_analysis/Chao_Tau_V7/4genotypes",
       width=4,height=4,units="in")


Idents(Tau_integrated) <- "orig.ident"
Tau_integrated <- RenameIdents(Tau_integrated,
                               `IKKbCtrl_1` = "ikbkb+/+_1", `IKKbCtrl_2`="ikbkb+/+_2", `IKKbWT_1`="ikbkb+/+_3", `IKKbWT_2`="ikbkb+/+_4",
                               `IKKbCtrl;P301S+_1`="ikbkb+/+;P301S+_1", `IKKbCtrl;P301S+_2`="ikbkb+/+;P301S+_2",`IKKbWT;P301S+_1`="ikbkb+/+;P301S+_3",
                               `IKKbKO;P301S+_1`="ikbkb-/-;P301S+_1", `IKKbKO;P301S+_2`="ikbkb-/-;P301S+_2",
                               `IKKbCA;P301S+_1`="ikbkbCA;P301S+_1", `IKKbCA;P301S+_2`="ikbkbCA;P301S+_2"
)

Tau_integrated$Sample <- Idents(Tau_integrated)

# QC
pdf("Tau_integrated_QC.pdf", width=10, height=4.5)
VlnPlot(Tau_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "Sample", ncol = 3, pt.size=0, idents=NULL)
dev.off()
pdf("Tau_integrated_UMI_gene.pdf", width=6.5, height=4.5)
FeatureScatter(Tau_integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample", pt.size=0.5)
dev.off()
pdf("Tau_integrated_UMI_mt.pdf", width=6.5, height=4.5)
FeatureScatter(Tau_integrated, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample", pt.size=0.5)
dev.off()

# calculate ratio of each cell type in each sample
a<-as.data.frame(table(Tau_integrated$Sample,Tau_integrated$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
        geom_bar(stat="identity")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        xlab("Sample")+
        ylab("Cell type ratio per sample") + RotatedAxis()

ggsave("sample_celltype_distribution.pdf",plot=last_plot(),path="~/Desktop/data_analysis/Chao_Tau_V7/4genotypes",
       width=5,height=4,units="in")

Idents(Tau_integrated) <- "celltype"
pdf("Tau_DotPlot_annotation.pdf", width=9, height=3.5)
DotPlot(Tau_integrated, features = c("Plp1", "Mbp", "Mobp", "Gad1", "Gad2","Slc17a7", "Nrgn", "Pla2g7","Gli3","Clu", 
                           "Cx3cr1", "P2ry12", "Csf1r","Vcan", "Pdgfra", "Cped1", "Flt1")) + RotatedAxis()
dev.off()


pdf("Tau_umap_annotation_split_Condition.pdf", width=16, height=4)
DimPlot(Tau_integrated, reduction = 'umap', split.by = "Condition", label = TRUE)
dev.off()


markers.to.plot <- c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Gad1", "Gad2", "Clu", "Aldoc", 
                     "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Scrg1", "Pdgfra", "Vtn", "Igfbp7",
                     "Bnc2", "Slc47a1", "Ttr")
pdf("DotPlot_top2-3.pdf", width=12, height=8)
DotPlot(Tau_integrated, features = rev(markers.to.plot), cols = c("orange", "plum", "pink", "grey"), dot.scale = 8, 
        split.by = "Condition") + RotatedAxis()
dev.off()



