# FIgure.4 and Suppl Figure.3E

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(RColorBrewer)
Set2 <- c("#66C2A5" ,"#FC8D62" ,"#8DA0CB", "#E78AC3" ,"#A6D854" ,"#FFD92F" ,"#E5C494" ,"#B3B3B3")

# data input - total
data_dir <- '/Volumes/Labo_SSD/#2023/#23-h5(Obesity_Visium)/23_31b_C3_epi/outs'
data.Cepi <- Load10X_Spatial(data.dir = data_dir, filename = "filtered_feature_bc_matrix.h5")
data.Cepi <- SCTransform(data.Cepi, assay = "Spatial", verbose = FALSE)
data.Cepi$orig.ident <- "Control"

data_dir <- '/Volumes/Labo_SSD/#2023/#23-h5(Obesity_Visium)/23_31b_O3_epi/outs'
data.Oepi <- Load10X_Spatial(data.dir = data_dir, filename = "filtered_feature_bc_matrix.h5")
data.Oepi <- SCTransform(data.Oepi, assay = "Spatial", verbose = FALSE)
data.Oepi$orig.ident <- "Obese"

DefaultAssay(data.merge) <- "SCT"
VariableFeatures(data.merge) <- c(VariableFeatures(data.Cepi), VariableFeatures(data.Oepi))
data <- RunPCA(data.merge, verbose = FALSE)

ElbowPlot(data, ndims = 30)

data <- FindNeighbors(data, reduction = "pca", dims = 1:25)

data <- FindClusters(data, verbose = FALSE, resolution = 0.5)
data <- RunUMAP(data, reduction = "pca", dims = 1:25)

p1 <- DimPlot(data, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(data, label = TRUE, label.size = 3, alpha = c(0.5, 1))
p1 + p2

DimPlot(data, reduction = "umap", label = TRUE)
SpatialDimPlot(data, label = TRUE, label.size = 3, alpha = c(0.5, 1))
DimPlot(data, reduction = "umap", split.by = "group", label = TRUE, pt.size = 0.8, label.size = 6, cols = "Set3", na.value = "grey50") + NoLegend() 
DimPlot(data, reduction = "umap",  label = TRUE, pt.size = 0.8, label.size = 6, cols = "Set3", na.value = "grey50")

data$group <- NA
data$group[data$orig.ident == "Control"] <- "Control"
data$group[data$orig.ident == "Obese"] <- "Obese"

setwd("/Volumes/Labo_SSD/#2023/#23-h5(Obesity_Visium)/240523_Visium_analysis")
saveRDS(data, file = "240524_data.merge_dims25_reso0.5.rds")  # 論文使用データとして保管するのはこのファイルとする

# Load data
data <- readRDS("/Volumes/Labo_SSD/#2023/#23-h5(Obesity_Visium)/240523_Visium_analysis/240524_data.merge_dims25_reso0.5.rds")

# Fig. 4A - Dimplot
DimPlot(data, reduction = "umap", split.by = "group", label = TRUE, pt.size = 0.8, label.size = 6, cols = "Set2", na.value = "grey50") + NoLegend() 

# Fig. 4B - SpatialDimplot
SpatialDimPlot(data, label = F, label.size = 3,  cols = "Set2", alpha = c(0.5, 1), image.alpha = 1)

# Fig. 4C - Stack violinplot
features <- c("Ptprc", "Cd68", "Itgax", "Col1a1", "Tgfb1", "Pecam1",  "Tek", "Vwf", "Cd34")
cols.features <- c("#8DA0CB", "#8DA0CB", "#8DA0CB", "#8DA0CB", "#8DA0CB", "#E5C494" ,"#E5C494" ,"#E5C494" ,"#B3B3B3")
VlnPlot(data, features, stack = TRUE, sort = F, flip = TRUE, cols = cols.features) + theme(legend.position = "none") + ggtitle("Identity on x-axis")

# Fig. 4D - Dimplot excluding vessels
data2 <- subset(data, idents = c(0,1,2,3,4,5))
SpatialDimPlot(data2, label = F, label.size = 3,  cols = "Set2", image.alpha = 1, alpha = c(0.5, 1))
DimPlot(data2, reduction = "umap", split.by = "group", label = TRUE, pt.size = 0.8, label.size = 6, cols = "Set2", na.value = "grey50") + NoLegend() 

# Fig. 4E - SpatialDimplot excluding vessels
SpatialDimPlot(data2, label = F, label.size = 3,  cols = "Set2", alpha = c(0.5, 1), image.alpha = 1)

# Fig. 4F - Dotplots
DotPlot(data2, features = c("Clec4e", "Wnt4", "Osm"), cols = c("gray","#B53471"), col.max = 1, dot.scale = 6)
DotPlot(data2, features = c("Itgb2"), cols = c("gray","#B53471"), col.min = 1, col.max = 1, dot.scale = 10)
DotPlot(data2, features = c("Fzd4", "Fzd1", "Lrp6", "Osmr", "Lifr", "Il6st", "Icam1"),
        col.max = 1, cols = c("#d1ccc0", "#009432"), dot.scale = 10)

# Suppl Fig.3E - Clec4e, Trem2, and Cd9, and Lyve1
DotPlot(data, features = c("Clec4e"), cols = c("gray","#B53471"), dot.scale = 4) 
DotPlot(data, features = c("Lyve1"), cols = c("gray","#B53471"), dot.scale = 4) 
DotPlot(data, features = c("Trem2", "Cd9"), cols = c("gray","#B53471"), dot.scale = 7) 
