# Figure.2 and 3, Suppl Fig.3

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(CellChat)

# data input 
mBC <- load("Seurat.rda")

# add metadata
mBC$group <- NA
mBC$group[mBC$orig.ident == "A0301"] <- "1.SD"
mBC$group[mBC$orig.ident == "A0302"] <- "1.SD"
mBC$group[mBC$orig.ident == "A0303"] <- "1.SD"
mBC$group[mBC$orig.ident == "A0304"] <- "2.HFD8"
mBC$group[mBC$orig.ident == "A0305"] <- "2.HFD8"
mBC$group[mBC$orig.ident == "A0307"] <- "2.HFD8"
mBC$group[mBC$orig.ident == "A0308"] <- "3.HFD16"
mBC$group[mBC$orig.ident == "A0309"] <- "3.HFD16"
mBC$group[mBC$orig.ident == "A0310"] <- "3.HFD16"


# Fig. 2C - total Dimplot
DimPlot(mBC, reduction = "FItSNE", label = T, pt.size = 0.8, 
        cols = c(
          '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
          '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
          '#EE5A24','#009432','#0652DD','#9980FA','#833471',
          '#EA2027','#006266','#1B1464','#5758BB','#6F1E51',
          '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
          '#2c2c54','#474787','#aaa69d','#227093','#218c74',
          '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
          '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
        )) 

# Fig. 2D - total annotation Dimplot
new.cluster.ids <- c("Fibroblast", "Macrophage", "Macrophage", "Fibroblast", "Fibroblast", "Fibroblast",
                     "Macrophage", "Others", "Others", "Macrophage", "Macrophage", "Macrophage", 
                     "Fibroblast", "Macrophage", "Others", "Macrophage", "Fibroblast", "Fibroblast", 
                     "Fibroblast", "Others", "Macrophage", "Others", "Others", "Others", 
                     "Others", "Others", "Others", "Macrophage", "Macrophage", "Others", "Others", "Others", "Others")
names(new.cluster.ids) <- levels(mBC)
mBC_anno <- RenameIdents(mBC, new.cluster.ids)
DimPlot(mBC_anno, reduction = "FItSNE", label = FALSE, pt.size = 0.8, cols = c("#009432", "#B53471", "gray")) + NoLegend()

# Fig, 2E - total stack violinplot for annotation
features<- c("Ptprc","Adgre1", "Itgam" ,"Col1a1", "Pdgfra", "Acta2")
VlnPlot(mBC, features, stack = TRUE, sort = TRUE, flip = TRUE, 
        cols =  c("#B53471", "#B53471", "#B53471", "#009432", "#009432", "#009432"))  + theme(legend.position = "none")

# Fig. 2G - macrophage Dimplot
DimPlot(mBC, reduction = "FItSNE", split.by = "group", label = TRUE, pt.size = 0.8, label.size = 6, 
        cols = c("#a75629" ,"#54b0e4" ,"#222f75", "#1c9e78" ,"#b2df8b" ,"#e3be00" ,"#fc9a9a" ,"#e7298a", "#910141", "#04cdd1", "#a6cee3"), , na.value = "grey50") + NoLegend() 
DimPlot(mBC, reduction = "FItSNE",  label = TRUE, pt.size = 0.8, label.size = 6, 
        cols = c("#a75629" ,"#54b0e4" ,"#222f75", "#1c9e78" ,"#b2df8b" ,"#e3be00" ,"#fc9a9a" ,"#e7298a", "#910141", "#04cdd1", "#a6cee3"), na.value = "grey50")

# Fig. 2I - fibroblast Dimplot
DimPlot(mBC, reduction = "FItSNE", split.by = "group", label = TRUE, pt.size = 0.8, label.size = 6, 
        cols = c("#e41a1b" ,"#377eb8" ,"#4eaf4a", "#984ea3" ,"#f29401" ,"#f781c1" ,"#bc9dcc"), , na.value = "grey50") + NoLegend() 
DimPlot(mBC, reduction = "FItSNE", label = TRUE, pt.size = 0.8, label.size = 6, 
        cols = c("#e41a1b" ,"#377eb8" ,"#4eaf4a", "#984ea3" ,"#f29401" ,"#f781c1" ,"#bc9dcc"), , na.value = "grey50") 

# Fig. 2K - cluster marker genes heatmap (Fibro)
cols <- c("#e41a1b" ,"#377eb8" ,"#4eaf4a", "#984ea3" ,"#f29401" ,"#f781c1" ,"#bc9dcc")
top2 <- c("Abca8a","Abca8b","Efhd1","Sbsn","mt-Rnr2","mt-Atp8",
          "Mdk","C4b","Duoxa1","Aif1l","Tgfbi","Cpe","Kcnq1ot1","Gm29055")
DoHeatmap(mBC, features = top2, group.colors = cols, label = F) + NoLegend() 
+   scale_fill_gradient2(low = "#227093",high = "#EE5A24", midpoint=0)

# Fig. 2L - fibroblast cluster marker genes (FALSE)
ALL.markers <- FindAllMarkers(mBC, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
ALL.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(ALL.markers, file = "WTA2_Fibro_ALL.markers_FALSE.txt", row.names = T, col.names = T, sep = "\t")　

# CellChat analysis - load data
# SD = cellchat.list[[2]]
# HF8W = cellchat.list[[1]]
# HF16W = cellchat.list[[3]]
# 1 = Fibro_0, 2 = Fibro_1, 3 = Fibro_2, 4 = Fibro_3, 5 = Fibro_4, 6 = Fibro_5, 7 = Fibro_6
# 8 = Mac_0,  9 = Mac_1, 10 = Mac_2, 11 = Mac_3,  12 = Mac_4, 13 = Mac_5, 
# 14 = Mac_6,  15 = Mac_7, 16 = Mac_8, 17 = Mac_9,  18 = Mac_ten

cellchat.list <- load("cellchat.rda")

# Fig.3A - CellChat circleplot  --> 松島研出力データ 

# Fig.3B - macrophage Clec4e Violin plot
VlnPlot(mBC, features = c("Clec4e","Osm"), split.by = "group",cols = c( "#FEEFBD", '#FF90BC', "#E82C9B"))

# Fig.3C - Macrophage subclusters #4 and #10 markers
features = c("Ly6c2", "Ceacam1", "Cx3cr1", "Ccr2", "Fn1","Cd226")
DotPlot(mBC, features = features, cols = c( "gray", "#B53471"), dot.scale = 6) + coord_flip()

# Fig.3E - CellChat bubble plot
netVisual_bubble(cellchat.list[[3]], sources.use = c(12, 18),　 
                 targets.use = c(4),　　　　　　　　　　 
                 signaling =  unique(c("ITGAL-ITGB2", "OSM", "WNT")),　　　　　
                 remove.isolate = F, return.data=T, font.size = 8)

# Fig.3F - Cellchat selected genes FeaturePlot
FeaturePlot(mBC, features = c("Osm","Wnt4","Itgb2"), cols = c( "gray", "#B53471"))
FeaturePlot(mBC, features = c("Osmr","Fzd4","Icam1"), cols = c( "gray", "#009432"))


# Suppl Fig. 3B -  macrophage LAM score
mBC <- AddModuleScore(obj = mBC, 
                      features = list(c("Trem2","Cd9", 
                                        "Cd36", "Fabp4", "Fabp5",
                                        "Plin2", "Lipa", "Lpl",
                                        "Lgals1", "Lgals3", "Ctsb", "Ctsl",
                                        "Mmp12", "Gpnmb", "Spp1")), 
                      name = "Lipid.associated_macrophage_score")
mBC <- AddModuleScore(obj = mBC, 
                      features = list(c("Lyve1","Cd163","Cd209f")), 
                      name = "Perivascular_macrophage_score")
DotPlot(mBC,features = c("Clec4e", "Lipid.associated_macrophage_score1","Perivascular_macrophage_score1"), 
        cols = c( "gray", "#B53471")) + coord_flip()

# Suppl Fig. 3C -  macrophage Clec4e, Lyve1 and Trem2
FeaturePlot(mBC, features = c("Clec4e","Trem2"), cols = c( "gray", "#C80000", "#3BC83B"), blend =  T)
FeaturePlot(mBC, features = c("Clec4e","Lyve1"), cols = c( "gray", "#C80000", "#3BC83B"), blend =  T)

# Suppl Fig. 3D -  LAM to Fibro#3 signaling plot
cols <- c("gray", "gray", "gray", "#984ea3", "gray", 
          "gray", "gray", "gray", "#54b0e4", "gray", 
          "#1c9e78", "gray", "gray", "gray", "gray",
          "gray", "gray", "gray")　　　　　　　　　　# Mac1とMac3のハイライト

pathway.name <- "PDGF"    # Pdgfb-Pdgfra, Pdgfb-Pdgfrb
pathway.name <- "TWEAK"   
pathway.name <- "VCAM"    # Itga4-receptor, Itga9-receptor

netAnalysis_contribution(cellchat.list[[3]], signaling = pathway.name)
pairLR <- extractEnrichedLR(cellchat.list[[3]], signaling = pathway.name, geneLR.return = FALSE)
LR.show <-  pairLR[4,] # show one ligand-receptor pair

netVisual_individual(cellchat.list[[3]],    # 全てのクラスターを対象とする方が誤解がなくて良い
                     color.use = cols,
                     signaling = pathway.name, pairLR.use = LR.show, layout = "circle")

