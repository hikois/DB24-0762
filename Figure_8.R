# Figure.8

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(DESeq2)
library(tidyverse)

# data input - total
mBC <- readRDS("WTA3total")

# data input - Mac label transferred data --> Mac_WTA3
Mac_WTA3 <- load("MacSubclusteringLabelAdd.rda")

# data input - Fibro label transferred data --> seu1
seu1 <- load("fibro_WTA3LabelTransfer.rda")

# Fig. 8B - genotype annotation
DimPlot(mBC, label = F, pt.size = 0.8, group.by = "group", cols = c('#aaa69d', '#EE5A24')) 

# Fig. 8C - Mac Umap plot
DimPlot(Mac_WTA3, label = TRUE, pt.size = 0.8, label.size = 6, 
        cols = c("#a75629" ,"#54b0e4" ,"#222f75", "#1c9e78" ,"#b2df8b" ,"#fc9a9a" ,"#e7298a", "#910141", "#04cdd1", "#a6cee3"), , na.value = "grey50")  

# Fig. 8C - Mac cluster composition
cluster_ids <- Idents(Mac_WTA3)
conditions <- Mac_WTA3@meta.data$orig.ident
cluster_condition_counts <- table(cluster_ids, conditions)
cluster_condition_counts_df <- as.data.frame(cluster_condition_counts)
colnames(cluster_condition_counts_df) <- c("Cluster", "Condition", "Cell_Count")
print(cluster_condition_counts_df)
write.csv(cluster_condition_counts_df, "Mac_label_cluster_condition_counts.csv", row.names = FALSE)

# Fig. 8D - Fibro Umap plot
DimPlot(seu1, label = TRUE, pt.size = 0.8, label.size = 6, 
        cols = c("#e41a1b" ,"#377eb8" ,"#4eaf4a", "#984ea3" ,"#f29401" ,"#f781c1" ,"#bc9dcc"), , na.value = "grey50")  

# Fig. 8D - Fibro cluster composition
cluster_ids <- Idents(seu1)
conditions <- seu1@meta.data$orig.ident
cluster_condition_counts <- table(cluster_ids, conditions)
cluster_condition_counts_df <- as.data.frame(cluster_condition_counts)
colnames(cluster_condition_counts_df) <- c("Cluster", "Condition", "Cell_Count")
print(cluster_condition_counts_df)
write.csv(cluster_condition_counts_df, "Fibro_label_cluster_condition_counts.csv", row.names = FALSE)

# Fig. 8E Fibro cluster#3 pseudobulk analysis
cell <- gsub('*._', '', colnames(seu1))
seu1[["cell"]] <- cell

seu1$genotype <- NA  
seu1$genotype[seu1$orig.ident == "A0309"] <- "WT"
seu1$genotype[seu1$orig.ident == "A0310"] <- "WT"
seu1$genotype[seu1$orig.ident == "A0311"] <- "KO"
seu1$genotype[seu1$orig.ident == "A0312"] <- "KO"

seu1$sample <- paste0(seu1$genotype, seu1$cell)

cts <- AggregateExpression(seu1, 
                           group.by = c("WTA3LabelTransfer", "sample"),
                           assays = 'RNA',　　　　　　　　　
                           slot = "counts",
                           return.seurat = FALSE)
cts <- cts$RNA
cts[1:3,1:3]　

cts <- t(cts)　
cts <- as.data.frame(cts)

splitRows <- gsub('_.*', '', rownames(cts))
cts <- split.data.frame(cts,
                        f = factor(splitRows))

cts$'g3'[1:10, 1:10]  #Fibro #3

cts <- lapply(cts, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)})

cts$'g3'[1:3, 1:3]

counts_cluster <- cts$'g3'　　
counts_cluster[1:10, 1:10]

colData <- data.frame(samples = colnames(counts_cluster))
colData <- colData %>%
  mutate(condition = ifelse(grepl('KO', samples), 'KO', 'WT')) %>%
  column_to_rownames(var = 'samples')

dds <- DESeqDataSetFromMatrix(countData = counts_cluster,　
                              colData = colData,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10,]
dds <- DESeq(dds)　
resultsNames(dds)
res <- results(dds, name = "condition_WT_vs_KO")
write.table(res, file = "WT-vs-KO_FibroCluster3_result.txt", row.names = T, col.names = T, sep = "\t")

