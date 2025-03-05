library(SingleR)
library(Seurat)
library(dplyr)
library(pheatmap)
library(ggplot2)


seu = readRDS("WTA3total")
load("./Fibro_results_Seurat.rda")
mBC = UpdateSeuratObject(mBC)

seu1 = subset(seu, idents=c(1,4))

singler_res = SingleR(
  test = GetAssayData(seu1, assay = 'RNA', slot = 'data'),
  ref = mBC@assays$RNA@data,
  labels = mBC@meta.data$res.0.6)

singler_res = as.data.frame(singler_res)

hoge = singler_res$pruned.labels
hoge = as.numeric(hoge)
names(hoge)=rownames(singler_res)
hoge = factor(hoge, levels=c(0:6))
seu1 = AddMetaData(seu1, metadata=hoge, col.name="WTA3LabelTransfer")
Idents(seu1)="WTA3LabelTransfer"
p = DimPlot(seu1, pt.size=0.3, label=T)

save(seu1, file="fibro_WTA3LabelTransfer.rda")

mBC = seu1

tmp = table(mBC@active.ident, mBC@meta.data$orig.ident)

tmp = as.data.frame(tmp)

tmp1 = tmp %>% tidyr::spread(key = Var2, value=Freq)
colnames(tmp1)[1]="Seurat_Clusters"
hoge = c("doublet", "not-detected", "not")
tmp1 = tmp1[,!colnames(tmp1) %in% hoge]
tmp_cellcount = tmp1

rownames(tmp1) = tmp1[,1]
tmp1 = tmp1[,2:ncol(tmp1)]

#cell composition analysis, calculate % of total
rownames(tmp_cellcount)=tmp_cellcount[,1]
tmp_cellcount = tmp_cellcount[,2:ncol(tmp_cellcount), drop=F]
nf = 1/colSums(tmp_cellcount)
tmp_cellcount = sweep(tmp_cellcount, 2, nf, "*")

temp_labels <- mBC@meta.data %>%
  group_by(orig.ident) %>%
  tally()


if (mBC@version == "2.3.4"){
  tmp2 = tmp_cellcount %>% tibble::rownames_to_column() %>%
    reshape2::melt(id.vars = 'rowname') %>%
    mutate(rowname = factor(rowname, levels = levels(mBC@ident)))
  
} else if (mBC@version >= "3.0"){
  tmp2 = tmp_cellcount %>% tibble::rownames_to_column() %>%
    reshape2::melt(id.vars = 'rowname') %>%
    mutate(rowname = factor(rowname, levels = levels(mBC@active.ident)))
}


colnames(tmp2)[1:2]=c("Seurat_clusters", "Sample")

p_PercentOfTotal = tmp2 %>%
  ggplot(aes(Sample, value)) +
  geom_bar(aes(fill = Seurat_clusters), position = 'fill', stat = 'identity') +
  scale_fill_manual(name = 'Seurat_clusters', values = custom_colors$discrete) +
  scale_y_continuous(name = '% of total cells', labels = scales::percent_format(), expand = c(0.01,0)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
file.name_PercentOfTotal="WTA3Fibro_cell_PercentOfTotal.png"
ggsave(file = file.name_PercentOfTotal, plot = p_PercentOfTotal, device="png", units="in", dpi = 300,
       width = 4, height = 4, limitsize=FALSE, bg="white")

