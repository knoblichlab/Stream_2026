##

# ------ load packages ----- #
library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(forcats)
library(patchwork)
library(scales)
library(scattermore)
library(qs)

# Colors 
cols.clusters <- c(brewer.pal(12, name = 'Set3'), '#689689', '#f7d486','#49beaa', '#f27a7d' )
cols.3 <- c('#f7d486','#49beaa', '#f27a7d')

# ------ format data ------ # 

### list files, sample names, and future object names
samples.df <- data.frame(sample = list.files(path = 'classification_cr710/classification/scv2/mtx')) %>%
  mutate(path = paste0('classification_cr710/classification/scv2/mtx/', sample, '/filtered_feature_bc_matrix'),
         out.object = paste0('cts.', sample))


### make each into Seurat object 
for(x in samples.df$sample){
  path = paste0('classification_cr710/classification/scv2/mtx/', x, '/filtered_feature_bc_matrix')
  out.ob <- paste0('cts.', x)
  cts.data <- Read10X(data.dir = path)
  cts <- CreateSeuratObject(counts = cts.data, project = x, min.cells = 3, min.features = 200)
  assign(out.ob, cts)
}

### Merge matrices together 
s.merge <- merge(x = `cts.176-241500`, y = c(`cts.176-241501`, `cts.176-241502`, `cts.TSC.p24.4-241500`, `cts.TSC.p24.4-241501`, `cts.TSC.p24.4-241502`, `cts.TSC.p5.32s-241500`, `cts.TSC.p5.32s-241501`, `cts.TSC.p5.32s-241502`))

length(Cells(s.merge)) # 19827 cells 

qsave(s.merge, file = 'Wong_MergedTimepoints_Unfiltered.qs')


# ----- QC Filtering ---- # 

## standard percent mt & feature threshold
s.merge[['percent.mt']] <- PercentageFeatureSet(s.merge, pattern = '^MT-')

metric.spread <- FetchData(s.merge, vars = c('orig.ident', 'percent.mt', 'nCount_RNA', 'nFeature_RNA'))

s.sub <- subset(s.merge, subset = percent.mt < 5 & nFeature_RNA < (mean(s.merge$nFeature_RNA) + (1.5 * sd(s.merge$nFeature_RNA))))

# Cluster and visualize 
s.sub <- NormalizeData(s.sub)
s.sub <- FindVariableFeatures(s.sub)
s.sub <- ScaleData(s.sub) 
s.sub <- RunPCA(s.sub) 
s.sub <- FindNeighbors(s.sub, dims = 1:20)

s.sub <- FindClusters(s.sub, resolution = 0.5, cluster.name = 'clust_0.5')
s.sub <- FindClusters(s.sub, resolution = 0.3, cluster.name = 'clust_0.3')

s.sub <- RunUMAP(s.sub, reduction = 'pca', dims = 1:20, reduction.name = 'umap_default')
s.sub <- RunUMAP(s.sub, reduction = 'pca', dims = 1:20, n.neighbors = 7L, min.dist = 0.05, reduction.name = 'umap_k7_m005')

# Timepoint metadata 
meta.tp <- FetchData(s.sub, vars = c('orig.ident')) %>%
  rownames_to_column(var = 'cell') %>% 
  rowwise() %>%
  mutate(key = str_split_1(orig.ident, '-')[[1]][1]) %>%
  ungroup() %>% 
  mutate(timepoint = case_when(key == '176' ~ 'D270', 
                               key == 'TSC.p5.32s' ~ 'D160',
                               key == 'TSC.p24.4' ~ 'D60')) %>%
  dplyr::select(cell, timepoint) %>%
  column_to_rownames(var = 'cell')
  
s.sub <- AddMetaData(s.sub, meta.tp)

# ------- integration ------- # 

# run harmony on sample 
ss.har <- RunHarmony(s.sub, 'orig.ident')
ss.har <- RunPCA(ss.har) 
ss.har <- FindNeighbors(ss.har, dims = 1:20)

ss.har <- FindClusters(ss.har, resolution = 0.7, cluster.name = 'hclust_0.7')
ss.har <- FindClusters(ss.har, resolution = 0.5, cluster.name = 'hclust_0.5')
ss.har <- FindClusters(ss.har, resolution = 0.3, cluster.name = 'hclust_0.3')

ss.har <- RunUMAP(ss.har, reduction = 'pca', dims = 1:20, reduction.name = 'har_default')
ss.har <- RunUMAP(ss.har, reduction = 'pca', dims = 1:20, n.neighbors = 10L, min.dist = 0.1, reduction.name = 'har_umap_k10_m01')
ss.har <- RunUMAP(ss.har, reduction = 'pca', dims = 1:20, n.neighbors = 7L, min.dist = 0.05, reduction.name = 'har_umap_k7_m005')

# examine markers 

## neural progenitors 
FeaturePlot(ss.har, reduction = 'har_umap_k10_m01', toupper(c("ascl1", "sox2", "sox1", "egfr", "vim")))

## apcs 
FeaturePlot(ss.har, reduction = 'har_umap_k10_m01', toupper(c("aqp4", "gfap", "s100b", "sparcl1", "f3")))

## opcs 
FeaturePlot(ss.har, reduction = 'har_umap_k10_m01', toupper(c("olig1", "olig2", "sox10", "nkx2-2", "pdgfra")))

## inh progenitors 
FeaturePlot(ss.har, reduction = 'har_umap_k10_m01', toupper(c("dlx2", "dlx5", "dlx6-as1", "scgn", "calb2")))

### cge lineage 
FeaturePlot(ss.har, reduction = 'har_umap_k10_m01', toupper(c("nr2f2", "prox1", "sp8", "reln", "vip")))

### mge lineage 
FeaturePlot(ss.har, reduction = 'har_umap_k10_m01', toupper(c("nkx2-1", "lhx6", "sox6", "dlx6-as1", "gsx6")))

### lge lineage 
FeaturePlot(ss.har, reduction = 'har_umap_k10_m01', toupper(c("meis2", "pax6", "gsx2", "npy", "ebf1")))

### excitatory neurons 
FeaturePlot(ss.har, reduction = 'har_umap_k10_m01', toupper(c("eomes", "tbr1", "cux1", "cux2", "emx1", "satb2","BCL11B" )))

## proliferating 
FeaturePlot(ss.har, reduction = 'har_umap_k10_m01', toupper(c("top2a", "mki67","hopx", "foxg1", "SLC1A3")))

## mature neurons 
FeaturePlot(ss.har, reduction = 'har_umap_k10_m01', toupper(c("NEUROD1", "GAD6","RBFOX3")))

## stream top genes 
FeaturePlot(ss.har, reduction = 'har_umap_k10_m01', toupper(c("DAB1", "FGF14","KAZN", "CNTN5", "RBFOX1")))

### other embedding options 
ss.har <- RunUMAP(ss.har, reduction = 'pca', dims = 1:20, n.neighbors = 15L, min.dist = 0.1, reduction.name = 'har_umap_k15_m01')
ss.har <- RunUMAP(ss.har, reduction = 'pca', dims = 1:20, n.neighbors = 20L, min.dist = 0.1, reduction.name = 'har_umap_k20_m01')

ss.har <- RunUMAP(ss.har, reduction = 'pca', dims = 1:20, n.neighbors = 20L, min.dist = 0.2, reduction.name = 'har_umap_k20_m02')
ss.har <- RunUMAP(ss.har, reduction = 'pca', dims = 1:20, n.neighbors = 30L, min.dist = 0.3, reduction.name = 'har_umap_k30_m03')

ss.har <- RunUMAP(ss.har, reduction = 'pca', dims = 1:20, n.neighbors = 40L, min.dist = 0.1, reduction.name = 'har_umap_k40_m01')

### add annotation (level 1) 
annot1.df <- FetchData(ss.har, vars = c('hclust_0.3')) %>% 
  rownames_to_column(var = 'cell') %>% 
  mutate(annot_1 = case_when(hclust_0.3 %in% c('0', '1', '10') ~ 'IN', 
                             hclust_0.3 %in% c('2', '5', '8') ~ 'IN_Progenitor', 
                             hclust_0.3 %in% c('3', '4', '9', '11') ~ 'IN_Immature', 
                             hclust_0.3 %in% c('6') ~ 'Astrocyte_Progenitor', 
                             hclust_0.3 %in% c('7') ~ 'OPC', 
                             .default = 'unknown')) %>% 
  dplyr::select(cell, annot_1) %>%
  column_to_rownames(var = 'cell')

ss.har <-AddMetaData(ss.har, annot1.df)

### add annotation (level 2)
annot2.df <- FetchData(ss.har, vars = c('hclust_0.3')) %>% 
  rownames_to_column(var = 'cell') %>% 
  mutate(annot_2 = case_when(hclust_0.3 %in% c('0') ~ 'IN_MGE', 
                             hclust_0.3 %in% c('1') ~ 'IN_CGE', 
                             hclust_0.3 %in% c('2', '5', '8') ~ 'IN_Progenitor', 
                             hclust_0.3 %in% c('3') ~ 'IN_Immature_CGE', 
                             hclust_0.3 %in% c('4') ~ 'IN_Immature_MGE', 
                             hclust_0.3 %in% c('9', '11') ~ 'IN_Immature', 
                             hclust_0.3 %in% c('10') ~ 'IN_MGE', 
                             hclust_0.3 %in% c('6') ~ 'Astrocyte_Progenitor', 
                             hclust_0.3 %in% c('7') ~ 'OPC', 
                             .default = 'unknown')) %>% 
  dplyr::select(cell, annot_2) %>%
  column_to_rownames(var = 'cell')

ss.har <-AddMetaData(ss.har, annot2.df)

# save 
ss.har <- JoinLayers(ss.har)
qsave(ss.har, file = 'Harmonized_Wong_AllTimepoints.qs')

### ----- Annotation: Find Markers ----- ### 
# reload 
cts <- qread(file = '../Integrate_Timepoints/Harmonized_Wong_AllTimepoints.qs')

# call markers
cts.markers <- FindAllMarkers(cts, group.by = 'hclust_0.3', logfc.threshold = 0.25, random.seed = 1)

## heatmap
cts.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 4) %>%
    ungroup() -> top4.markers

DoHeatmap(cts, group.by = 'hclust_0.3', features = top4.markers$gene) + NoLegend() + scale_fill_viridis() + theme(aspect.ratio = 0.35)

## annotation level 3 

md.annot3 <- FetchData(cts, vars = c('hclust_0.3')) %>% 
  rownames_to_column(var = 'cell') %>%
  mutate(annot_3 = case_when(hclust_0.3 == '0' ~ '0_IN_Mixed', 
                           hclust_0.3 == '1' ~ '1_IN_CGE',
                           hclust_0.3 == '2' ~ '2_Progenitor_Dividing', 
                           hclust_0.3 == '3' ~ '3_IN_CGE_Immature', 
                           hclust_0.3 == '4' ~ '4_IN_Immature', 
                           hclust_0.3 == '5' ~ '5_IN_Progenitor', 
                           hclust_0.3 == '6' ~ '6_Progenitor_APC', 
                           hclust_0.3 == '7' ~ '7_Progenitor_OPC',
                           hclust_0.3 == '8' ~ '8_Progenitor_OPC',
                           hclust_0.3 == '9' ~ '9_IN_Immature', 
                           hclust_0.3 == '10' ~ '10_IN_Mixed', 
                           hclust_0.3 == '11' ~ '11_Stress', 
                           hclust_0.3 == '12' ~ '12_CP')) %>%
    dplyr::select(cell, annot_3) %>% 
    column_to_rownames(var = 'cell') 

cts <- AddMetaData(cts, md.annot3)

### dimplot 
DimPlot(cts, reduction = 'har_umap_k10_m01', group.by = 'annot_3', label = T, pt.size = .1) + 
  scale_color_manual(values = cols.clusters) + 
  theme(aspect.ratio = 1)

# ------ Identify Stream ----- # 

# subset D270 timepoint 
cts.270 <- subset(cts, subset = timepoint == 'D270')

# harmonize samples & re-embed 
cts.270 <- NormalizeData(cts.270)
cts.270 <- FindVariableFeatures(cts.270)
cts.270 <- ScaleData(cts.270) 
cts.270 <- RunPCA(cts.270)
cts.270 <- RunHarmony(cts.270, 'orig.ident', reduction.save = 'harmony_late') 
cts.270 <- FindNeighbors(cts.270, dims = 1:20)
cts.270 <- FindClusters(cts.270, resolution = 0.3, cluster.name = 'har_late_0.3')
cts.270 <- FindClusters(cts.270, resolution = 0.4, cluster.name = 'har_late_0.4')
cts.270 <- FindClusters(cts.270, resolution = 0.5, cluster.name = 'har_late_0.5')
cts.270 <- RunUMAP(cts.270, reduction = 'harmony_late', dims = 1:20, min.dist = 0, reduction.name = 'harm_ex_d0_umap')
cts.270 <- RunUMAP(cts.270, reduction = 'harmony_late', dims = 1:20, n.neighbors = 70L, min.dist = .1, reduction.name = 'harm_ex_k70d1_umap')
cts.270 <- RunUMAP(cts.270, reduction = 'harmony_late', dims = 1:20, n.neighbors = 20L, min.dist = .9, reduction.name = 'harm_ex_k20d9_umap')
cts.270 <- RunUMAP(cts.270, reduction = 'harmony_late', dims = 1:20, n.neighbors = 20L, min.dist = .8, reduction.name = 'harm_ex_k20d8_umap')
cts.270 <- RunUMAP(cts.270, reduction = 'harmony_late', dims = 1:20, n.neighbors = 40L, min.dist = .8, reduction.name = 'harm_ex_k40d8_umap')
cts.270 <- RunUMAP(cts.270, reduction = 'harmony_late', dims = 1:20, n.neighbors = 10L, min.dist = .8, reduction.name = 'harm_ex_k10d8_umap')
cts.270 <- RunUMAP(cts.270, reduction = 'harmony_late', dims = 1:20, n.neighbors = 20L, min.dist = .7, reduction.name = 'harm_ex_k20d7_umap')
cts.270 <- RunUMAP(cts.270, reduction = 'harmony_late', dims = 1:20, reduction.name = 'harm_late_umap')

# visualize
DimPlot(cts.270, reduction = 'harm_late_umap', group.by = 'annot_3', label = T) + 
  scale_color_manual(values = cols.clusters) + 
  theme(aspect.ratio = 1)

## find markers
cts.markers <- FindAllMarkers(cts.270, group.by = 'har_late_0.5', only.pos = T, min.pct = 0.05, logfc.threshold = 0.25)

# find top markers
cts.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top15.markers

### cluster 10 is clearly stressed (HSP, DNAJ) -- remove 
cts.sub <- subset(cts.270, subset = ! har_late_0.5 == '10')

## module score using spatial stream top genes 

spa.stream.markers <- read_tsv(file = 'Step0b_Top200StreamMarkers.tsv') %>% pull(gene) 

mark.20 <- spa.stream.markers[1:20]

cts.sub <- AddModuleScore(cts.sub, features = list(spa.stream.markers), name = 'Stream_Top200_')
cts.sub <- AddModuleScore(cts.sub, features = list(mark.20), name = 'Stream_Top20_')

# visualize enrichment 
FeaturePlot(cts.sub, reduction = 'harm_ex_k20d8_umap', features = c('Stream_Top200_1')) + 
  scale_color_viridis() + 
  theme(aspect.ratio = 1)


### new annotation 
annot.df <- FetchData(cts.sub, vars = c('har_late_0.5', 'orig.ident')) %>% 
  mutate(annot_late_1 = case_when(har_late_0.5 %in% c(0, 1, 3, 5) ~ 'IN_Neuron_CGE', 
                             har_late_0.5 %in% c(4, 6, 8) ~ 'IN_Progenitor',
                             har_late_0.5 %in% c(2, 9, 13, 14) ~ 'OPC/APC',
                             har_late_0.5 %in% c(11) ~ 'APC/Astrocyte',
                             har_late_0.5 %in% c(7) ~ 'IN_Stream',
                             har_late_0.5 %in% c(12) ~ 'CP', .default = 'flag'), 
         annot_late_2 = case_when(har_late_0.5 == '0' ~ 'IN_CGE_1', 
                                  har_late_0.5 == '1' ~ 'IN_CGE_2', 
                                  har_late_0.5 == '2' ~ 'OPC_APC_1', 
                                  har_late_0.5 == '3' ~ 'IN_CGE_3', 
                                  har_late_0.5 == '4' ~ 'IN_Prog_1', 
                                  har_late_0.5 == '5' ~ 'IN_Mixed', 
                                  har_late_0.5 == '6' ~ 'IN_Prog_2', 
                                  har_late_0.5 == '7' ~ 'IN_Stream', 
                                  har_late_0.5 == '8' ~ 'IN_Prog_3', 
                                  har_late_0.5 == '9' ~ 'OPC_APC_2', 
                                  har_late_0.5 == '11' ~ 'APC_Astro', 
                                  har_late_0.5 == '12' ~ 'CP',
                                  har_late_0.5 == '13' ~ 'OPC_APC_3', 
                                  har_late_0.5 == '14' ~ 'OPC_APC_4', )
         ) 

annot1.md <- annot.df %>% dplyr::select(annot_late_1)
cts.sub <- AddMetaData(cts.sub, annot1.md)
annot2.md <- annot.df %>% dplyr::select(annot_late_2)
cts.sub <- AddMetaData(cts.sub, annot2.md)

# cell type composition across samples 
ggplot(annot.df, aes(x = orig.ident, fill = annot_late_1)) + 
  geom_bar(position = 'fill', width = 0.6, alpha = 0.8) + 
  scale_fill_manual(values = annot.pal) + 
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) + 
  theme_cowplot() + 
  xlab('Sample') + 
  ylab('Cell Proportion') + 
  guides(fill = guide_legend(title = 'Cell Type')) + 
  theme(aspect.ratio = 1)

## plot and compare enrichments to confirm cell type annotation 
mod.df <- FetchData(cts.sub, vars = c('annot_late_1', 'annot_late_2', 'Stream_Top200_1', 'orig.ident'))

# level 1
ggplot(mod.df, aes(x = fct_rev(fct_reorder(annot_late_1, Stream_Top200_1, .fun = mean)), y = Stream_Top200_1, 
                   fill = annot_late_1)) + 
  geom_jitter(width = 0.15, size = .4, alpha = 0.8) + 
  geom_violin(alpha = 0.8) + 
  scale_fill_manual(values = annot.pal) + 
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = 'black') + 
  theme_cowplot() + 
  xlab('Organoid D270 Cluster') +
  ylab('Stream Marker Enrichment') + 
  theme(legend.position = 'none', aspect.ratio = 0.75, axis.text.x = element_text(angle = 45, hjust = 1))

## level 2 
ggplot(mod.df, aes(x = fct_rev(fct_reorder(annot_late_2, Stream_Top200_1, .fun = mean)), y = Stream_Top200_1, 
                   fill = annot_late_2)) + 
  geom_jitter(width = 0.15, size = .2, alpha = 0.8, color = 'grey50') + 
  geom_violin(alpha = 0.8) + 
  scale_fill_manual(values = pal.step) + 
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = 'black') + 
  theme_cowplot() + 
  xlab('Organoid D270 Cluster') +
  ylab('Stream Marker Enrichment') + 
  theme(legend.position = 'none', aspect.ratio = 0.75, axis.text.x = element_text(angle = 45, hjust = 1))

# Repeat module score for full timecourse object 

# Module score for full timecourse object 
cts <- AddModuleScore(cts, features = list(spa.stream.markers), name = 'Stream_Top200_')

# adjust plot order
fp2 <- FeaturePlot(cts, reduction = 'har_umap_k10_m01', features = c('Stream_Top200_1')) 
fp2.data <- fp2$data
fp2.data <- fp2.data[order(fp2.data$Stream_Top200_1), ]

ggplot(fp2.data, aes(x = harumapk10m01_1, y = harumapk10m01_2, color = Stream_Top200_1)) + 
  geom_point(size = 0.3) + 
  scale_color_viridis(option = 'viridis') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1)










