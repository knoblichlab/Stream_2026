##

# ------ load packages ----- #
#
require(Seurat)
require(tidyverse)
require(ggplot2)
require(cowplot)
require(forcats)
require(patchwork)
require(scales)
require(scattermore)




# ------ format data ------ # 

# list files, sample names, and future object names
samples.df <- data.frame(sample = list.files(path = '/groups/knoblich/bioinfo/users/sakurako.wong/scRNAseq/R15864/cellranger/v2_GFP/classification_cr710/classification/scv2/mtx')) %>%
  mutate(path = paste0('/groups/knoblich/bioinfo/users/sakurako.wong/scRNAseq/R15864/cellranger/v2_GFP/classification_cr710/classification/scv2/mtx/', sample, '/filtered_feature_bc_matrix'),
         out.object = paste0('cts.', sample))


# Make each into Seurat object 
for(x in samples.df$sample){
  path = paste0('/groups/knoblich/bioinfo/users/sakurako.wong/scRNAseq/R15864/cellranger/v2_GFP/classification_cr710/classification/scv2/mtx/', x, '/filtered_feature_bc_matrix')
  out.ob <- paste0('cts.', x)
  cts.data <- Read10X(data.dir = path)
  cts <- CreateSeuratObject(counts = cts.data, project = x, min.cells = 3, min.features = 200)
  assign(out.ob, cts)
}

# Merge matrices together 
s.merge <- merge(x = `cts.176-241500`, y = c(`cts.176-241501`, `cts.176-241502`, `cts.TSC.p24.4-241500`, `cts.TSC.p24.4-241501`, `cts.TSC.p24.4-241502`, `cts.TSC.p5.32s-241500`, `cts.TSC.p5.32s-241501`, `cts.TSC.p5.32s-241502`))

length(Cells(s.merge)) # 19827 cells 

s.merge <- JoinLayers(s.merge)

qsave(s.merge, file = 'Wong_MergedTimepoints_Unfiltered.qs')


# ----- QC Filtering ---- # 

## standard percent mt & feature threshold
s.merge[['percent.mt']] <- PercentageFeatureSet(s.merge, pattern = '^MT-')

metric.spread <- FetchData(s.merge, vars = c('orig.ident', 'percent.mt', 'nCount_RNA', 'nFeature_RNA'))

ggplot(metric.spread, aes(x = orig.ident, y = percent.mt)) + 
  geom_scattermore(stat = 'identity', position = position_jitter(width = 0.1), size = 1.8) + 
  geom_violin() + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 12)) + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8))

s.sub <- subset(s.merge, subset = percent.mt < 5 & nFeature_RNA < (mean(s.merge$nFeature_RNA) + (1.5 * sd(s.merge$nFeature_RNA))))

length(Cells(s.sub))

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
  
# -------- init annotation -------- $ 



s.sub <- FindClusters(s.sub, resolution = 0.7, cluster.name = 'clust_0.7')
s.sub <- FindClusters(s.sub, resolution = 0.5, cluster.name = 'clust_0.5')
s.sub <- FindClusters(s.sub, resolution = 0.3, cluster.name = 'clust_0.3')
s.sub <- RunUMAP(s.sub, reduction = 'pca', dims = 1:20, n.neighbors = 10L, min.dist = 0.1, reduction.name = 'umap_k10_m01')
s.sub <- RunUMAP(s.sub, reduction = 'pca', dims = 1:20, n.neighbors = 7L, min.dist = 0.05, reduction.name = 'umap_k7_m005')

s.sub <- RunUMAP(s.sub, reduction = 'pca', dims = 1:20, n.neighbors = 15L, min.dist = 0.1, reduction.name = 'umap_k15_m01')
s.sub <- RunUMAP(s.sub, reduction = 'pca', dims = 1:20, n.neighbors = 20L, min.dist = 0.1, reduction.name = 'umap_k20_m01')
s.sub <- RunUMAP(s.sub, reduction = 'pca', dims = 1:20, n.neighbors = 20L, min.dist = 0.2, reduction.name = 'umap_k20_m02')
s.sub <- RunUMAP(s.sub, reduction = 'pca', dims = 1:20, n.neighbors = 30L, min.dist = 0.3, reduction.name = 'umap_k30_m03')
s.sub <- RunUMAP(s.sub, reduction = 'pca', dims = 1:20, n.neighbors = 40L, min.dist = 0.1, reduction.name = 'umap_k40_m01')

# annot level 1
annot1.df <- FetchData(s.sub, vars = c('clust_0.3')) %>% 
  rownames_to_column(var = 'cell') %>% 
  mutate(annot_1 = case_when(clust_0.3 %in% c('0', '1', '10') ~ 'IN', 
                             clust_0.3 %in% c('2', '5', '8') ~ 'IN_Progenitor', 
                             clust_0.3 %in% c('3', '4', '9', '11') ~ 'IN_Immature', 
                             clust_0.3 %in% c('6') ~ 'Astrocyte_Progenitor', 
                             clust_0.3 %in% c('7') ~ 'OPC', 
                             .default = 'unknown')) %>% 
  dplyr::select(cell, annot_1) %>%
  column_to_rownames(var = 'cell')

s.sub <-AddMetaData(s.sub, annot1.df)

# Annot_Level_2
annot2.df <- FetchData(s.sub, vars = c('clust_0.3')) %>% 
  rownames_to_column(var = 'cell') %>% 
  mutate(annot_2 = case_when(clust_0.3 %in% c('0') ~ 'IN_MGE', 
                             clust_0.3 %in% c('1') ~ 'IN_CGE', 
                             clust_0.3 %in% c('2', '5', '8') ~ 'IN_Progenitor', 
                             clust_0.3 %in% c('3') ~ 'IN_Immature_CGE', 
                             clust_0.3 %in% c('4') ~ 'IN_Immature_MGE', 
                             clust_0.3 %in% c('9', '11') ~ 'IN_Immature', 
                             clust_0.3 %in% c('10') ~ 'IN_MGE', 
                             clust_0.3 %in% c('6') ~ 'Astrocyte_Progenitor', 
                             clust_0.3 %in% c('7') ~ 'OPC', 
                             .default = 'unknown')) %>% 
  dplyr::select(cell, annot_2) %>%
  column_to_rownames(var = 'cell')

s.sub <- AddMetaData(s.sub, annot2.df)

DimPlot(s.sub, reduction = 'har_umap_k10_m01', group.by = 'annot_2', label = T)

# Top Markers 
cts.markers <- FindAllMarkers(s.sub, group.by = 'clust_0.3', logfc.threshold = 0.25, random.seed = 1)

write.table(cts.markers, file = 'Step0_FullObject_Markers.tsv', sep = '\t', quote = F, row.names = T, col.names = T)


cts.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 4) %>%
    ungroup() -> top4.markers


p0 <- DoHeatmap(cts, group.by = 'clust_0.3', features = top4.markers$gene) + NoLegend() + scale_fill_viridis() + theme(aspect.ratio = 0.35)

ggsave(p0, file = 'Step0_FullObject_MarkerHeatmap.pdf', width = 12, height = 6, units = 'in', device = 'pdf')

p0


md.annot3 <- FetchData(s.sub, vars = c('clust_0.3')) %>% 
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

s.sub <- AddMetaData(s.sub, md.annot3)

qsave(s.sub, file = '../Integrate_Timepoints/Wong_AllTimepoints.qs')

# ------- integration ------- # 

# run harmony by sample
s.sub <- NormalizeData(s.sub)
s.sub <- FindVariableFeatures(s.sub)
s.sub <- ScaleData(s.sub) 
s.sub <- RunPCA(s.sub) 
s.sub <- RunHarmony(s.sub, 'orig.ident', reduction.save = 'Harmony')

s.sub <- FindNeighbors(s.sub, reduction = 'Harmony', dims = 1:30)
s.sub <- FindClusters(s.sub, reduction = 'Harmony', resolution = 0.3, cluster.name = 'harm_0.3')
s.sub <- RunUMAP(s.sub, reduction = 'Harmony', dims = 1:20, reduction.name = 'harm_umap')
s.sub <- RunUMAP(s.sub, reduction = 'Harmony', dims = 1:30, n.neighbors = 10L, min.dist = 0.1, reduction.name = 'har_umap_k10_m01')
s.sub <- RunUMAP(s.sub, reduction = 'Harmony', dims = 1:30, n.neighbors = 20L, min.dist = 0.0, reduction.name = 'har_umap_k20_m0')
s.sub <- RunUMAP(s.sub, reduction = 'Harmony', dims = 1:30, n.neighbors = 20L, min.dist = 0.2, reduction.name = 'har_umap_k20_m02')

# markers
s.markers.df <- FindAllMarkers(s.sub, group_by = 'harm_0.3', only.pos = T)
s.markers.top <- s.markers.df %>% group_by(cluster) %>% slice_head(n = 20)

# h annotation
df.hannot <- FetchData(s.sub, vars = c('harm_0.3', 'timepoint')) %>%
  mutate(h_annot = case_when(harm_0.3 %in% c('0') ~ 'IN_CGE', 
                             harm_0.3 %in% c('1') ~ 'IN_Mixed', 
                             harm_0.3 %in% c('2') ~ 'Stress', 
                             harm_0.3 %in% c('3') ~  'IN_CGE_Imm', 
                             harm_0.3 %in% c('4') ~ 'IN_Prog', 
                             harm_0.3 %in% c('5') ~ 'APC_OPC', 
                             harm_0.3 %in% c('6', '9', '12', '14') ~ 'OPC_Olig', 
                             harm_0.3 %in% c('7') ~ 'IN_Prog_MGE', 
                             harm_0.3 %in% c('8') ~ 'IN_Mixed_Imm', 
                             harm_0.3 %in% c('10') ~ 'IN_Prog_Cycling', 
                             harm_0.3 %in% c('11') ~ 'CP', 
                             harm_0.3 %in% c('13') ~ 'NE_cells'))

s.sub <- AddMetaData(s.sub, df.hannot['h_annot'])

h.markers <- FindAllMarkers(s.sub, group_by = 'h_annot', only.pos = T)

h.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 4) %>%
    ungroup() -> top4.markers


p0 <- DoHeatmap(s.sub, group.by = 'h_annot', features = top4.markers$gene) + NoLegend() + scale_fill_viridis() + theme(aspect.ratio = 0.35)

ggsave(p0, file = 'Step1_H_Annot_MarkerHeatmap.pdf', width = 12, height = 6, units = 'in', device = 'pdf')

p0

# score stream 
spa.stream.markers <- read_tsv(file = '../Pipeline_Sept_2025/Step0c_Top200StreamMarkers_DupFilter_Res05.tsv') %>% pull(gene) 
s.sub <- AddModuleScore(s.sub, features = list(spa.stream.markers), name = 'Stream_Top200_')

## adjust plot order
fp2 <- FeaturePlot(s.sub, reduction = 'har_umap_k20_m02', features = c('Stream_Top200_1')) 
fp2.data <- fp2$data
fp2.data <- fp2.data[order(fp2.data$Stream_Top200_1), ]

p.d1 <- ggplot(fp2.data, aes(x = harumapk20m02_1, y = harumapk20m02_2, color = Stream_Top200_1)) + 
  geom_point(size = 0.3) + 
  scale_color_viridis(option = 'viridis') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1)

ggsave(p.d1, file = 'Step1_ReMerge_Feature_StreamModuleScore.pdf', width = 7, height = 6, units = 'in', device = 'pdf')

p.d1



qsave(s.sub, file = 'Full_Harmonized_ReMerged.qs')


# ------ Identify Stream ----- # 

# subset D270 timepoint 
cts.270 <- subset(s.sub, subset = timepoint == 'D270')

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










