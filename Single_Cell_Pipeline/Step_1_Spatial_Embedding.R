
# ------ packages ------- # 

library(qs)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(lemon)

# ------ embed ------- #

# load objects 

## full spatial embedding (all cells/bins) 
spa.full <- readRDS(file = 'objects/D02168B2.bin50.ndim20.600.filt.rds')

## subset of stream and adjacent cells (visually identified) 
spa <- qread(file = 'objects/nF50.nvar2000.ndim20.SCT_snn_res.1filt.100.qs')
spa$part <- 'stream_call'

# fetch coordinates & plot


spa.embed <- FetchData(spa, vars = c('x', 'y', 'cl', 'geschwind', 'EC_Stream')) %>% 
  rownames_to_column(var = 'cell') %>% 
  mutate(coords = paste0(x, '_', y))

p0 <- ggplot(spa.embed, aes(x = x, y = y, color = cl)) + 
  geom_point(size = 0.35) + 
  scale_color_manual(values = c('#c77dff','#38b000' ,'grey80'), na.value = 'grey80') +
  theme_cowplot() + 
  theme(aspect.ratio = .6)

## full

full.embed <- FetchData(spa.full, vars = c('x', 'y')) %>% 
  rownames_to_column(var = 'cell') %>% 
  mutate(coords = paste0(x, '_', y))

### QC Plots 
ggplot(full.embed, aes(x = x, y = y, color = nCount_Spatial)) + 
  scale_color_viridis() + 
  geom_point(size = .35) + 
  theme_cowplot() + 
  theme(aspect.ratio = .6)

ggplot(full.embed, aes(x = x, y = y, color = nFeature_Spatial)) + 
  scale_color_viridis() + 
  geom_point(size = .35) + 
  theme_cowplot() + 
  theme(aspect.ratio = .6)

ggplot(full.embed, aes(x = x, y = y, color = percent.mito)) + 
  scale_color_viridis() + 
  geom_point(size = .35) + 
  theme_cowplot() + 
  theme(aspect.ratio = .6)

### highlighted stream in spatial embedding 
ggplot(full.embed, aes(x = x, y = y)) + 
  geom_point(size = 0.35, color = 'grey') + 
  geom_point(data = spa.embed, aes(x = x, y = y, color = cl), size = 0.35) + 
  scale_color_manual(values = c('#c77dff','#38b000' ,'grey80'), na.value = 'grey80') + 
  theme_cowplot() + 
  theme(aspect.ratio = .6) 


# merge objects & metadata 

## add RNA layer to full spatial and merge

spa.full[['RNA']] <- spa.full[['Spatial']]
spa.merge <- merge(x = spa.full, y = spa)
spa.merge <- JoinLayers(spa.merge, assay = 'RNA')

# re-embed merged object 

DefaultAssay(spa.merge) <- "RNA"
spa.merge <- NormalizeData(spa.merge)
spa.merge <- FindVariableFeatures(spa.merge)
spa.merge <- ScaleData(spa.merge) 
spa.merge <- RunPCA(spa.merge)
spa.merge <- FindNeighbors(spa.merge, dims = 1:20)
### options for culstering resolutions
spa.merge <- FindClusters(spa.merge, resolution = 0.3, cluster.name = 'clust_0.3')
spa.merge <- FindClusters(spa.merge, resolution = 0.35, cluster.name = 'clust_0.35')
spa.merge <- FindClusters(spa.merge, resolution = 0.4, cluster.name = 'clust_0.4')
spa.merge <- FindClusters(spa.merge, resolution = 0.5, cluster.name = 'clust_0.5')
spa.merge <- RunUMAP(spa.merge, reduction = 'pca', dims = 1:20, min.dist = 0, reduction.name = 'd0_umap')

### plot stream versus neighbors 
DimPlot(spa.merge, group.by = 'cl', reduction = 'd0_umap') + 
  scale_color_manual(values = c('#c77dff','#38b000' ,'grey80'), na.value = 'grey80') +
  theme(aspect.ratio = 1)

### plot clusters 
DimPlot(spa.merge, group.by = 'clust_0.5', reduction = 'd0_umap', label = T, pt.size = 0.1) + 
  scale_color_manual(values = c(brewer.pal(12, 'Set3'), '#38b000')) + 
  theme(aspect.ratio = 1) 

## ------ call markers ------ # 

spa.full.markers <- FindAllMarkers(spa.merge, group.by = 'clust_0.5', only.pos = T)

# find top markers
full.cells <- clust.stream.df %>% 
  rownames_to_column(var = 'cell') %>% 
  dplyr::filter(is.na(part)) %>% 
  droplevels() %>%
  pull(cell)

spa.sub <- subset(spa.merge, cells = full.cells)

spa.full.markers.5 <- FindAllMarkers(spa.sub, group.by = 'clust_0.5', only.pos = T, max.cells.per.ident = 300)

# find top markers
spa.full.markers.5 %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 30) %>%
    ungroup() -> spa.top30.5

# label clusters
spa.merge.meta.df <- FetchData(spa.merge, vars = c('x', 'y', 'clust_0.5', 'cl')) %>%
  mutate(annot_spa = case_when(clust_0.5 == '0' ~ 'CGE_IN', 
                               clust_0.5 == '1' ~ 'Glia/NPC', 
                               clust_0.5 == '2' ~ 'Neuron_1', 
                               clust_0.5 == '3' ~ 'Neuron_2', 
                               clust_0.5 == '4' ~ 'Neuron_3', 
                               clust_0.5 == '5' ~ 'Mes', 
                               clust_0.5 == '6' ~ 'CGE_Stream', 
                               clust_0.5 == '7' ~ 'Unknown', 
                               clust_0.5 == '8' ~ 'Neuron_4', 
                               clust_0.5 == '9' ~ 'Neuron_5', 
                               clust_0.5 == '10' ~ 'Neuron_6', 
                               clust_0.5 == '11' ~ 'Cycling_Prog', 
                               clust_0.5 == '12' ~ 'Neuron_7', 
                               .default = 'flag'))

### add metadata
spa.merge <- AddMetaData(spa.merge, spa.an.meta.df[annot_spa])

### plot new clusters 
DimPlot(spa.merge, group.by = 'annot_spa', reduction = 'd0_umap', label = T, pt.size = 0.1, label.size = 3) + 
  scale_color_manual(values = c(brewer.pal(12, 'Set3'), '#38b000')) + 
  theme(aspect.ratio = 1) 

### ratio of stream per cluster 
clust.stream.df <- FetchData(spa.merge, vars = c('clust_0.5', 'cl', 'part', 'annot_spa'))

ggplot(clust.stream.df %>% dplyr::filter(!is.na(cl)), aes(x = annot_spa, fill = cl)) + 
  geom_bar(position = 'fill', color = 'black') + 
  scale_fill_manual(values = c('#c77dff','#38b000'), na.value = 'grey80') + 
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), aspect.ratio = 1)


### plot labels with spatial coordinates
ggplot(clust.embed, aes(x = x, y = y, color = annot_spa)) + 
  geom_point(size = 0.35) + 
  scale_color_manual(values = c(brewer.pal(8, 'Set2'), palette)) + 
  theme_cowplot() + 
  theme(aspect.ratio = .6) 

### ------ extract stream markers ------ ###

# Module Score
spa.full.markers.5 %>%
    dplyr::filter(avg_log2FC > 1 & cluster == '6') %>%
    slice_head(n = 200) %>%
    ungroup() -> stream.top200

## Plot Volcano 
### get positive & negative markers for volcano 
spa.full.bidir.markers <- FindAllMarkers(spa.sub, group.by = 'clust_0.5', only.pos = F, max.cells.per.ident = 300)

### find top markers
spa.full.bidir.markers %>%
    dplyr::filter(cluster == '6') %>% 
  mutate(cat = case_when(avg_log2FC > 1  & p_val_adj < 0.05  ~ 'up', 
                           avg_log2FC < -1  & p_val_adj < 0.05 ~ 'down', 
                         .default = 'ns')) -> 
    bidir.spa.stream.markers

### Top upregulated genes 
top.genes <- bidir.spa.stream.markers %>% group_by(cat) %>% arrange(p_val_adj) %>% slice_head(n = 20) %>% dplyr::filter(!cat == 'ns') %>%  pull(gene)

### make volcano frame
volc.df <- bidir.spa.stream.markers %>% 
  mutate(label = case_when(gene %in% top.genes ~ gene, .default = NA_character_)) %>% 
  mutate(p_val_adj = if_else(p_val_adj == 0, min(p_val_adj[p_val_adj > 0], na.rm = TRUE), p_val_adj))

ggplot(volc.df, aes(x = avg_log2FC, y = -log(p_val_adj), color = cat, label = label)) + 
  geom_point(size = 0.5) + 
  scale_color_manual(values = c('#7209b7','grey70','#38b000')) + 
  geom_hline(yintercept = -log(0.05, 10), color = 'grey50') + 
  geom_vline(xintercept = 1, color = 'grey50') + 
  geom_vline(xintercept = -1, color = 'grey50') + 
  xlab('Average Log2FC') + 
  ylab('Average P.Adjusted') + 
  ggtitle('Spatial Stream Markers') + 
  geom_label_repel(size = 5) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1, legend.position = 'none')






