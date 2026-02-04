

# ----- packages ------ # 

# packages
require(Seurat)
require(presto)

# Plotting
require(tidyverse)
require(ggplot2)
require(cowplot)
require(forcats)
require(patchwork)
require(scales)
require(scattermore)
require(viridis)


# Colors 
cols.clusters <- c(brewer.pal(12, name = 'Set3'), '#689689')
cols.clusters <- c(brewer.pal(12, name = 'Set3'), '#689689', '#247BA0',  '#F25F5C', '#ffafcc')
cols.3 <- c('#f7d486','#49beaa', '#f27a7d')

#  ----- subset ----- # 


# Load full object 
cts <- qread(file = '../Integrate_Timepoints/Harmonized_Wong_AllTimepoints.qs')

keep.list <- c('1_IN_CGE', '2_Progenitor_Dividing', '3_IN_CGE_Immature', '4_IN_Immature',
               '5_IN_Progenitor', '9_IN_Immature')

cts.cge <- subset(cts, subset = annot_3 %in% keep.list)

DimPlot(cts.cge, group.by = 'annot_3', label = T)

# re-embed 
cts.cge <- FindVariableFeatures(cts.cge)
cts.cge <- ScaleData(cts.cge) 
cts.cge <- RunPCA(cts.cge) 

cts.cge <- FindNeighbors(cts.cge, dims = 1:20)

cts.cge <- FindClusters(cts.cge, resolution = 0.5, cluster.name = 'cge_clust_0.5')
cts.cge <- FindClusters(cts.cge, resolution = 0.3, cluster.name = 'cge_clust_0.3')

cts.cge <- RunUMAP(cts.cge, reduction = 'pca', dims = 1:20, reduction.name = 'cge_umap_default')
cts.cge <- RunUMAP(cts.cge, reduction = 'pca', dims = 1:20, n.neighbors = 7L, min.dist = 0.05, reduction.name = 'cge_umap_k7_m005')

# harmonize 
cge.har <- RunHarmony(cts.cge, 'orig.ident')
cge.har <- RunPCA(cge.har, reduction.name = 'cge_har_pca') 
cge.har <- FindNeighbors(cge.har, reduction = 'cge_har_pca', dims = 1:20)
cge.har <- FindClusters(cge.har, reduction = 'cge_har_pca', resolution = 0.5, cluster.name = 'cge_hclust_0.5')

cge.har <- RunUMAP(cge.har, reduction = 'cge_har_pca', dims = 1:20, reduction.name = 'cge_har_default')

dp1 <- DimPlot(cge.har, reduction = 'cge_har_default', group.by = 'cge_hclust_0.5', label = T) + 
  scale_color_manual(values = cols.clusters) + 
  theme(aspect.ratio = 1)

dp2 <- DimPlot(cge.har, reduction = 'cge_har_default', group.by = 'timepoint', label = T) + 
  scale_color_manual(values = cols.3) + 
  theme(aspect.ratio = 1)

dp1 + dp2

### find markers 
cge.markers <- FindAllMarkers(cge.har, group.by = 'cge_hclust_0.5', logfc.threshold = 0.25, random.seed = 1)


write.table(cge.markers, file = 'Step1_CGE_Trajectory_Cluster_AllMarkers.tsv', sep = '\t', quote = F, row.names = T, col.names = T)

cge.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 6) %>%
    ungroup() -> top6.cge.markers

cge.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top15.cge.markers


DoHeatmap(cge.har, group.by = 'cge_hclust_0.5', features = top6.cge.markers$gene) + NoLegend() + scale_fill_viridis() 

## ---- second subset ----- ##
### there are clusters which do not belong: 15 is stessed; 3 is MGE-like; 7 is stressed 
cge2 <- subset(cge.har, subset = !(cge_hclust_0.5 %in% c('3', '7', '15')))
cge2 <- NormalizeData(cge2)
cge2 <- FindVariableFeatures(cge2)
cge2 <- ScaleData(cge2) 
cge2 <- RunPCA(cge2) 
cge2 <- FindNeighbors(cge2, dims = 1:20)
cge2 <- FindClusters(cge2, resolution = 0.5, cluster.name = 'cge2_clust_0.5')
cge2 <- FindClusters(cge2, resolution = 0.3, cluster.name = 'cge2_clust_0.3')
cge2 <- RunUMAP(cge2, reduction = 'pca', dims = 1:20, reduction.name = 'cge2_umap_default')
DimPlot(cge2, reduction = 'cge2_umap_default', group.by = 'cge2_clust_0.5', label = T) + 
  theme(aspect.ratio = 1)

cge2.markers <- FindAllMarkers(cge2, group.by = 'cge2_clust_0.5', logfc.threshold = 0.25, random.seed = 1)

cge2.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top15.cge2.markers

## ----- third subset ----- ##
### Cluster 3 contains CRABP1 and LHX6 neurons, non-CGE lineage. 

cge3 <- subset(cge2, subset = !(cge2_clust_0.5 == '3'))
cge3 <- NormalizeData(cge3)
cge3 <- FindVariableFeatures(cge3)
cge3 <- ScaleData(cge3) 
cge3 <- RunPCA(cge3) 
cge3 <- FindNeighbors(cge3, dims = 1:20)
cge3 <- FindClusters(cge3, resolution = 0.5, cluster.name = 'cge3_clust_0.5')
cge3 <- FindClusters(cge3, resolution = 0.3, cluster.name = 'cge3_clust_0.3')
cge3 <- RunUMAP(cge3, reduction = 'pca', dims = 1:20, reduction.name = 'cge3_umap_default')

DimPlot(cge3, reduction = 'cge3_umap_default', group.by = 'cge3_clust_0.5', label = T) + 
  theme(aspect.ratio = 1)

DimPlot(cge3, reduction = 'cge3_umap_default', group.by = 'timepoint', label = T) + 
  theme(aspect.ratio = 1)
### markers
cge3.markers <- FindAllMarkers(cge3, group.by = 'cge3_clust_0.5', logfc.threshold = 0.25, random.seed = 1)

# Top 5
cge3.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 6) %>%
    ungroup() -> top6.cge3.markers
# Top 15
cge3.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top15.cge3.markers

DoHeatmap(cge3, group.by = 'cge3_clust_0.5', features = top6.cge3.markers$gene) + NoLegend() + scale_fill_viridis() 

### harmonize 

cge3.har <- RunHarmony(cge3, 'orig.ident')
cge3.har <- RunPCA(cge3.har, reduction.name = 'cge3_har_pca') 
cge3.har <- FindNeighbors(cge3.har, reduction = 'cge3_har_pca', dims = 1:20)
cge3.har <- FindClusters(cge3.har, reduction = 'cge3_har_pca', resolution = 0.5, cluster.name = 'cge3_hclust_0.5')
cge3.har <- RunUMAP(cge3.har, reduction = 'cge3_har_pca', dims = 1:20, reduction.name = 'cge3_har_default')
DimPlot(cge3.har, reduction = 'cge3_har_default', group.by = 'cge3_hclust_0.5', label = T) + 
  scale_color_manual(values = cols.clusters) + 
  theme(aspect.ratio = 1)
DimPlot(cge3.har, reduction = 'cge3_har_default', group.by = 'timepoint', label = T) + 
  scale_color_manual(values = cols.3) + 
  theme(aspect.ratio = 1)

## ---- fourth subset ---- ##
### More off targets to filter - Cluster 10 is late APC progenitors; cluster 12 look like mispatterned melanocytes 
cge4 <- subset(cge3.har, subset = !(cge3_hclust_0.5 %in% c('10', '12')))
cge4 <- NormalizeData(cge4)
cge4 <- FindVariableFeatures(cge4)
cge4 <- ScaleData(cge4) 
cge4 <- RunPCA(cge4) 
cge4 <- FindNeighbors(cge4, dims = 1:20)
cge4 <- FindClusters(cge4, resolution = 0.5, cluster.name = 'cge4_clust_0.5')
cge4 <- FindClusters(cge4, resolution = 0.3, cluster.name = 'cge4_clust_0.3')
cge4 <- RunUMAP(cge4, reduction = 'pca', dims = 1:20, reduction.name = 'cge4_umap_default')
DimPlot(cge4, reduction = 'cge4_umap_default', group.by = 'cge4_clust_0.5', label = T) + 
  theme(aspect.ratio = 1)
cge4.har <- RunHarmony(cge4, 'orig.ident')
cge4.har <- RunPCA(cge4.har, reduction.name = 'cge4_har_pca') 
cge4.har <- FindNeighbors(cge4.har, reduction = 'cge4_har_pca', dims = 1:20)
cge4.har <- FindClusters(cge4.har, reduction = 'cge4_har_pca', resolution = 0.5, cluster.name = 'cge4_hclust_0.5')

cge4.har <- RunUMAP(cge4.har, reduction = 'cge4_har_pca', dims = 1:20, reduction.name = 'cge4_har_default')

dp6 <- DimPlot(cge4.har, reduction = 'cge4_har_default', group.by = 'cge4_hclust_0.5', label = T) + 
  scale_color_manual(values = cols.clusters) + 
  theme(aspect.ratio = 1)
dp7 <- DimPlot(cge4.har, reduction = 'cge4_har_default', group.by = 'timepoint', label = T) + 
  scale_color_manual(values = cols.3) + 
  theme(aspect.ratio = 1)

dp6 + dp7

### markers 
cge4.markers <- FindAllMarkers(cge4.har, group.by = 'cge4_hclust_0.5', logfc.threshold = 0.25, random.seed = 1)

# Top 5
cge4.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 6) %>%
    ungroup() -> top6.cge4.markers

# Top 15
cge4.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top15.cge4.markers

DoHeatmap(cge4.har, group.by = 'cge4_hclust_0.5', features = top6.cge4.markers$gene) + NoLegend() + scale_fill_viridis() 



### add annotation 
md.annot4 <- FetchData(cge4.har, vars = c('cge4_hclust_0.5')) %>%
  rownames_to_column(var = 'cell') %>% 
  mutate(annot_4 = case_when(cge4_hclust_0.5 %in% c('0', '3') ~ 'CGE_Interneurons', 
                             cge4_hclust_0.5 %in% c('10') ~ 'Stream_INs', 
                             cge4_hclust_0.5 %in% c('5') ~ 'Early1_IN_Immature', 
                             cge4_hclust_0.5 %in% c('1', '8') ~ 'Early2_IN_Immature',
                             cge4_hclust_0.5 %in% c('6', '2') ~ 'Late_IN_Immature',
                             cge4_hclust_0.5 %in% c('7', '9') ~ 'IN_Progenitor',
                             cge4_hclust_0.5 %in% c('4') ~ 'Late_CGE_Progenitor', 
                             .default = 'flag') ) %>%
  dplyr::select(cell, annot_4) %>%
  column_to_rownames(var = 'cell')

cge4.har <- AddMetaData(cge4.har, md.annot4)

### another level by cluster 
md.annot_cluster <- FetchData(cge4.har, vars = c('cge4_hclust_0.5', 'annot_4')) %>% 
  mutate(annot_cluster = paste0(cge4_hclust_0.5, '_', annot_4)) %>%
  dplyr::select(annot_cluster)

cge4.har <- AddMetaData(cge4.har, md.annot_cluster)

qsave(cge4.har, file = 'Step1_CGE_Trajectory_LastSubset.qs')





