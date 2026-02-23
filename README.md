# Stream_2026

Scripts for analyzing and generating visualizations for the CGE Organoid Streams, in Wong et al., 2026. 

## Single Cell Pipeline 

* Step 1 - Spatial Embedding  - Code to analyze and visualize the Day 270 CGE-Enriched organoids from the spatial single-cell RNA seq, and extract the stream signature. 

* Step 2 - Integrate Timepoints  - To QC filter, integrate various timepoints of CGE-Enriched organoids via Harmony, and identify the stream population using the stream signature; as well as call cluster markers and make cell type assignments. 

* Step 3 - Subset CGE  - Subsetting of the CGE lineages from Step 2, with additional cluster annotation. 

* Step 4 - CGE Trajectory  - Application of Monocle3 to the lineages from Step 3 and computation of pseudo-time.

* Step 5 - Call Markers  - Calling of Stream-specific markers within the CGE lineage and intersection with the Nascimento et al., 2024 markers. 

* Step 6 - WGCNA  - Identification of co-expression modules (via hdWGCNA) and their expression across pseudotime. 

* Step 7 - singleCellNet  - Random forest classification via singleCellNet of cell types in the CGE lineages using training data from the human EC Stream (Nascimento et al, 2024). 
