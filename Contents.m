% scGEAToolbox - Single-Cell Gene Expression Analysis Toolbox
% Version 25.8.3 23-Aug-2025
%
% Functions
%   cdgea                     - CDGEA - Change working directory to the scGEAToolbox folder
%   demo_Getting_Started      - Learn to use SCGEATOOL - Single-cell Gene Expression Analysis Tool
%   online_landing            - online_landing is a script.
%   sc_cellcyclescore         - Score cell cycle phases
%   sc_cellscore_admdl        - SC_CELLSCORE_ADMDL   Cell-level gene signature scoring (Seurat/AddModuleScore).
%   sc_cellscore_ucell        - SC_CELLSCORE_UCELL   Cell-level gene signature scoring (UCell-style).
%   sc_celltypeanno           - sc_celltypeanno is a function.
%   sc_cluster_s              - sc_cluster_s - cluster cells using cell embeding s
%   sc_cluster_x              - sc_cluster_x - cluster cells using UMI matrix X
%   sc_csubtypeanno           - sc_csubtypeanno is a function.
%   sc_deg                    - SC_DEG - DEG analysis using Mann–Whitney U test or t-test
%   sc_diffumap               - X = sc_transform(X);
%   sc_distcorr               - sc_distcorr is a function.
%   sc_dpg                    - X = log1p(sc_norm(X));
%   sc_filterc                - sc_filterc is a function.
%   sc_filterg                - sc_filterg is a function.
%   sc_fullembed              - SC_FULLEMBED  Ensure all desired embeddings exist in single-cell object.
%   sc_genestat               - SC_GENESTAT  Compute per-gene statistics and optionally filter invalid values
%   sc_genestats              - SC_GENESTATS  Compute per-gene statistics into a tidy table
%   sc_grn                    - SC_GRN Construct single‑cell gene regulatory network (scGRN)
%   sc_grnview                - SC_GRNVIEW  Display a gene regulatory network as a graph GUI.
%   sc_grnview2               - SC_GRNVIEW2  Display two gene regulatory networks side‑by‑side.
%   sc_hvg                    - Identify HVGs
%   sc_ifft                   - https://github.com/Sanofi-Public/PMCB-scGFT/blob/master/R/Core.R
%   sc_ifft_dev               - https://github.com/Sanofi-Public/PMCB-scGFT/blob/master/R/Core.R
%   sc_impute                 - Imputation
%   sc_knngraph               - Generate KNN group network from cell embeddings
%   sc_mergedata              - sc_mergedata is a function.
%   sc_mergesces              - Merges two SCE objects
%   sc_mtratio                - sc_mtratio is a function.
%   sc_multiembed             - sc_multiembed is a function.
%   sc_nnmf                   - Expression program signatures analysis
%   sc_norm                   - sc_norm is a function.
%   sc_pcnet                  - Construct GRN network A using PC regression (pcnet)
%   sc_pcnetbat               - Construct network using PC regression
%   sc_pcnetpar               - [A]=sc_pcnetpar(X,ncom)
%   sc_phate                  - PHATE embedding of cells
%   sc_pickmarkers            - sc_pickmarkers is a function.
%   sc_potency                - Estimate differentiation potency of cells
%   sc_qcfilter               - Basic QC filter
%   sc_qcmetrics              - SC_QCMETRICS Calculate quality control metrics for single-cell data
%   sc_read10xdir             - Read 10x folder
%   sc_read10xh5file          - Read 10x Genomics H5 file
%   sc_readfile               - SC_READFILE  Load single‑cell expression data files
%   sc_readgeoaccess          - SC_READGEOACCESS  Download and parse GEO single-cell dataset by accession
%   sc_readh5adfile           - Read H5AD file
%   sc_readhdf5file           - Read HDF5 file
%   sc_readloomfile           - Read LOOM file
%   sc_readmtxfile            - Read MTX file
%   sc_readparsebio           - SC_READPARSEBIO  Read ParseBio matrix + feature + cell metadata
%   sc_readrdsfile            - Read Seurat/RDS file
%   sc_readtsvfile            - Read TSV/TXT file
%   sc_rmdugenes              - sc_rmdugenes is a function.
%   sc_rmmtcells              - Remove cells with high mtDNA ratio
%   sc_rmmtgenes              - Remove mt-genes
%   sc_scatter                - SC_SCATTER Single-cell Gene Expression Analysis Toolbox Graphical User Interface.
%   sc_scattermarker          - SC_SCATTERMARKER(X,genelist,g,s,methodid)
%   sc_sce2h5ad               - Write SCE to H5AD file
%   sc_sce2hdf5               - Write SCE to HDF5 file
%   sc_sce2rds                - Write SCE to Seurat/RDS file
%   sc_selectc                - Select cells by library size and number of genes
%   sc_selectg                - Select genes by expression levels
%   sc_simudata               - SC_SIMUDATA  Simulate single-cell RNA‑seq count data
%   sc_snndpc                 - Clustering cell embeddings using SNNDPC - a SNN clustering algorithm
%   sc_splinefit              - SC_SPLINEFIT identify genes with a profile deviated from normal
%   sc_splinefit2             - Compare genes' expression profiles with spline fit regression
%   sc_stemness               - SC_STEMNESS   Compute stemness score for single-cell data
%   sc_tenifoldnet            - scTenifoldNet pipeline for constructing scGRNs from scRNAseq data X0 and
%   sc_tfactivity             - The activity level of a transcription factor (TF) in a given cell is the
%   sc_trajectory             - SC_TRAJECTORY Estimate pseudotime trajectories from single-cell data
%   sc_transform              - SC_TRANSFORM Transformations for single-cell data
%   sc_tsne                   - tSNE embedding of cells
%   sc_umap                   - UMAP embedding of cells
%   sc_writefile              - sc_writefile is a function.
%   scgeatool                 - SCGEATOOL - Launch either original figure GUI or App Designer version
%   scgeatool_legacy          - SCGEATOOL_LEGACY - Single-cell Gene Expression Analysis Toolbox Graphical User Interface.
