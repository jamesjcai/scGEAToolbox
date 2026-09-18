.. _function_reference:

Function Reference
==================

A map of what is available, by task. Every entry has help text in MATLAB --
``help sc_grn``, ``help run.py_scimilarity`` -- which is the authoritative
description of arguments and defaults. This page is for finding the name.

.. contents:: On this page
   :local:
   :depth: 2

Data import
-----------

Readers return ``[X, genelist, celllist]`` unless noted.

=========================  ===============================================
Function                   Reads
=========================  ===============================================
``sc_read10xdir``          a 10x Genomics folder
``sc_read10xh5file``       a 10x Genomics ``.h5`` file
``sc_readfile``            delimited text (csv/tsv/tab/txt), mtx or h5
``sc_readtsvfile``         TSV/TXT expression matrix
``sc_readmtxfile``         Matrix Market ``.mtx`` plus features/barcodes
``sc_readh5adfile``        AnnData ``.h5ad``
``sc_readh5mufile``        MuData ``.h5mu`` (multimodal)
``sc_h5mumodalities``      lists the modalities inside a ``.h5mu``
``sc_readhdf5file``        generic HDF5
``sc_readloomfile``        Loom
``sc_readrdsfile``         Seurat / RDS (requires R)
``sc_readparsebio``        ParseBio matrix + features + cell metadata
``sc_readgeoaccess``       a GEO accession, downloaded -- returns an SCE
=========================  ===============================================

Data export
-----------

=====================  ===================================================
Function               Writes
=====================  ===================================================
``sc_writefile``       delimited text
``sc_sce2h5ad``        an SCE to ``.h5ad``
``sc_sce2hdf5``        an SCE to HDF5
``sc_sce2rds``         an SCE to Seurat / RDS
``sc_sce2mp4``         a rotating 3-D embedding to video
=====================  ===================================================

Quality control and filtering
-----------------------------

=====================  ===================================================
Function               Purpose
=====================  ===================================================
``sc_qcfilter``        combined QC filter, cells and genes
``sc_qcmetrics``       robust scale of a per-cell quality statistic
``sc_filterc``         filter cells
``sc_filterg``         filter genes
``sc_selectc``         select cells by library size and gene count
``sc_selectg``         select genes by expression level
``sc_rmmtcells``       remove cells with a high mtDNA ratio
``sc_rmmtgenes``       remove mitochondrial genes
``sc_rmdugenes``       remove duplicate genes
``sc_mtratio``         per-cell mitochondrial read fraction
``sc_scrublet``        doublet detection against simulated doublets
``sc_geosketch``       geometric-sketching subsample of cells
=====================  ===================================================

``sc_selectg(X, g, min_cells_nonzero, isexpressed_cutoff)`` takes the cell
count first and the per-cell read threshold second. A value below 1 for
``min_cells_nonzero`` is read as a fraction of all cells.

Normalization and transformation
--------------------------------

=========================  ===============================================
Function                   Purpose
=========================  ===============================================
``sc_norm``                ``'libsize'`` (default), ``'deseq'``, ``'shiftedclr'``
``sc_transform``           Pearson residuals, kNN smoothing, SCTransform, Freeman-Tukey
``sc_sctransformv2``       pure-MATLAB port of ``sctransform::vst(v2)``
``sc_regressout``          regress out unwanted per-cell covariates
``sc_impute``              imputation
``sc_mergesces``           merge SCE objects, ``'union'`` or ``'intersect'`` genes
``sc_mergedata``           merge count matrices
=========================  ===============================================

Feature selection
-----------------

=========================  ===============================================
Function                   Purpose
=========================  ===============================================
``sc_hvg``                 highly variable genes, squared-CV method (Brennecke)
``sc_splinefit``           genes deviating from the fitted 3-D curve
``sc_splinefit2``          the same, comparing two data sets (DD genes)
``sc_analyticfit``         closed-form replacement for the spline curve
``sc_analyticfit2``        the same, comparing two conditions
``sc_genestat``            per-gene statistics as arrays
``sc_genestats``           per-gene statistics as a tidy table
=========================  ===============================================

Dimensionality reduction
------------------------

=====================  ===================================================
Function               Purpose
=====================  ===================================================
``sc_tsne``            t-SNE embedding
``sc_umap``            UMAP embedding
``sc_phate``           PHATE embedding
``sc_multiembed``      several embeddings in one call
``sc_fullembed``       embedding over the full gene set
``sc_diffumap``        differential UMAP
``sc_nnmf``            non-negative matrix factorization
``sc_nmfpattern``      NMF pattern discovery with pathway enrichment
=====================  ===================================================

Clustering
----------

=====================  ===================================================
Function               Purpose
=====================  ===================================================
``sc_cluster_s``       cluster on an embedding: kmeans, kmedoids, dbscan,
                       spectclust, snndpc, mbkmeans
``sc_cluster_x``       cluster on the expression matrix: sc3, simlr,
                       soptsc, sinnlrr
``sc_snndpc``          SNN density-peak clustering
``sc_knngraph``        kNN group network from an embedding
=====================  ===================================================

Differential analysis
---------------------

=========================  ===============================================
Function                   Purpose
=========================  ===============================================
``sc_deg``                 Mann-Whitney U (default) or t-test
``sc_degnb``               negative-binomial test on counts
``sc_degmast``             MAST hurdle model
``sc_dvg``                 differential variability between two SCE objects
``sc_dpg``                 differential program analysis
``sc_diffabundance``       differential cell-type abundance between conditions
``sc_pickmarkers``         marker genes per cluster
=========================  ===============================================

Gene sets, pathways and scoring
-------------------------------

=========================  ===============================================
Function                   Purpose
=========================  ===============================================
``sc_cellscore``           signature scoring: UCell, AddModuleScore, AUCell
``sc_cellcyclescore``      cell cycle phase scores
``sc_stemness``            stemness score
``sc_potency``             differentiation potency
``sc_competitionscore``    cell competition between conditions or regions
``sc_malignscore``         malignancy score from a CNV profile
``sc_fgsea``               preranked GSEA against Enrichr libraries
``sc_gsettest``            competitive gene-set tests on a ranked list
``sc_pathwayactivity``     per-cell signaling pathway activity
``sc_tfactivity``          per-cell transcription factor activity
=========================  ===============================================

Cell type annotation
--------------------

=========================  ===============================================
Function                   Purpose
=========================  ===============================================
``sc_annotatecells``       one front door over every method below
``sc_celltypeanno``        marker matching against PanglaoDB
``sc_csubtypeanno``        subtypes within one primary cell type
``sc_singler``             reference-based, SingleR algorithm
``sc_scibettrain``         train a SciBet reference classifier
``sc_scibetpredict``       assign types with a trained SciBet model
=========================  ===============================================

Trajectory and dynamics
-----------------------

=====================  ===================================================
Function               Purpose
=====================  ===================================================
``sc_trajectory``      pseudotime: ``'splinefit'`` (default) or ``'tscan'``
``sc_meld``            relative likelihood each cell came from each sample
``sc_ifft``            scGFT synthetic cell generation
``sc_simudata``        simulated data sets
=====================  ===================================================

Gene regulatory networks
------------------------

=====================  ===================================================
Function               Purpose
=====================  ===================================================
``sc_grn``             network construction -- see the methods below
``sc_grnview``         display a network as an interactive graph
``sc_grnview2``        display two networks side by side
``sc_prs``             perturbation response scanning on an adjacency matrix
``sc_tenifoldnet``     scTenifoldNet pipeline
``sc_distcorr``        distance correlation for a gene set
``sc_resnet_causal2``  reservoir-based causal inference
``sc_causalcccnet``    local approximation of a causalCCC / MIIC network
=====================  ===================================================

``sc_grn`` methods: ``pcrnet`` (default), ``pcrnet_batch``,
``pcrnet_denoised``, ``genie3``, ``pearson``, ``xicor``, ``distcorr``, ``mi``,
``grnformer``, ``tn``. The older name ``pcnet`` is gone; use ``pcrnet``.
``sc_grn`` expects data that is already normalized and transformed -- no branch
normalizes for you. The implementations live in the ``+net/`` package.

Copy number and malignancy
--------------------------

=====================  ===================================================
Function               Purpose
=====================  ===================================================
``sc_infercnv``        large-scale CNV inference from expression
``sc_malignscore``     score and label cells as malignant from a CNV profile
=====================  ===================================================

Glycobiology
------------

=====================  ===================================================
Function               Purpose
=====================  ===================================================
``sc_glycostate``      per-cell glycobiological state
``sc_glycoccc``        glyco-lectin cell-cell communication
``sc_glycoshield``     cell-type-specific N-glycan shielding
``sc_glycoweight``     re-weight communication edges by glyco-state
=====================  ===================================================

Visualization
-------------

================================  ========================================
Function                          Purpose
================================  ========================================
``sc_scattermarker``              expression of one gene over an embedding
``gui.sc_scattergenes``           gene scatter: ``mean_cv``,
                                  ``meanlg_varlg``, ``mean_dropr``
``gui.sc_stem3``                  stem plot of selected genes across cells
``gui.i_hvgsplinefitplot``        3-D spline-fit scatter
``gui.i_plot_pseudotimeseries``   expression against pseudotime
``gui.sc_celltypeexplorer_auto``  cluster and explore cell types
================================  ========================================

The ``SingleCellExperiment`` class
----------------------------------

=========================  ===============================================
Method                     Purpose
=========================  ===============================================
``qcfilter``               QC filtering in place
``qcfilterwhitelist``      QC filtering that protects named genes
``embedcells``             compute an embedding
``clustercells``           cluster cells
``assigncelltype``         label clusters from PanglaoDB markers
``estimatecellcycle``      assign cell cycle phase
``estimatepotency``        assign differentiation potency
``sortcells``              reorder cells
``onestepraw2anno``        raw counts through to annotation in one call
``exportToJsonl``          export for external tools
``toSCE2``                 convert to ``SingleCellExperiment2``
=========================  ===============================================

Wrappers for external tools
---------------------------

``run.ml_*`` are MATLAB implementations bundled with the toolbox and need
nothing extra. ``run.py_*`` and ``run.r_*`` shell out to Python and R, which
must be installed and configured separately. ``run.web_*`` call web services.

**MATLAB** -- ``ml_alona``, ``ml_alona_new``, ``ml_cogaps``, ``ml_ComBat``,
``ml_diffuse``, ``ml_Enrichr``, ``ml_geneagent``, ``ml_GENIE3``,
``ml_Harmony``, ``ml_Harmony2``, ``ml_MAGIC``, ``ml_metaviz``, ``ml_PHATE``,
``ml_PickMarkers``, ``ml_SC3``, ``ml_scDock``, ``ml_SCEVAN``,
``ml_scGeneFit``, ``ml_SIMLR``, ``ml_SinNLRR``, ``ml_SnnDpc``, ``ml_SoptSC``,
``ml_talklr``, ``ml_TENET``, ``ml_TSCAN``, ``ml_UMAP``

**Python** -- ``py_cellbender``, ``py_GenKI``, ``py_geosketch``,
``py_GSEApy_enr``, ``py_harmonypy``, ``py_MELD``, ``py_memento``,
``py_panhumanpy``, ``py_scimilarity``, ``py_scrublet``,
``py_scTenifoldCko_gene``, ``py_scTenifoldCko_path``, ``py_scTenifoldXct``,
``py_SERGIO``, ``py_writeh5ad``

**R** -- ``r_clustermole``, ``r_cogaps``, ``r_copykat``, ``r_decontX``,
``r_DESeq2``, ``r_fgsea``, ``r_harmony``, ``r_infercnv``, ``r_MAST``,
``r_monocle3``, ``r_readSeuratRds``, ``r_saveSeuratRds``, ``r_SCEVAN``,
``r_seurat``, ``r_SeuratCellCycle``, ``r_SeuratSctransform``

**Web** -- ``web_Enrichr``, ``web_Enrichr_bkg``, ``web_STRING``

Packages
--------

=================  =======================================================
Package            Contents
=================  =======================================================
``+cli/``          the ``scgea`` command-line dispatcher
``+gui/``          GUI callbacks and plotting helpers
``+llm/``          LLM and MCP integration, including the GEOcellar agent
``+net/``          GRN inference algorithms behind ``sc_grn``
``+pkg/``          core computational and utility functions
``+run/``          wrappers for external tools
``+ten/``          scTenifold and tensor methods
``+gly/``          glycobiology methods
``+qtm/``          quantum and quantum-inspired methods
=================  =======================================================
