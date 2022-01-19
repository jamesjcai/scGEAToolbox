Code Formulas
=============

Example codes for common tasks.

Import 10x Genomics files
-------------------------
In the 10x Genomics folder, there are three files, namely, matrix.mtx, features.tsv (or genes.tsv) and barcodes.tsv. Here is how to import them:

.. code-block:: matlab

  mtxf='GSM3535276_AXLN1_matrix.mtx';
  genf='GSM3535276_AXLN1_genes.tsv';
  bcdf='GSM3535276_AXLN1_barcodes.tsv';
  [X,genelist,barcodelist]=sc_readmtxfile(mtxf,genf,bcdf,2);

If the barcodees.tsv is not available, then use the following

.. code-block:: matlab

  mtxf='GSM3535276_AXLN1_matrix.mtx';
  genf='GSM3535276_AXLN1_genes.tsv';
  [X,g]=sc_readmtxfile(mtxf,genf,[],2);


Process expression matrix, `X` and gene list, `g`
-------------------------------------------------------
Here is an example of raw data processing.

.. code-block:: matlab
  
  [X,g,b]=sc_readmtxfile('matrix.mtx','features.tsv','barcodes.tsv',2);
  [X,g]=sc_qcfilter(X,g);
  [X,g]=sc_selectg(X,g,1,0.05);
  [s]=sc_tsne(X);
  scgeatool(X,g,s)

t-SNE embedding of cells using highly varible genes (HVGs)
----------------------------------------------------------

.. code-block:: matlab
  
  [~,Xhvg]=sc_hvg(X,g);
  [s]=sc_tsne(Xhvg(1:2000,:));
  scgeatool(X,g,s)
  
An example pipeline for raw data processing
-------------------------------------------

.. code-block:: matlab

  [X,g]=sc_readmtxfile('matrix.mtx','features.tsv');
  [X,g]=sc_qcfilter(X,g);                   % basic QC
  [X,g]=sc_selectg(X,g,1,0.05);             % select genes expressed in at least 5% of cells
  [~,Xhvg]=sc_hvg(X,g);                     % identify highly variable genes (HVGs) 
  [s]=sc_tsne(Xhvg(1:2000,:));              % using expression of top 2000 HVGs for tSNE
  sce=SingleCellExperiment(X,g,s);          % make SCE class
  sce=sce.estimatepotency(2);               % estimate differentiation potency (1-human; 2-mouse)
  sce=sce.estimatecellcycle;                % estimate cell cycle phase
  id=sc_cluster_s(s,10);                    % clustering on tSNE coordinates using k-means
  sce.c_cluster_id=id;                      % assigning cluster Ids to SCE class
  scgeatool(sce)                            % visualize cells  

An example pipeline for processing 10x data folder
--------------------------------------------------
Assuming the .m file containing the following code is in the folder ./filtered_feature_bc_matrix. In this folder, three files: matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz, are present.

.. code-block:: matlab

  [X,genelist,celllist]=sc_read10xdir(pwd);
  sce=SingleCellExperiment(X,genelist);
  sce.c_cell_id=celllist;
  sce=sce.qcfilter;
  sce=sce.estimatecellcycle;
  sce=sce.estimatepotency("mouse");
  sce=sce.embedcells('tSNE',true);
  save clean_data sce -v7.3
  scgeatool(sce)

Merge two data sets (WT and KO)
-------------------------------

.. code-block:: matlab

  load WT/clean_data.mat sce
  sce_wt=sce;
  load KO/clean_data.mat sce
  sce_ko=sce;
  sce=sc_mergesces({sce_wt,sce_ko},'union');    % use parameter 'union' or 'intersect' to merge genes
  sce.c=sce.c_batch_id;
  scgeatool(sce)                                % blue - WT and red - KO  
  
You may want to re-compute tSNE coordinates after merging.
