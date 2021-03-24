Code Formulas
=============

Example codes for common tasks.

Import files downloaded from GEO 
--------------------------------

.. code-block::

  mtxf='GSM3535276_AXLN1_matrix.mtx';
  genf='GSM3535276_AXLN1_genes.tsv';
  bcdf='GSM3535276_AXLN1_barcodes.tsv';
  [X,genelist,barcodelist]=sc_readmtxfile(mtxf,genf,bcdf,2);

If the barcodees.tsv is not available, then use the following

.. code-block::

  mtxf='GSM3535276_AXLN1_matrix.mtx';
  genf='GSM3535276_AXLN1_genes.tsv';
  [X,genelist]=sc_readmtxfile(mtxf,genf,[],2);


Process 10x Genomics expression matrix
--------------------------------------
In the 10x Genomics folder, there are three files, namely, matrix.mtx, features.tsv and barcodes.tsv. Here is the best practice of raw data processing.

.. code-block::
  
  [X,genelist,barcodelist]=sc_readmtxfile('matrix.mtx','features.tsv','barcodes.tsv',2);
  [X,genelist]=sc_qcfilter(X,genelist);
  [X,genelist]=sc_selectg(X,genelist,1,0.05);
  [s_tsne]=sc_tsne(X,3,false,false);
  sc_scatter(X,genelist,s_tsne)

t-SNE using top 2000 highly varible genes (HVGs)
------------------------------------------------

.. code-block::
  
  [Xnorm]=sc_norm(X,'type','deseq');
  [~,Xhvg]=sc_hvg(Xnorm,genelist,true);
  [s_tsne]=sc_tsne(Xhvg(1:2000,:),3,false,false);
  sc_scatter(X,genelist,s_tsne)
  
A complete pipeline of raw data processing
------------------------------------------

.. code-block::

  [X,genelist,barcodelist]=sc_readmtxfile('matrix.mtx','features.tsv','barcodes.tsv',2);
  [X,genelist]=sc_qcfilter(X,genelist);         % basic QC
  [X,genelist]=sc_selectg(X,genelist,1,0.05);   % select genes expressed in at least 5% of cells
  [Xnorm]=sc_norm(X,'type','deseq');            % normalize using DeSeq method
  [T,Xhvg]=sc_hvg(Xnorm,genelist,true);         % identify highly variable genes (HVGs) 
  [s_tsne]=sc_tsne(Xhvg(1:2000,:),3,false,false);   % using expression of top 2000 HVGs to t-SNE for cells
  sce=SingleCellExperiment(X,genelist,s_tsne);      % make SCE class
  sce=sce.estimatepotency(2);                   % estimate differentiation potency (1-human; 2-mouse)
  % sce=sce.estimatecellcycle;                  % estimate cell cycle phase. Need R/Seurat to be installed.
  id=sc_cluster_s(s_tsne,10);                   % clustering on the t-SNE coordinates using k-means
  sce.c_cluster_id=id;                          % assigning cluster Ids to SCE class
  sc_scatter(sce)                               % visualize cells  

