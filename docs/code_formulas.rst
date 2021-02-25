Code Formulas
=============

Example codes for common tasks.

Process 10x Genomics raw data
-----------------------------
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
  [X,genelist]=sc_qcfilter(X,genelist);
  [X,genelist]=sc_selectg(X,genelist,1,0.05);
  [Xnorm]=sc_norm(X,'type','deseq');
  [~,Xhvg]=sc_hvg(Xnorm,genelist,true);
  [s_tsne]=sc_tsne(Xhvg(1:2000,:),3,false,false);
  sce=SingleCellExperiment(X,genelist,s_tsne);
  sce=sce.estimatepotency(2);
  sce=sce.estimatecellcycle;
  id=sc_cluster_s(s_tsne,10);
  sce.c_cluster_id=id;
  sc_scatter_sce(sce)
