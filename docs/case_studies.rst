.. _case_studies:

Case Studies and Tutorials
==========================
`sc_scatter_sce` can be used to visualize `SCE` class/object. Below are links to several case studies and examples using sc_scatter_sce to explore high-dimensional data. All examples are below are publically available through GitHub.

Process 10x Genomics raw data
-----------------------------
In the 10x Genomics folder, there are three files, namely, matrix.mtx, features.tsv and barcodes.tsv. Here is the best practice of raw data processing.

.. code-block::
  
  [X,genelist,barcodelist]=sc_readmtxfile('matrix.mtx','features.tsv','barcodes.tsv',2);
  [X,genelist]=sc_qcfilter(X,genelist);
  [X,genelist]=sc_selectg(X,genelist,1,0.05);
  [Xnorm]=sc_norm(X,'type','deseq');
  [~,Xhvg]=sc_hvg(Xnorm,genelist,true);
  [s_tsne]=sc_tsne(Xhvg(1:2000,:),3,false,false);
  sc_scatter(X,genelist,s_tsne)

Overview of `SC_SCATTER`
-----------------------

|Overview of sc_scatter|

Using SC_SCATTER_SCE to explore scRNA-seq data stored as SCE class
------------------------------------------------------------------
.. raw:: html

         <div style="position: relative; padding-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">

         <iframe width="560" height="315" src="https://www.youtube.com/embed/6IYl35dTwGw" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
         </div>


Label cell type interactively with SC_SCATTER_SCE
-------------------------------------------------
.. raw:: html

         <div style="position: relative; padding-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">

         <iframe width="560" height="315" src="https://www.youtube.com/embed/HRQiXX3Jwpg" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
         </div>

.. We used :ref:`clustergrammer2`, the plotting library `bqplot`_, the Jupyter dashboard library `voila`_, and the Jupyter notebook hosting service `Binder`_ to build an interactive data exploration dashboard for `Visium`_ data from the mouse brain from `10X Genomics`_ (try dashboard: `Visium-Clustergrammer2 Dashboard`_, code: `https://github.com/ismms-himc/visium-clustergrammer2`_). This dashboard generates linked views of spatial tissue data and high-dimensional gene expression data - see GitHub repo `https://github.com/ismms-himc/visium-clustergrammer2`_ for more information.


.. |Overview of sc_scatter| image:: https://pbs.twimg.com/media/ErAPaabW4AEwAvf?format=png&name=large
   :target: https://youtu.be/6IYl35dTwGw
   

