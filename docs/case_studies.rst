.. _case_studies:

Case Studies and Tutorials
--------------------------
sc_scatter_sce was developed to visualize high-dimensional biological data (e.g. genome-wide expression data), but it can also generally be applied to any high-dimensional data. Below are links to several case studies and examples using sc_scatter_sce to explore high-dimensional data. All examples are below are publically available through GitHub.

Process 10x Genomics raw data
=============================
In the 10x Genomics folder, there are three files, namely, matrix.mtx, features.tsv and barcodes.tsv. Here is the best practice of raw data processing.

.. code-block:: matlab
  
  [X,genelist,barcodelist]=sc_readmtxfile('matrix.mtx','features.tsv','barcodes.tsv',2);
  [X,genelist]=sc_qcfilter(X,genelist);
  [X,genelist]=sc_selectg(X,genelist,1,0.05);
  [Xnorm]=sc_norm(X,'type','deseq');
  [~,Xhvg]=sc_hvg(Xnorm,genelist,false);
  [s_tsne]=sc_tsne(Xhvg(1:2000,:),3,false,false);
  sc_scatter(X,genelist,s_tsne)

Label cell type interactively
=============================
.. raw:: html

         <div style="position: relative; padding-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">

         <iframe width="560" height="315" src="https://www.youtube.com/embed/HRQiXX3Jwpg" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
         </div>

We used :ref:`clustergrammer2`, the plotting library `bqplot`_, the Jupyter dashboard library `voila`_, and the Jupyter notebook hosting service `Binder`_ to build an interactive data exploration dashboard for `Visium`_ data from the mouse brain from `10X Genomics`_ (try dashboard: `Visium-Clustergrammer2 Dashboard`_, code: `https://github.com/ismms-himc/visium-clustergrammer2`_). This dashboard generates linked views of spatial tissue data and high-dimensional gene expression data - see GitHub repo `https://github.com/ismms-himc/visium-clustergrammer2`_ for more information.


.. _athero_plaques:

Single-cell Gene Expression and Proteomics from Human Atherosclerotic Plaques
==============================================================================

.. figure:: _static/chiara_citeseq_adt.gif
  :width: 450px
  :align: left
  :alt: CyTOF Screenshot
  :target: http://nbviewer.jupyter.org/github/MaayanLab/Cytof_Plasma_PMA/blob/master/notebooks/Plasma_vs_PMA_Phosphorylation.ipynb

  Immune Profiling of Atherosclerotic Plaques Identifies Innate and Adaptive Dysregulations Associated with Ischemic Cerebrovascular Events (`Fernandez et al.`_).



Our collaborators in the `Giannarelli Lab`_ used single-cell proteomics and transcriptomics to investigate the immune landscape of atheroscelerotic plaques (`Fernandez et al.`_) and identify features of T cells and macrophages that were associated with clinical symptomatic disease state. We used Clustergrammer2 to analyze scRNA-seq and CITE-seq data as well as infer cell-cell communication pathways. Interactive notebooks can be found in the Giannarelli lab GitHub repo: `Single-Cell-Immune-Profiling-of-Atherosclerotic-Plaques`_.


.. _clustergrammer2_2700:

scRNA-seq Gene Expression 2,700 PBMC
=======================================

|MyBinder-scRNA-seq| |NBViewer-scRNA-seq|

.. raw:: html

         <div style="position: relative; padding-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">

         <iframe width="560" height="315" src="https://www.youtube.com/embed/BEPspcC7vIY" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
         </div>

Single cell RNA-seq (scRNA-seq) is a powerful method to interrogate gene expression across thousands of single cells. This method provides thousands of measurements (single cells) across thousands of dimensions (genes). Clustergrammer2 is used to interactively explore an example dataset of 2,700 PBMCs obtained from `10X Genomics`_. Bulk gene expression signatures of cell types obtained from `CIBERSORT`_ were used to obtain a tentative cell type for each cell. The data and code can be found on GitHub at `clustergrammer2-notebooks`_ and the notebook can be viewed and re-run on the cloud - see below.

.. _clustergrammer2_citeseq_7800:

CITE-seq 7,800 PBMC
=======================================

|MyBinder-CITE-seq| |NBViewer-CITE-seq|

.. raw:: html

         <div style="position: relative; padding-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">

         <iframe width="560" height="315" src="https://www.youtube.com/embed/oG9TunM1Bvw" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
         </div>

CITE-seq (a.k.a feature barcoding from 10X genomics) is a new method that enabels researchers to simultaneously measure gene expression and protein levels in single cells. This notebook uses Clustergrammer2 to interactively explore an example dataset measuring the gene expression and surface marker proteins of 7,800 PBMCs obtained from 10X Genomics. Cell type was assigned based on unbiased hierarchical clustering of cells in surface marker space (ADTs) and transferred to cells in gene expression space. Please see the video tutorial above for more information.


Mouse Organogenesis Cell Atlas 2 Million Cells
==================================================
|MyBinder-Mouse-Atlas| |NBViewer-Mouse-Atlas|

.. raw:: html

         <div style="position: relative; padding-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">

         <iframe width="560" height="315" src="https://www.youtube.com/embed/pNPxuPDavO8" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
         </div>

`Cao, J and Spielmann, M et al`_ profiled gene expression from ~2 million mouse cells between 9.5 and 13.5 days of gestation. They identified 38 major cell types and measured ~25,000 genes. We generated a downsampled view of this data representing the ~1.3 million single cells (excluding ~600K suspected doublets) in the dataset by averaging expression for each cell type in each embryo, resulting in ~2,000 cell-type and embryo representative clusters. We demonstrate how Clustergrammer2 can be used to explore cell type clustering, find genes associated with cell type clusters, as well as identify genes that are differentially regulated across developmental stage. For more information, see the video tutorial above and launch or view the notebook using the badges.

CODEX Single Cell Multiplexed Imaging Dashboard
=================================================
|MyBinder-Codex|

.. raw:: html

         <div style="position: relative; padding-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">

         <iframe width="560" height="315" src="https://www.youtube.com/embed/JlUvt4rpF-s" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
         </div>


`Goltsev et al`_ used a highly multiplexed cytometric approach called CODEX to measure ~30 surface markers in spatially resolved single cells from mouse spleens. We utilized Clustergrammer2 to hierarchically cluster ~5,000 sinlge cells (from a subset of a segmented spleen image). We also used the Jupyter Widget `bqplot`_ to visualize single cell location data using voronoi plots. We then built a dasnboard using the library `voila`_, which converts Jupyter notebooks to dashboards/web-apps, and linked our heatmap to the spatial map. This allows to interact with the Clustergrammer2 heatmap and highlight cells in the spatially resolved map. These kind of linked views are crucial for exploration of spatially resolved high-dimensional single cell data. Finally, we are running this dashboard using MyBinder. See `CODEX Dashboard`_ for code.


.. _clustergrammer2_CCLE:

Cancer Cell Line Encyclopedia Gene Expression Data
==================================================
|MyBinder-CCLE| |NBViewer-CCLE|

.. raw:: html

         <div style="position: relative; padding-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">

         <iframe width="560" height="315" src="https://www.youtube.com/embed/6wZ0E6Veod0" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
         </div>

The Cancer Cell Line Encyclopedia (`CCLE`_) is a publicly available project that has characterized (e.g. genetic characterization) over 1,000 cancer cell lines. We used Clustergrammer to re-analyze and visualize CCLE's gene expression data in the `CCLE Explorer`_. The CCLE Explorer allows users to explore the CCLE by tissue type and visualize the most commonly differentially expressed genes for each tissue type as an interactive heatmap. The `CCLE Jupyter Notebook`_ generates an overview of the CCLE gene expression data, investigates specific tissues, and explains how to use :ref:`Enrichrgram <enrichrgram>` to understand the biological functions of differentially expressed genes.


Lung Cancer Post-Translational Modification and Gene Expression Regulation
==========================================================================

.. figure:: _static/CST_screenshot.png
  :width: 450px
  :align: left
  :alt: CST Screenshot
  :target: http://nbviewer.jupyter.org/github/MaayanLab/CST_Lung_Cancer_Viz/blob/master/notebooks/CST_Data_Viz.ipynb

  Screenshot from the `CST_Data_Viz.ipynb`_ Jupyter notebook showing hierarchical clustering of differential phosphorylation, methylation, acetylation, and gene expression data across 37 lung cancer cell lines. See the interactive Jupyter notebook `CST_Data_Viz.ipynb`_ for more information.

Lung cancer is a complex disease that is known to be regulated at the post-translational modification (PTM) level, e.g. phosphorylation driven by kinases. Our collaborators at `Cell Signaling Technology Inc`_ used Tandem Mass Tag (TMT) mass spectrometry to measure differential phosphorylation, acetylation, and methylation in a panel of 42 lung cancer cell lines compared to non-cancerous lung tissue. Gene expression data from 37 of these lung cancer cell lines was also independently obtained from the publicly available Cancer Cell Line Encyclopedia (`CCLE`_). In the Jupyter notebook `CST_Data_Viz.ipynb`_ we:

- Visualize PTM data, gene expression data, and merged PTM/gene-expression data
- Identify co-regulated clusters of PTMs/genes in distinct lung cancer cell line subtypes
- Perform enrichment analysis to understand the biological processes involved in PTM/expression clusters



Large Network: Kinase Substrate Similarity Network
==================================================
.. figure:: _static/kinase_network_screenshot.png
  :width: 450px
  :align: left
  :alt: Kinase Network Screenshot
  :target: https://maayanlab.github.io/kinase_substrate_similarity_network/

  Screenshot from the `Kinase Substrate Similarity Network`_ example that demonstrates how Clustergrammer can be used to visualize a large network of kinases based on shared substrates.

Clustergrammer can be used to visualize large networks without the formation of 'hairballs'. In the `Kinase Substrate Similarity Network`_ example we use Clustergrammer to visualize a network kinases based on shared substrate that includes 404 kinases and 163,216 links. Kinases are shown as rows and columns. For more information see the `Kinase Substrate Similarity Network`_ example.

Machine Learning and Miscellaneous Datasets
===========================================
.. figure:: _static/MNIST_screenshot.png
  :width: 450px
  :align: left
  :alt: MNIST Screenshot
  :target: http://nbviewer.jupyter.org/github/MaayanLab/MNIST_heatmaps/blob/master/notebooks/MNIST_Notebook.ipynb#Visualize-Downsampled-Version-of-MNIST

  Screenshot from the `MNIST Notebook`_ that demonstrates how the :ref:`clustergrammer_widget` can be used to visualize the `MNIST Data`_. Downsampled handwritten digits (K-means downsampled from 70,0000 handwritten digits to 300 digit-clusters) are shown as columns with digit-type categories and pixels are shown as rows. For more information see the `MNIST Notebook`_.

Clustergrammer was used to visualize several widely used machine learning Datasets and other miscellaneous Datasets:

- `MNIST Handwritten Digit Dataset`_
- `Iris Flower Dataset`_
- `USDA Nutrient Dataset`_

These examples demonstrate the generality of heatmap visualizations and enable users to interactively explore familiar Datasets.


.. _`Kinase Substrate Similarity Network`: https://maayanlab.github.io/kinase_substrate_similarity_network/
.. _`MNIST Data`: http://yann.lecun.com/exdb/mnist/

.. _`Giannarelli Lab`: http://labs.icahn.mssm.edu/giannarellilab/
.. _`Fernandez et al.`: https://www.nature.com/articles/s41591-019-0590-4
.. _`Single-Cell-Immune-Profiling-of-Atherosclerotic-Plaques`: https://github.com/giannarelli-lab/Single-Cell-Immune-Profiling-of-Atherosclerotic-Plaques

.. _`bqplot`: https://github.com/bloomberg/bqplot
.. _`Binder`: https://mybinder.org/
.. _`https://github.com/ismms-himc/visium-clustergrammer2`: https://github.com/ismms-himc/visium-clustergrammer2
.. _`Visium`: https://www.10xgenomics.com/spatial-transcriptomics/
.. _`Visium-Clustergrammer2 Dashboard`: http://bit.ly/visium-clustergrammer2

.. _`Icahn School of Medicine Human Immune Monitoring Core`: http://icahn.mssm.edu/research/portal/resources/deans-cores/human-immune-monitoring-core
.. _`CST_Data_Viz.ipynb`: http://nbviewer.jupyter.org/github/MaayanLab/CST_Lung_Cancer_Viz/blob/master/notebooks/CST_Data_Viz.ipynb
.. _`Cell Signaling Technology Inc`: https://www.cellsignal.com/
.. _`CCLE Explorer`: http://amp.pharm.mssm.edu/clustergrammer/CCLE/
.. _`CCLE Jupyter Notebook`: http://nbviewer.jupyter.org/github/MaayanLab/CCLE_Clustergrammer/blob/master/notebooks/Clustergrammer_CCLE_Notebook.ipynb
.. _`Iris Flower Dataset`: http://nbviewer.jupyter.org/github/MaayanLab/iris_clustergrammer_visualization/blob/master/Iris%20Dataset.ipynb
.. _`MNIST Notebook`: http://nbviewer.jupyter.org/github/MaayanLab/MNIST_heatmaps/blob/master/notebooks/MNIST_Notebook.ipynb
.. _`MNIST Handwritten Digit Dataset`: http://nbviewer.jupyter.org/github/MaayanLab/MNIST_heatmaps/blob/master/notebooks/MNIST_Notebook.ipynb
.. _`CCLE`: https://portals.broadinstitute.org/ccle/home
.. _`USDA Nutrient Dataset`: http://nbviewer.jupyter.org/github/MaayanLab/USDA_Nutrients_Viz/blob/master/USDA_Nutrients.ipynb
.. _`10X Genomics`: https://www.10xgenomics.com/resources/datasets/
.. _`CIBERSORT`: https://cibersort.stanford.edu/
.. _`clustergrammer2-notebooks`: https://github.com/ismms-himc/clustergrammer2-notebooks
.. _`MyBinder`: https://gke.mybinder.org/

.. _`Cao, J and Spielmann, M et al`: https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/landing:

.. _`Goltsev et al`: https://linkinghub.elsevier.com/retrieve/pii/S0092867418309048

.. _`bqplot`: https://github.com/bloomberg/bqplot
.. _`voila`: https://github.com/QuantStack/voila
.. _`CODEX Dashboard`: https://github.com/ismms-himc/codex_dashboard



.. |visium-clustergrammer2| image:: https://mybinder.org/badge_logo.svg?style=flat
    :alt: visium-clustergrammer2
    :scale: 100%
    :target: http://bit.ly/visium-clustergrammer2
