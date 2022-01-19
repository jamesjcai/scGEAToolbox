`SCGEATOOL`
===========

Overview of `SCGEATOOL`
----------------------
`SCGEATOOL` can be used to visualize `SCE` class/object. Below are links to several case studies and examples using sc_scatter_sce to explore high-dimensional data. All examples are below are publically available through GitHub.

|Overview of scgeatool|

.. |Overview of scgeatool| image:: https://github.com/jamesjcai/scGEAToolbox/raw/main/resources/Tooltips.png
   :target: https://github.com/jamesjcai/scGEAToolbox/raw/main/resources/Tooltips.png
  
Using `SCGEATOOL` to explore
-----------------------------
For a quick exploratry data analysis using `SCGEATOOL` function

.. code-block::

  cdgea;
  load example_data\testXgs.mat
  scgeatool(X,g,s)

If everything goes right, you will see the figure like this:

|gui|

Making scRNA-seq data into `SCE`
--------------------------------
`SingleCellExperiment` defines a Single-cell Experiment (SCE) class in order to store scRNAseq data and variables. To make an SCE class, you need two variables: :math:`X` and :math:`g`, which are gene expression matrix and gene list, respectively. 

.. code-block::

  cdgea;
  load example_data\testXgs.mat
  sce=SingleCellExperiment(X,g,s);
  scgeatool(sce)
  
.. |gui| image:: https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/main/resources/sc_scatter.png
   :width: 250
   :target: https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/main/resources/sc_scatter.png

