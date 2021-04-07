`SC_SCATTER`
============

Overview of `SC_SCATTER`
-----------------------
`SC_SCATTER` can be used to visualize `SCE` class/object. Below are links to several case studies and examples using sc_scatter_sce to explore high-dimensional data. All examples are below are publically available through GitHub.

|Overview of sc_scatter|

.. |Overview of sc_scatter| image:: https://github.com/jamesjcai/scGEAToolbox/raw/master/resources/Tooltips.png
   :target: https://github.com/jamesjcai/scGEAToolbox/raw/master/resources/Tooltips.png
  
Using `SC_SCATTER` to explore
-----------------------------
For a quick exploratry data analysis using `SC_SCATTER` function

.. code-block::

  cdgea;
  load example_data\testXgs.mat
  sc_scatter(X,g,s)

If everything goes right, you will see the figure like this:

|gui|

Making scRNA-seq data into `SCE`
--------------------------------
`scGEAToolbox` defines a Single-cell Experiment (SCE) class in order to store scRNA-seq data and variables. To make an SCE class, you need two variables: :math:`X` and :math:`g`, which are gene expression matrix and gene list, respectively. 

.. code-block::

  cdgea;
  load example_data\testXgs.mat
  sce=SingleCellExperiment(X,g,s);
  sc_scatter(sce)
  
.. |gui| image:: https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/master/resources/sc_scatter.png
   :width: 250
   :target: https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/master/resources/sc_scatter.png

