.. _getting_started:

Getting Started
===============

Run `DEMO_GETTING_STARTED`
------------------------

.. code-block::

 demo_Getting_Started


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
  sce=SingleCellExperiment(X,g);
  sc_scatter(sce)
  
.. |gui| image:: https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/master/resources/sc_scatter.png
   :width: 250
   :target: https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/master/resources/sc_scatter.png

