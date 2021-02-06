.. _getting_started:

Getting Started with `SC_SCATTER`
================================

This section shows how to quickly install `scGEAToolbox` and run `SC_SCATTER`--an interactive, explorary analysis tool for scRNA-seq data.

Qick installation
-----------------
Run the following code in `MATLAB`:

.. code-block::

  tic
  disp('Installing scGEAToolbox...')
  unzip('https://github.com/jamesjcai/scGEAToolbox/archive/master.zip');
  addpath('./scGEAToolbox-master');
  toc
  if exist('cdgea.m','file')
      disp('scGEAToolbox installed!')
  end

Run DEMO_GETTING_STARTED
------------------------

.. code-block::

 demo_Getting_Started


Using `SC_SCATTER`
------------------
For a quick exploratry data analysis using `SC_SCATTER` function

.. code-block::

  cdgea;
  load example_data\testXgs.mat
  sc_scatter(X,g,s)

If everything goes right, you will see the figure like this:
|gui|

Making SCE
----------
`scGEAToolbox` defines a Single-cell Experiment (SCE) class in order to store scRNA-seq data and variables. To make an SCE class, you need two variables: X and g, which are for gene expression matrix and gene list, respectively. 

.. code-block::

  cdgea;
  load example_data\testXgs.mat
  sce=SingleCellExperiment(X,g);
  sc_scatter_sce(sce)
  
.. |gui| image:: https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/master/resources/sc_scatter.png
   :target: https://twitter.com/hashtag/scGEAToolbox?src=hashtag_click

