SCGEATOOL
===========

SCGEATOOL is a lightweight and blazing fast MATLAB application that provides interactive visualization functionality to analyze single-cell transcriptomic data. SCGEATOOL allows you to easily interrogate different views of your scRNA-seq data to quickly gain insights into the underlying biology.

Overview
--------
In MATLAB, scgeatool function can be used to start SCGEATOOL to visualize `SCE` class/object. Below are links to several case studies and examples using the scgeatool function to explore scRNA-seq data. All examples are below are publically available through GitHub.

|Overview of scgeatool|

.. |Overview of scgeatool| image:: https://github.com/jamesjcai/scGEAToolbox/raw/main/resources/Images/Tooltips.png
   :target: https://github.com/jamesjcai/scGEAToolbox/raw/main/resources/Images/Tooltips.png
  
Using SCGEATOOL to explore
-----------------------------
For a quick exploratry data analysis using `scgeatool` function

.. code-block::

  cdgea;
  load example_data\testXgs.mat
  sce=SingleCellExperiment(X,g,s);
  scgeatool(X,g,s)
  
where X is the expression matrix, g is the list of genes, and s is the coordinates of embedding.

You can also load an example SCE (`SingleCellExperiment` object) variable using the following code:

.. code-block::

  cdgea;
  load example_data\testSce.mat
  scgeatool(sce)

If everything goes right, you will see the main inferface of SCGEATOOL like this:

|gui|

Making scRNA-seq data into `SCE`
--------------------------------
`SingleCellExperiment` defines a Single-cell Experiment (SCE) class in order to store scRNAseq data and variables. To make an SCE class, you need two variables: :math:`X` and :math:`g`, which are gene expression matrix and gene list, respectively. 

.. code-block::

  cdgea;
  load example_data\testXgs.mat
  sce=SingleCellExperiment(X,g,s);
  scgeatool(sce)
  
.. |gui| image:: https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/main/resources/Images/scgeatool.png
   :width: 250
   :target: https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/main/resources/Images/scgeatool.png

SCGEATOOL standalone for Windows
--------------------------------
`SCGEATOOL standalone <https://scgeatool.github.io/>`__ is a lightweight and blazing fast desktop application that provides interactive visualization functionality to analyze single-cell transcriptomic data. SCGEATOOL allows you to easily interrogate different views of your scRNA-seq data to quickly gain insights into the underlying biology. SCGEATOOL is a pre-compiled standalone application developed in MATLAB. Pre-compiled standalone releases are meant for those environments without access to MATLAB licenses. Standalone releases provide access to all of the functionality of the SCGEATOOL standard MATLAB release encapsulated in a single application. SCGEATOOL is open-sourced to allow you to experience the added flexibility and speed of the MATLAB environment when needed.
