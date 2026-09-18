SCGEATOOL
===========

SCGEATOOL is a lightweight and blazing fast MATLAB application that provides interactive visualization functionality to analyze single-cell transcriptomic data. SCGEATOOL allows you to easily interrogate different views of your scRNA-seq data to quickly gain insights into the underlying biology.

Overview
--------
In MATLAB, the ``scgeatool`` function starts SCGEATOOL on a ``SingleCellExperiment`` (SCE) object. See :ref:`case_studies` for worked examples. All the examples below use data bundled with the toolbox.

|Overview of scgeatool|

.. |Overview of scgeatool| image:: https://github.com/jamesjcai/scGEAToolbox/raw/main/assets/Images/Tooltips.png
   :target: https://github.com/jamesjcai/scGEAToolbox/raw/main/assets/Images/Tooltips.png
  
Using SCGEATOOL to explore
-----------------------------
For a quick exploratory data analysis using `scgeatool` function

.. code-block:: matlab

  cdgea;
  load example_data/testXgs.mat
  sce=SingleCellExperiment(X,g,s);
  scgeatool(sce);
  
where X is the expression matrix, g is the list of genes, and s is the coordinates of embedding.

You can also load an example SCE (`SingleCellExperiment` object) variable using the following code:

.. code-block:: matlab

  cdgea;
  load example_data/testSce.mat
  scgeatool(sce);

If everything goes right, you will see the main inferface of SCGEATOOL like this:

|gui|

Making scRNA-seq data into `SCE`
--------------------------------
`SingleCellExperiment` stores the expression matrix and everything derived from
it in one object. Only two inputs are required: :math:`X`, the gene expression
matrix, and :math:`g`, the gene list. An embedding :math:`s` and a grouping
:math:`c` are optional.

.. code-block:: matlab

  sce = SingleCellExperiment(X, g);          % minimum
  sce = SingleCellExperiment(X, g, s);       % with a precomputed embedding
  sce = SingleCellExperiment(X, g, s, c);    % and a grouping

Any of the readers listed in :doc:`function_reference` can supply :math:`X` and
:math:`g`. If you omit :math:`s`, the constructor fills it with random
coordinates as a placeholder -- run ``sce.embedcells(...)`` before relying on
it, and before clustering on it.


.. |gui| image:: https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/main/assets/Images/scgeatool.png
   :width: 250
   :target: https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/main/assets/Images/scgeatool.png

SCGEATOOL standalone for Windows
--------------------------------
`SCGEATOOL standalone <https://scgeatool.github.io/>`__ is a lightweight and blazing fast desktop application that provides interactive visualization functionality to analyze single-cell transcriptomic data. SCGEATOOL allows you to easily interrogate different views of your scRNA-seq data to quickly gain insights into the underlying biology. SCGEATOOL is a pre-compiled standalone application developed in MATLAB. Pre-compiled standalone releases are meant for those environments without access to MATLAB licenses. Standalone releases provide access to all of the functionality of the SCGEATOOL standard MATLAB release encapsulated in a single application. SCGEATOOL is open-sourced to allow you to experience the added flexibility and speed of the MATLAB environment when needed.
