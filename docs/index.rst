.. scGEAToolbox documentation master file, created by
   sphinx-quickstart on Sun Dec 20 17:44:55 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to scGEAToolbox's documentation!
========================================

Single-cell RNA sequencing (scRNA-seq) technology has revolutionized the way research is done in biomedical sciences. It provides an unprecedented level of resolution across individual cells for studying cell heterogeneity and gene expression variability. Analyzing scRNA-seq data is challenging though, due to the sparsity and high dimensionality of the data. scGEAToolbox is a MATLAB toolbox for scRNA-seq data analysis. It contains a comprehensive set of functions for quality control, normalization, feature selection, batch correction, imputation, cell clustering, trajectory/pseudotime analysis, cell type annotation, differential expression and variability testing, gene set scoring, and gene regulatory network construction, which can be combined into custom workflows. Most functions are implemented in native MATLAB; wrapper functions let you call third-party tools written in MATLAB, Python or R from the same session. scGEAToolbox is equipped with sophisticated graphical user interfaces (GUIs), making it an easy-to-use application for quick data processing.

Three ways to use it
--------------------

- **SCGEATOOL**, the App Designer GUI: run ``scgeatool`` in MATLAB. See
  :doc:`scgeatool`.
- **The function API**, at the matrix level or through the
  ``SingleCellExperiment`` class. See :doc:`getting_started` and
  :doc:`function_reference`.
- **The** ``scgea`` **command line**, for batch and pipeline work without
  opening MATLAB. See :doc:`cli`.

Official Websites and Social Networks
-------------------------------------

Please, visit the official website |View scGEAToolbox on File Exchange| of scGEAToolbox for further information.

.. |View scGEAToolbox on File Exchange| image:: https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg
   :target: https://www.mathworks.com/matlabcentral/fileexchange/72917-scgeatoolbox


.. toctree::
   :maxdepth: 2
   :caption: Main Documents

   quick_installation
   getting_started
   scgeatool
   cli
   code_formulas
   case_studies
   function_reference

.. toctree::
   :maxdepth: 2
   :caption: Social Networks
   
   social_networks

.. toctree::
   :maxdepth: 2
   :caption: Appendices
   
   publications
   license
