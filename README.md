scGEAToolbox - a Matlab toolbox for single-cell RNA-seq data analyses
---------------------------------------------------------------------

## Introduction
**Motivation**: Single-cell RNA sequencing (scRNA-seq) technology has revolutionized the way research is done in biomedical sciences. It provides an unprecedented level of resolution across individual cells for studying cell heterogeneity and gene expression variability. Analyzing scRNA-seq data is challenging though, due to the sparsity and high dimensionality of the data.
**Results**: I developed scGEAToolbox—a Matlab toolbox for scRNA-seq data analysis. It contains a comprehensive set of functions for data normalization, feature selection, batch correction, imputation, cell clustering, trajectory/pseudotime analysis, and network construction, which can be combined and integrated to building custom workflow. While most of the functions are implemented in native Matlab, wrapper functions are provided to allow users to call the “third-party” tools developed in Matlab or other languages. Furthermore, scGEAToolbox is equipped with sophisticated graphical user interfaces (GUIs) generated with App Designer, making it an easy-to-use application for quick data processing.
**Availability**: [https://github.com/jamesjcai/scGEAToolbox](https://github.com/jamesjcai/scGEAToolbox)
**Contact**: jcai@tamu.edu

## To get started immediately, check out demo scripts:

* [Demo script 1](http://htmlpreview.github.io/?https://github.com/jamesjcai/scGEAToolbox/blob/master/demo_script1.html)
* [Demo script 2](http://htmlpreview.github.io/?https://github.com/jamesjcai/scGEAToolbox/blob/master/demo_script2.html)
* [Demo script 3](http://htmlpreview.github.io/?https://github.com/jamesjcai/scGEAToolbox/blob/master/demo_script3.html)
* [Demo script 4](http://htmlpreview.github.io/?https://github.com/jamesjcai/scGEAToolbox/blob/master/demo_script4.html)
* [Demo script 5](http://htmlpreview.github.io/?https://github.com/jamesjcai/scGEAToolbox/blob/master/demo_script5.html)

## GUI interface

After installing the toolbox, the main GUI can be run by calling `scGEApp`.

![](https://github.com/jamesjcai/scGEAToolbox/blob/master/example_data/Fig_2.png?raw=true)
**Fig. 1. Screenshots of an execution of scGEApp -- the app interface of scGEAToolbox.** (a) Two example panels of the main GUI scGEApp; (b) A 3-D scatter plot showing genes whose position is determined by expression mean, CV and dropout rate; (c) A stem plot showing expression level of 50 selected genes across 2,000 cells: 1,000 in one state (blue) and the other 1,000 in the other state (red).

## Analytical workflow built with scGEAToolbox

![](https://github.com/jamesjcai/scGEAToolbox/blob/master/example_data/Fig_1.png?raw=true)
**Fig. 2. A software workflow built with scGEAToolbox for single-cell gene regulatory network (scGRN) analyses.** High-dimensional scRNA-seq data is filtered, normalized, and used as input for two paths. The first is a combination of (A) dimension-ality reduction and (B) trajectory/psedotime analysis to provide pseudotime-series data. The second is using network inference algorithms to generate (C) a global, coarse GRN structure. The integration of results from the two paths produces (D) pseudotime-series scGRNs, which can be further analyzed through regulatory modeling using parameter estimation algorithms to infer (E) a refined dynamic scGRN.


## Help

If you have any questions or require assistance using scGEAToolbox, please contact me at [scgeatoolbox.slack.com](https://join.slack.com/t/scgeatoolbox/shared_invite/enQtNTQ3NjQ0MjIwNjc4LTRiODU5MDI4N2RkNzNkYmZlYWViY2FmNDdhM2EwYzEwY2VjNzk0Y2YwOWUyNzYwNzNkYjEyZmI3M2Y3MjEwNWE)