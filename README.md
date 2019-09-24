# scGEAToolbox - a Matlab toolbox for single-cell RNA-seq data analyses


![](https://github.com/jamesjcai/scGEAToolbox/blob/master/example_data/Fig_1.png?raw=true)
Fig. 1. A software workflow built with scGEAToolbox for single-cell gene regulatory network (scGRN) analyses. High-dimensional scRNA-seq data is filtered, normal-ized, and used as input for two paths. The first is a combination of (A) dimension-ality reduction and (B) trajectory/psedotime analysis to provide pseudotime-series data. The second is using network inference algorithms to generate (C) a global, coarse GRN structure. The integration of results from the two paths produces (D) pseudotime-series scGRNs, which can be further analyzed through regulatory modeling using parameter estimation algorithms to infer (E) a refined dynamic scGRN.


![](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/b282868e-8263-496e-b52a-29588685a3c6/b590d970-4886-4bc6-87f2-c07467dea006/images/screenshot.png)
Fig. 2. Screenshots of an execution of scGEApp -- the app interface of scGEAToolbox. (a) Two example panels of the main GUI scGEApp; (b) A 3-D scatter plot showing genes whose position is determined by µ, CV and rdrop; (c) A stem plot showing expression level of 50 selected genes across 2000 cells: 1000 in one state (blue) and the other 1000 in the other state (red).
