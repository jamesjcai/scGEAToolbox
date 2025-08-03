scGEAToolbox - a Matlab toolbox for single-cell RNA-seq data analyses
---------------------------------------------------------------------

[![View scGEAToolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/72917-scgeatoolbox-single-cell-gene-expression-analysis-toolbox)

## Installation Guide

There are three ways to install `scGEAToolbox`:  

### :: Method 1: Install as a MATLAB Add-On (Recommended)  
1. **Open the Add-On Explorer**
   - In MATLAB, go to the **Home** tab.  
   - Click on the **Add-Ons** icon to open the Add-On Explorer.

2. **Search for scGEAToolbox**
   - In the search bar, type **"scGEAToolbox"** and press Enter.

3. **Select and Install**
   - Find **"scGEAToolbox (single-cell Gene Expression Analysis Toolbox)"** in the results.  
   - Click the **Add** button to install the toolbox.

5. **Launch scGEAToolbox**  
   - To start using scGEAToolbox, enter:
     ```matlab
     scgeatool
     ```
   - This will open the graphical user interface for single-cell transcriptomic data analysis.

### :: Method 2: Install via MATLAB Add-On (`.mltbx`)  
1. **Download** `scGEAToolbox.mltbx` from:  
   - 🔗 [GitHub Releases](https://github.com/jamesjcai/scGEAToolbox/releases)  
   - 🔗 [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/72917-scgeatoolbox-single-cell-gene-expression-analysis-toolbox)  

2. **Install** the toolbox:  
   - **Double-click** `scGEAToolbox.mltbx`  
   - MATLAB will open the **Add-On Manager**  
   - Click **Install**  

3. **Verify the installation** by running:  
   ```matlab
   scgeatool
   ```

### :: Method 3: Install via GitHub (For Developers)   
1. **Run** the following code in MATLAB:  
    ```matlab
    unzip('https://github.com/jamesjcai/scGEAToolbox/archive/main.zip');  
    ```
2. **Add** the toolbox to the MATLAB path
     ```matlab
     addpath('./scGEAToolbox-main');
     savepath(fullfile(userpath,'pathdef.m'));
     % savepath;
     ```

## Getting Started Quickly Using the SCGEATOOL App
Run the following code in MATLAB to start
```matlab
scgeatool
```

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=jamesjcai/scGEAToolbox&file=online_landing.m)
 
## Read Documentation

[scGEAToolbox Documentation](https://scgeatoolbox.readthedocs.io/)

<!---
## To get started immediately, check out demo scripts:

* [Demo script 1](http://htmlpreview.github.io/?https://github.com/jamesjcai/scGEAToolbox/blob/main/demo_script1.html) Filter, Normalization and Batch Correction of Data
* [Demo script 2](http://htmlpreview.github.io/?https://github.com/jamesjcai/scGEAToolbox/blob/main/demo_script2.html) Feature Selection
* [Demo script 3](http://htmlpreview.github.io/?https://github.com/jamesjcai/scGEAToolbox/blob/main/demo_script3.html) Visualization
* [Demo script 4](http://htmlpreview.github.io/?https://github.com/jamesjcai/scGEAToolbox/blob/main/demo_script4.html) Clustering
* [Demo script 5](http://htmlpreview.github.io/?https://github.com/jamesjcai/scGEAToolbox/blob/main/demo_script5.html) Pseudotime Analysis and Gene Network 
* [Demo script 6](http://htmlpreview.github.io/?https://github.com/jamesjcai/scGEAToolbox/blob/main/demo_script6.html) DE Analysis and Marker Gene Identification

## GUI interface

After installing the toolbox, the main GUI can be run by calling `scGEApp`. 
![](https://github.com/jamesjcai/scGEAToolbox/blob/main/example_data/Fig_2.png?raw=true)
**Fig. 1. Screenshots of an execution of scGEApp -- the app interface of scGEAToolbox.** (a) Two example panels of the main GUI scGEApp; (b) A 3-D scatter plot showing genes whose position is determined by expression mean, CV and dropout rate; (c) A stem plot showing expression level of 50 selected genes across 2,000 cells: 1,000 in one state (blue) and the other 1,000 in the other state (red).

## Analytical workflow built with scGEAToolbox

![](https://github.com/jamesjcai/scGEAToolbox/blob/main/example_data/Fig_1.png?raw=true)  
**Fig. 2. A software workflow built with scGEAToolbox for single-cell gene regulatory network (scGRN) analyses.** High-dimensional scRNA-seq data is filtered, normalized, and used as input for two paths. The first is a combination of (A) dimension-ality reduction and (B) trajectory/psedotime analysis to provide pseudotime-series data. The second is using network inference algorithms to generate (C) a global, coarse GRN structure. The integration of results from the two paths produces (D) pseudotime-series scGRNs, which can be further analyzed through regulatory modeling using parameter estimation algorithms to infer (E) a refined dynamic scGRN.

## Interactive cell type annotation with scGEAToolbox (sc_celltype explorer)

[![scGEAToolbox sc_celltypeexplorer - interactive cell type annotation](https://img.youtube.com/vi/HRQiXX3Jwpg/0.jpg)](https://youtu.be/HRQiXX3Jwpg)
-->

## SEE ALSO - Standalone Application

[SCGEATOOL.exe](https://scgeatool.github.io/) is a standalone application running on Windows machines that do not have MATLAB installed. SCGEATOOL is a lightweight and blazing fast desktop application that provides interactive visualization functionality to analyze single-cell transcriptomic data. SCGEATOOL allows you to easily interrogate different views of your scRNA-seq data to quickly gain insights into the underlying biology. SCGEATOOL is a pre-compiled standalone application developed in MATLAB. Pre-compiled standalone releases are meant for those environments without access to MATLAB licenses. Standalone releases provide access to all of the functionality of the SCGEATOOL standard MATLAB release encapsulated in a single application. SCGEATOOL is open-sourced to allow you to experience the added flexibility and speed of the MATLAB environment when needed.

## Help

If you have any questions or require assistance using scGEAToolbox, please [contact me](https://scgeatool.github.io/#contact).

## Citation

**Cai JJ**. scGEAToolbox: a Matlab toolbox for single-cell RNA sequencing data analysis. *Bioinformatics*. 2019;btz830. [doi:10.1093/bioinformatics/btz830](https://doi.org/10.1093/bioinformatics/btz830)
