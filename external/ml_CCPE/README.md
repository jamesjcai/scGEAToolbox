# CCPE

## Introduction
CCPE is a cell cycle pseudotime estimation method characterizing cell cycle timing from single-cell RNA-seq data. CCPE maps high-dimensional scRNA-seq data onto a helix in three-dimensional space, where 2D space is used to capture the cycle information in scRNA-seq data, and one dimension to predict the chronological orders of cells along the cycle, which we called cell cycle pseudotime. ScRNA-seq data is repeatedly transformed from high dimensional to low dimensional and then mapped back to high dimensional. At the same time, CCPE iteratively optimizes the discriminative dimensionality reduction via learning a helix until convergence (Figure 1). CCPE is applied to several downstream analyses and applications to demonstrate its ability to accurately estimate the cell cycle pseudotime and stages.<br/>
![image](https://github.com/LiuJJ0327/CCPE/blob/main/images/figure1.PNG)<br/>

## Reference
Our paper was published in Nucleic Acids Research, available at [Jiajia Liu, Mengyuan Yang, Weiling Zhao, Xiaobo Zhou, CCPE: cell cycle pseudotime estimation for single cell RNA-seq data, Nucleic Acids Research, 2021](https://doi.org/10.1093/nar/gkab1236)

## Quick Start<br/>
drtoolbox was downloaded from https://lvdmaaten.github.io/drtoolbox/<br/>
```bash
wget https://github.com/LiuJJ0327/CCPE/archive/refs/heads/main.zip
unzip CCPE-main.zip
cd CCPE-main/
matlab run.m
```
### Tansfer pseudotime to discrete stages<br/>
```bash
Rscript pseudotime_to_label.r 
```

## Applications and analyses of CCPE<br/>
The reproduction of the applications and analyses in our paper can be available at `reproduction/` <br/>
### 1. Pseudotime analysis<br/>
```bash
Rscript reproduction/pseudotime_analysis/4_analysis_R/plot_pseudotime.r
```
![image](https://github.com/LiuJJ0327/CCPE/blob/main/images/1_pseudotime.PNG)<br/>

### 2. Gene expression analysis (example: Aurka)<br/>
```bash
Rscript reproduction/pseudotime_analysis/4_analysis_R/correlation.r
```
![image](https://github.com/LiuJJ0327/CCPE/blob/main/images/2_gene_expression_1.PNG)<br/>

```bash
library(Seurat)
mesc<-read.table("./pseudotime_analysis/1_mesc_Quartz_data/mESC_Quartz_preprocessed.txt",header=T,row.names = 1)
marrow <- CreateSeuratObject(raw.data  = mesc)
marrow@ident<-factor(marrow@ident,levels = c("G1","S","G2M"))
RidgePlot(object =marrow, features = 'Aurka')
```
![image](https://github.com/LiuJJ0327/CCPE/blob/main/images/2_gene_expression_2.PNG)<br/>

### 3. Evaluation of the performance using then clustering metrics (example: E-MTAB-2805 mESCs data)<br/>
```bash
matlab estimate_mtrics/Evaluate.m
```
plot evaluation result of CCPE
```bash
library(ggplot2)
metrics<-read.csv("./comparison/2_mesc_288/evaulation_mesc_avg.csv",header = T)
CCPE_metrics<-metrics[metrics$method=='CCPE',]
ggplot(CCPE_metrics,aes(metrics,value,fill=metrics))+geom_bar(stat="identity",position="dodge")
```
![image](https://github.com/LiuJJ0327/CCPE/blob/main/images/3_evaluation.PNG)<br/>

### 4. Differential gene expression analysis <br/>
Differentially expressed genes were identified using [Deseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)<br/>
Gene set enrichment analysis was using [Enrichr](https://maayanlab.cloud/Enrichr/)<br/>
![image](https://github.com/LiuJJ0327/CCPE/blob/main/images/4_DEG.PNG)<br/>

### 5. UMAP plot of simulated dataset with dropouts<br/>
```bash
Rscript UMAP_cancer_cellline.r
```
![image](https://github.com/LiuJJ0327/CCPE/blob/main/images/5_UMAP.PNG)<br/>

### 6. Cell cycle effect removal<br/>
```bash
Rscript reproduction/cell_cycle_effect_removal/cellcycle_removal.r
```
![image](https://github.com/LiuJJ0327/CCPE/blob/main/images/6_removal.PNG)<br/>

