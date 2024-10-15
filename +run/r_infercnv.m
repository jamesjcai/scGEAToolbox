%{
https://academic.oup.com/bib/article/22/3/bbaa127/5867555

Cell malignancy estimation
In cancer research, estimating cell malignancy and distinguishing malignant and non-malignant cells is also a critical issue. Generally, copy number alterations (CNV) inferred from scRNA-seq data are potential to identify malignant tumor cells, which was widely used in the research of various cancer types [1, 3, 5, 6]. Here, we used the algorithm of R package infercnv [26] to estimate single-cell CNVs. The algorithm calculated the moving average of expression values across each chromosome and then compared them with some normal reference cells to estimate the CNVs. For the convenience of users, we curated 3366 precompiled normal microenvironmental cells (including immune cells, fibroblasts and endothelial cells) from a few high-quality samples as the default reference data set. Users can also set their own reference cells. Furthermore, considering the impact of dropout, we proposed to use cells neighbor information to smooth CNV values to estimate malignancy better. In detail, we took advantage of the similarity matrix of clustering to identify the neighbors 
 of each cell 
⁠, and used the normalized similarities as weights 
⁠. By calculating the weighted average of each cell 
 and its neighbors initial CNV values 
 and 
⁠, we got the final CNV profiles 
⁠:
 
 
 
Then, we defined the malignancy score as the mean of the squares of final CNV values across chromosome. By comparing the distribution of malignancy scores with reference and detecting its bimodality, we assigned malignancy labels for cells (Figure 4A and the supplementary materials section 6).

%}