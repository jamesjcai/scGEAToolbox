Test 1
======

# Demonstration of Filter, Normalization and Batch Correction of Data in scGEAToolbox
# Read scRNA-seq data, X and Y

```matlab:Code
cdgea; % set working directory
[X,genelistx]=sc_readfile('example_data/GSM3204304_P_P_Expr.csv');
```


```text:Output
Reading example_data/GSM3204304_P_P_Expr.csv ...... done.
```


```matlab:Code
[Y,genelisty]=sc_readfile('example_data/GSM3204305_P_N_Expr.csv');
```


```text:Output
Reading example_data/GSM3204305_P_N_Expr.csv ...... done.
```

# Select genes with at least 3 cells having more than 5 reads per cell.

```matlab:Code
[X,genelistx]=sc_selectg(X,genelistx,5,3);
[Y,genelisty]=sc_selectg(Y,genelisty,5,3);
```

# Obtain gene intersection of X and Y

```matlab:Code
[genelist,i,j]=intersect(genelistx,genelisty,'stable');
X=X(i,:);
Y=Y(j,:);
% libsizex=sum(X);
% libsizey=sum(Y);
% X=X(:,libsizex>quantile(libsizex,0.3)&libsizex<quantile(libsizex,0.95));
% Y=Y(:,libsizey>quantile(libsizey,0.3)&libsizey<quantile(libsizey,0.95));
clearvars -except X Y genelist
```

# Show raw data

```matlab:Code
figure; imagesc([X(1:100,1:500) Y(1:100,1:500)]); title('Raw Counts'); colorbar; xline(500,'y-');
```


![figure_0.eps](demo_script1_images/figure_0.eps)

# Show DESeq normalized data

```matlab:Code
[Xs]=sc_norm(X,'type','deseq');
[Ys]=sc_norm(Y,'type','deseq');
figure; imagesc([Xs(1:100,1:500) Ys(1:100,1:500)]); title('DESeq Normalized'); colorbar; xline(500,'y-');
```


![figure_1.eps](demo_script1_images/figure_1.eps)

# Show library-size normalized data

```matlab:Code
[X]=sc_norm(X,'type','libsize');
[Y]=sc_norm(Y,'type','libsize');
figure; imagesc([X(1:100,1:500) Y(1:100,1:500)]); title('Library Size Normalized'); colorbar; xline(500,'y-');
```


![figure_2.eps](demo_script1_images/figure_2.eps)

# Log(x+1) transformed normalized data

```matlab:Code
X=log(X+1);
Y=log(Y+1);
figure; imagesc([X(1:100,1:500) Y(1:100,1:500)]); title('Log(x+1) Transformed'); colorbar; xline(500,'y-');
```


![figure_3.eps](demo_script1_images/figure_3.eps)

# Show data subject to MAGIC imputation

```matlab:Code
Xo=run.MAGIC(X);
```


```text:Output
doing PCA
computing kernel
Computing alpha decay kernel:
Number of samples = 835
First iteration: k = 300
Number of samples below the threshold from 1st iter: 790
Using radius based search for the rest
   Symmetrize affinities
   Done computing kernel
imputing using optimal t
t = 1
t = 2
t = 3
t = 4
t = 5
t = 6
t = 7
t = 8
t = 9
optimal t = 9
done.
```


```matlab:Code
Yo=run.MAGIC(Y);
```


```text:Output
doing PCA
computing kernel
Computing alpha decay kernel:
Number of samples = 644
First iteration: k = 300
Number of samples below the threshold from 1st iter: 620
Using radius based search for the rest
   Symmetrize affinities
   Done computing kernel
imputing using optimal t
t = 1
t = 2
t = 3
t = 4
t = 5
t = 6
t = 7
t = 8
t = 9
optimal t = 9
done.
```


```matlab:Code
figure; imagesc([Xo(1:100,1:500) Yo(1:100,1:500)]); title('MAGIC Imputated'); colorbar; xline(500,'y-');
```


![figure_4.eps](demo_script1_images/figure_4.eps)

# Show HCP normalized data

```matlab:Code
[Xm,Ym]=run.HCP(X,Y);
figure; imagesc([Xm(1:100,1:500) Ym(1:100,1:500)]); title('HCP Normalized'); colorbar; xline(500,'y-');
```


![figure_5.eps](demo_script1_images/figure_5.eps)

# Show data with ComBat batch correction

```matlab:Code
[Xn,Yn]=run.ComBat(X,Y);
```


```text:Output
[combat] Found 2 batches
[combat] Adjusting for 0 covariate(s) of covariate level(s)
[combat] Standardizing Data across features
[combat] Fitting L/S model and finding priors
[combat] Finding parametric adjustments
[combat] Adjusting the Data
```


```matlab:Code
figure; imagesc([Xn(1:100,1:500) Yn(1:100,1:500)]); title('ComBat Corrected'); colorbar; xline(500,'y-');
```


![figure_6.eps](demo_script1_images/figure_6.eps)

# Visulize cells before and after ComBat batch correction

```matlab:Code
batchidx=[1*ones(size(X,2),1); 2*ones(size(Y,2),1)];

figure;
subplot(2,2,1)
[s]=sc_tsne([X Y]);
gscatter(s(:,1),s(:,2),batchidx,'','',5);
subplot(2,2,2)
[s]=sc_tsne([Xn Yn]);
gscatter(s(:,1),s(:,2),batchidx,'','',5);
```


![figure_7.eps](demo_script1_images/figure_7.eps)

# The End
