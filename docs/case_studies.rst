.. _case_studies:

Case Studies and Tutorials
==========================

This page contains demo scripts illustrating the main features of scGEAToolbox. Each demo is self-contained and uses example data bundled with the toolbox.

.. contents:: On this page
   :local:
   :depth: 2

Demo 1: Filter, Normalization and Batch Correction of Data
-----------------------------------------------------------

Read two scRNA-seq data sets and apply gene selection, normalization, imputation, and batch correction.

**Read scRNA-seq data, X and Y**

.. code-block:: matlab

  cdgea;
  [X, genelistx] = sc_readfile('example_data/GSM3204304_P_P_Expr.csv');
  [Y, genelisty] = sc_readfile('example_data/GSM3204305_P_N_Expr.csv');

**Select genes with at least 3 cells having more than 5 reads per cell**

.. code-block:: matlab

  [X, genelistx] = sc_selectg(X, genelistx, 5, 3);
  [Y, genelisty] = sc_selectg(Y, genelisty, 5, 3);

**Obtain gene intersection of X and Y**

.. code-block:: matlab

  [genelist, i, j] = intersect(genelistx, genelisty, 'stable');
  X = X(i, :);
  Y = Y(j, :);
  clearvars -except X Y genelist

**Normalization methods**

.. code-block:: matlab

  % DESeq normalization
  [Xs] = sc_norm(X, 'type', 'deseq');
  [Ys] = sc_norm(Y, 'type', 'deseq');

  % Library-size normalization
  [X] = sc_norm(X, 'type', 'libsize');
  [Y] = sc_norm(Y, 'type', 'libsize');

  % Log(x+1) transformation
  X = log1p(X);
  Y = log1p(Y);

**MAGIC imputation**

.. code-block:: matlab

  Xo = run.ml_MAGIC(X);
  Yo = run.ml_MAGIC(Y);

**ComBat batch correction**

.. code-block:: matlab

  [Xn, Yn] = run.ml_ComBat(X, Y);

**Visualize cells before and after batch correction**

.. code-block:: matlab

  batchidx = [ones(size(X, 2), 1); 2*ones(size(Y, 2), 1)];

  figure;
  subplot(2, 2, 1)
  [s] = sc_tsne([X Y], 2, false, false);
  gscatter(s(:,1), s(:,2), batchidx, '', '', 5);

  subplot(2, 2, 2)
  [s] = sc_tsne([Xn Yn], 2, false, false);
  gscatter(s(:,1), s(:,2), batchidx, '', '', 5);


Demo 2: Feature Selection Functions
------------------------------------

Identify highly variable genes (HVGs) and differentially deviated (DD) genes using HVG analysis and spline-fit methods.

**HVG analysis with single data X**

.. code-block:: matlab

  cdgea;
  [X, genelist] = sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
  [X, genelist] = sc_selectg(X, genelist, 3, 1);

  % Normalize data with DESeq method
  Xn = sc_norm(X, 'type', 'deseq');
  [T] = sc_hvg(Xn, genelist, true, true);

  % Highly variable genes (HVGenes), FDR < 0.05
  HVGenes = T.genes(T.fdr < 0.05);
  disp(HVGenes(1:10))

**Spline-fit feature selection with single data X**

.. code-block:: matlab

  [X, genelist] = sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
  [X, genelist] = sc_selectg(X, genelist, 3, 1);

  sortit = true;
  [T1] = sc_splinefit(X, genelist, sortit);

  % Top 10 featured genes with highest deviation (D) values
  T1.genes(1:10)

  % Show data points and the spline-fit curve
  dofit = true;
  showdata = true;
  gui.i_hvgsplinefitplot(X, genelist, dofit, showdata);

**Analysis of differentially deviated (DD) genes with data X and Y**

.. code-block:: matlab

  [X, genelistx] = sc_readfile('example_data/GSM3204304_P_P_Expr.csv');
  [Y, genelisty] = sc_readfile('example_data/GSM3204305_P_N_Expr.csv');
  [X, genelistx] = sc_selectg(X, genelistx, 3, 1);
  [Y, genelisty] = sc_selectg(Y, genelisty, 3, 1);

  % Spline-fit plots for each data set
  gui.i_hvgsplinefitplot(X, genelistx, true, true);
  title('Data 1')

  gui.i_hvgsplinefitplot(Y, genelisty, true, true);
  title('Data 2')

  % Fit X and Y separately and obtain DD value of each gene
  [T2] = sc_splinefit2(X, Y, genelistx, genelisty, true);

  % Top 10 genes with highest DD value
  T2.genes(1:10)


Demo 3: Visualization Functions
--------------------------------

Dimensionality reduction and gene-level scatter plots.

**Load and pre-process three data sets**

.. code-block:: matlab

  cdgea;
  [X, genelistx] = sc_readfile('example_data/GSM3204304_P_P_Expr.csv');
  [Y, genelisty] = sc_readfile('example_data/GSM3204305_P_N_Expr.csv');
  [Z, genelistz] = sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
  [X, genelistx] = sc_selectg(X, genelistx, 3, 1);
  [Y, genelisty] = sc_selectg(Y, genelisty, 3, 1);
  [Z, genelistz] = sc_selectg(Z, genelistz, 3, 1);

**Intersection of common genes and remove mitochondrial genes**

.. code-block:: matlab

  [genelist] = intersect(intersect(genelistx, genelisty, 'stable'), genelistz, 'stable');
  i = startsWith(genelist, 'MT-');
  genelist(i) = [];
  [~, i1] = ismember(genelist, genelistx);
  [~, i2] = ismember(genelist, genelisty);
  [~, i3] = ismember(genelist, genelistz);
  X = X(i1, :);
  Y = Y(i2, :);
  Z = Z(i3, :);

  cellidx = [ones(size(X, 2), 1); 2*ones(size(Y, 2), 1); 3*ones(size(Z, 2), 1)];

**PCA, t-SNE, and Diffusion Map embeddings**

.. code-block:: matlab

  % PCA
  [~, s] = pca([X Y Z]');
  gscatter(s(:,1), s(:,2), cellidx, '', '', 8);

  % t-SNE
  [s] = sc_tsne([X Y Z], 2, false, true);
  gscatter(s(:,1), s(:,2), cellidx, '', '', 8);

  % Diffusion Map
  [s] = run.ml_diffuse([X Y Z]);
  gscatter(s(:,1), s(:,2), cellidx, '', '', 8);

**Gene-level scatter plots**

.. code-block:: matlab

  figure;
  gui.sc_scattergenes(X, genelistx, 'mean_cv');

  figure;
  gui.sc_scattergenes(X, genelistx, 'meanlg_varlg');

  figure;
  gui.sc_scattergenes(X, genelistx, 'mean_dropr');

  % 3D scatter plot with spline fit
  figure;
  gui.i_hvgsplinefitplot(X, genelistx, true, true);

**Feature selection: top 50 DD genes**

.. code-block:: matlab

  T = sc_splinefit2(X, Y, genelistx, genelisty);
  T = sortrows(T, size(T, 2), 'descend');
  [~, idx1] = ismember(table2array(T(:,1)), genelistx);
  [~, idx2] = ismember(table2array(T(:,1)), genelisty);
  figure;
  gui.sc_stem3(X(idx1, :), Y(idx2, :), genelistx(idx1), 50);


Demo 4: Clustering Functions
------------------------------

Cluster cells using SIMLR, SC3, and SoptSC, and compare results.

**Load and merge two data sets**

.. code-block:: matlab

  cdgea;
  [X, genelistx] = sc_readfile('example_data/GSM3204304_P_P_Expr.csv');
  [Y, genelisty] = sc_readfile('example_data/GSM3204305_P_N_Expr.csv');
  [X, genelistx] = sc_selectg(X, genelistx, 3, 1);
  [Y, genelisty] = sc_selectg(Y, genelisty, 3, 1);

  [genelist] = intersect(genelistx, genelisty, 'stable');
  i = startsWith(genelist, 'MT-');
  genelist(i) = [];
  [~, i1] = ismember(genelist, genelistx);
  [~, i2] = ismember(genelist, genelisty);
  X = X(i1, :);
  Y = Y(i2, :);
  cellidx = [ones(size(X, 2), 1); 2*ones(size(Y, 2), 1)];

**Cluster cells using SIMLR**

.. code-block:: matlab

  C = sc_cluster_x([X Y], 5, 'type', 'simlr');
  s = sc_tsne([X Y], 2, true, false);

  figure;
  scatter(s(:,1), s(:,2), 20, C, 'filled')

  figure;
  scatter(s(:,1), s(:,2), 20, cellidx, 'filled')

**Comparing SC3, SIMLR, and SoptSC on yan.csv data**

.. code-block:: matlab

  [X, genelist] = sc_readtsvfile('example_data/yan.csv');
  t = readtable('example_data/yan_celltype.txt');
  celltypelist = string(t.cell_type1);
  rng(235);

  s = sc_tsne(X, 2);
  c1 = run.ml_SC3(X, 6);
  c2 = run.ml_SIMLR(X, 6);
  c3 = run.ml_SoptSC(X, 'k', 6);

  % Compare with SC3/R reference results
  load example_data/sc3_results.txt
  c0 = sc3_results;

  Cal_NMI(c0, c1)
  Cal_NMI(c0, c2)
  Cal_NMI(c0, c3)

  figure;
  subplot(2, 2, 1)
  gscatter(s(:,1), s(:,2), celltypelist)
  legend('Location', 'northwest');
  title('Cell type (Biological Truth)')

  subplot(2, 2, 2)
  gscatter(s(:,1), s(:,2), c1)
  title('SC3')

  subplot(2, 2, 3)
  gscatter(s(:,1), s(:,2), c2)
  title('SIMLR')

  subplot(2, 2, 4)
  gscatter(s(:,1), s(:,2), c3)
  title('SoptSC')


Demo 5: Pseudotime Analysis and Gene Network Functions
-------------------------------------------------------

Trajectory analysis and single-cell gene regulatory network (scGRN) construction.

**Load example data and select genes**

.. code-block:: matlab

  cdgea;
  [X, genelist] = sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
  [X, genelist] = sc_selectg(X, genelist, 5, 3);

**Trajectory analysis using splinefit**

.. code-block:: matlab

  figure;
  t = sc_trajectory(X, "type", "splinefit", "plotit", true);

  % Correlate gene expression with pseudotime
  r = corr(t, X', 'type', 'spearman');
  [~, idxp] = maxk(r, 4);    % top 4 positively correlated genes
  [~, idxn] = mink(r, 3);    % top 3 negatively correlated genes
  selectedg = genelist([idxp idxn]);

  figure;
  gui.i_plot_pseudotimeseries(log1p(X), genelist, t, selectedg)

**Trajectory analysis using TSCAN**

.. code-block:: matlab

  figure;
  t = sc_trajectory(X, "type", "tscan", "plotit", true);

  r = corr(t, X', 'type', 'spearman');
  [~, idxp] = maxk(r, 4);
  [~, idxn] = mink(r, 3);
  selectedg = genelist([idxp idxn]);

  figure;
  gui.i_plot_pseudotimeseries(log1p(X), genelist, t, selectedg)

**Construct scGRN using PCNet**

.. code-block:: matlab

  X50 = X(1:50, :);
  genelist50 = genelist(1:50);
  A = sc_grn(X50, 'pcnet');

  A = A .* (abs(A) > quantile(abs(A(:)), 0.9));
  G = digraph(A, genelist50);
  LWidths = abs(5 * G.Edges.Weight / max(G.Edges.Weight));
  LWidths(LWidths == 0) = 1e-5;
  figure;
  plot(G, 'LineWidth', LWidths);

**Construct scGRN using GENIE3**

.. code-block:: matlab

  X20 = X(1:20, :);
  genelist20 = genelist(1:20);
  A = run.ml_GENIE3(X20);

  A = A .* (abs(A) > quantile(abs(A(:)), 0.9));
  G = digraph(A, genelist20);
  LWidths = abs(5 * G.Edges.Weight / max(G.Edges.Weight));
  LWidths(LWidths == 0) = 1e-5;
  figure;
  plot(G, 'LineWidth', LWidths);


Demo 6: DE Analysis and Marker Gene Identification Functions
-------------------------------------------------------------

Differential expression analysis, marker gene identification, and cell type annotation.

**Load example data**

.. code-block:: matlab

  cdgea;
  load example_data/markergeneident_demo X genelist s_tsne

**Automatically cluster cells and explore cell type**

.. code-block:: matlab

  figure;
  gui.sc_celltypeexplorer_auto(X, genelist, s_tsne, "species", "mouse");

**Group cells into clusters (k=6)**

.. code-block:: matlab

  figure;
  rng(1234)
  cluster_kmedoids = sc_cluster_s(s_tsne, 6, 'type', 'kmedoids', 'plotit', true);

**Identify marker genes for cluster #4**

.. code-block:: matlab

  gmarkers = sc_pickmarkers(X, genelist, cluster_kmedoids, 9);
  gmarkers = gmarkers{4};

**Show expression level of top 9 marker genes**

.. code-block:: matlab

  figure;
  for k = 1:9
      g = gmarkers(k);
      subplot(3, 3, k)
      grp = log1p(X(genelist == g, :));
      scatter3(s_tsne(:,1), s_tsne(:,2), s_tsne(:,3), 10, grp, 'filled');
      title(g)
  end

**Show expression level of a single marker gene**

.. code-block:: matlab

  figure;
  g = gmarkers(5);
  sc_scattermarker(X, genelist, s_tsne, g);
  title(g)

**Determine cell type for each cluster using PanglaoDB marker genes**

.. code-block:: matlab

  Tct = run.ml_alona(X, genelist, cluster_kmedoids, 'species', 'mouse');

**Use the scgeatool interactive tool**

.. code-block:: matlab

  sce = SingleCellExperiment(X, genelist, s_tsne);
  scgeatool(sce);
