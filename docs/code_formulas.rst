Code Formulas
=============

Example codes for common tasks.

Import 10x Genomics files
-------------------------
In the 10x Genomics folder, there are three files, namely, matrix.mtx, features.tsv (or genes.tsv) and barcodes.tsv. Here is how to import them:

.. code-block:: matlab

  mtxf = 'GSM3535276_AXLN1_matrix.mtx';
  genf = 'GSM3535276_AXLN1_genes.tsv';
  bcdf = 'GSM3535276_AXLN1_barcodes.tsv';
  [X, genelist, barcodelist] = sc_readmtxfile(mtxf, genf, bcdf, 2);

If the barcodes.tsv is not available, then use the following

.. code-block:: matlab

  mtxf = 'GSM3535276_AXLN1_matrix.mtx';
  genf = 'GSM3535276_AXLN1_genes.tsv';
  [X,g] = sc_readmtxfile(mtxf, genf, [], 2);


Process expression matrix, `X` and gene list, `g`
-------------------------------------------------
Here is an example of raw data processing.

.. code-block:: matlab
  
  [X, g, b] = sc_readmtxfile('matrix.mtx', 'features.tsv', 'barcodes.tsv', 2);
  [X, g] = sc_qcfilter(X, g);
  [X, g] = sc_selectg(X, g, 0.05);
  [s] = sc_tsne(X);
  sce=SingleCellExperiment(X,g,s);
  scgeatool(sce)

t-SNE embedding of cells using highly variable genes (HVGs)
-----------------------------------------------------------

.. code-block:: matlab

  [~, Xsorted, gsorted] = sc_hvg(X, g);     % T, sorted matrix, sorted gene list
  [s] = sc_tsne(Xsorted(1:2000, :));
  sce = SingleCellExperiment(X, g, s);
  scgeatool(sce)

``sc_splinefit`` returns the same three outputs and is the method the GUI and
``embedcells`` prefer; swap it in where you want the spline-fit ranking instead
of the squared-CV one.


An example pipeline for raw data processing
-------------------------------------------

.. code-block:: matlab

  [X,g] = sc_readmtxfile('matrix.mtx', 'features.tsv');
  [X,g] = sc_qcfilter(X, g);                % basic QC
  [X,g] = sc_selectg(X, g, 0.05);           % genes expressed in at least 5% of cells
  [~,Xsorted] = sc_hvg(X, g);               % identify highly variable genes (HVGs)
  [s] = sc_tsne(Xsorted(1:2000, :));        % tSNE on the top 2000 HVGs
  sce = SingleCellExperiment(X, g, s);      % make SCE object
  sce = sce.estimatepotency("mouse");       % differentiation potency ("human" / "mouse")
  sce = sce.estimatecellcycle();            % cell cycle phase
  id = sc_cluster_s(s, 10);                 % k-means on the tSNE coordinates
  sce.c_cluster_id = id;                    % assign cluster ids
  scgeatool(sce)                            % visualize cells

Clustering through ``sc_cluster_s`` and assigning ``c_cluster_id`` yourself,
as above, is one way to cluster on an embedding computed outside the object.
``sce.clustercells`` is not a drop-in for those two lines: its embedding-based
methods check that ``embedcells`` produced the embedding, because the
constructor fills ``s`` with random coordinates when you do not supply any and
an ``isempty(sce.s)`` test would pass on those. Pass the coordinates as the
fifth argument if you want the method to use them::

  sce = sce.clustercells(10, 'kmeans', true, s);

An example pipeline for processing 10x data folder
--------------------------------------------------
Assuming the .m file containing the following code is in the folder ./filtered_feature_bc_matrix. In this folder, three files: matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz, are present.

.. code-block:: matlab

  [X, genelist, celllist] = sc_read10xdir(pwd);
  sce = SingleCellExperiment(X, genelist);
  sce.c_cell_id = celllist;
  sce = sce.qcfilter;
  sce = sce.estimatecellcycle;
  sce = sce.estimatepotency("mouse");
  sce = sce.embedcells('tsne3d', true);
  save clean_data sce -v7.3
  scgeatool(sce)

Merge two data sets (WT and KO)
-------------------------------

.. code-block:: matlab

  load WT/clean_data.mat sce
  sce_wt = sce;
  load KO/clean_data.mat sce
  sce_ko = sce;
  sce = sc_mergesces({sce_wt, sce_ko}, 'union');    % use parameter 'union' or 'intersect' to merge genes
  sce.c = sce.c_batch_id;
  scgeatool(sce)                                    % blue - WT and red - KO  
  
You may want to re-compute tSNE coordinates after merging.

Import other file formats
-------------------------

.. code-block:: matlab

  [X, g, b] = sc_read10xdir('path/to/filtered_feature_bc_matrix');  % 10x folder
  [X, g, b] = sc_read10xh5file('sample.h5');                        % 10x HDF5
  [X, g, b] = sc_readh5adfile('sample.h5ad');                       % AnnData
  [X, g, b] = sc_readloomfile('sample.loom');                       % Loom
  [X, g, b] = sc_readparsebio('path/to/parsebio_output');           % ParseBio

Two readers return a ready-made ``SingleCellExperiment`` rather than a matrix:

.. code-block:: matlab

  sce = sc_readrdsfile('sample.rds');      % Seurat / RDS, needs R
  sce = sc_readgeoaccess('GSE123456');     % downloads the accession first

A ``.h5mu`` file can carry several assays; list them before reading:

.. code-block:: matlab

  [names, sizes] = sc_h5mumodalities('sample.h5mu');
  [X, g, b, batchid, celltype, adt, adtlist] = sc_readh5mufile('sample.h5mu');

Export an SCE for other tools
-----------------------------

.. code-block:: matlab

  sc_sce2h5ad(sce, 'out.h5ad');    % AnnData, for scanpy (needs Python)
  sc_sce2rds(sce, 'out.rds');      % Seurat (needs R)
  sc_sce2hdf5(sce, 'out.h5');      % generic HDF5
  save('out.mat', 'sce', '-v7.3'); % native, and the format the CLI reads

Variance-stabilizing transformation
-----------------------------------

``sc_norm`` scales; ``sc_transform`` stabilizes variance. The SCTransform v2
algorithm is available both through R and as a native MATLAB port that needs
no R installation:

.. code-block:: matlab

  X = sc_transform(X, 'type', 'PearsonResiduals');
  X = sc_transform(X, 'type', 'SCTransformMATLAB');   % native port
  X = sc_transform(X, 'type', 'SCTransform');         % via R/Seurat

Doublet detection
-----------------

.. code-block:: matlab

  [isdoublet, score] = sc_scrublet(X);
  X = X(:, ~isdoublet);

Subsample a large data set
--------------------------

Geometric sketching keeps rare populations that uniform subsampling loses:

.. code-block:: matlab

  idx = sc_geosketch(X, 5000);     % indices of 5000 representative cells
  Xs = X(:, idx);

Score cells against a gene signature
------------------------------------

.. code-block:: matlab

  signature = ["CD3D" "CD3E" "CD2"];
  score = sc_cellscore(X, g, signature);        % AddModuleScore (default)
  score = sc_cellscore(X, g, signature, [], 1); % 1 = UCell, 3 = AUCell
  sce.setCellAttribute('tcell_score', score);

Build and inspect a gene regulatory network
-------------------------------------------

``sc_grn`` assumes its input is already normalized and transformed -- none of
its branches does that for you:

.. code-block:: matlab

  Xn = log1p(sc_norm(X));
  idx = 1:200;                       % keep the network small enough to read
  A = sc_grn(Xn(idx, :), 'pcrnet');  % also: genie3, pearson, mi, xicor, distcorr
  sc_grnview(A, g(idx));

Run a step from the shell
-------------------------

Every recipe above has a command-line equivalent that needs no MATLAB desktop:

.. code-block:: bash

  scgea run --input data.csv --norm deseq --embed umap --cluster --k 8 --output result.mat
  scgea deg --input result.mat --group1 1 --group2 2 --output deg.csv

See :doc:`cli` for the full set.
