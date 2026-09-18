.. _cli:

Command-Line Interface
======================

``scgea`` runs the toolbox from a shell, without opening the MATLAB desktop.
It is meant for batch processing, remote jobs and pipeline steps: each
subcommand loads a file, does one thing, and writes a file.

.. contents:: On this page
   :local:
   :depth: 2

Setup
-----

MATLAB must be on ``PATH`` (or set the ``MATLAB_BIN`` environment variable to
the executable). The wrappers live in the toolbox root and work out the
toolbox location from their own path, so no other configuration is needed.

.. code-block:: bash

   # Linux / macOS / WSL
   ./scgea.sh help

.. code-block:: bash

   :: Windows CMD / PowerShell
   scgea help

Both wrappers start MATLAB with ``-batch`` and hand off to ``cli.main``. You
can call that directly instead, which is useful inside an existing MATLAB
session or a job script:

.. code-block:: matlab

   cli.main('deg', '--input', 'data.h5ad', '--group1', '1', '--group2', '2')

Subcommands
-----------

==============  ===================================================
Subcommand      What it does
==============  ===================================================
``run``         full pipeline: filter, normalize, embed, cluster
``filter``      QC filtering of cells and genes
``norm``        normalization
``embed``       dimensionality reduction
``cluster``     cell clustering
``deg``         differential expression
``dvg``         differential variability
``grn``         gene regulatory network construction
``trajectory``  pseudotime estimation
``markers``     marker gene identification
``celltypes``   cell type annotation
``tfactivity``  transcription factor activity
``score``       cell scoring
``help``        print usage
==============  ===================================================

Input and output formats
------------------------

Input: ``.h5ad``, ``.mat``, ``.csv`` / ``.tsv`` / ``.txt``, ``.mtx``,
``.loom``, ``.h5`` / ``.hdf5``, ``.rds``, or a GEO accession written as
``GEO:GSMxxxxxx``.

Output: ``.mat`` (default), ``.h5ad`` (requires Python), ``.csv``, ``.txt``.
Several subcommands choose what to write from the output extension -- ``grn``
writes the full adjacency matrix to ``.mat`` and an edge list to ``.csv``, and
``tfactivity`` writes a score matrix to ``.mat`` and per-cluster averages to
``.csv``.

Preprocessing
-------------

**run** -- the whole pipeline in one call:

.. code-block:: bash

   scgea run --input data.csv --norm deseq --embed umap --cluster --k 8 --output result.mat

Options: ``--norm libsize|deseq``, ``--log1p``, ``--hvg``,
``--embed umap|tsne|phate``, ``--ndim <int>`` (default 3), ``--cluster``,
``--k <int>`` (default 6), ``--cluster-type <method>``, ``--min-cells <int>``
(default 3), ``--min-genes <int>`` (default 200).

**filter**, **norm**, **embed**, **cluster** -- the same steps individually:

.. code-block:: bash

   scgea filter --input raw.mtx --min-cells 3 --min-genes 200 --max-mt-ratio 0.25 --output qc.mat
   scgea norm   --input qc.mat --type deseq --log1p --output norm.mat
   scgea embed  --input norm.mat --method umap --ndim 2 --output emb.mat
   scgea cluster --input emb.mat --k 8 --type kmeans --output clustered.mat

``cluster`` works on the embedding by default; ``--use-expression`` clusters on
the expression matrix instead, which is what the ``sc3``, ``simlr``, ``soptsc``
and ``sinnlrr`` methods need.

Analysis
--------

**deg** -- differential expression between two labelled groups:

.. code-block:: bash

   scgea deg --input result.mat --group1 1 --group2 2 --output deg.csv

Options: ``--method mwu|ttest`` (default ``mwu``),
``--label-by cluster|celltype``. Add ``--batch --output-dir <dir>`` to run every
pairwise comparison. Output columns: ``gene, log2FC, p_val, p_val_adj, avg_1,
avg_2``.

**dvg** -- differential variability, either between two files or between two
groups inside one file:

.. code-block:: bash

   scgea dvg --input1 group1.mat --input2 group2.mat --output dvg.csv
   scgea dvg --input result.mat --group1 1 --group2 2 --output dvg.csv

Options: ``--method splinefit|analytic|brennecke``, ``--direction mean|deviation``.

**grn** -- gene regulatory network:

.. code-block:: bash

   scgea grn --input result.mat --method pcrnet --top-k 500 --output grn.csv

Options: ``--method pcrnet|genie3|pearson|mi|xicor|distcorr``,
``--genes <genelist.txt>`` to restrict the network to named genes,
``--top-k <int>`` to write only the strongest edges, ``--norm``, ``--log1p``.

**trajectory** -- pseudotime:

.. code-block:: bash

   scgea trajectory --input result.mat --method splinefit --output pseudotime.csv

Annotation
----------

**markers** -- marker genes per cluster:

.. code-block:: bash

   scgea markers --input result.mat --top-n 10 --method fast --output markers.csv

``--method fast`` and ``lasso`` return a genuine ranking, best marker first.
``scgenefit`` is different on both counts: it picks one global set by linear
programming and then assigns each gene to the group expressing it most, so a
group can end up with more markers than ``--top-n``, fewer, or none, and the
order within a group is selection order rather than strength.

**celltypes** -- cell type annotation against PanglaoDB, or against your own
marker list:

.. code-block:: bash

   scgea celltypes --input result.mat --species mouse --output celltypes.csv
   scgea celltypes --input result.mat --markers my_markers.csv --output celltypes.csv

The custom marker file is a CSV with ``cell_type`` and ``gene`` columns. Output
is one row per cluster: ``cluster_id, cell_type, score``, where a cluster whose
markers matched nothing comes back as ``Unknown`` with a score of 0.

**tfactivity** -- transcription factor activity:

.. code-block:: bash

   scgea tfactivity --input result.mat --species human --output tf.csv

**score** -- per-cell scores:

.. code-block:: bash

   scgea score --input result.mat --mode cellcycle --output scores.csv
   scgea score --input result.mat --mode signature --signature-genes genes.txt --output scores.csv

Modes are ``cellcycle``, ``stemness`` and ``signature``. The ``cellcycle``
output has ``cell_id, phase, S_score, G2M_score`` -- there is no G1 score,
because a cell is called G1 when neither of the other two is positive.
``--species`` is accepted but has no effect on any of the three modes.

Adding a subcommand
-------------------

The dispatcher in ``+cli/main.m`` holds a registry table. To add a subcommand,
create ``+cli/cmd_<name>.m`` with the signature ``cmd_<name>(varargin)`` and
append one row to that table. No other file needs to change.
