.. _getting_started:

Getting Started
===============

.. contents:: On this page
   :local:
   :depth: 2

Launch the GUI
--------------

Run the following in MATLAB to start SCGEATOOL, the App Designer interface:

.. code-block:: matlab

   scgeatool

Called with no argument it opens empty, ready to import data. Called with a
``SingleCellExperiment`` object it opens on that object:

.. code-block:: matlab

   cdgea;
   load example_data/testSce.mat sce
   scgeatool(sce);

``cdgea`` changes the working folder to the toolbox root, which is what makes
the ``example_data/`` paths used throughout this documentation resolve.

Requirements
------------

MATLAB R2026a, with the Statistics and Machine Learning Toolbox. Individual
functions may need more: the Bioinformatics Toolbox for some file readers, the
Parallel Computing Toolbox for the parallel GRN methods. Wrapper functions
under ``run.py_*`` and ``run.r_*`` call out to Python and R respectively and
need those installed separately.

Two ways to work
----------------

The toolbox offers a matrix-level API and an object-level API. Both are
supported; they differ in what they carry between steps.

**Matrix level.** The readers return ``[X, genelist, celllist]`` and the
``sc_*`` functions thread that matrix through. Most root-level ``sc_*``
functions operate on the count matrix, not on an SCE object -- passing an
``sce`` where ``X`` is expected fails the ``mustBeNumeric`` validator.

.. code-block:: matlab

   [X, genelist] = sc_read10xdir('path/to/10x_data');
   [X, genelist] = sc_qcfilter(X, genelist);   % combined QC filter
   X = sc_norm(X);                             % 'libsize' (default), 'deseq', 'shiftedclr'
   T = sc_hvg(X, genelist);                    % table of HVG statistics
   s = sc_umap(X, 2);                          % embedding, cells x ndim

**Object level.** ``SingleCellExperiment`` wraps the matrix, gene list,
embedding and every derived annotation in one object and exposes the same
steps as methods:

.. code-block:: matlab

   sce = SingleCellExperiment(X, genelist);
   sce = sce.qcfilter();
   sce = sce.embedcells('umap');
   sce = sce.clustercells(8, 'kmeans', true);
   scgeatool(sce);

``SingleCellExperiment`` is a ``handle`` class with a ``Copyable`` mixin.
Plain assignment aliases the object rather than copying it, so the methods
above mutate ``sce`` in place and the ``sce =`` on the left is a convention,
not a copy. To keep an independent snapshot, use ``copy``:

.. code-block:: matlab

   sce_before = copy(sce);

What is in an SCE
-----------------

===========================  ===============================================
Property                     Contents
===========================  ===============================================
``X``                        counts, genes x cells (sparse)
``g``                        gene list, string array
``s``                        cell embedding, cells x ndim
``c``                        the currently active grouping
``c_cluster_id``             clustering result, one entry per cell
``c_cell_type_tx``           cell type labels
``c_cell_cycle_tx``          cell cycle phase labels
``c_batch_id``               batch identifier
``c_cell_id``                cell barcodes
``list_cell_attributes``     derived per-cell values (potency, scores, ...)
``struct_cell_embeddings``   every embedding computed, by method
``struct_cell_clusterings``  every clustering computed, by method
``struct_modalities``        additional assays on the same cells (ADT, ATAC)
``metadata``                 provenance, appended to as steps run
===========================  ===============================================

``NumCells`` and ``NumGenes`` are dependent properties, computed from ``X``.

Where to go next
----------------

- :ref:`case_studies` -- worked demos, from filtering through to annotation.
- :doc:`code_formulas` -- short recipes for common tasks.
- :doc:`cli` -- running the toolbox from a shell, without opening MATLAB.
- :doc:`function_reference` -- what is available, by category.
