# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Goal

Reimplement [scDock](https://github.com/Andrewneteye4343/scDock) in MATLAB, leveraging functions from **scGEAToolbox_dev** (`C:\Users\jingc\Documents\GitHub\scGEAToolbox_dev`) as much as possible. Reference paper: `scDock_Bioinfo_btag103.pdf` in this repository.

## scDock Architecture (Three Modules)

### Module 1: scRNA-seq Analysis
- Quality control, normalization, scaling, dimensionality reduction, clustering, integration
- **Geometric elbow method** for PC selection (max perpendicular distance from line connecting first/last PCs on variance-explained plot)
- Batch correction: CCA, RPCA, or Harmony
- Cell-type annotation via marker genes

### Module 2: Cell-Cell Communication Inference
- Ligand-receptor expression pattern analysis (analogous to CellChat v2 in R)
- Permutation-based communication probability scoring (100 bootstraps, P < 0.05, ≥10 cells/type)
- Single-group: rank by interaction probability
- Multi-group: rank by absolute difference in communication probability between conditions

### Module 3: Molecular Docking
- Interface to AutoDock Vina for binding affinity estimation
- Protein preprocessing: remove non-standard residues/heteroatoms/water, add hydrogens, assign Gasteiger charges, convert to PDBQT
- Compound library options: FDA-approved drugs, PubChem CAS lookup, user-supplied structures
- Output ranked drug candidates by binding affinity

## scGEAToolbox_dev Integration

Toolbox location: `C:\Users\jingc\Documents\GitHub\scGEAToolbox_dev`

### Key functions to reuse for Module 1 (scRNA-seq)
- `sc_qcfilter`, `sc_qcmetrics` — quality control
- `sc_norm`, `sc_transform` — normalization
- `sc_hvg` — highly variable genes
- `sc_umap`, `sc_tsne` — dimensionality reduction/visualization
- `sc_cluster_x`, `sc_cluster_s` — clustering
- `sc_celltypeanno`, `sc_csubtypeanno` — cell type annotation
- `sc_deg` — differential expression
- `+run/ml_Harmony`, `+run/ml_Harmony2` — batch correction
- `+pkg/i_grpmean.m` — group mean computation
- `@SingleCellExperiment` — primary data structure

### Key functions to reuse for Module 2 (Cell-Cell Communication)
- `+pkg/e_cellscores.m` — cell scoring
- `+pkg/e_determinecelltype.m` — cell type utilities
- `sc_grn`, `sc_pcnet` — network construction patterns (adapt for LR networks)
- `+run/web_STRING.m` — external network data retrieval pattern

### Key functions for supporting infrastructure
- `+gui/i_waitbar.m` — progress bars
- `+pkg/e_uint2sparse.m` — sparse matrix handling
- `+pkg/e_makeembedstruct.m`, `+pkg/e_makecluststruct.m` — struct patterns
- `+pkg/e_fdr_bh.m` — FDR correction for permutation tests
- `+pkg/e_grptest.m` — group statistical testing

## MATLAB File Naming Conventions (follow scGEAToolbox_dev patterns)

- Root-level pipeline functions: `sc_dock_*.m` (e.g., `sc_dock_rna.m`, `sc_dock_ccc.m`, `sc_dock_vina.m`)
- Internal helper functions: `+pkg/` with `i_` prefix (internal) or `e_` prefix (exported utility)
- External tool wrappers: `+run/` with `ml_` or `py_` prefix
- GUI callbacks: `+gui/` with `callback_` prefix

## Implementation Notes

- AutoDock Vina integration should use `system()` calls (pattern: see `+run/` wrappers in scGEAToolbox_dev)
- Protein/compound file handling uses PDBQT format; use OpenBabel via `system()` for format conversion
- `@SingleCellExperiment` is the primary data container — store docking results as custom attributes via `list_cell_attributes`
- Ligand-receptor database: bundle CellChatDB human/mouse tables as `.mat` files
