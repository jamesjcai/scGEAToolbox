# scDock — MATLAB Reimplementation

MATLAB reimplementation of [scDock](https://github.com/Andrewneteye4343/scDock) (Bioinformatics, btag103).
Reference paper: `scDock_Bioinfo_btag103.pdf`.

Depends on [scGEAToolbox_dev](https://github.com/jamesjcai/scGEAToolbox) — must be on the MATLAB path.

---

## Pipeline Overview

```
Raw count matrix (X, g)
        │
        ▼
┌─────────────────────┐
│  Module 1           │  sc_dock_rna
│  scRNA-seq analysis │  QC → Norm → HVG → PCA (elbow) → Harmony → UMAP → Cluster → Annotate
└────────┬────────────┘
         │  cell-type labels + normalized expression
         ▼
┌─────────────────────┐
│  Module 2           │  sc_dock_ccc
│  Cell-cell comms    │  LR scoring → permutation test → ranked interactions
└────────┬────────────┘
         │  top receptor genes (auto) or user-supplied PDB IDs
         ▼
┌─────────────────────┐
│  Module 3           │  sc_dock_vina
│  Molecular docking  │  Protein prep → Compound prep → AutoDock Vina → ranked affinities
└─────────────────────┘
```

---

## Quick Start

```matlab
% Add scGEAToolbox_dev to path
addpath(genpath('C:\path\to\scGEAToolbox_dev'));

% Load data
[X, g] = sc_readmtxfile('path/to/10x/');

% Run full pipeline
result = sc_dock(X, g, ...
    'compounds', {'aspirin', 'metformin', 'dexamethasone'}, ...
    'species',   'human', ...
    'outdir',    './dock_output');

% Inspect results
disp(result.ccc.T_interactions(1:10, :))   % top LR interactions
disp(result.vina.T_results(1:10, :))       % top drug candidates
```

---

## sc_dock — Top-Level Orchestrator

### Inputs

| Argument | Type | Description |
|---|---|---|
| `X` | `genes × cells` numeric matrix (sparse OK) | Raw UMI count matrix |
| `g` | string array / cellstr, length = nGenes | Gene symbols |

### Parameters — Module 1 (scRNA-seq)

| Parameter | Default | Description |
|---|---|---|
| `'batch_id'` | `[]` | Per-cell batch label vector; `[]` skips batch correction |
| `'species'` | `'human'` | `'human'` or `'mouse'` |
| `'libsz_cutoff'` | `500` | Minimum UMIs per cell |
| `'mt_ratio'` | `0.15` | Maximum mitochondrial read fraction |
| `'min_cells'` | `10` | Minimum cells a gene must be detected in |
| `'n_hvg'` | `2000` | Number of highly variable genes for PCA |
| `'n_pcs'` | `0` | PCs to retain; `0` = auto via geometric elbow |
| `'n_clusters'` | `10` | Number of clusters |
| `'clust_method'` | `'kmeans'` | `'kmeans'`, `'kmedoids'`, `'spectclust'`, `'snndpc'`, `'mbkmeans'` |
| `'batch_method'` | `'harmony'` | `'harmony'` or `'none'` |

### Parameters — Module 2 (cell-cell communication)

| Parameter | Default | Description |
|---|---|---|
| `'lr_db'` | `[]` | MATLAB table with columns `ligand`, `receptor` (and optionally `pathway`). Default: auto-loads from `scGEAToolbox_dev/assets/Ligand_Receptor/` |
| `'condition'` | `[]` | Per-cell condition label; triggers multi-group mode when provided |
| `'cond_labels'` | `{}` | `{cond1, cond2}` — required for multi-group |
| `'n_perm'` | `100` | Permutation test iterations |
| `'pval_cutoff'` | `0.05` | FDR-adjusted p-value threshold |
| `'ccc_min_cells'` | `10` | Minimum cells per cell type to include in CCC |
| `'sender'` | `[]` | Restrict to these sender cell types (`[]` = all) |
| `'receiver'` | `[]` | Restrict to these receiver cell types (`[]` = all) |

### Parameters — Module 3 (docking)

| Parameter | Default | Description |
|---|---|---|
| `'compounds'` | `[]` | PubChem CIDs, compound names, or struct array (see below). `[]` skips docking |
| `'proteins'` | `[]` | PDB IDs or struct array (see below). `[]` = auto-infer from top CCC receptors |
| `'outdir'` | `'./dock_output'` | Directory for all docking intermediate and output files |
| `'vina_exe'` | `'vina'` | Path to AutoDock Vina executable |
| `'obabel_exe'` | `'obabel'` | Path to OpenBabel executable |
| `'exhaustiveness'` | `8` | Search exhaustiveness (1–32; higher = slower but more thorough) |
| `'n_modes'` | `9` | Number of binding poses to generate per run |
| `'global_docking'` | `true` | Auto-define docking box from full protein bounding box + 10 Å padding |

### Output

`result` struct:

| Field | Type | Description |
|---|---|---|
| `result.rna` | struct | Output of `sc_dock_rna` (see Module 1 below) |
| `result.ccc` | struct | Output of `sc_dock_ccc` (see Module 2 below); `[]` if skipped |
| `result.vina` | struct | Output of `sc_dock_vina` (see Module 3 below); `[]` if no compounds |

---

## Module 1 — sc_dock_rna

### Inputs

| Argument | Description |
|---|---|
| `X` | `genes × cells` raw count matrix |
| `g` | Gene list |

Same optional parameters as the Module 1 section of `sc_dock` above.

### Output fields

| Field | Size | Description |
|---|---|---|
| `X_norm` | `nHVG × nCells` | Normalized log1p expression (HVGs only) |
| `g_hvg` | `nHVG × 1` string | Gene names for rows of `X_norm` |
| `pca_scores` | `nCells × nPCs` | PCA cell embeddings |
| `pca_var` | `nPCs × 1` | Variance explained (%) per PC |
| `n_pcs` | scalar | Number of PCs selected by geometric elbow |
| `s_umap` | `nCells × 3` | UMAP 3-D coordinates |
| `c_cluster_id` | `nCells × 1` int | Cluster assignment (1-based) |
| `c_cell_type` | `nCells × 1` string | Cell-type annotation label |
| `cell_kept` | `nCellsOrig × 1` logical | Cells surviving QC (index into original `X`) |

---

## Module 2 — sc_dock_ccc

### Inputs

| Argument | Description |
|---|---|
| `X_norm` | `genes × cells` normalized log1p matrix (typically `result.rna.X_norm`) |
| `g` | Gene list matching rows of `X_norm` |
| `c_cell_type` | Per-cell cell-type label string array |

Same optional parameters as the Module 2 section of `sc_dock` above.

### Output fields

| Field | Description |
|---|---|
| `T_interactions` | Table of significant LR pairs, sorted by communication probability (single-group) or \|ΔP\| (multi-group). Columns: `ligand`, `receptor`, `sender`, `receiver`, `prob`, `pval`, `adj_pval` [+ `prob_cond1`, `prob_cond2`, `delta_prob` in multi-group] [+ `pathway` if in lr_db] |
| `T_pathways` | Pathway-level aggregation (only when `lr_db` has a `pathway` column). Columns: `pathway`, `sender`, `receiver`, `prob_sum` |
| `cell_types` | String array of cell-type labels used (those with ≥ `min_cells` cells) |
| `mode` | `'single'` or `'multi'` |

### LR database format

The `'lr_db'` table must have at minimum:

| Column | Type | Description |
|---|---|---|
| `ligand` | string | HGNC gene symbol of ligand |
| `receptor` | string | HGNC gene symbol of receptor |
| `pathway` | string (optional) | Signaling pathway name (enables `T_pathways` output) |

Default search order for auto-loaded database (from `scGEAToolbox_dev/assets/Ligand_Receptor/`):
1. `CellChatDB_human_LR.txt` — exported via `export_CellChatDB.R` (richest; includes pathway)
2. `Ligand_Receptor_more.txt` — minimal ligand/receptor columns
3. `Ligand_Receptor.txt` — HPMR/HPRD curated pairs
4. `connectomeDB2020.txt` — ConnectomeDB 2020

---

## Module 3 — sc_dock_vina

### Inputs

**`proteins`** — specify target proteins in one of two ways:
- `string array` of PDB IDs, e.g. `["1ABC", "2XYZ"]` — structures are fetched from RCSB automatically
- `struct array` with fields `.id` (string) and `.pdbqt` (path to pre-prepared PDBQT file), or `.pdb` (path to PDB file, converted automatically)

**`compounds`** — specify small molecules in one of two ways:
- `string array` of PubChem CIDs or compound names, e.g. `["2244", "aspirin"]` — SDFs fetched from PubChem automatically
- `struct array` with fields `.id` (string) and `.pdbqt` (prepared PDBQT), or `.sdf` (SDF file, converted automatically)

Same optional parameters as the Module 3 section of `sc_dock` above.

### Output fields

| Field | Description |
|---|---|
| `T_results` | Table sorted by `affinity_kcal_mol` ascending (most negative = tightest binding). Columns: `protein_id`, `compound_id`, `affinity_kcal_mol`, `pose_file`, `log_file` |
| `outdir` | Output directory containing all PDBQT, pose, and log files |
| `n_docked` | Number of successfully completed docking runs |

### Output files (in `outdir`)

| File | Description |
|---|---|
| `<protein>__<compound>_out.pdbqt` | Predicted binding poses (all modes) |
| `<protein>__<compound>_log.txt` | AutoDock Vina log with per-mode affinities |
| `<protein>.pdbqt` | Prepared receptor structure |
| `<compound>.pdbqt` | Prepared ligand structure |

---

## External Requirements (Module 3 only)

| Tool | Install |
|---|---|
| [AutoDock Vina](https://vina.scripps.edu/) | Download binary; add to PATH or pass via `'vina_exe'` |
| [OpenBabel](https://openbabel.org/) | `conda install -c conda-forge openbabel` |

---

## Supported Input Data Formats

Use scGEAToolbox_dev readers to load data before calling `sc_dock`:

| Format | Function |
|---|---|
| 10x Genomics Cell Ranger (directory) | `sc_read10xdir` |
| 10x HDF5 (`.h5`) | `sc_read10xh5file` |
| HDF5 / h5ad (AnnData) | `sc_readh5adfile` |
| Plain text (genes × cells `.txt`) | `sc_readtsvfile` |
| MEX (`.mtx` + barcodes + features) | `sc_readmtxfile` |
