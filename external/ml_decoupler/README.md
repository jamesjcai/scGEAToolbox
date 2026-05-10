# ml_decoupler — Enrichment Analysis for scGEAToolbox

MATLAB port of the Python [decoupler](https://github.com/scverse/decoupler-py) enrichment framework (Badia-i-Mompel et al., *Bioinformatics* 2022). Provides five complementary methods for computing per-cell biological activity scores from single-cell RNA-seq data using predefined gene networks.

All functions accept a `SingleCellExperiment` object and a **network table**, and return a `n_cells × n_terms` score matrix that can be passed to `sc_rankby_group` to identify the top active terms per cluster.

---

## Setup

Add the folder to the MATLAB path once per session (or add to `startup.m`):

```matlab
addpath('path/to/scGEAToolbox_dev/external/ml_decoupler');
```

---

## Quick Start

```matlab
% 1. Load a gene network
net = sc_net_panglaodb("human");      % cell-type markers (178 types)
% or
net = sc_net_progeny("human");        % 14 signaling pathways

% 2. Score every cell against every term
[scores, pvals, terms] = sc_ulm(sce, net, MinGenes=5);

% 3. Find top-ranked terms per cluster
T = sc_rankby_group(scores, terms, sce.c_cluster_id, TopN=3);
disp(T)
```

---

## Enrichment Methods

All scoring functions share the same interface:

```matlab
[scores, pvals, terms] = sc_METHOD(sce, net)
[scores, pvals, terms] = sc_METHOD(sce, net, MinGenes=5, ...)
```

| Argument | Type | Description |
|----------|------|-------------|
| `sce` | `SingleCellExperiment` | Log-normalised counts in `sce.X` (genes × cells) |
| `net` | `table` | Network with columns `source`, `target`, `weight` |
| `scores` | `double n_cells × n_terms` | Enrichment score per cell per term |
| `pvals` | `double n_cells × n_terms` | Two-sided p-value (NaN where not applicable) |
| `terms` | `string 1 × n_terms` | Term names in score column order |

> **Note on weights.** Methods that use weights (ULM, MLM, ZSCORE) require non-uniform weights to produce meaningful scores. Databases like PROGENy and PanglaoDB already supply non-uniform weights. For plain gene lists where all weights = 1, use AUCell or ORA instead.

---

### `sc_ulm` — Univariate Linear Model

For each cell and each gene set, fits `expression_g = β·weight_g + intercept` via OLS. Returns the t-statistic of β. Each gene set is fit independently, making this fast and interpretable.

```matlab
[scores, pvals, terms] = sc_ulm(sce, net)
[scores, pvals, terms] = sc_ulm(sce, net, MinGenes=5, Verbose=true)
```

| Option | Default | Description |
|--------|---------|-------------|
| `MinGenes` | `3` | Skip terms with fewer overlapping genes |
| `Verbose` | `true` | Print summary line |

**Best for:** Weighted networks (PROGENy pathways, DoRothEA TF regulons). Scores are t-statistics — larger magnitude = stronger enrichment; sign indicates over- (positive) or under-representation.

**Reference:** Badia-i-Mompel et al., *Bioinformatics* 2022. https://doi.org/10.1093/bioinformatics/btac832

---

### `sc_mlm` — Multivariate Linear Model

Fits all gene sets **simultaneously** in a single multivariate OLS per cell:

```
expression_g = β₀ + β₁·w_{g,1} + β₂·w_{g,2} + … + βₜ·w_{g,T} + ε
```

Returns the t-statistic of each βₖ. Unlike ULM, shared target genes across gene sets are handled jointly, reducing cross-talk artifacts in dense networks.

```matlab
[scores, pvals, terms] = sc_mlm(sce, net)
[scores, pvals, terms] = sc_mlm(sce, net, MinGenes=5, ReturnCoef=false)
```

| Option | Default | Description |
|--------|---------|-------------|
| `MinGenes` | `3` | Skip terms with fewer overlapping genes |
| `ReturnCoef` | `false` | Return OLS coefficients instead of t-statistics |

**Best for:** Dense TF networks (DoRothEA) where many regulators share targets. Requires `n_genes >> n_terms`.

**Reference:** Badia-i-Mompel et al., *Bioinformatics* 2022.

---

### `sc_aucell` — Area Under the recovery Curve

Genes are ranked by expression per cell. For each gene set the AUC of the recovery curve measures how many set members fall within the top-ranked genes. Scores are in [0, 1]; no weights are used.

```matlab
[scores, terms] = sc_aucell(sce, net)
[scores, terms] = sc_aucell(sce, net, TopGenes=500, MinGenes=3)
```

| Option | Default | Description |
|--------|---------|-------------|
| `TopGenes` | `ceil(0.05 × n_genes)` | Size of the top-gene pool (5% rule) |
| `MinGenes` | `3` | Skip terms with fewer overlapping genes |

> AUCell does not return p-values. The second output `terms` is returned directly (no `pvals`).

**Best for:** Unweighted gene sets, robustness to outlier cells, low-quality or sparse data.

**Reference:** Aibar et al., *Nature Methods* 2017. https://doi.org/10.1038/nmeth.4463

---

### `sc_ora` — Over-Representation Analysis

For each cell, defines the top `TopGenes` expressed genes as "significant". Tests each gene set for over-representation in this list using a one-sided Fisher's exact test (hypergeometric). Returns the log odds ratio (Haldane-Anscombe corrected) as the score.

```matlab
[scores, pvals, terms] = sc_ora(sce, net)
[scores, pvals, terms] = sc_ora(sce, net, TopGenes=300, Background=15000, MinGenes=3)
```

| Option | Default | Description |
|--------|---------|-------------|
| `TopGenes` | `ceil(0.05 × n_genes)` | Top expressed genes treated as "significant" |
| `Background` | `n_genes` | Background gene pool size |
| `MinGenes` | `3` | Skip terms with fewer overlapping genes |

**Best for:** Discrete gene lists, fast exploratory enrichment, situations where a binary active/inactive threshold is appropriate.

---

### `sc_zscore` — Weighted Z-score

Computes the weighted mean expression of target genes, normalised by the within-cell standard deviation and the square root of the number of targets. Two flavors:

- **RoKAI** (default): subtracts the per-cell global mean before scaling
- **KSEA**: uses raw weighted mean without global-mean correction

```
ES = (µ_set − µ_global) × √m / σ_global     [RoKAI]
ES =  µ_set              × √m / σ_global     [KSEA]
```

```matlab
[scores, pvals, terms] = sc_zscore(sce, net)
[scores, pvals, terms] = sc_zscore(sce, net, Flavor="KSEA", MinGenes=3)
```

| Option | Default | Description |
|--------|---------|-------------|
| `Flavor` | `"RoKAI"` | `"RoKAI"` or `"KSEA"` |
| `MinGenes` | `3` | Skip terms with fewer overlapping genes |

**Best for:** Fast first-pass scoring, kinase substrate enrichment (KSEA flavor), phosphoproteomics.

**References:** Hernandez-Armenta et al., *Bioinformatics* 2017 (KSEA). Yılmaz et al., *Cell Systems* 2021 (RoKAI).

---

## Method Comparison

| Method | Uses weights | Output range | p-values | Best use case |
|--------|:---:|---------|:---:|--------------|
| `sc_ulm` | Required | (−∞, +∞) t-stat | Yes | Weighted networks; fast, per-term |
| `sc_mlm` | Required | (−∞, +∞) t-stat | Yes | Dense networks with shared genes |
| `sc_aucell` | Ignored | [0, 1] AUC | No | Any gene set; robust to outliers |
| `sc_ora` | Ignored | (−∞, +∞) log OR | Yes | Discrete top-N list enrichment |
| `sc_zscore` | Required | (−∞, +∞) z-score | Yes | Fast; kinase/phospho data |

---

## Network Loaders

### `sc_net_panglaodb` — PanglaoDB cell-type markers

Returns a network table of canonical cell-type marker genes for 178 human or mouse cell types. Weights are log₂-fold-enrichment scores from the PanglaoDB database.

```matlab
net = sc_net_panglaodb()                       % human, all markers
net = sc_net_panglaodb("mouse")
net = sc_net_panglaodb("human", MinGeneScore=1.5)
```

| Option | Default | Description |
|--------|---------|-------------|
| `MinGeneScore` | `0` | Keep markers with weight ≥ this value |

**Reference:** Franzén et al., *Database* 2019. https://doi.org/10.1093/database/baz046

---

### `sc_net_progeny` — PROGENy signaling pathways

Returns a weighted network for 14 curated signaling pathways (Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFκB, p53, PI3K, TGFβ, TNFα, Trail, VEGF, WNT). Weights are from the PROGENy perturbation-response regression model.

Fetches from OmniPath on first use (internet required); result is cached in `tempdir()`. Tries Python decoupler first if installed.

```matlab
net = sc_net_progeny()                         % human, top 500 genes/pathway
net = sc_net_progeny("mouse")
net = sc_net_progeny("human", TopGenes=100)
```

| Option | Default | Description |
|--------|---------|-------------|
| `TopGenes` | `500` | Top N genes per pathway ranked by \|weight\| |

**Reference:** Schubert et al., *Nature Communications* 2018. https://doi.org/10.1038/s41467-017-02391-6

---

## Post-scoring Utilities

### `sc_rankby_group` — Top terms per cluster

Ranks gene sets within each cell cluster by mean enrichment score. Equivalent to `dc.tl.rankby_group()`.

```matlab
T = sc_rankby_group(scores, terms, labels)
T = sc_rankby_group(scores, terms, sce.c_cluster_id, TopN=5)
```

**Inputs:**

| Argument | Description |
|----------|-------------|
| `scores` | `n_cells × n_terms` score matrix from any scoring function |
| `terms` | `1 × n_terms` string array from the same scoring function |
| `labels` | `n_cells × 1` cluster/cell-type labels |

**Output table columns:**

| Column | Description |
|--------|-------------|
| `Group` | Cluster / cell-type label |
| `Rank` | Rank within group (1 = top) |
| `Term` | Gene set name |
| `MeanScore` | Mean enrichment score across cells in this group |
| `PctPositive` | Fraction of cells in group with score > 0 |
| `MeanVsRest` | Mean score in group minus mean score in all other cells |

---

## Full Workflow Example

```matlab
addpath('external/ml_decoupler');

% --- Cell type annotation with PanglaoDB ---
net  = sc_net_panglaodb("human");
[scores, ~, terms] = sc_ulm(sce, net, MinGenes=5);
T    = sc_rankby_group(scores, terms, sce.c_cluster_id, TopN=3);
disp(T(T.Group == "1", :))   % top cell types for cluster 1

% --- Pathway activity with PROGENy ---
net_pg = sc_net_progeny("human", TopGenes=100);
[sc_pg, pv_pg, tr_pg] = sc_ulm(sce, net_pg, MinGenes=5);
T_pg   = sc_rankby_group(sc_pg, tr_pg, sce.c_cluster_id);
disp(T_pg)

% --- Unweighted gene sets: use AUCell ---
[auc, tr_ac] = sc_aucell(sce, net, TopGenes=500);
T_ac = sc_rankby_group(auc, tr_ac, sce.c_cluster_id, TopN=3);

% --- Discrete enrichment: use ORA ---
[lor, pv_or, tr_or] = sc_ora(sce, net);
T_or = sc_rankby_group(lor, tr_or, sce.c_cluster_id, TopN=3);

% --- Compare methods on the same network ---
[s_ul, ~, t_ul] = sc_ulm(sce,    net_pg, MinGenes=5);
[s_ml, ~, t_ml] = sc_mlm(sce,    net_pg, MinGenes=5);
[s_zs, ~, t_zs] = sc_zscore(sce, net_pg, MinGenes=5);
```

---

## File Reference

| File | Description |
|------|-------------|
| `sc_ulm.m` | Univariate Linear Model scoring |
| `sc_mlm.m` | Multivariate Linear Model scoring |
| `sc_aucell.m` | AUCell rank-based scoring |
| `sc_ora.m` | Over-Representation Analysis |
| `sc_zscore.m` | Weighted z-score scoring |
| `sc_net_panglaodb.m` | Load PanglaoDB marker network |
| `sc_net_progeny.m` | Load PROGENy pathway network |
| `sc_rankby_group.m` | Rank terms by cell group |
| `i_build_adj.m` | Internal: build genes × terms adjacency matrix |

---

## References

1. Badia-i-Mompel P et al. "decoupleR: ensemble of methods that infer biological activities from omics data." *Bioinformatics* 2022. https://doi.org/10.1093/bioinformatics/btac832
2. Aibar S et al. "SCENIC: single-cell regulatory network inference and clustering." *Nature Methods* 2017. https://doi.org/10.1038/nmeth.4463
3. Franzén O et al. "PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data." *Database* 2019. https://doi.org/10.1093/database/baz046
4. Schubert M et al. "Perturbation-response genes reveal signaling footprints in cancer gene expression." *Nature Communications* 2018. https://doi.org/10.1038/s41467-017-02391-6
5. Hernandez-Armenta C et al. "Benchmarking substrate-based kinase activity inference." *Bioinformatics* 2017. https://doi.org/10.1093/bioinformatics/btx082
6. Yılmaz S et al. "Robust inference of kinase activity using functional networks." *Cell Systems* 2021. https://doi.org/10.1016/j.cels.2021.08.012
