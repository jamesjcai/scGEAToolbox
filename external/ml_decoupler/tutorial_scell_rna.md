# Single-Cell RNA-seq Enrichment Analysis

> **MATLAB equivalent of the [decoupler single-cell tutorial](https://decoupler.readthedocs.io/en/latest/notebooks/scell/rna_sc.html)**
>
> This tutorial demonstrates how to compute per-cell enrichment scores that
> summarize molecular activity from predefined gene networks, using the same
> PBMC 3 k dataset as the original decoupler vignette.

---

## Contents

1. [Setup](#1-setup)
2. [Load and Prepare Data](#2-load-and-prepare-data)
3. [Cell Type Annotation — PanglaoDB](#3-cell-type-annotation--panglaodb)
4. [Transcription Factor Activity — DoRothEA](#4-transcription-factor-activity--dorothea)
5. [Pathway Activity — PROGENy](#5-pathway-activity--progeny)
6. [Hallmark Gene Sets — MSigDB](#6-hallmark-gene-sets--msigdb)
7. [Method Comparison](#7-method-comparison)
8. [Tips and Caveats](#8-tips-and-caveats)
9. [References](#9-references)

---

## 1. Setup

```matlab
% Add the ml_decoupler folder to the MATLAB path
addpath('/path/to/scGEAToolbox_dev');
addpath('/path/to/scGEAToolbox_dev/external/ml_decoupler');
```

**Dependencies** (all bundled with scGEAToolbox):

| Resource | Location |
|----------|----------|
| PanglaoDB markers | `external/fun_alona_panglaodb/` |
| DoRothEA TF network | `assets/DoRothEA_TF_Target_DB/` |
| MSigDB Hallmark | `assets/MSigDB/msigdb_H.mat` |

PROGENy is fetched from the OmniPath web service on first use and cached locally.

---

## 2. Load and Prepare Data

### 2.1 Download PBMC 3k

Download the filtered 10x Genomics matrix from the Cell Ranger website:

```matlab
url = 'https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz';
websave('pbmc3k.tar.gz', url);
untar('pbmc3k.tar.gz');   % extracts to filtered_gene_bc_matrices/hg19/
```

Alternatively, load a pre-processed H5AD file (e.g. from Scanpy) using
`sc_readh5adfile`, or any other compatible format supported by scGEAToolbox.

### 2.2 Load and QC

```matlab
% Load raw counts (genes x cells)
sce = sc_read10xdir('filtered_gene_bc_matrices/hg19/');
fprintf('Loaded: %d genes x %d cells\n', sce.NumGenes, sce.NumCells);
% Expected: ~32,738 genes x 2,700 cells

% Filter low-quality cells (< 200 or > 2,500 genes expressed)
[Xf, keptc] = sc_filterc(sce.X, 200);
[Xf]        = sc_filterc(Xf',  2500)';   % upper bound via transpose trick
sce.X = Xf;
sce.g = sce.g;   % gene names unchanged

% Filter genes expressed in fewer than 3 cells
[Xf, gf] = sc_filterg(sce.X, sce.g, 3);
sce = SingleCellExperiment(Xf, gf);
fprintf('After QC: %d genes x %d cells\n', sce.NumGenes, sce.NumCells);
```

### 2.3 Normalize and Embed

```matlab
% Library-size normalization and log1p transform
sce.X = sparse(log1p(sc_norm(sce.X)));

% Highly variable genes (top 2,000) → PCA → UMAP
[T_hvg] = sc_hvg(full(sce.X), sce.g);
hvg_idx = ismember(sce.g, T_hvg.gene(1:2000));

X_hvg = full(sce.X(hvg_idx, :));           % HVG subset for embedding
s_umap = sc_umap(X_hvg, 2);                % 2-D UMAP coordinates
sce.s  = [s_umap, zeros(sce.NumCells,1)];  % store in sce.s

% Leiden / k-NN clustering on UMAP coordinates
c_ids = sc_cluster_s(s_umap, 8);           % k=8 neighbours
sce.c_cluster_id = c_ids;

% Quick UMAP overview
figure;
gscatter(sce.s(:,1), sce.s(:,2), categorical(c_ids));
xlabel('UMAP 1'); ylabel('UMAP 2');
title('PBMC 3k — Leiden clusters');
legend('Location','bestoutside');
```

---

## 3. Cell Type Annotation — PanglaoDB

Enrichment analysis requires a gene network that maps biological concepts
(cell types, TFs, pathways) to gene sets with associated weights. Here we use
**PanglaoDB** [1] canonical marker genes, weighted by their log₂-enrichment
scores, and the **ULM** method [2] to score every cell.

> **Why ULM for cell types?**  
> Unlike a simple mean, ULM detects whether highly-weighted markers are
> co-expressed at proportionally higher levels than low-weighted markers.
> This makes it sensitive to the marker hierarchy and robust to ambient RNA.

### 3.1 Load the network

```matlab
net_ct = sc_net_panglaodb("human");

% Inspect
fprintf('%d cell types, %d marker entries\n', ...
    numel(unique(net_ct.source)), height(net_ct));
disp(head(net_ct, 5))
```

```
178 cell types, 7843 marker entries
         source          target    weight
    ______________    __________   ______
    "Acinar cells"    "CTRB1"      1.970
    "Acinar cells"    "KLK1"       1.804
    "Acinar cells"    "RBPJL"      2.000
    "Acinar cells"    "PTF1A"      1.970
    "Acinar cells"    "CELA3A"     2.000
```

### 3.2 Score cells with ULM

```matlab
[ct_scores, ct_pvals, ct_terms] = sc_ulm(sce, net_ct, MinGenes=5);
% SC_ULM: 22,699 genes x 2,638 cells
% SC_ULM: scored 101/178 terms (77 skipped with <5 genes)
```

### 3.3 Visualize on UMAP

Color the UMAP by the ULM score for two canonical cell types:

```matlab
cell_types_to_plot = ["T cells", "B cells", "NK cells", "Monocytes"];

figure('Position', [100 100 1200 300]);
for k = 1:numel(cell_types_to_plot)
    ct = cell_types_to_plot(k);
    idx = strcmp(ct_terms, ct);
    if ~any(idx), continue; end

    sv = ct_scores(:, idx);
    sv_norm = (sv - min(sv)) / (max(sv) - min(sv) + eps);  % [0,1] for display

    subplot(1, numel(cell_types_to_plot), k);
    scatter(sce.s(:,1), sce.s(:,2), 5, sv_norm, 'filled');
    colormap(gca, flipud(colormap('hot')));
    clim([0 1]);  colorbar;
    title(ct);  axis off;
end
sgtitle('PanglaoDB ULM scores');
```

### 3.4 Top predicted cell type per cluster

```matlab
T_ct = sc_rankby_group(ct_scores, ct_terms, sce.c_cluster_id, TopN=3);
disp(T_ct(T_ct.Rank == 1, :))   % single best prediction per cluster
```

```
 Group  Rank        Term         MeanScore  PctPositive  MeanVsRest
 _____  ____  _______________    _________  ___________  __________
   "0"    1   "T cells"            12.41       0.94         10.83
   "1"    1   "Monocytes"          18.77       0.97         15.02
   "2"    1   "B cells"            20.55       0.99         19.88
   "3"    1   "NK cells"           14.22       0.91          9.64
   ...
```

### 3.5 Matrix plot — top 3 cell types per cluster

```matlab
% Collect top-3 unique terms across all clusters
top_terms = unique(T_ct.Term(T_ct.Rank <= 3), 'stable');
clusters  = string(unique(sce.c_cluster_id));

% Build mean-score matrix: clusters x terms
hm_data = zeros(numel(clusters), numel(top_terms));
for g = 1:numel(clusters)
    mask = string(sce.c_cluster_id) == clusters(g);
    for t = 1:numel(top_terms)
        tidx = strcmp(ct_terms, top_terms(t));
        hm_data(g, t) = mean(ct_scores(mask, tidx), 'omitnan');
    end
end

% Z-score across clusters for display
hm_z = (hm_data - mean(hm_data, 1)) ./ (std(hm_data, 0, 1) + eps);

figure;
heatmap(top_terms, clusters, hm_z, ...
    'Colormap', redbluecmap(256), ...
    'ColorLimits', [-2 2], ...
    'Title', 'Cell type scores per cluster (z-scored)');
```

---

## 4. Transcription Factor Activity — DoRothEA

Transcription factor transcript levels often correlate poorly with their
activity. By instead scoring each cell's gene expression against the TF's
known regulon, we get a more informative protein-activity proxy.

The **DoRothEA** database [3] provides signed TF→target relationships with
five confidence levels (A = highest, E = lowest). We use levels A–C for a
balance of coverage and reliability.

### 4.1 Build the network from the bundled database

```matlab
pw_root = fileparts(fileparts(which('sc_ulm')));  % toolbox root
load(fullfile(pw_root, 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat'), 'T');

% Keep confidence levels A, B, C
conf  = string(T.confidence);
T_abc = T(ismember(conf, ["A","B","C"]), :);

net_tf = table( ...
    string(T_abc.tf),     ...   % source (TF name)
    string(T_abc.target), ...   % target (gene symbol)
    double(T_abc.mor),    ...   % weight: +1 activation, -1 repression
    'VariableNames', {'source','target','weight'});

fprintf('DoRothEA (A-C): %d TFs, %d interactions\n', ...
    numel(unique(net_tf.source)), height(net_tf));
% DoRothEA (A-C): 271 TFs, 39,820 interactions
```

### 4.2 Score cells with ULM

```matlab
[tf_scores, tf_pvals, tf_terms] = sc_ulm(sce, net_tf, MinGenes=5);
```

> **Biological insight**: PAX5 drives B-cell identity but its mRNA is often
> low in scRNA-seq. ULM detects it as highly active in B cells because its
> target genes are co-ordinately expressed — even when PAX5 itself is absent
> from the count matrix.

### 4.3 Top TFs per cluster

```matlab
T_tf = sc_rankby_group(tf_scores, tf_terms, sce.c_cluster_id, TopN=5);
disp(T_tf(T_tf.Rank <= 3, :))
```

### 4.4 Highlight a key TF — PAX5 in B cells

```matlab
pax5_idx = strcmp(tf_terms, "PAX5");
figure('Position', [100 100 900 350]);

% UMAP colored by PAX5 activity
subplot(1,2,1);
sv = tf_scores(:, pax5_idx);
scatter(sce.s(:,1), sce.s(:,2), 6, sv, 'filled');
colormap(gca, redbluecmap(256));
clim([-max(abs(sv)), max(abs(sv))]);
colorbar;  title('PAX5 TF activity (ULM)');  axis off;

% Violin / box plot per cluster
subplot(1,2,2);
boxchart(categorical(sce.c_cluster_id), sv, ...
    'BoxFaceColor', [0.2 0.5 0.8], 'MarkerStyle', '.');
yline(0, '--k');
xlabel('Cluster');  ylabel('ULM t-statistic');
title('PAX5 activity by cluster');
```

### 4.5 TF network visualization

Use the `sc_grnview` function for an interactive network diagram of a TF
and its top-scoring targets in a cluster of interest:

```matlab
% Get mean target-gene expression for B-cell cluster
b_mask  = strcmp(string(T_tf.Term(T_tf.Group == "2" & T_tf.Rank == 1)), "PAX5");
% (Adjust cluster label to match your data's B-cell cluster)

% Subset the DoRothEA PAX5 regulon
pax5_net = net_tf(net_tf.source == "PAX5", :);
pax5_net = pax5_net(ismember(pax5_net.target, sce.g), :);

% sc_grnview(sce.X, sce.g, pax5_net, ...)  % see sc_grnview documentation
```

---

## 5. Pathway Activity — PROGENy

**PROGENy** [4] provides 14 signaling pathways whose gene-weight sets were
derived from large-scale perturbation experiments. The signed weights capture
which genes respond most strongly (and in which direction) when each pathway
is activated, making ULM a natural fit.

### 5.1 Load PROGENy

```matlab
% Fetches from OmniPath on first call; cached in tempdir() afterwards
net_pg = sc_net_progeny("human", TopGenes=100);
fprintf('PROGENy: %d pathways\n', numel(unique(net_pg.source)));
% PROGENy: 14 pathways
```

Pathways covered:

| | | | |
|---|---|---|---|
| Androgen | EGFR | Estrogen | Hypoxia |
| JAK-STAT | MAPK | NFkB | p53 |
| PI3K | TGFb | TNFa | Trail |
| VEGF | WNT | | |

### 5.2 Score cells with ULM

```matlab
[pg_scores, pg_pvals, pg_terms] = sc_ulm(sce, net_pg, MinGenes=5);
```

### 5.3 All 14 pathways — heatmap by cluster

```matlab
clusters = string(unique(sce.c_cluster_id));
hm_pg = zeros(numel(clusters), numel(pg_terms));
for g = 1:numel(clusters)
    mask = string(sce.c_cluster_id) == clusters(g);
    hm_pg(g, :) = mean(pg_scores(mask, :), 1, 'omitnan');
end
% Z-score across clusters
hm_pg_z = (hm_pg - mean(hm_pg,1)) ./ (std(hm_pg,0,1) + eps);

figure;
heatmap(pg_terms, clusters, hm_pg_z, ...
    'Colormap', redbluecmap(256), ...
    'ColorLimits', [-2 2], ...
    'Title', 'PROGENy pathway activity per cluster (z-scored)');
```

Expected patterns consistent with the literature:
- **Trail** — elevated in B cells (apoptosis resistance)
- **TNFα** — highest in FCGR3A+ non-classical monocytes (inflammatory)
- **PI3K** — highest in megakaryocytes (platelet-related)

### 5.4 Single-pathway UMAP and distribution

```matlab
pw_name = "TNFa";
pw_idx  = strcmp(pg_terms, pw_name);

figure('Position', [100 100 900 350]);
subplot(1,2,1);
sv = pg_scores(:, pw_idx);
scatter(sce.s(:,1), sce.s(:,2), 6, sv, 'filled');
colormap(gca, redbluecmap(256));
clim([-max(abs(sv)), max(abs(sv))]);
colorbar;  title(pw_name + " (PROGENy ULM)");  axis off;

subplot(1,2,2);
boxchart(categorical(sce.c_cluster_id), sv, 'MarkerStyle', '.');
yline(0, '--k');
xlabel('Cluster');  ylabel('ULM t-statistic');
title(pw_name + " by cluster");
```

---

## 6. Hallmark Gene Sets — MSigDB

The **MSigDB Hallmark** collection [5] contains 50 curated biological state
gene sets designed to minimise redundancy. Unlike PROGENy, Hallmark sets have
**uniform weights** (each gene counts equally).

> **Method choice**: ULM requires non-uniform weights to produce meaningful
> scores (it models co-variation of expression with weight). For uniform-weight
> gene sets, use **AUCell** (rank-based AUC) or **ZSCORE** instead.

### 6.1 Load Hallmark from bundled MSigDB

```matlab
pw_root = fileparts(fileparts(which('sc_ulm')));
load(fullfile(pw_root, 'assets', 'MSigDB', 'msigdb_H.mat'), ...
    'setmatrx', 'setnames', 'setgenes');

% Convert logical membership matrix (50 sets x 4,384 genes) to network table
[set_idx, gene_idx] = find(setmatrx);
net_hm = table( ...
    setnames(set_idx)',         ...  % source: Hallmark set name
    string(setgenes(gene_idx))', ... % target: gene symbol
    ones(numel(set_idx), 1),    ...  % weight: 1 (uniform)
    'VariableNames', {'source','target','weight'});

fprintf('Hallmark: %d gene sets, %d entries\n', ...
    numel(unique(net_hm.source)), height(net_hm));
% Hallmark: 50 gene sets, 22,478 entries
```

### 6.2 Score cells with AUCell

AUCell ranks all genes by expression per cell and measures how many Hallmark
set members fall within the top-ranked genes, without requiring gene weights.

```matlab
[hm_scores, hm_terms] = sc_aucell(sce, net_hm, TopGenes=500, MinGenes=5);
% SC_AUCELL: scored 50/50 terms (TopGenes=500, skipped 0 with <5 genes)
```

Alternatively, use ZSCORE (produces signed scores and p-values):

```matlab
[hm_scores_z, hm_pvals_z, hm_terms_z] = sc_zscore(sce, net_hm, MinGenes=5);
```

### 6.3 Top Hallmark sets per cluster

```matlab
T_hm = sc_rankby_group(hm_scores, hm_terms, sce.c_cluster_id, TopN=3);
disp(T_hm)
```

Expected top markers (consistent with the literature):

| Cell type | Top Hallmark sets |
|-----------|-------------------|
| CD4 T cells | `MYC_TARGETS_V1`, `MYC_TARGETS_V2`, `WNT_BETA_CATENIN_SIGNALING` |
| B cells | `KRAS_SIGNALING_DN`, `INTERFERON_GAMMA_RESPONSE` |
| Monocytes | `COMPLEMENT`, `ANGIOGENESIS`, `EPITHELIAL_MESENCHYMAL_TRANSITION` |
| NK cells | `INTERFERON_ALPHA_RESPONSE`, `INTERFERON_GAMMA_RESPONSE` |

### 6.4 Matrix heatmap of all 50 Hallmark sets

```matlab
clusters = string(unique(sce.c_cluster_id));
hm_mat = zeros(numel(clusters), numel(hm_terms));
for g = 1:numel(clusters)
    mask = string(sce.c_cluster_id) == clusters(g);
    hm_mat(g, :) = mean(hm_scores(mask, :), 1, 'omitnan');
end
hm_mat_z = (hm_mat - mean(hm_mat,1)) ./ (std(hm_mat,0,1) + eps);

% Shorten term names for display
short_names = regexprep(hm_terms, "^HALLMARK_", "");

figure('Position', [100 100 1400 450]);
heatmap(short_names, clusters, hm_mat_z, ...
    'Colormap', redbluecmap(256), ...
    'ColorLimits', [-2 2], ...
    'FontSize', 7, ...
    'Title', 'MSigDB Hallmark AUCell scores per cluster (z-scored)');
```

---

## 7. Method Comparison

The five available methods differ in their assumptions and sensitivity.
Run them all on the same PROGENy network to compare:

```matlab
[s_ulm, ~,  t_ulm]  = sc_ulm(sce,    net_pg, MinGenes=5);
[s_mlm, ~,  t_mlm]  = sc_mlm(sce,    net_pg, MinGenes=5);
[s_zs,  ~,  t_zs]   = sc_zscore(sce, net_pg, MinGenes=5);
[s_ac,      t_ac]   = sc_aucell(sce,  net_pg, MinGenes=5);
[s_ora, ~,  t_ora]  = sc_ora(sce,    net_pg, MinGenes=5);

methods  = {"ULM", "MLM", "ZSCORE", "AUCell", "ORA"};
score_sets = {s_ulm, s_mlm, s_zs, s_ac, s_ora};
term_sets  = {t_ulm, t_mlm, t_zs, t_ac, t_ora};

% Top pathway per cluster for each method
fprintf('%-10s  %-20s  %-20s\n', 'Method', 'Cluster 0 top', 'Cluster 1 top');
for m = 1:numel(methods)
    T_m = sc_rankby_group(score_sets{m}, term_sets{m}, sce.c_cluster_id, TopN=1);
    c0  = T_m.Term(T_m.Group == "0");
    c1  = T_m.Term(T_m.Group == "1");
    if isempty(c0), c0 = "—"; end
    if isempty(c1), c1 = "—"; end
    fprintf('%-10s  %-20s  %-20s\n', methods{m}, c0, c1);
end
```

### Summary of methods

| Method | Weights needed | Output | p-values | Best for |
|--------|:-:|---------|:-:|---|
| `sc_ulm` | Yes (non-uniform) | t-statistic | Yes | Weighted networks; per-term |
| `sc_mlm` | Yes (non-uniform) | t-statistic | Yes | Dense networks with shared targets |
| `sc_zscore` | Yes | z-score | Yes | Fast first-pass; kinase/phospho data |
| `sc_aucell` | Ignored | AUC [0,1] | No | Any gene set; robust to outliers |
| `sc_ora` | Ignored | log odds ratio | Yes | Discrete top-N gene lists |

For most applications, **ULM** (weighted networks) and **AUCell** (unweighted
gene sets) cover complementary use cases.

---

## 8. Tips and Caveats

**Input normalization**  
All scoring methods expect log-normalized counts in `sce.X` (i.e.
`log1p(sc_norm(raw_counts))`). Feeding raw counts will inflate scores for
high-library-size cells and distort p-values.

**Minimum genes per term**  
The `MinGenes` threshold (default 3) controls the minimum gene overlap
between your dataset and the network. Increase to 5–10 to exclude sparsely
covered terms.

**ULM requires non-uniform weights**  
If all weights in a gene set are identical (e.g. a plain gene list), the OLS
regression is degenerate and returns a score of zero. Use `sc_aucell` or
`sc_zscore` for uniform-weight gene sets.

**DoRothEA confidence levels**  
Level A interactions are supported by four independent evidence types;
level E by literature only. For a focused TF analysis, use levels A–B.
For broader coverage, extend to A–C.

**PROGENy TopGenes**  
`sc_net_progeny` returns up to `TopGenes` genes per pathway ranked by
|weight|. The default (500) is appropriate for most datasets; reduce to
100–200 for small panels or sparse data.

**Scores vs. p-values**  
Scores (t-statistics, z-scores, AUC) reflect enrichment magnitude.
P-values from `sc_ulm`, `sc_mlm`, `sc_zscore`, and `sc_ora` are
cell-level and are not corrected for multiple testing. Use them for
filtering highly significant cells, not for genome-wide inference.

---

## 9. References

1. Franzén O, Gan L-M, Björkegren JLM. "PanglaoDB: a web server for
   exploration of mouse and human single-cell RNA sequencing data."
   *Database* 2019. https://doi.org/10.1093/database/baz046

2. Badia-i-Mompel P et al. "decoupleR: ensemble of methods that infer
   biological activities from omics data within a unified framework."
   *Bioinformatics* 2022. https://doi.org/10.1093/bioinformatics/btac832

3. Garcia-Alonso L et al. "Benchmark and integration of resources for
   the estimation of human transcription factor activities."
   *Genome Research* 2019. https://doi.org/10.1101/gr.240663.118

4. Schubert M et al. "Perturbation-response genes reveal signaling
   footprints in cancer gene expression."
   *Nature Communications* 2018.
   https://doi.org/10.1038/s41467-017-02391-6

5. Liberzon A et al. "The Molecular Signatures Database Hallmark Gene
   Set Collection." *Cell Systems* 2015.
   https://doi.org/10.1016/j.cels.2015.12.004

6. Aibar S et al. "SCENIC: single-cell regulatory network inference
   and clustering." *Nature Methods* 2017.
   https://doi.org/10.1038/nmeth.4463
