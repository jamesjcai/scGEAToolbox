# Glycobiology Gene Set Collection

A curated collection of **27 gene sets (448 unique genes)** covering the major
functional modules of the glycome. Scoring cells against these sets yields a
per-cell, per-module readout that collectively defines a cell's
**glycobiological state** — the relative activity of its glycan biosynthesis,
degradation, transport, and recognition programs.

Gene symbols are human (HGNC) and are matched **case-insensitively**, so mouse
orthologues that share a spelling still match. A few genes are listed under both
their current and legacy symbols (e.g. `OGA`/`MGEA5`, `GFUS`/`TSTA3`,
`CRPPA`/`ISPD`, `RXYLT1`/`TMEM5`) so the sets match regardless of annotation.

Curated from KEGG glycan pathways, GO glycosylation terms, the CAZy / GlycoGene
classification, and the Consortium for Functional Glycomics (CFG) glycogene
lists.

## The 27 gene sets

| # | Name | Genes | Description | Reference |
|---|------|:-----:|-------------|-----------|
| 1 | `Glyco_N_glycan_biosynthesis_ER` | 31 | N-glycan lipid-linked oligosaccharide assembly and en-bloc transfer (OST) | KEGG hsa00510 |
| 2 | `Glyco_N_glycan_processing_Golgi` | 19 | N-glycan trimming and branching in the Golgi | KEGG hsa00510 |
| 3 | `Glyco_O_GalNAc_mucin_type` | 26 | Mucin-type O-GalNAc glycan initiation and core extension | KEGG hsa00512 |
| 4 | `Glyco_O_GlcNAc_cycling` | 9 | Nucleocytoplasmic O-GlcNAc addition and removal (nutrient sensing) | GO:0006493 |
| 5 | `Glyco_hexosamine_biosynthesis` | 14 | Hexosamine biosynthetic pathway (UDP-GlcNAc / CMP-sialic acid supply) | GO:0006047 |
| 6 | `Glyco_sialylation` | 25 | Sialyltransferases and CMP-sialic acid biosynthesis | CAZy GT29 / GO:0097503 |
| 7 | `Glyco_sialidases` | 4 | Sialidases / neuraminidases (desialylation) | GO:0004308 |
| 8 | `Glyco_fucosylation` | 18 | Fucosyltransferases and GDP-fucose supply | CAZy GT10/GT23 / GO:0036065 |
| 9 | `Glyco_galactosylation` | 17 | Galactosyltransferases and galactose (Leloir) metabolism | CAZy GT7/GT31 / GO:0006486 |
| 10 | `Glyco_heparan_sulfate_biosynthesis` | 28 | Heparan sulfate chain polymerization and modification | KEGG hsa00534 |
| 11 | `Glyco_chondroitin_dermatan_sulfate` | 21 | Chondroitin and dermatan sulfate biosynthesis | KEGG hsa00532 |
| 12 | `Glyco_keratan_sulfate` | 8 | Keratan sulfate biosynthesis | KEGG hsa00533 |
| 13 | `Glyco_hyaluronan_metabolism` | 10 | Hyaluronan synthesis and turnover | GO:0030212 |
| 14 | `Glyco_glycosphingolipid_biosynthesis` | 15 | Glycosphingolipid (ganglio/globo/lacto) biosynthesis | KEGG hsa00600/00601/00603/00604 |
| 15 | `Glyco_glycan_lysosomal_degradation` | 31 | Lysosomal glycosidases and sulfatases (glycan/GAG/GSL degradation) | KEGG hsa00511/00531 |
| 16 | `Glyco_GPI_anchor_biosynthesis` | 29 | Glycosylphosphatidylinositol (GPI) anchor biosynthesis and remodeling | KEGG hsa00563 |
| 17 | `Glyco_ER_glycoprotein_quality_control` | 18 | Calnexin/calreticulin cycle and glycan-dependent ERAD | GO:0006986 / R-HSA-901042 |
| 18 | `Glyco_nucleotide_sugar_metabolism` | 27 | Biosynthesis and interconversion of nucleotide sugar donors | GO:0009225 |
| 19 | `Glyco_nucleotide_sugar_transporters` | 25 | SLC35 family nucleotide sugar transporters | SLC35 / GO:0015780 |
| 20 | `Glyco_O_mannosylation_dystroglycan` | 20 | O-mannosylation and the alpha-dystroglycan matriglycan pathway | KEGG hsa00515 |
| 21 | `Glyco_Notch_EGF_O_glycosylation` | 14 | O-fucose/O-glucose glycosylation of EGF repeats and Fringe modification | GO:0018345 / R-HSA-1912399 |
| 22 | `Glyco_carbohydrate_sulfotransferases` | 31 | Carbohydrate sulfotransferases and endosulfatases | GO:0008146 |
| 23 | `Glyco_galectins` | 17 | Galectins (S-type beta-galactoside-binding lectins) | GO:0005534 |
| 24 | `Glyco_siglecs` | 15 | Siglecs (sialic acid-binding immunoglobulin-type lectins) | GO:0033691 |
| 25 | `Glyco_selectins_and_ligands` | 15 | Selectins and the enzymes that build their sialyl-Lewis-x ligands | GO:0070492 |
| 26 | `Glyco_C_type_lectins` | 27 | C-type lectins and collectins (glycan recognition) | GO:0005537 |
| 27 | `Glyco_proteoglycan_core_proteins` | 29 | Proteoglycan core proteins (GAG scaffolds) | GO:0005539 |

The sets group into six themes:

- **Glycan biosynthesis** (1–14): N-glycan (ER assembly/OST + Golgi processing),
  O-GalNAc mucin-type, O-GlcNAc cycling, hexosamine pathway, sialylation,
  fucosylation, galactosylation, glycosaminoglycans (heparan/chondroitin/keratan
  sulfate, hyaluronan), and glycosphingolipids.
- **Glycan degradation** (15): lysosomal glycosidases and sulfatases.
- **Support machinery** (16–19): GPI anchoring, ER glycoprotein quality control,
  nucleotide-sugar metabolism, and SLC35 transporters.
- **Specialized O-glycosylation** (20–21): dystroglycan O-mannosylation and
  Notch/EGF O-fucose/O-glucose glycosylation.
- **Glycan modification** (22): carbohydrate sulfotransferases.
- **Glycan recognition & scaffolds** (23–27): galectins, siglecs, selectins,
  C-type lectins, and proteoglycan core proteins.

Gene sets intentionally overlap where pathways share enzymes (e.g. `XYLT1`,
`B4GALT7`, and `B3GAT3` appear in both the heparan and chondroitin/dermatan
sulfate sets, since they build the shared GAG–protein linker tetrasaccharide).

## How the modules relate to glycobiology

"Glyco" here is **glycosylation** — the enzymatic attachment of sugar chains
(glycans) to proteins and lipids, plus the machinery that supplies, builds,
remodels, recognizes, and degrades them. Unlike transcription and translation,
glycan synthesis is **non-templated and combinatorial**: the final structure
emerges from the balance of many competing enzymes in the ER and Golgi. That is
why glycobiology is best summarized as **functional modules** rather than single
genes — a module score (the mean expression of a pathway's enzymes) tracks
pathway *capacity* far better than any one transcript.

The 27 modules tile the complete glycan lifecycle:

1. **Supply the building blocks.** Sugars must first be activated as
   nucleotide sugars (UDP-GlcNAc, GDP-fucose, CMP-sialic acid, UDP-Gal…) and
   imported into the secretory lumen. → `hexosamine_biosynthesis`,
   `nucleotide_sugar_metabolism`, `nucleotide_sugar_transporters` (SLC35).
2. **Attach glycans to proteins.**
   - *N-linked* — an oligosaccharide is preassembled on a dolichol lipid in the
     ER, transferred *en bloc* to asparagine, then trimmed/branched in the
     Golgi. → `N_glycan_biosynthesis_ER`, `N_glycan_processing_Golgi`.
   - *O-linked* — sugars added one at a time to Ser/Thr: mucin-type O-GalNAc on
     secreted/surface mucins, plus specialized forms for α-dystroglycan and for
     Notch/EGF repeats. → `O_GalNAc_mucin_type`, `O_mannosylation_dystroglycan`,
     `Notch_EGF_O_glycosylation`.
   - *O-GlcNAc* — a distinct **intracellular** single-sugar mark that cycles on
     and off nucleocytoplasmic proteins like phosphorylation, a nutrient/stress
     sensor. → `O_GlcNAc_cycling`.
3. **Glycosylate lipids and anchors.** Glycosphingolipids (on ceramide) and GPI
   anchors (which tether proteins to the membrane). →
   `glycosphingolipid_biosynthesis`, `GPI_anchor_biosynthesis`.
4. **Build matrix sugars (glycosaminoglycans / proteoglycans).** Long sulfated
   GAG chains — heparan sulfate, chondroitin/dermatan sulfate, keratan sulfate —
   assembled on proteoglycan core proteins, plus non-sulfated, membrane-made
   hyaluronan. → `heparan_sulfate_biosynthesis`, `chondroitin_dermatan_sulfate`,
   `keratan_sulfate`, `hyaluronan_metabolism`, `proteoglycan_core_proteins`.
5. **Decorate and finish (terminal modifications).** The caps that set
   recognition and half-life — sialylation, fucosylation, galactosylation,
   sulfation — with sialidases providing the reverse. → `sialylation`,
   `fucosylation`, `galactosylation`, `carbohydrate_sulfotransferases`,
   `sialidases`.
6. **Quality control.** The calnexin/calreticulin cycle reads N-glycan tags to
   gate folding and route misfolded glycoproteins to ERAD. →
   `ER_glycoprotein_quality_control`.
7. **Degrade and recycle.** Lysosomal glycosidases and sulfatases dismantle
   glycans, GAGs, and glycolipids (defects cause lysosomal storage diseases). →
   `glycan_lysosomal_degradation`.
8. **Read the glyco-code.** Lectins interpret surface sugars: galectins
   (β-galactosides), siglecs (sialic acids; largely immune-inhibitory),
   selectins with their sialyl-Lewis-x ligand enzymes (leukocyte rolling), and
   C-type lectins/collectins (pattern recognition). → `galectins`, `siglecs`,
   `selectins_and_ligands`, `C_type_lectins`.

Read together, a cell's 27 module scores are a compact fingerprint of its
**glycobiological state** — how much it is transcriptionally investing in donor
supply, protein/lipid/matrix glycosylation, terminal decoration, quality
control, turnover, and glycan recognition.

## Provenance and caveats

The collection is defined in
[`../../+pkg/e_glycogenesets.m`](../../+pkg/e_glycogenesets.m); each module
carries a `Reference` naming its primary source:

- **KEGG PATHWAY** glycan maps — `hsa00510` (N-glycan), `hsa00512` (mucin
  O-glycan), `hsa00515` (mannose-type O-glycan), `hsa00532/00533/00534` (GAG
  biosynthesis), `hsa00563` (GPI), `hsa00600/00601/00603/00604`
  (glycosphingolipid), `hsa00511/00531` (glycan degradation).
- **Gene Ontology** process/function terms — e.g. `GO:0006493` (O-GlcNAc),
  `GO:0097503` (sialylation), `GO:0008146` (sulfotransferase),
  `GO:0005534`/`GO:0033691` (glycan binding).
- **Reactome** — `R-HSA-901042` (ER quality control), `R-HSA-1912399` (Notch
  glycosylation).
- **CAZy** glycosyltransferase families (GT29 sialyl-, GT10/GT23 fucosyl-,
  GT7/GT31 galactosyl-) and the **GlycoGene / CFG** (Consortium for Functional
  Glycomics) glycogene lists.

These are **hand-curated representative panels** — the well-established genes of
each module, chosen for balanced coverage of biosynthesis, degradation,
transport, and recognition — not exhaustive database dumps, and not a formally
benchmarked signature set. Symbols are human HGNC, matched case-insensitively
(so shared-spelling mouse orthologues match), with a few legacy aliases
included. For database-exact membership (e.g. for a publication), cross-check
each list against its cited KEGG/GO/Reactome source.

## Files

| File | Purpose |
|------|---------|
| [`../../+pkg/e_glycogenesets.m`](../../+pkg/e_glycogenesets.m) | **Source of truth.** Returns the collection as `[setmatrx, setnames, setgenes, T]` — the same native form as `pkg.e_getgenesets`. Optionally exports a GMT file. |
| [`../../sc_glycostate.m`](../../sc_glycostate.m) | Top-level scoring function. Returns a modules × cells glyco-state matrix. |
| [`glycobiology.mat`](glycobiology.mat) | Precomputed `setmatrx` / `setnames` / `setgenes` cache (loaded by `pkg.e_getgenesets`). |
| [`glycobiology.gmt`](glycobiology.gmt) | Portable MSigDB-style GMT export for use outside MATLAB. |

## Usage

### Programmatic — score every module per cell

```matlab
[cs, modules] = sc_glycostate(sce.X, sce.g);   % modules × cells score matrix
% cluster or embed cells on cs, or attach one module as a cell attribute:
sce = sce.addcellattr('Sialylation', cs(modules == "Glyco_sialylation", :)');
```

`sc_glycostate(X, genelist, methodid, minGenes)`:
- `methodid` — scoring method forwarded to `sc_cellscore`: `1` = UCell,
  `2` = AddModuleScore (default), `3` = AUCell.
- `minGenes` — a module is scored only if at least this many of its genes are
  present in the data (default `3`); modules below the threshold are dropped.

### Programmatic — get the raw collection

```matlab
[setmatrx, setnames, setgenes, T] = pkg.e_glycogenesets();   % native form + table
[setmatrx, setnames, setgenes]    = pkg.e_getgenesets(4);    % via the standard loader
pkg.e_glycogenesets("my_glyco.gmt");                          % export a GMT copy
```

### GUI

The collection appears as **"Glycobiology Gene Sets"** in the gene-set collection
pickers:

- **Analyze → Compare Cell Score Between Clusters** — pick one or more modules and
  score cells directly (violin / heatmap / stem-scatter output).
- The DP-analysis gene-set collection picker (`e_getgenesets` option `4`).

## Regenerating the cached files

The `.mat` and `.gmt` are snapshots of `pkg.e_glycogenesets()`. After editing the
gene lists in `+pkg/e_glycogenesets.m`, rebuild them:

```matlab
[setmatrx, setnames, setgenes] = pkg.e_glycogenesets();
save(fullfile('assets','GeneSets','glycobiology.mat'), ...
    'setmatrx', 'setnames', 'setgenes');
pkg.e_glycogenesets(fullfile('assets','GeneSets','glycobiology.gmt'));
```
