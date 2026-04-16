# sc_dock_vina — Module 3: AutoDock Vina Molecular Docking

MATLAB refactoring of the original R `Vina_Docking` function from
[scDock](https://github.com/Andrewneteye4343/scDock) (Bioinformatics btag103).

---

## Design rationale (vs original R function)

The original R `Vina_Docking` took 13 parameters, many of which were unnecessary:

| Original R parameter | Disposition in MATLAB |
|---|---|
| `Run_CellChat_output_path` | Removed — caller passes proteins directly (or uses `sc_dock` pipeline) |
| `Vina_Docking_output_path` | Defaults to `'./dock_output'` |
| `Vina_Docking_ligand_ref_file` | Removed — bundled as `receptor_reference.csv` co-located with `sc_dock.m` |
| `Vina_Docking_receptor_ref_file` | Removed — same; UniProt REST API used as fallback |
| `Vina_Docking_cas_txt_file` | Merged into `compounds` argument (pass a string array of CAS numbers) |
| `Vina_Docking_use_fda` | Removed — pass `"fda"` as `compounds` shorthand instead |
| `Vina_Docking_fda_txt` | Removed — bundled `fda.txt` is loaded automatically when `compounds = "fda"` |
| `Vina_Docking_docking_ligand_dir` | Removed — pass pre-built structs with `.pdbqt` field instead |
| `Vina_Docking_docking_receptor_dir` | Removed — same |
| `Vina_Docking_vina_exhaustiveness` | `'exhaustiveness'` name-value, default `8` |
| `Vina_Docking_vina_num_modes` | `'n_modes'` name-value, default `9` |
| `Vina_Docking_vina_seed` | Removed — non-essential; Vina uses its own default |
| `Vina_Docking_vina_cpu` | Removed — Vina auto-detects CPU count |

Result: **2 required inputs** instead of 13, with sensible defaults for everything else.

The mouse→human gene conversion (`homologene`) is also eliminated: inputs are
protein/compound identifiers, not gene names. Gene-to-PDB mapping happens one
level up in `sc_dock.m` via `i_extract_top_proteins`.

---

## Function signature

```matlab
result = sc_dock_vina(proteins, compounds)
result = sc_dock_vina(proteins, compounds, Name, Value, ...)
```

### Required inputs

| Argument | Type | Description |
|---|---|---|
| `proteins` | string array of PDB IDs **or** struct array with `.id`/`.pdbqt` fields | Target proteins to dock against. If a PDB ID string is given, the `.pdb` file is downloaded from RCSB and converted to PDBQT automatically. |
| `compounds` | string array of PubChem CIDs / CAS numbers **or** `"fda"` **or** struct array with `.id`/`.sdf` fields | Small-molecule ligands. Pass `"fda"` to use the bundled FDA-approved drug list (`fda.txt`). SDF files are downloaded from PubChem and converted to PDBQT automatically. |

### Optional name-value pairs

| Name | Default | Description |
|---|---|---|
| `'outdir'` | `'./dock_output'` | Root directory for all output files. Created if absent. |
| `'vina_exe'` | `'vina'` | Path to AutoDock Vina executable. Must be on PATH or specified here. |
| `'obabel_exe'` | `'obabel'` | Path to OpenBabel executable. |
| `'exhaustiveness'` | `8` | Vina search exhaustiveness (1–32). Higher = more thorough, slower. |
| `'n_modes'` | `9` | Maximum binding modes to generate per run. |
| `'box_size'` | `[30 30 30]` | Docking box dimensions in Ångströms. Used only when `global_docking = false`. |
| `'global_docking'` | `true` | Auto-compute bounding box from protein ATOM records (+ 10 Å padding). Recommended for blind docking. |
| `'ph'` | `7.4` | pH for hydrogen assignment by OpenBabel during PDBQT conversion. |

### Output

```matlab
result.T_results   % table: protein_id, compound_id, affinity_kcal_mol,
                   %        pose_file, log_file  (sorted ascending by affinity)
result.outdir      % output directory path (string)
result.n_docked    % number of successfully completed runs (integer)
```

Binding affinity is in kcal/mol; **more negative = tighter binding**.

---

## Usage examples

### Minimal — FDA drugs against two receptors

```matlab
result = sc_dock_vina(["1IYT", "2GS6"], "fda");
disp(result.T_results(1:10, :))   % top 10 hits
```

### User-supplied CAS compound list

```matlab
result = sc_dock_vina(["1IYT", "2GS6"], ...
    ["50-78-2", "15687-27-1", "103-90-2"], ...   % aspirin, ibuprofen, paracetamol
    'exhaustiveness', 12, ...
    'outdir', './my_docking');
```

### Pre-built PDBQT files (skip download/conversion)

```matlab
proteins = struct('id', {'VEGFR2', 'EGFR'}, ...
                  'pdbqt', {'vegfr2_prepared.pdbqt', 'egfr_prepared.pdbqt'});
result = sc_dock_vina(proteins, "fda");
```

### Full pipeline (recommended) — from raw counts to docking hits

```matlab
% X: genes×cells raw count matrix; g: gene list
result = sc_dock(X, g, 'compounds', "fda", 'outdir', './run1');
% result.vina.T_results contains the ranked docking hits
```

`sc_dock` (the top-level pipeline) calls `sc_dock_vina` automatically after
running scRNA-seq analysis (Module 1) and cell-cell communication inference
(Module 2). Protein targets are inferred from the top-ranked CCC interactions
via `receptor_reference.csv` (bundled), with UniProt REST API as fallback.

---

## File outputs

All files are written under `outdir`:

```
dock_output/
  <protein_id>.pdb          # downloaded receptor structure
  <protein_id>.pdbqt        # preprocessed receptor
  <compound_id>.sdf         # downloaded compound structure
  <compound_id>.pdbqt       # preprocessed compound
  <protein>__<compound>_out.pdbqt   # docked pose
  <protein>__<compound>_log.txt     # full Vina scoring log
```

The result table (`result.T_results`) contains the path to each pose and log
file for downstream inspection or re-scoring.

---

## Dependencies

| Tool | Role | Install |
|---|---|---|
| AutoDock Vina | Docking engine | [vina.scripps.edu](https://vina.scripps.edu) |
| OpenBabel | PDB/SDF → PDBQT conversion | `conda install -c conda-forge openbabel` |
| MATLAB `websave` / `webread` | PDB and PubChem download | built-in (MATLAB R2014b+) |

MATLAB toolboxes required: none beyond base MATLAB.
scGEAToolbox_dev does not need to be on the path for `sc_dock_vina` alone,
but is required for `sc_dock_rna` and `sc_dock_ccc` (Modules 1–2).

---

## Comparison with original R `Vina_Docking`

| Aspect | R (`Vina_Docking.R`) | MATLAB (`sc_dock_vina.m`) |
|---|---|---|
| Required inputs | 13 (2 truly required) | 2 |
| Gene→PDB mapping | Inside function via ref CSV | In `sc_dock.m`; UniProt API fallback |
| Mouse→human conversion | `homologene` R package | Not needed (caller provides protein IDs) |
| Structure download | Python subprocess calls | Native `websave` (MATLAB built-in) |
| PDBQT conversion | `prepare_receptor.py` (AutoDockTools) | `obabel` system call |
| Docking box | Pre-computed `_grid.txt` required | Auto-computed from protein ATOM records |
| CAS → SDF → PDBQT | `download_cas_pubchem.py` | `i_fetch_pubchem` + `obabel` (native) |
| FDA drug shorthand | `use_fda = TRUE` boolean + separate txt path | `compounds = "fda"` string |
| Output format | Per-receptor CSV + raw PDBQT files | Unified MATLAB table + PDBQT files |
