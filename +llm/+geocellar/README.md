# GEOcellar — Native MATLAB Multi-Agent Pipeline

GEOcellar is a four-agent LLM pipeline that automates hypothesis generation and single-cell RNA-seq analysis from GEO datasets. This is a native MATLAB reimplementation of the Python/AutoGen-based `GEOcellar_autogen` pipeline, integrated into scGEAToolbox as the `+llm/+geocellar/` package.

---

## Architecture

```
Phase 1a  Retrieval Agent     →  queries ChromaDB Cloud for relevant GEO studies
Phase 1b  Hypothesis Agent    →  generates 2–4 testable hypotheses (3-pass critique loop)
Phase 2a  Data Analysis Agent →  downloads GEO samples, runs llm.run_* analysis functions
Phase 2b  Interpret Agent     →  synthesises results into a markdown report
```

The orchestrator (`geocellar_main`) runs the phases in sequence, checkpointing after each successful analysis to prevent re-running work across sessions.

---

## Files

| File | Role |
|------|------|
| `geocellar_main.m` | Pipeline orchestrator — entry point for programmatic use |
| `geocellar_app.m` | uifigure-based GUI (launched from scGEAToolbox main menu) |
| `geocellar_config.m` | Reads API keys and paths from environment / MATLAB prefs |
| `retrieval_agent.m` | Phase 1a: calls `chroma_query`, formats study context as JSON |
| `hypothesis_agent.m` | Phase 1b: 3-pass draft/critique/revise loop via `openAIChat` |
| `data_analysis_agent.m` | Phase 2a: downloads samples, calls `llm.run_*` functions directly |
| `interpret_agent.m` | Phase 2b: single LLM call → markdown report |
| `chroma_query.m` | ChromaDB Cloud query via Python subprocess (see below) |
| `download_sample.m` | Downloads `{gse}/{gsm}/cleandata.mat` from GCS via `websave` |
| `load_checkpoint.m` | Reads `checkpoint.json` → struct array of past analyses |
| `save_checkpoint.m` | Writes struct array → `checkpoint.json` |
| `i_chroma_query.py` | Python helper script called by `chroma_query.m` |
| `hypothesis_system.txt` | System prompt for the hypothesis agent (reused from Python pipeline) |
| `interpret_system.txt` | System prompt for the interpret agent (reused from Python pipeline) |

---

## Usage

### Programmatic
```matlab
% Basic run
llm.geocellar.geocellar_main("lactation and breast cancer")

% With seed hypothesis
llm.geocellar.geocellar_main("T cell exhaustion", ...
    SeedHypothesis="PDCD1 drives CD8 exhaustion")

% Custom working directory
llm.geocellar.geocellar_main("cancer immunotherapy", DataDir="C:/mydata")
```

### GUI
Launch from the scGEAToolbox main menu: **Tools → GEOcellar Agent**

The GUI callback (`+gui/callback_LaunchGEOcellar.m`) validates the LLM Add-On and API key, prompts for a working directory via `gui.gui_setprgmwkdir`, then opens `geocellar_app`.

---

## Configuration (`geocellar_config`)

All settings are read from environment variables (loaded via `loadenv` from the file registered in the `scgeatoolbox` → `llapikeyenvfile` preference).

| Variable | Default | Purpose |
|----------|---------|---------|
| `OPENAI_API_KEY` | — | TAMU AI API key (required) |
| `OPENAI_API_BASE` | `https://chat-ai.tamu.ai/api` | LLM base URL |
| `OPENAI_MODEL` | `protected.Claude Sonnet 4.5` | Model for analysis/interpretation |
| `HYPOTHESIS_MODEL` | same as `OPENAI_MODEL` | Model override for hypothesis agent |
| `CHROMA_API_KEY` | — | ChromaDB Cloud API key (required) |
| `DATA_DIR` | `../GEOcellar_MATLAB/data` | Where GEO sample `.mat` files are saved |
| `OUT_DIR` | `../GEOcellar_MATLAB/output` | Where reports and checkpoints are written (overridden to `<DataDir>/output` when `DataDir` is set) |

---

## ChromaDB Integration (`chroma_query`)

ChromaDB Cloud requires client-side embeddings — the REST API does not embed `query_texts` server-side (returns 405 if attempted). The solution delegates to Python:

1. `chroma_query.m` calls `i_chroma_query.py` via `uv run` inside the `GEOcellar_autogen` project directory.
2. The Python script uses `chromadb.HttpClient` + `SentenceTransformerEmbeddingFunction` (`all-MiniLM-L6-v2`) to embed the query locally and send `query_embeddings` to ChromaDB Cloud.
3. Results are printed as JSON and parsed back in MATLAB.

**Prerequisites:** `GEOcellar_autogen/` must be a sibling of `scGEAToolbox_dev/` and its `uv` environment must have `chromadb` and `sentence-transformers` installed (they are listed in `pyproject.toml`).

---

## Hypothesis Agent — 3-Pass Critique Loop

```
Pass 1  Draft    openAIChat generates initial JSON hypothesis array
Pass 2  Critique LLM reviews against: GSM ID validity, analysis type fit,
                 diversity of analysis types, scientific motivation
Pass 3  Revise   LLM incorporates critique, returns final JSON array
```

After parsing, hypotheses are filtered: only those whose GSM IDs are present in the retrieved ChromaDB text are kept.

Each hypothesis struct has fields: `id`, `statement`, `gse`, `sample1`, `sample2`, `sample1_description`, `sample2_description`, `analysis`, `kogene`, `celltype1`, `celltype2`, `implications`, `analysis_plan`.

---

## Data Analysis Agent — Supported Analysis Types

| `analysis` value | Function called | Required fields |
|-----------------|-----------------|-----------------|
| `DEG` | `llm.run_de_analysis` + `llm.run_enrichr` | `sample1`, `sample2` |
| `DV` | `llm.run_dv_analysis` | `sample1`, `sample2` |
| `scTenifoldNet` | `llm.run_sctenifoldnet` | `sample1`, `sample2` |
| `scTenifoldKnk` | `llm.run_sctenifoldknk` | `sample1`, `kogene` |
| `scTenifoldXct` | `llm.run_sctenifoldxct` | `sample1`, `celltype1`, `celltype2` |

Results from each analysis are JSON-encoded into the checkpoint record.

---

## Outputs (in `OUT_DIR`)

| File | Contents |
|------|----------|
| `retrieved_data.json` | Raw ChromaDB query results |
| `hypotheses.json` | Generated hypothesis structs |
| `checkpoint.json` | Persistent record of all completed analyses |
| `report.md` | Final markdown interpretation report |
| `H{n}_{analysis}/` | Per-hypothesis output (DEG tables, Enrichr results, etc.) |
| `{timestamp}/` | Archived copy of JSON and report from each run |

---

## GUI (`geocellar_app`)

Follows the `LassoAnalysisApp` pattern used across scGEAToolbox:

- Signature: `geocellar_app(wrkdir, parentfig)` — `parentfig` last, default `[]`
- Shows `gui.myWaitbar(parentfig)` while building the figure
- Creates figure with `Visible = "off"`, calls `gui.i_movegui2parent`, then shows it
- Left panel: topic input, seed hypothesis input, scrolling log textarea, Run / Open Output buttons
- Right panel: report textarea (populated after pipeline completes)

---

## Key Design Decisions

**No MCP bridge needed.** The Python pipeline calls MATLAB via MCP. In native MATLAB the Data Analysis Agent calls `llm.run_*` functions directly — simpler and faster.

**One function per file.** MATLAB only exports the first function in a file and requires the filename to match. `checkpoint_io.m` was split into `load_checkpoint.m` and `save_checkpoint.m` for this reason.

**`import llm.geocellar.*` in orchestrators.** Functions that call sibling package functions (`geocellar_main`, `retrieval_agent`, `data_analysis_agent`) must declare `import llm.geocellar.*` after the `arguments` block (or after `nargin` checks). Functions that only call `llm.*` or `openAIChat` do not need the import.

**Python only for ChromaDB.** Everything else (LLM chat, hypothesis generation, interpretation, analysis, download) is pure MATLAB. The Python subprocess is limited to `chroma_query` because ChromaDB Cloud mandates client-side embeddings and the sentence-transformers model is already available in the `GEOcellar_autogen` uv environment.
