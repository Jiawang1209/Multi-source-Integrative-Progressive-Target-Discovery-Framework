# Single-CAS Figure 1 Pipeline

## Purpose

This document defines the current Figure 1 analysis pipeline as a **single-CAS project workflow**.

The basic unit is:

- one `CAS` number
- one `case_dir`
- one dated project directory named `CAS_xxxxx_Date`
- one complete `a -> b -> c -> d -> e` analysis chain
- one final result set

This is the constraint that should be preserved when the workflow is later packaged into a reusable tool or package.

## Core Principle

The pipeline must not mix multiple compounds in one case-level result.

A correct project structure is:

- `CAS_491-71-4_2026-03-21` is one independent case
- `CAS_117-02-2_2026-03-21` is another independent case
- each case has its own source files, intermediate files, figures, and model outputs

If a downstream table or prediction file contains multiple compounds, that file is not a valid case-level final deliverable. It may be acceptable as a temporary development artifact, but not as the final output of a single-CAS workflow.

## Figure 1 Definition

For one CAS molecule, Figure 1 should be constructed as:

1. `a`: workflow schematic of target collection from `ChEMBL / PPB2 / Swiss / SEA`
2. `b`: Venn diagram of the four human target sets
3. `c`: GO + KEGG enrichment of the union target set, followed by keyword-based filtering and visualization
4. `d`: KEGG re-analysis of key target genes, visualized as a `circlize` chord diagram
5. `e`: Chemprop multi-task model built strictly from the `d`-derived key KEGG genes, followed by target ranking and filtering

## Current Script Map

The current workflow is implemented by these scripts in [`scripts`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts):

- [`build_case_sources.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/build_case_sources.py)
- [`fetch_chembl_targets.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/fetch_chembl_targets.py)
- [`fetch_ppb2_targets.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/fetch_ppb2_targets.py)
- [`run_figure1_bcd.R`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/run_figure1_bcd.R)
- [`build_chemprop_targets_from_kegg.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/build_chemprop_targets_from_kegg.py)
- [`prepare_chemprop_data.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/prepare_chemprop_data.py)
- [`train_chemprop_multitask.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/train_chemprop_multitask.py)
- [`evaluate_chemprop_model.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/evaluate_chemprop_model.py)
- [`analyze_chemprop_predictions.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/analyze_chemprop_predictions.py)
- [`plot_figure1e.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/plot_figure1e.py)

## Recommended Case Directory Layout

Each CAS should have its own dated project directory, for example:

- `CAS_491-71-4_2026-03-21`

Recommended structure:

```text
CAS_491-71-4_2026-03-21/
  a_source_collection/
  b_venn/
  c_go_kegg/
  d_kegg_circlize/
  e_chemprop/
```

Each step must write only into its own folder.

### Folder A: `a_source_collection`

- raw `Swiss` result
- raw `SEA` result
- raw `ChEMBL` target result
- raw `PPB2` target result
- `idmapping.tsv`
- source summary files

### Folder B: `b_venn`

- `venn.rds`
- `venn_inputs.json`
- `1.p_venn.pdf`

### Folder C: `c_go_kegg`

- `2.GO_result_2.xlsx`
- `2.KEGG_result_2.xlsx`
- `GO_KEGG_plot.xlsx`
- `2.GO_KEGG.pdf`

### Folder D: `d_kegg_circlize`

- `key_targets_from_kegg.csv`
- `3.KEGG_circos.pdf`

### Folder E: `e_chemprop`

- `*_targets_from_d.csv`
- `chemprop_long.csv`
- `chemprop_multitask.csv`
- `task_summary.csv`
- `inference_template.csv`
- `training/model_0/best.pt`
- `predictions.csv`
- `test_metrics_by_target.csv`
- `ranked_targets.csv`
- filtered `Figure 1e`

## Directory Rule

The directory convention should be treated as mandatory:

- one case directory per CAS
- directory name format: `CAS_xxxxx_Date`
- one subfolder per panel step
- folder order must follow `a -> b -> c -> d -> e`
- files from different CAS runs must never be mixed

## Step A: Target Source Definition

### Objective

Collect compound identity and human target proteins for one real CAS molecule from four online sources:

- `ChEMBL`
- `PPB2`
- `SwissTargetPrediction`
- `SEA`

### Input

- one real `CAS`
- automatic online resolution of compound identity
- automatic retrieval of:
  - compound name
  - molecular formula
  - SMILES
  - InChIKey

### Human-only constraint

All target lists must be limited to `Homo sapiens`.

### Online-first rule

The intended package behavior is:

- input a real CAS
- resolve compound identity online
- submit the resolved structure to the four target sources online
- write all source outputs into `a_source_collection`

Local cached files may exist in the repository for development and debugging, but they must not define the primary package behavior.

### Current online identity source

Current compound identity resolution is based on `PubChem PUG REST`.

The resolver should preferentially derive:

- `Title` or `IUPACName`
- `MolecularFormula`
- `CanonicalSMILES` or `ConnectivitySMILES`
- `InChIKey`

Current implementation:

- `ChEMBL`: keep `Organism == "Homo sapiens"`
- `SEA`: keep rows whose target identifier ends with `_HUMAN`
- `PPB2`: keep standardized human gene symbols after mapping
- `Swiss`: use the returned target names and normalize them to gene-like symbols where available

### Relevant scripts

- [`src/miptd/discovery.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/discovery.py)
- [`fetch_chembl_targets.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/fetch_chembl_targets.py)
- [`fetch_ppb2_targets.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/fetch_ppb2_targets.py)
- [`fetch_swiss_targets.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/fetch_swiss_targets.py)
- [`fetch_sea_targets.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/fetch_sea_targets.py)
- [`build_case_sources.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/build_case_sources.py)

### Output

Per-case standardized source tables and a summary of counts for the four source lists. These outputs belong in folder `a_source_collection`.

## Step B: Venn Diagram

### Objective

Build the four-set Venn diagram from the human target proteins returned by:

- `Swiss`
- `SEA`
- `ChEMBL`
- `PPB2`

### Key rule

The Venn is based on the **per-source target sets**, not on downstream filtered targets.

### Relevant script

- [`run_figure1_bcd.R`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/run_figure1_bcd.R)

### Output

- `venn.rds`
- `venn_inputs.json`
- `1.p_venn.pdf`

### Interpretation

The Venn is used to show overlap across the four databases.
The downstream enrichment is not restricted to overlap only.
Instead, the next step uses the **union** of all targets.

## Step C and D Filtering Rule

The filtering used for panels `c` and `d` is driven by the command-line disease keywords.

The workflow is:

1. run full GO/KEGG enrichment on the union target set
2. convert `--disease-keywords` into a regex-style matching rule
3. match that rule against the enrichment `Description` field
4. preferentially keep matched rows for panel `c`
5. derive panel `d` genes from the selected KEGG rows

Example:

- `--disease-keywords "NAFLD,liver,hepatic"`

This behaves approximately like:

- `NAFLD|liver|hepatic`

Important fallback:

- if no enrichment description matches the user keywords, the workflow falls back to top-ranked terms by significance and count so panels `c` and `d` are still produced

## Step C: GO and KEGG Enrichment

### Objective

Use the union of all targets from the four databases for enrichment analysis.

### Input

Union of:

- `Swiss`
- `SEA`
- `ChEMBL`
- `PPB2`

### Method

Run:

- `GO enrichment`
- `KEGG enrichment`

Then select terms related to liver biology and liver disease using a keyword-based filter.

Current keyword logic includes patterns such as:

- `liver`
- `hepatic`
- `lipid`
- `fatty`
- `cholesterol`
- `bile`
- `oxidative`
- `apopt`
- `inflamm`
- `insulin`
- `AMPK`
- `PPAR`
- `TNF`
- `NF-kappa`
- `FoxO`
- `MAPK`
- `PI3K`
- `AKT`
- `ABC`

### Relevant script

- [`run_figure1_bcd.R`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/run_figure1_bcd.R)

### Output

- `2.GO_result_2.xlsx`
- `2.KEGG_result_2.xlsx`
- `GO_KEGG_plot.xlsx`
- `2.GO_KEGG.pdf`

### Interpretation

This step reduces the union target space into biologically relevant pathways and functions with a liver-focused interpretation.

## Step D: KEGG Re-analysis and Circlize Chord Diagram

### Objective

Extract key target genes from the selected liver-relevant KEGG pathways and visualize their pathway-gene relationships using `circlize`.

### Input

Selected KEGG pathways from Step C.

### Method

For each selected KEGG pathway:

- extract gene members
- merge all genes from the selected pathways
- build a pathway-gene relationship table
- plot it with a chord diagram

### Relevant script

- [`run_figure1_bcd.R`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/run_figure1_bcd.R)

### Output

- `key_targets_from_kegg.csv`
- `3.KEGG_circos.pdf`

### Interpretation

This file is the correct upstream input for panel `e`.

Important rule:

- panel `e` should start from `key_targets_from_kegg.csv`
- panel `e` should not be built directly from raw Venn overlap or from an earlier shortcut target list

## Step E: Chemprop Multi-task Modeling

### Objective

Build a Chemprop multi-task model from the key KEGG genes derived from panel `d`, then score the single CAS query molecule against those trainable targets.

### Required rule

Panel `e` must be built **strictly from panel `d`**.

The correct chain is:

- `key_targets_from_kegg.csv`
- `build_chemprop_targets_from_kegg.py`
- `prepare_chemprop_data.py`
- `train_chemprop_multitask.py`
- `evaluate_chemprop_model.py`
- `analyze_chemprop_predictions.py`
- `plot_figure1e.py`

### Step E1: Convert d-derived genes into Chemprop candidate targets

Relevant script:

- [`build_chemprop_targets_from_kegg.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/build_chemprop_targets_from_kegg.py)

Output:

- `*_targets_from_d.csv`

This table represents the target candidates that come from the `d` panel KEGG genes.
It should be written into folder `e_chemprop`.

### Step E2: Prepare Chemprop training data

Relevant script:

- [`prepare_chemprop_data.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/prepare_chemprop_data.py)

Key logic:

- use local target catalog and `idmapping`
- map genes to `target_chembl_id`
- prefer `SINGLE PROTEIN`
- fetch ChEMBL bioactivity data
- keep only `Homo sapiens`
- keep standard activity types such as `IC50`, `Ki`, `Kd`, `EC50`
- aggregate measurements to one `pChEMBL` value per molecule-target pair
- keep only targets with enough training samples

Outputs:

- `chemprop_long.csv`
- `chemprop_multitask.csv`
- `task_summary.csv`
- `inference_template.csv`

These outputs should be written into folder `e_chemprop`.

### Step E3: Train Chemprop model

Relevant script:

- [`train_chemprop_multitask.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/train_chemprop_multitask.py)

Current implementation:

- task type: regression
- split: `SCAFFOLD_BALANCED`
- save train/validation/test SMILES splits
- save best checkpoint and final predictions

Outputs:

- `training/model_0/best.pt`
- `predictions.csv`
- `run_metadata.json`

These outputs should be written into folder `e_chemprop`.

### Step E4: Evaluate per-target model quality

Relevant script:

- [`evaluate_chemprop_model.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/evaluate_chemprop_model.py)

This script re-predicts the saved test split and exports per-target metrics.

Outputs:

- `test_predictions.csv`
- `test_metrics_by_target.csv`
- `test_metrics_summary.json`

These outputs should be written into folder `e_chemprop`.

Current required test metrics:

- `RMSE`
- `MAE`
- `R2`
- `n_test_molecules`

### Step E5: Rank targets and apply model-quality filtering

Relevant script:

- [`analyze_chemprop_predictions.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/analyze_chemprop_predictions.py)

Current exported columns include:

- `predicted_pchembl`
- `n_training_molecules`
- `test_rmse`
- `test_mae`
- `test_r2`
- `n_test_molecules`
- `filter_basis`
- `filter_status`

Current filtering rule:

- `n_training_molecules >= 100`
- `n_test_molecules >= 30`
- `test_r2 >= 0.30`

This rule is used to separate:

- raw predicted targets
- high-confidence targets

### Step E6: Plot Figure 1e

Relevant script:

- [`plot_figure1e.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/plot_figure1e.py)

Two figure styles are currently possible:

- unfiltered `Figure 1e`
- filtered `Figure 1e`, showing only `filter_status = keep`

Recommended final output for publication:

- filtered `Figure 1e`

This output should be written into folder `e_chemprop`.

## Current Development Lessons

### 1. Single-CAS isolation is mandatory

The most important workflow requirement is:

- one CAS
- one project
- one result space

The package should prevent accidental mixing of multiple compounds in one inference or ranking table.

### 2. Panel e must start from panel d

This is now a fixed workflow rule.

The correct upstream source for Chemprop is:

- `key_targets_from_kegg.csv`

It should not be replaced by:

- pure overlap targets from panel `b`
- an ad hoc consensus target shortcut

### 3. High prediction score alone is not enough

A target with a high predicted score but poor task-level model quality should not be treated as a final high-confidence target.

Examples of rejection conditions:

- `R2 < 0.30`
- `n_test_molecules < 30`

### 4. Human-only filtering must be enforced twice

Human-only restriction should exist at both levels:

- target catalog level
- activity record level

This is already implemented in the current Chemprop data preparation script.

## Suggested Package Architecture

When this workflow is packaged, the functions should be separated roughly as:

1. `resolve_case(cas)`
2. `fetch_targets(case)`
3. `standardize_targets(case)`
4. `build_venn(case)`
5. `run_enrichment(case)`
6. `extract_kegg_key_genes(case)`
7. `prepare_chemprop(case)`
8. `train_chemprop(case)`
9. `evaluate_chemprop(case)`
10. `rank_targets(case)`
11. `plot_figure1e(case)`
12. `assemble_figure1(case)`

Each function should accept a single case object or case directory and return case-scoped outputs only.

## Suggested Final Deliverables Per CAS

For one CAS case, the final package-level outputs should be:

- source summary
- `1.p_venn.pdf`
- `2.GO_KEGG.pdf`
- `3.KEGG_circos.pdf`
- `key_targets_from_kegg.csv`
- `chemprop_multitask.csv`
- `test_metrics_by_target.csv`
- `ranked_targets.csv`
- filtered `Figure 1e`
- one assembled `Figure 1`

These deliverables should remain inside the same `CAS_xxxxx_Date` project directory and be separated by the `a-e` panel-step folders.

## Final Rule

For future development, the workflow should be treated as:

- not a multi-compound batch report
- not a shared target ranking table across compounds
- but a **single-CAS analysis package**

That package may support repeated runs across many CAS numbers, but each run must remain isolated.
