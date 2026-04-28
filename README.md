# MIPTD

`MIPTD` stands for `Multi-source Integrative Progressive Target Discovery`.

`MIPTD` is a single-CAS, online-first analysis package for reproducing the Figure 1 target-discovery workflow.  
The package starts from one real `CAS` number, resolves compound identity online, queries four target sources, and runs the complete `a -> b -> c -> d -> e` pipeline in a structured case directory.

## Project Title

**MIPTD: Multi-source Integrative Progressive Target Discovery**

License: `MIT`

## Project Overview

This project is designed to turn a single compound, identified by its `CAS` number, into a reproducible target-discovery workflow.

The pipeline does the following:

1. Resolves compound identity online through `PubChem`
2. Retrieves compound name, formula, `SMILES`, and `InChIKey`
3. Queries four target sources:
   - `SwissTargetPrediction`
   - `SEA`
   - `ChEMBL`
   - `PPB2`
4. Standardizes human target information to the gene level
5. Builds the Figure 1 analysis chain:
   - `a`: source collection workflow
   - `b`: Venn diagram
   - `c`: GO/KEGG enrichment and disease-related filtering
   - `d`: KEGG circlize chord diagram
   - `e`: Chemprop multi-task target prioritization

The basic unit of the package is:

- one `CAS`
- one dated case directory
- one complete result set

## Technology Stack

The project uses a mixed Python + R stack.

### Python

- `Python 3.11`
- `requests`
- `beautifulsoup4`
- `playwright`
- `chemprop`
- `matplotlib` (declared in `environment.yml` as `matplotlib-base`, used by `scripts/plot_figure1e.py`)
- `pillow` (used by `scripts/build_figure1.py` for image composition)
- standard library modules for API access, HTML parsing, and pipeline orchestration

### R

- `r-base`
- `tidyverse`
- `readxl`
- `writexl`
- `ggvenn`
- `clusterProfiler`
- `org.Hs.eg.db`
- `ggnewscale`
- `circlize`
- `ComplexHeatmap`
- `legendry` (installed from CRAN by `scripts/install_miptd.sh`; not on conda-forge / bioconda)

### External Online Sources

- `PubChem PUG REST`
- `SwissTargetPrediction`
- `SEA`
- `ChEMBL API`
- `PPB2`

## Repository Structure

- [`src/miptd`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd)
  Python package source code
- [`scripts`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts)
  analysis, crawling, and resource-generation scripts
- [`src/miptd/resources`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/resources)
  bundled formal runtime resource files
- [`doc`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/doc)
  workflow and package documentation

## Bundled Resources

The package currently includes two formal resource files:

- [`src/miptd/resources/idmapping.tsv`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/resources/idmapping.tsv)
  `UNIPROT -> SYMBOL` mapping for human genes
- [`src/miptd/resources/ChEMBL_target_catalog.csv`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/resources/ChEMBL_target_catalog.csv)
  human `ChEMBL target` catalog used for target normalization and Chemprop task preparation

These files are part of the package runtime requirements and must be shipped together with the code.

Resource generation scripts are also preserved in [`scripts`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts):

- [`scripts/build_idmapping_tsv.R`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/build_idmapping_tsv.R)
- [`scripts/build_chembl_target_catalog.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/build_chembl_target_catalog.py)

## Installation

Recommended environment setup:

### Recommended Install Script

Use the bundled installer script. It creates a minimal environment first, then installs `mamba` inside that environment, then installs the remaining dependencies in the correct order.

```bash
bash scripts/install_miptd.sh
```

Custom environment name:

```bash
bash scripts/install_miptd.sh my_miptd_env
```

### Manual Install

Use this if you want to run each step yourself.

1. Create the isolated conda environment.
2. Install `mamba` inside that environment.
3. Install the conda dependencies with `mamba`.
4. Install the missing R package `legendry`.
5. Install the package in editable mode.

```bash
conda create -n miptd -y python=3.11 pip
conda install -n miptd -y -c conda-forge mamba
conda run -n miptd conda install -y -c pytorch -c bioconda -c conda-forge \
  nodejs pytorch rdkit r-base r-tidyverse r-readxl r-writexl r-ggvenn \
  r-jsonlite bioconductor-clusterprofiler bioconductor-org.hs.eg.db \
  r-ggnewscale r-circlize bioconductor-complexheatmap \
  matplotlib-base pillow
conda run -n miptd Rscript -e "install.packages('legendry', repos='https://cloud.r-project.org')"
conda run -n miptd pip install -e . chemprop==2.2.2
```

If browser automation is needed for `SwissTargetPrediction`, make sure the local browser dependency is available.

`legendry` is installed as an additional R package after environment creation because it is not bundled in the conda environment definition.
This two-stage install is preferred over `conda env create -f environment.yml` because it ensures the heavy dependency installation happens after `mamba` is available inside the target environment.

### Fully Locked Install

Use this when you want the exact verified package versions from the validated `miptd` environment.

```bash
conda env create -f environment.lock.yml
conda activate miptd
Rscript -e "install.packages('legendry', repos='https://cloud.r-project.org')"
pip install -e .
```

Locked environment file:

- [`environment.lock.yml`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/environment.lock.yml)

## Usage

### Run the Full Pipeline

```bash
MIPTD \
  --cas 491-71-4 \
  --disase-keywords "NAFLD,liver,hepatic,steatosis"
```

Main arguments:

- `--cas`
  input CAS number
- `--disease-keywords` / `--disase-keywords`
  comma-separated disease keywords
- `--output-root`
  output root directory
- `--run-date`
  optional date suffix for the case directory
- `--project-root`
  optional repository root
- `--dry-run`
  print the planned workflow without executing

### Keyword-Based Filtering

The filtering used in `Step C` and `Step D` is driven by the keywords passed through:

```bash
--disease-keywords "NAFLD,liver,hepatic"
```

The package first runs full GO/KEGG enrichment on the union target set, then converts the comma-separated keywords into a regex-style matching rule, and finally applies that rule to the enrichment `Description` field.

The actual logic is:

1. run full GO/KEGG enrichment on the union target set
2. build a keyword-matching rule from `--disease-keywords`
3. match that rule against GO/KEGG `Description`
4. preferentially keep the matched rows for panel `c`
5. derive panel `d` genes from the selected KEGG rows

Fallback rule:

- if none of the descriptions match the user-provided keywords, the workflow falls back to top-ranked rows by significance and count so `c` and `d` are still produced

### Validate a Case Directory

```bash
MIPTD-validate \
  --case-dir CAS_491-71-4_2026-03-21
```

## Tests

Run the local unit and lightweight integration tests with:

```bash
python3 -m unittest discover -s tests -v
```

## Output

Each run creates one independent case directory:

```text
CAS_491-71-4_YYYY-MM-DD/
  compound_identity.json
  compound_identity.csv
  status.json
  run.log
  a_source_collection/
  b_venn/
  c_go_kegg/
  d_kegg_circlize/
  e_chemprop/
```

This structure enforces the core rule:

- one `CAS`
- one case directory
- one result set

## Failure Handling

The package uses an explicit degradation strategy for online source collection:

- failure of a single source database does not immediately abort the run
- failed sources are recorded in `status.json`
- a placeholder empty source file is written so the case remains structurally consistent
- the pipeline continues only if at least 2 source databases produce non-empty human target sets after normalization
- if fewer than 2 non-empty sources remain, the run stops with a clear error

## Documentation

Detailed documentation is available in [`doc`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/doc):

- [`doc/README.md`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/doc/README.md)
- [`doc/包功能与使用说明.md`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/doc/包功能与使用说明.md)
- [`doc/单CAS_Figure1分析流程说明.md`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/doc/单CAS_Figure1分析流程说明.md)
- [`doc/Single_CAS_Figure1_Pipeline.md`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/doc/Single_CAS_Figure1_Pipeline.md)
- [`TROUBLESHOOTING.md`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/TROUBLESHOOTING.md)
- [`TROUBLESHOOTING_CN.md`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/TROUBLESHOOTING_CN.md)

## Current Status

The package currently has:

- single-CAS workflow design
- online compound identity resolution
- four-source target collection
- structured `a -> e` pipeline orchestration
- bundled `idmapping.tsv` and `ChEMBL target catalog`
- case validation support
- pipeline logging support

## Release Notes

Before packaging or publishing:

- keep `src/miptd/resources/`
- keep `src/miptd/resources/idmapping.tsv`
- keep `src/miptd/resources/ChEMBL_target_catalog.csv`
- verify `MIPTD --help`
- verify `MIPTD-validate --help`
- verify the environment can be created from [`environment.yml`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/environment.yml)
