# MIPTD

## Package Mode

This repository is being refactored into a **single-CAS, online-first pipeline** for reproducing the Figure 1 analysis workflow.

`MIPTD` stands for `Multi-source Integrative Progressive Target Discovery`.

Core package files:

- [`pyproject.toml`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/pyproject.toml)
- [`environment.yml`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/environment.yml)
- [`src/miptd/cli.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/cli.py)
- [`src/miptd/pipeline.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/pipeline.py)
- [`src/miptd/discovery.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/discovery.py)

## Environment

Recommended setup:

1. Create the isolated conda environment.
2. Activate `miptd`.
3. The environment includes `mamba` for follow-up package management.
4. Install `legendry`.
5. Install the package itself.

```bash
conda env create -f environment.yml
conda activate miptd
Rscript -e "install.packages('legendry', repos='https://cloud.r-project.org')"
pip install -e .
```

The release bundle must include the package-level [`src/miptd/resources`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/resources) directory.

Required bundled resources:

- [`src/miptd/resources/idmapping.tsv`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/resources/idmapping.tsv)
- [`src/miptd/resources/ChEMBL_target_catalog.csv`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/resources/ChEMBL_target_catalog.csv)

These two files are treated as formal package resources rather than temporary development artifacts.

## Online-First Design

The intended production workflow is:

1. input one real `CAS`
2. resolve compound identity online
3. obtain compound name, molecular formula, SMILES, and InChIKey automatically
4. query the four target sources online:
   - `SwissTargetPrediction`
   - `SEA`
   - `ChEMBL`
   - `PPB2`
5. run the full `a -> b -> c -> d -> e` workflow

Current compound identity resolution is handled through **PubChem PUG REST** from the CAS number.
In production, `CAS` is only the entry key. The package first resolves the compound identity, then submits the resolved `SMILES` to `SwissTargetPrediction`, `SEA`, and `PPB2`. `ChEMBL` matching uses the resolved name, InChIKey, and SMILES together.

Local fallback tables may exist in development environments, but they are not the intended primary path for the package.

## Command Line

Package-style entry:

```bash
MIPTD \
  --cas 491-71-4 \
  --disase-keywords "NAFLD,liver,hepatic,steatosis"
```

Validation entry:

```bash
MIPTD-validate \
  --case-dir CAS_491-71-4_2026-03-21
```

## Release Notes

For installation and release packaging, keep these constraints:

- ship the `src/miptd/resources/` directory together with the code
- do not remove `src/miptd/resources/idmapping.tsv`
- do not remove `src/miptd/resources/ChEMBL_target_catalog.csv`
- verify that the runtime can resolve both resource files before publishing a release

## Output Layout

Expected output structure:

```text
CAS_491-71-4_YYYY-MM-DD/
  compound_identity.json
  compound_identity.csv
  a_source_collection/
    compound_identity.json
    compound_identity.csv
  b_venn/
  c_go_kegg/
  d_kegg_circlize/
  e_chemprop/
```

Each case directory represents one CAS-only project.
The identity files are written before source fetching so the resolved compound metadata remains traceable even if one of the online sources times out.

## Current Status

- the package entry for `a -> e` is defined and runnable
- CAS-based compound identity can be resolved online through PubChem
- `Swiss`, `SEA`, `ChEMBL`, and `PPB2` all have online fetch scripts
- a real online single-CAS run has been verified through `e_chemprop`
- `e` is now constrained to a single query compound rather than replaying legacy demo compounds
- resolved compound identity is saved as standalone JSON/CSV outputs in the case root and `a_source_collection`
- case outputs are organized as `CAS_xxxxx_Date`
- case outputs can be checked with `MIPTD-validate` for directory completeness and single-CAS consistency

## Documentation

Detailed workflow and directory rules:

- [`doc/Single_CAS_Figure1_Pipeline.md`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/doc/Single_CAS_Figure1_Pipeline.md)
- [`doc/单CAS_Figure1分析流程说明.md`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/doc/单CAS_Figure1分析流程说明.md)
