# Changelog

## 0.1.0

- Rebranded the package and CLI to `MIPTD` / `MIPTD-validate`.
- Standardized the package layout under [`src/miptd`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd).
- Switched the workflow to a single-CAS, online-first pipeline.
- Added bundled formal runtime resources under [`src/miptd/resources`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/src/miptd/resources):
  - `idmapping.tsv`
  - `ChEMBL_target_catalog.csv`
- Added resource generation scripts:
  - [`scripts/build_idmapping_tsv.R`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/build_idmapping_tsv.R)
  - [`scripts/build_chembl_target_catalog.py`](/Users/liuyue/Desktop/workspace/Cadisum_MM_Project/Nature_Communication_Figure1/scripts/build_chembl_target_catalog.py)
- Added explicit failure handling with `status.json` and degraded-source continuation rules.
- Added local tests for:
  - package resource discovery
  - case directory creation
  - CLI argument parsing
  - dry-run pipeline execution
  - case validation
- Verified clean installation in an isolated `miptd` conda environment.
