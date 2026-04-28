# Troubleshooting

## Installation Problems

### `conda env create -f environment.yml` is slow

This environment mixes:

- Python
- PyTorch
- RDKit
- Node.js
- R
- Bioconductor

That combination is heavy for dependency solving. This is expected.

Recommended approach:

1. Create the `miptd` environment.
2. Activate it.
3. Use the bundled `mamba` inside `miptd` for any follow-up package operations.

```bash
conda env create -f environment.yml
conda activate miptd
```

### `legendry` is missing in R

`legendry` is not currently included in the conda environment definition.

Install it manually inside `miptd`:

```bash
Rscript -e "install.packages('legendry', repos='https://cloud.r-project.org')"
```

### `MIPTD` command is not found

Make sure both conditions are true:

1. the `miptd` environment is activated
2. the package has been installed in editable mode

```bash
conda activate miptd
pip install -e .
```

Then check:

```bash
MIPTD --help
MIPTD-validate --help
```

## Runtime Problems

### `SwissTargetPrediction` times out

This is one of the most likely external failure points.

Current package behavior:

- the Swiss source is marked as failed in `status.json`
- an empty placeholder file is created
- the pipeline continues only if at least 2 source databases still produce non-empty human target sets

Check:

- `run.log`
- `status.json`
- `a_source_collection/swiss_fetch/`

### One source fails but the run continues

This is expected behavior.

The package uses a degraded-source strategy:

- single-source failure does not immediately abort the run
- the failure is recorded in `status.json`
- downstream steps continue only if evidence is still sufficient

### The run stops after Step A

If fewer than 2 source databases produce non-empty human target sets after normalization, the pipeline stops on purpose.

This prevents a false-success case with insufficient evidence.

Check:

- `status.json`
- `a_source_collection/source_summary.json`

## Resource Problems

### `idmapping.tsv` or `ChEMBL_target_catalog.csv` cannot be found

The package expects these formal resources inside:

- [`src/miptd/resources/idmapping.tsv`](src/miptd/resources/idmapping.tsv)
- [`src/miptd/resources/ChEMBL_target_catalog.csv`](src/miptd/resources/ChEMBL_target_catalog.csv)

If they are missing, rebuild them with:

```bash
Rscript scripts/build_idmapping_tsv.R --output src/miptd/resources/idmapping.tsv
python3 scripts/build_chembl_target_catalog.py \
  --output-csv src/miptd/resources/ChEMBL_target_catalog.csv \
  --summary-json src/miptd/resources/ChEMBL_target_catalog.summary.json
```

## Validation Problems

### `MIPTD-validate` fails

The validator checks:

- required step directories
- required output files
- that `inference_template.csv` contains exactly one compound
- that `prediction_summary.json` contains exactly one compound
- that CAS values match the case manifest

If validation fails, inspect:

- `case_manifest.json`
- `e_chemprop/chemprop_data/inference_template.csv`
- `e_chemprop/chemprop_model/analysis/prediction_summary.json`

## Recommended Debug Files

For most failures, inspect these files first:

- `run.log`
- `status.json`
- `case_manifest.json`
- `a_source_collection/source_summary.json`

