#!/usr/bin/env python3

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


DEFAULT_MPLCONFIGDIR = "tmp/mplconfig"
DEFAULT_XDG_CACHE_HOME = "tmp/fontconfig"
DEFAULT_METADATA_FILENAME = "run_metadata.json"
DEFAULT_TRAINING_METRICS = ["rmse", "mae"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train a Chemprop multi-task model and run inference on query compounds."
    )
    parser.add_argument("--data-csv", required=True, help="Chemprop multi-task training CSV.")
    parser.add_argument("--inference-csv", required=True, help="CSV containing compounds for prediction.")
    parser.add_argument("--output-dir", required=True, help="Directory for model outputs and predictions.")
    parser.add_argument("--chemprop-bin", default=None, help="Path to the chemprop executable.")
    parser.add_argument("--smiles-column", default="canonical_smiles")
    parser.add_argument("--ignore-columns", nargs="*", default=["molecule_chembl_id"])
    parser.add_argument("--target-columns", nargs="*", default=None)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--split", default="SCAFFOLD_BALANCED")
    parser.add_argument("--split-sizes", nargs=3, type=float, default=[0.8, 0.1, 0.1])
    parser.add_argument("--data-seed", type=int, default=2026)
    parser.add_argument("--pytorch-seed", type=int, default=2026)
    parser.add_argument("--accelerator", default="cpu")
    parser.add_argument("--devices", default="auto")
    parser.add_argument("--message-hidden-dim", type=int, default=300)
    parser.add_argument("--depth", type=int, default=3)
    parser.add_argument("--dropout", type=float, default=0.0)
    parser.add_argument("--ffn-hidden-dim", type=int, default=300)
    parser.add_argument("--ffn-num-layers", type=int, default=2)
    parser.add_argument("--ensemble-size", type=int, default=1)
    parser.add_argument("--metrics", nargs="*", default=DEFAULT_TRAINING_METRICS)
    parser.add_argument("--molecule-featurizers", nargs="*", default=None)
    return parser.parse_args()


def detect_chemprop_bin(explicit_path):
    if explicit_path:
        return explicit_path

    candidates = [
        os.environ.get("CHEMPROP_BIN"),
        "/Users/liuyue/miniconda3/envs/Chemprop/bin/chemprop",
        shutil.which("chemprop"),
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return candidate

    raise FileNotFoundError("Unable to locate `chemprop` executable. Pass --chemprop-bin explicitly.")


def detect_target_columns(data_csv, smiles_column, ignore_columns, target_columns):
    if target_columns:
        return target_columns

    with open(data_csv, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)

    ignored = set(ignore_columns) | {smiles_column}
    inferred = [column for column in header if column not in ignored]
    if not inferred:
        raise ValueError("No target columns detected in training CSV.")

    return inferred


def ensure_env(project_root):
    env = os.environ.copy()
    mpl_dir = project_root / DEFAULT_MPLCONFIGDIR
    xdg_dir = project_root / DEFAULT_XDG_CACHE_HOME
    mpl_dir.mkdir(parents=True, exist_ok=True)
    xdg_dir.mkdir(parents=True, exist_ok=True)
    env["MPLCONFIGDIR"] = str(mpl_dir)
    env["XDG_CACHE_HOME"] = str(xdg_dir)
    return env


def run_command(cmd, env, cwd):
    print("$", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, env=env, cwd=cwd)


def write_metadata(path, payload):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)


def main():
    args = parse_args()

    project_root = Path(__file__).resolve().parents[1]
    output_dir = (project_root / args.output_dir).resolve() if not Path(args.output_dir).is_absolute() else Path(args.output_dir)
    data_csv = (project_root / args.data_csv).resolve() if not Path(args.data_csv).is_absolute() else Path(args.data_csv)
    inference_csv = (project_root / args.inference_csv).resolve() if not Path(args.inference_csv).is_absolute() else Path(args.inference_csv)
    output_dir.mkdir(parents=True, exist_ok=True)

    chemprop_bin = detect_chemprop_bin(args.chemprop_bin)
    target_columns = detect_target_columns(data_csv, args.smiles_column, args.ignore_columns, args.target_columns)
    env = ensure_env(project_root)

    train_cmd = [
        chemprop_bin,
        "train",
        "--data-path",
        str(data_csv),
        "--output-dir",
        str(output_dir / "training"),
        "--smiles-columns",
        args.smiles_column,
        "--target-columns",
        *target_columns,
        "--ignore-columns",
        *args.ignore_columns,
        "--task-type",
        "regression",
        "--epochs",
        str(args.epochs),
        "--patience",
        str(args.patience),
        "--batch-size",
        str(args.batch_size),
        "--num-workers",
        str(args.num_workers),
        "--split",
        args.split,
        "--split-sizes",
        *(str(x) for x in args.split_sizes),
        "--data-seed",
        str(args.data_seed),
        "--pytorch-seed",
        str(args.pytorch_seed),
        "--accelerator",
        args.accelerator,
        "--devices",
        args.devices,
        "--message-hidden-dim",
        str(args.message_hidden_dim),
        "--depth",
        str(args.depth),
        "--dropout",
        str(args.dropout),
        "--ffn-hidden-dim",
        str(args.ffn_hidden_dim),
        "--ffn-num-layers",
        str(args.ffn_num_layers),
        "--ensemble-size",
        str(args.ensemble_size),
        "--metrics",
        *args.metrics,
        "--show-individual-scores",
        "--save-smiles-splits",
    ]
    if args.molecule_featurizers:
        train_cmd.extend(["--molecule-featurizers", *args.molecule_featurizers])

    predict_cmd = [
        chemprop_bin,
        "predict",
        "--test-path",
        str(inference_csv),
        "--model-paths",
        str(output_dir / "training"),
        "--smiles-columns",
        args.smiles_column,
        "--output",
        str(output_dir / "predictions.csv"),
    ]

    metadata = {
        "chemprop_bin": chemprop_bin,
        "data_csv": str(data_csv),
        "inference_csv": str(inference_csv),
        "output_dir": str(output_dir),
        "smiles_column": args.smiles_column,
        "target_columns": target_columns,
        "ignore_columns": args.ignore_columns,
        "epochs": args.epochs,
        "patience": args.patience,
        "batch_size": args.batch_size,
        "num_workers": args.num_workers,
        "split": args.split,
        "split_sizes": args.split_sizes,
        "data_seed": args.data_seed,
        "pytorch_seed": args.pytorch_seed,
        "accelerator": args.accelerator,
        "devices": args.devices,
        "message_hidden_dim": args.message_hidden_dim,
        "depth": args.depth,
        "dropout": args.dropout,
        "ffn_hidden_dim": args.ffn_hidden_dim,
        "ffn_num_layers": args.ffn_num_layers,
        "ensemble_size": args.ensemble_size,
        "metrics": args.metrics,
        "molecule_featurizers": args.molecule_featurizers,
    }
    write_metadata(output_dir / DEFAULT_METADATA_FILENAME, metadata)

    run_command(train_cmd, env=env, cwd=project_root)
    run_command(predict_cmd, env=env, cwd=project_root)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        print(f"Chemprop command failed with exit code {exc.returncode}.", file=sys.stderr)
        raise
