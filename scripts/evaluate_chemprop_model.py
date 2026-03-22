#!/usr/bin/env python3

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
from pathlib import Path


DEFAULT_MPLCONFIGDIR = "tmp/mplconfig"
DEFAULT_XDG_CACHE_HOME = "tmp/fontconfig"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate a trained Chemprop model on the saved test split and export per-target metrics."
    )
    parser.add_argument("--model-dir", required=True, help="Chemprop model output dir containing training/ and run_metadata.json.")
    parser.add_argument("--output-dir", required=True, help="Directory for evaluation outputs.")
    parser.add_argument("--chemprop-bin", default=None, help="Path to the chemprop executable.")
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


def ensure_env(project_root):
    env = os.environ.copy()
    mpl_dir = project_root / DEFAULT_MPLCONFIGDIR
    xdg_dir = project_root / DEFAULT_XDG_CACHE_HOME
    mpl_dir.mkdir(parents=True, exist_ok=True)
    xdg_dir.mkdir(parents=True, exist_ok=True)
    env["MPLCONFIGDIR"] = str(mpl_dir)
    env["XDG_CACHE_HOME"] = str(xdg_dir)
    return env


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def safe_float(value):
    if value in (None, ""):
        return None
    return float(value)


def rmse(y_true, y_pred):
    if not y_true:
        return None
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(y_true, y_pred)) / len(y_true))


def mae(y_true, y_pred):
    if not y_true:
        return None
    return sum(abs(a - b) for a, b in zip(y_true, y_pred)) / len(y_true)


def r2(y_true, y_pred):
    if not y_true:
        return None
    mean_true = sum(y_true) / len(y_true)
    ss_res = sum((a - b) ** 2 for a, b in zip(y_true, y_pred))
    ss_tot = sum((a - mean_true) ** 2 for a in y_true)
    if ss_tot == 0:
        return None
    return 1 - (ss_res / ss_tot)


def main():
    args = parse_args()
    project_root = Path(__file__).resolve().parents[1]
    model_dir = (project_root / args.model_dir).resolve() if not Path(args.model_dir).is_absolute() else Path(args.model_dir)
    output_dir = (project_root / args.output_dir).resolve() if not Path(args.output_dir).is_absolute() else Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = json.loads((model_dir / "run_metadata.json").read_text(encoding="utf-8"))
    chemprop_bin = detect_chemprop_bin(args.chemprop_bin)
    env = ensure_env(project_root)

    train_dir = model_dir / "training"
    data_csv = Path(metadata["data_csv"])
    smiles_column = metadata["smiles_column"]
    target_columns = metadata["target_columns"]
    test_smiles_path = train_dir / "test_smiles.csv"
    test_predictions_path = output_dir / "test_predictions.csv"

    predict_cmd = [
        chemprop_bin,
        "predict",
        "--test-path",
        str(test_smiles_path),
        "--model-paths",
        str(train_dir),
        "--smiles-columns",
        smiles_column,
        "--output",
        str(test_predictions_path),
    ]
    subprocess.run(predict_cmd, check=True, env=env, cwd=project_root)

    data_rows = read_csv(data_csv)
    truth_by_smiles = {row[smiles_column]: row for row in data_rows}
    pred_rows = read_csv(test_predictions_path)

    metrics_rows = []
    for target in target_columns:
        y_true = []
        y_pred = []
        for row in pred_rows:
            smiles = row[smiles_column]
            truth_row = truth_by_smiles.get(smiles)
            if not truth_row:
                continue
            true_value = safe_float(truth_row.get(target))
            pred_value = safe_float(row.get(target))
            if true_value is None or pred_value is None:
                continue
            y_true.append(true_value)
            y_pred.append(pred_value)

        metrics_rows.append(
            {
                "target": target,
                "rmse": rmse(y_true, y_pred),
                "mae": mae(y_true, y_pred),
                "r2": r2(y_true, y_pred),
                "n_test_molecules": len(y_true),
            }
        )

    write_csv(
        output_dir / "test_metrics_by_target.csv",
        metrics_rows,
        ["target", "rmse", "mae", "r2", "n_test_molecules"],
    )

    valid_rmse = [row["rmse"] for row in metrics_rows if row["rmse"] is not None]
    valid_mae = [row["mae"] for row in metrics_rows if row["mae"] is not None]
    valid_r2 = [row["r2"] for row in metrics_rows if row["r2"] is not None]
    summary = {
        "n_targets": len(metrics_rows),
        "mean_rmse": sum(valid_rmse) / len(valid_rmse) if valid_rmse else None,
        "mean_mae": sum(valid_mae) / len(valid_mae) if valid_mae else None,
        "mean_r2": sum(valid_r2) / len(valid_r2) if valid_r2 else None,
    }
    (output_dir / "test_metrics_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
