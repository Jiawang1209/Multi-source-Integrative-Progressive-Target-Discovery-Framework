#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize Chemprop predictions into ranked target tables."
    )
    parser.add_argument("--predictions-csv", required=True)
    parser.add_argument("--task-summary-csv", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--metrics-csv", default=None)
    parser.add_argument("--compound-column", default="compound")
    parser.add_argument("--cas-column", default="cas")
    parser.add_argument("--name-column", default="name")
    parser.add_argument("--smiles-column", default="canonical_smiles")
    parser.add_argument("--min-test-r2", type=float, default=0.30)
    parser.add_argument("--min-test-molecules", type=int, default=30)
    parser.add_argument("--min-training-molecules", type=int, default=100)
    return parser.parse_args()


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def safe_float(value):
    if value in (None, ""):
        return None
    return float(value)


def read_metrics(path):
    if not path:
        return {}
    rows = read_csv(path)
    return {
        row["target"]: {
            "test_rmse": safe_float(row.get("rmse")),
            "test_mae": safe_float(row.get("mae")),
            "test_r2": safe_float(row.get("r2")),
            "n_test_molecules": int(row["n_test_molecules"]) if row.get("n_test_molecules") not in (None, "") else None,
        }
        for row in rows
    }


def build_filter_status(metrics, n_training_molecules, min_test_r2, min_test_molecules, min_training_molecules):
    reasons = []
    keep = True

    test_r2 = metrics.get("test_r2")
    n_test_molecules = metrics.get("n_test_molecules")

    if n_training_molecules is None or n_training_molecules < min_training_molecules:
        keep = False
        reasons.append(f"train_n<{min_training_molecules}")

    if n_test_molecules is None or n_test_molecules < min_test_molecules:
        keep = False
        reasons.append(f"test_n<{min_test_molecules}")

    if test_r2 is None or test_r2 < min_test_r2:
        keep = False
        reasons.append(f"R2<{min_test_r2:.2f}")

    if keep:
        return "keep", f"pass: train_n>={min_training_molecules}; test_n>={min_test_molecules}; R2>={min_test_r2:.2f}"
    return "filtered", "; ".join(reasons)


def main():
    args = parse_args()
    project_root = Path(__file__).resolve().parents[1]
    predictions_csv = (project_root / args.predictions_csv).resolve() if not Path(args.predictions_csv).is_absolute() else Path(args.predictions_csv)
    task_summary_csv = (project_root / args.task_summary_csv).resolve() if not Path(args.task_summary_csv).is_absolute() else Path(args.task_summary_csv)
    output_dir = (project_root / args.output_dir).resolve() if not Path(args.output_dir).is_absolute() else Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics_csv = None
    if args.metrics_csv:
        metrics_csv = (project_root / args.metrics_csv).resolve() if not Path(args.metrics_csv).is_absolute() else Path(args.metrics_csv)

    prediction_rows = read_csv(predictions_csv)
    task_rows = read_csv(task_summary_csv)
    metric_meta = read_metrics(metrics_csv)
    task_meta = {
        row["target"]: {
            "consensus_score": safe_float(row.get("consensus_score")),
            "platform_vote": int(row["platform_vote"]),
            "pathway_n": int(row["pathway_n"]),
            "sources": row["sources"],
            "target_chembl_ids": row["target_chembl_ids"],
            "n_unique_molecules": int(row["n_unique_molecules"]),
        }
        for row in task_rows
    }

    target_columns = [
        column for column in prediction_rows[0].keys()
        if column not in {args.compound_column, args.cas_column, args.name_column, args.smiles_column}
    ]

    ranked_rows = []
    summary = {"compounds": []}
    for row in prediction_rows:
        compound_name = row.get(args.compound_column, "")
        cas = row.get(args.cas_column, "")
        display_name = row.get(args.name_column, "")
        compound_rows = []
        for target in target_columns:
            pred = safe_float(row.get(target))
            if pred is None:
                continue
            meta = task_meta.get(target, {})
            metrics = metric_meta.get(target, {})
            n_training_molecules = meta.get("n_unique_molecules")
            filter_status, filter_basis = build_filter_status(
                metrics,
                n_training_molecules,
                args.min_test_r2,
                args.min_test_molecules,
                args.min_training_molecules,
            )
            compound_rows.append(
                {
                    "compound": compound_name,
                    "cas": cas,
                    "name": display_name,
                    "target": target,
                    "predicted_pchembl": pred,
                    "consensus_score": meta.get("consensus_score"),
                    "platform_vote": meta.get("platform_vote"),
                    "pathway_n": meta.get("pathway_n"),
                    "n_training_molecules": meta.get("n_unique_molecules"),
                    "sources": meta.get("sources"),
                    "target_chembl_ids": meta.get("target_chembl_ids"),
                    "test_rmse": metrics.get("test_rmse"),
                    "test_mae": metrics.get("test_mae"),
                    "test_r2": metrics.get("test_r2"),
                    "n_test_molecules": metrics.get("n_test_molecules"),
                    "filter_basis": filter_basis,
                    "filter_status": filter_status,
                }
            )
        compound_rows.sort(
            key=lambda x: (
                -(x["predicted_pchembl"] if x["predicted_pchembl"] is not None else -999),
                -(x["consensus_score"] if x["consensus_score"] is not None else -999),
                x["target"],
            )
        )
        ranked_rows.extend(compound_rows)
        summary["compounds"].append(
            {
                "compound": compound_name,
                "cas": cas,
                "name": display_name,
                "top_targets": compound_rows[:5],
            }
        )

    ranked_path = output_dir / "ranked_targets.csv"
    with open(ranked_path, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "compound",
            "cas",
            "name",
            "target",
            "predicted_pchembl",
            "consensus_score",
            "platform_vote",
            "pathway_n",
            "n_training_molecules",
            "sources",
            "target_chembl_ids",
            "test_rmse",
            "test_mae",
            "test_r2",
            "n_test_molecules",
            "filter_basis",
            "filter_status",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(ranked_rows)

    with open(output_dir / "prediction_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
