#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
except ModuleNotFoundError as exc:  # pragma: no cover - environment guard
    sys.stderr.write(
        f"[plot_figure1e] Missing required dependency: {exc.name}.\n"
        "Install it into the active MIPTD conda environment with:\n"
        "    conda install -c conda-forge matplotlib-base\n"
        "or, if conda is unavailable, fall back to:\n"
        "    pip install matplotlib\n"
    )
    sys.exit(2)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot Figure 1e from Chemprop ranked targets and test metrics."
    )
    parser.add_argument("--ranked-targets-csv", required=True)
    parser.add_argument("--metrics-csv", required=True)
    parser.add_argument("--cas", required=True)
    parser.add_argument("--compound-name", required=True)
    parser.add_argument("--output-png", required=True)
    parser.add_argument("--output-pdf", default=None)
    parser.add_argument("--top-n", type=int, default=8)
    parser.add_argument("--filter-status", default=None, help="Optional filter_status value to keep, e.g. keep")
    return parser.parse_args()


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def safe_float(value):
    if value in (None, ""):
        return None
    return float(value)


def wrap_text(text, max_len=28):
    words = text.split()
    lines = []
    current = []
    current_len = 0
    for word in words:
        extra = len(word) + (1 if current else 0)
        if current and current_len + extra > max_len:
            lines.append(" ".join(current))
            current = [word]
            current_len = len(word)
        else:
            current.append(word)
            current_len += extra
    if current:
        lines.append(" ".join(current))
    return "\n".join(lines)


def main():
    args = parse_args()
    project_root = Path(__file__).resolve().parents[1]
    ranked_path = (project_root / args.ranked_targets_csv).resolve() if not Path(args.ranked_targets_csv).is_absolute() else Path(args.ranked_targets_csv)
    metrics_path = (project_root / args.metrics_csv).resolve() if not Path(args.metrics_csv).is_absolute() else Path(args.metrics_csv)
    output_png = (project_root / args.output_png).resolve() if not Path(args.output_png).is_absolute() else Path(args.output_png)
    output_pdf = None
    if args.output_pdf:
        output_pdf = (project_root / args.output_pdf).resolve() if not Path(args.output_pdf).is_absolute() else Path(args.output_pdf)

    ranked_rows = [row for row in read_csv(ranked_path) if row["cas"] == args.cas]
    if args.filter_status:
        ranked_rows = [row for row in ranked_rows if row.get("filter_status") == args.filter_status]
    ranked_rows = sorted(ranked_rows, key=lambda row: safe_float(row["predicted_pchembl"]) or -999, reverse=True)[: args.top_n]
    metric_rows = read_csv(metrics_path)
    metric_map = {row["target"]: row for row in metric_rows}

    targets = [row["target"] for row in ranked_rows]
    scores = [safe_float(row["predicted_pchembl"]) for row in ranked_rows]
    train_n = [int(row["n_training_molecules"]) for row in ranked_rows]
    test_r2 = [safe_float(metric_map.get(row["target"], {}).get("r2")) for row in ranked_rows]

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig = plt.figure(figsize=(14, 7.5), dpi=300)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 1.45], wspace=0.12)
    ax_left = fig.add_subplot(gs[0, 0])
    ax_right = fig.add_subplot(gs[0, 1])

    ax_left.set_xlim(0, 1)
    ax_left.set_ylim(0, 1)
    ax_left.axis("off")

    box_specs = [
        (0.08, 0.72, 0.34, 0.16, "#f6e7c1", "d-derived KEGG genes", "Key genes from panel d"),
        (0.56, 0.72, 0.34, 0.16, "#dcefdc", "Human ChEMBL activities", "Trainable targets with enough bioactivity records"),
        (0.08, 0.42, 0.34, 0.16, "#d9e8f6", "Chemprop multi-task", "Shared graph neural network"),
        (0.56, 0.42, 0.34, 0.16, "#f6d9d4", "Target ranking", "Predicted pChEMBL for the query compound"),
    ]
    for x, y, w, h, color, title, subtitle in box_specs:
        patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.03",
                               facecolor=color, edgecolor="#2f2f2f", linewidth=1.2)
        ax_left.add_patch(patch)
        ax_left.text(x + w / 2, y + h * 0.62, title, ha="center", va="center", fontsize=12, weight="bold")
        ax_left.text(x + w / 2, y + h * 0.26, wrap_text(subtitle, 24), ha="center", va="center", fontsize=9)

    arrows = [
        ((0.42, 0.80), (0.56, 0.80)),
        ((0.25, 0.72), (0.25, 0.58)),
        ((0.73, 0.72), (0.73, 0.58)),
        ((0.42, 0.50), (0.56, 0.50)),
    ]
    for start, end in arrows:
        ax_left.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=14, linewidth=1.5, color="#333333"))

    ax_left.text(0.08, 0.95, "e", fontsize=22, weight="bold")
    ax_left.text(0.08, 0.08, f"{args.compound_name}\nCAS {args.cas}", fontsize=12, weight="bold", ha="left", va="bottom")
    footer = "Panel e rebuilt strictly from panel d KEGG genes"
    if args.filter_status:
        footer += f"\nOnly targets with filter_status = {args.filter_status}"
    ax_left.text(0.08, 0.02, footer, fontsize=9, color="#555555", ha="left", va="bottom")

    colors = []
    for value in test_r2:
        if value is None:
            colors.append("#bdbdbd")
        elif value >= 0.5:
            colors.append("#c95f3c")
        elif value >= 0.3:
            colors.append("#e39f5b")
        else:
            colors.append("#8aa6c8")

    y_pos = list(range(len(targets)))
    bars = ax_right.barh(y_pos, scores, color=colors, edgecolor="#2f2f2f", linewidth=0.8)
    ax_right.set_yticks(y_pos, targets)
    ax_right.invert_yaxis()
    ax_right.set_xlabel("Predicted pChEMBL")
    ax_right.set_title(f"Top predicted targets for {args.compound_name}", fontsize=13, pad=10)
    ax_right.grid(axis="x", color="#d8d8d8", linewidth=0.7, alpha=0.8)
    ax_right.set_axisbelow(True)

    max_score = max(scores) if scores else 7
    ax_right.set_xlim(0, max_score + 1.2)
    for idx, bar in enumerate(bars):
        score = scores[idx]
        r2_value = test_r2[idx]
        train_value = train_n[idx]
        label = f"{score:.2f}"
        if r2_value is not None:
            label += f" | R2 {r2_value:.2f}"
        label += f" | n={train_value}"
        ax_right.text(bar.get_width() + 0.06, bar.get_y() + bar.get_height() / 2, label, va="center", fontsize=9)

    ax_right.text(
        0.99,
        -0.12,
        "Color by test R2: >=0.50 high, 0.30-0.49 moderate, <0.30 weak",
        transform=ax_right.transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
        color="#555555",
    )

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, bbox_inches="tight")
    if output_pdf:
        fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
