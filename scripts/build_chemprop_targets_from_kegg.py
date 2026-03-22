#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate d-panel KEGG gene results into a Chemprop target candidate table."
    )
    parser.add_argument("--kegg-genes-csv", required=True)
    parser.add_argument("--output-csv", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.kegg_genes_csv)
    output_path = Path(args.output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    pathways_by_gene = defaultdict(set)
    pathway_ids_by_gene = defaultdict(set)
    with input_path.open(newline="", encoding="utf-8") as f:
      reader = csv.DictReader(f)
      for row in reader:
        gene = (row.get("geneID") or "").strip().upper()
        if not gene:
            continue
        pathways_by_gene[gene].add((row.get("Description") or "").strip())
        pathway_ids_by_gene[gene].add((row.get("ID") or "").strip())

    rows = []
    for gene, pathways in pathways_by_gene.items():
        pathway_names = sorted(x for x in pathways if x)
        pathway_ids = sorted(x for x in pathway_ids_by_gene[gene] if x)
        pathway_n = len(pathway_names)
        rows.append(
            {
                "target": gene,
                "consensus_score": float(pathway_n),
                "platform_vote": 2,
                "pathway_n": pathway_n,
                "best_pathway_padj": "",
                "pathways": " | ".join(pathway_names),
                "sources": "KEGG_from_d",
                "pathway_ids": "|".join(pathway_ids),
            }
        )

    rows.sort(key=lambda x: (-x["pathway_n"], x["target"]))
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "target",
                "consensus_score",
                "platform_vote",
                "pathway_n",
                "best_pathway_padj",
                "pathways",
                "sources",
                "pathway_ids",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
