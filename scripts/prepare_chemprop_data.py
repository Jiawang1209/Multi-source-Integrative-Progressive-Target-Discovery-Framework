#!/usr/bin/env python3

import argparse
import csv
import json
import statistics
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold


ACTIVITY_TYPES = ("IC50", "Ki", "Kd", "EC50")
USER_AGENT = "ChempropDataPrep/1.0"

# Silence RDKit's stderr spam during the sanitize_smiles_rows pass — we deliberately
# probe each SMILES and want to count failures, not log them one by one.
RDLogger.DisableLog("rdApp.*")


def smiles_passes_scaffold_check(smiles: str) -> bool:
    """Whether a SMILES survives the same RDKit / MurckoScaffold path that
    chemprop's astartes splitter exercises during training.

    Recent RDKit versions (>=2024.03) raise hard `Pre-condition Violation:
    bad bond stereo` errors on malformed double-bond stereo annotations that
    older versions silently corrected. Such rows must be dropped before they
    reach `chemprop train`, otherwise the entire training run aborts.
    """
    if not smiles:
        return False
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:
        return False
    return True


def sanitize_smiles_rows(rows: list[dict]) -> tuple[list[dict], int]:
    """Drop rows whose `canonical_smiles` would crash RDKit downstream.

    Returns the kept rows and the number of dropped rows. SMILES strings are
    cached so that the same molecule is only validated once across many target
    rows, which keeps the cost negligible relative to the ChEMBL API fetches.
    """
    cache: dict[str, bool] = {}
    kept: list[dict] = []
    dropped = 0
    for row in rows:
        smiles = row.get("canonical_smiles") or ""
        if smiles not in cache:
            cache[smiles] = smiles_passes_scaffold_check(smiles)
        if cache[smiles]:
            kept.append(row)
        else:
            dropped += 1
    return kept, dropped


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare Chemprop multi-task training tables from local target predictions and ChEMBL activities."
    )
    parser.add_argument("--key-targets-csv", required=True)
    parser.add_argument("--chembl-target-catalog", required=True)
    parser.add_argument("--idmapping-tsv", required=True)
    parser.add_argument("--query-compound", default=None)
    parser.add_argument("--query-cas", default=None)
    parser.add_argument("--query-name", default=None)
    parser.add_argument("--query-smiles", default=None)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--min-pathway-n", type=int, default=2)
    parser.add_argument("--min-platform-vote", type=int, default=2)
    parser.add_argument("--min-task-samples", type=int, default=100)
    parser.add_argument("--sleep-seconds", type=float, default=0.2)
    parser.add_argument("--targets", nargs="*", default=None)
    parser.add_argument("--max-retries", type=int, default=8)
    parser.add_argument("--request-timeout", type=float, default=60.0)
    return parser.parse_args()


def read_key_targets(path, top_n, min_pathway_n, min_platform_vote):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            row["consensus_score"] = float(row["consensus_score"])
            row["platform_vote"] = int(row["platform_vote"])
            row["pathway_n"] = int(row["pathway_n"])
            if row["platform_vote"] < min_platform_vote:
                continue
            if row["pathway_n"] < min_pathway_n:
                continue
            rows.append(row)
    rows.sort(key=lambda x: (-x["consensus_score"], -x["pathway_n"], -x["platform_vote"], x["target"]))
    return rows[:top_n]


def read_id_mapping(path):
    mapping = defaultdict(set)
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mapping[row["From"]].add(row["To"].upper())
    return mapping


def read_chembl_target_catalog(path, accession_to_gene):
    priority = {
        "SINGLE PROTEIN": 0,
        "PROTEIN COMPLEX": 1,
        "PROTEIN COMPLEX GROUP": 2,
        "PROTEIN FAMILY": 3,
        "CHIMERIC PROTEIN": 4,
    }
    gene_to_candidates = defaultdict(list)
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=";")
        for row in reader:
            if row.get("Organism") != "Homo sapiens":
                continue
            accessions = [x for x in (row.get("Accessions") or "").split("|") if x]
            target_id = row.get("ChEMBL ID")
            if not target_id:
                continue
            target_type = row.get("Type") or ""
            genes = set()
            for accession in accessions:
                genes.update(accession_to_gene.get(accession, set()))
            for gene in genes:
                gene_to_candidates[gene].append(
                    {
                        "target_chembl_id": target_id,
                        "target_type": target_type,
                        "priority": priority.get(target_type, 99),
                    }
                )
    gene_to_target_ids = {}
    for gene, rows in gene_to_candidates.items():
        best_priority = min(x["priority"] for x in rows)
        gene_to_target_ids[gene] = {
            x["target_chembl_id"]
            for x in rows
            if x["priority"] == best_priority
        }
    return gene_to_target_ids


def log(message):
    print(message, flush=True)


def write_progress(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def fetch_json(url, max_retries=8, request_timeout=60.0):
    last_error = None
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=request_timeout) as response:
                return json.load(response)
        except Exception as exc:
            last_error = exc
            log(f"[retry] {attempt + 1}/{max_retries} failed for {url}: {exc}")
            time.sleep(min(10, 1.5 * (attempt + 1)))
    raise last_error


def fetch_activities_for_target(target_chembl_id, sleep_seconds, max_retries, request_timeout):
    rows = []
    for standard_type in ACTIVITY_TYPES:
        offset = 0
        page_i = 0
        while True:
            page_i += 1
            query = urllib.parse.urlencode(
                {
                    "target_chembl_id": target_chembl_id,
                    "standard_relation": "=",
                    "standard_type": standard_type,
                    "limit": 1000,
                    "offset": offset,
                }
            )
            url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?{query}"
            log(f"    [activity-page] target={target_chembl_id} type={standard_type} page={page_i} offset={offset}")
            try:
                payload = fetch_json(url, max_retries=max_retries, request_timeout=request_timeout)
            except Exception as exc:
                log(
                    f"    [skip] target={target_chembl_id} type={standard_type} "
                    f"page={page_i} offset={offset}: {exc}"
                )
                break
            batch = payload.get("activities", [])
            if not batch:
                break
            for row in batch:
                pchembl = row.get("pchembl_value")
                smiles = row.get("canonical_smiles")
                units = row.get("standard_units")
                relation = row.get("standard_relation")
                target_organism = row.get("target_organism")
                if pchembl in (None, "") or smiles in (None, ""):
                    continue
                if units != "nM" or relation != "=":
                    continue
                if target_organism != "Homo sapiens":
                    continue
                rows.append(
                    {
                        "molecule_chembl_id": row.get("molecule_chembl_id"),
                        "canonical_smiles": smiles,
                        "target_chembl_id": row.get("target_chembl_id"),
                        "standard_type": row.get("standard_type"),
                        "standard_value": row.get("standard_value"),
                        "standard_units": units,
                        "pchembl_value": float(pchembl),
                        "document_chembl_id": row.get("document_chembl_id"),
                        "assay_chembl_id": row.get("assay_chembl_id"),
                    }
                )
            page_meta = payload.get("page_meta", {})
            if not page_meta.get("next"):
                break
            offset += page_meta.get("limit", 1000)
            time.sleep(sleep_seconds)
    return rows


def median(values):
    return statistics.median(values) if values else None


def aggregate_long_rows(gene, activity_rows):
    grouped = defaultdict(list)
    types = defaultdict(set)
    metadata = defaultdict(set)
    for row in activity_rows:
        key = (
            row["molecule_chembl_id"],
            row["canonical_smiles"],
            gene,
        )
        grouped[key].append(row["pchembl_value"])
        types[key].add(row["standard_type"])
        metadata[key].add(row["target_chembl_id"])

    aggregated = []
    for key, values in grouped.items():
        molecule_chembl_id, canonical_smiles, gene_symbol = key
        aggregated.append(
            {
                "molecule_chembl_id": molecule_chembl_id,
                "canonical_smiles": canonical_smiles,
                "target_gene": gene_symbol,
                "pchembl_value": median(values),
                "n_measurements": len(values),
                "standard_types": "|".join(sorted(types[key])),
                "target_chembl_ids": "|".join(sorted(metadata[key])),
            }
        )
    return aggregated


def build_wide_table(long_rows, selected_genes):
    rows_by_molecule = {}
    for row in long_rows:
        molecule_id = row["molecule_chembl_id"]
        if molecule_id not in rows_by_molecule:
            rows_by_molecule[molecule_id] = {
                "molecule_chembl_id": molecule_id,
                "canonical_smiles": row["canonical_smiles"],
            }
        rows_by_molecule[molecule_id][row["target_gene"]] = row["pchembl_value"]

    wide_rows = []
    for molecule_id, row in rows_by_molecule.items():
        for gene in selected_genes:
            row.setdefault(gene, "")
        wide_rows.append(row)
    wide_rows.sort(key=lambda x: x["molecule_chembl_id"])
    return wide_rows


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_task_summary(selected_targets, gene_to_target_ids, long_rows):
    rows_by_gene = defaultdict(list)
    for row in long_rows:
        rows_by_gene[row["target_gene"]].append(row)

    summary = []
    for row in selected_targets:
        gene = row["target"]
        gene_rows = rows_by_gene.get(gene, [])
        summary.append(
            {
                "target": gene,
                "consensus_score": row["consensus_score"],
                "platform_vote": row["platform_vote"],
                "pathway_n": row["pathway_n"],
                "sources": row["sources"],
                "target_chembl_ids": "|".join(sorted(gene_to_target_ids.get(gene, set()))),
                "n_unique_molecules": len({x["molecule_chembl_id"] for x in gene_rows}),
                "n_aggregated_rows": len(gene_rows),
            }
        )
    return summary


def build_single_inference_template(compound, cas, name, canonical_smiles, selected_genes):
    if not canonical_smiles:
        return []
    out = {
        "compound": compound or name or cas or "",
        "cas": cas or "",
        "name": name or compound or "",
        "canonical_smiles": canonical_smiles,
    }
    for gene in selected_genes:
        out[gene] = ""
    return [out]


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    key_targets = read_key_targets(
        args.key_targets_csv,
        top_n=args.top_n,
        min_pathway_n=args.min_pathway_n,
        min_platform_vote=args.min_platform_vote,
    )
    if args.targets:
        selected = {x.upper() for x in args.targets}
        key_targets = [row for row in key_targets if row["target"].upper() in selected]

    accession_to_gene = read_id_mapping(args.idmapping_tsv)
    gene_to_target_ids = read_chembl_target_catalog(args.chembl_target_catalog, accession_to_gene)
    progress_path = output_dir / "progress.json"

    selected_genes = [row["target"] for row in key_targets]
    log(f"[selected-targets] {', '.join(selected_genes)}")
    target_mapping_rows = []
    long_rows = []

    for gene_index, row in enumerate(key_targets, start=1):
        gene = row["target"]
        chembl_ids = sorted(gene_to_target_ids.get(gene, set()))
        log(f"[target {gene_index}/{len(key_targets)}] {gene}: {len(chembl_ids)} target_chembl_ids")
        target_mapping_rows.append(
            {
                "target": gene,
                "consensus_score": row["consensus_score"],
                "platform_vote": row["platform_vote"],
                "pathway_n": row["pathway_n"],
                "sources": row["sources"],
                "target_chembl_ids": "|".join(chembl_ids),
                "n_target_chembl_ids": len(chembl_ids),
            }
        )
        write_progress(
            progress_path,
            {
                "stage": "fetching_activities",
                "current_target": gene,
                "target_index": gene_index,
                "target_total": len(key_targets),
                "target_chembl_ids": chembl_ids,
                "long_rows_so_far": len(long_rows),
            },
        )
        for chembl_index, target_chembl_id in enumerate(chembl_ids, start=1):
            log(f"  [fetch {chembl_index}/{len(chembl_ids)}] {gene} <- {target_chembl_id}")
            raw_rows = fetch_activities_for_target(
                target_chembl_id,
                args.sleep_seconds,
                max_retries=args.max_retries,
                request_timeout=args.request_timeout,
            )
            log(f"  [done] {gene} <- {target_chembl_id}: {len(raw_rows)} filtered activity rows")
            long_rows.extend(aggregate_long_rows(gene, raw_rows))
            write_progress(
                progress_path,
                {
                    "stage": "fetching_activities",
                    "current_target": gene,
                    "current_target_chembl_id": target_chembl_id,
                    "target_index": gene_index,
                    "target_total": len(key_targets),
                    "target_chembl_ids": chembl_ids,
                    "long_rows_so_far": len(long_rows),
                },
            )

    merged_by_gene_molecule = defaultdict(list)
    details = {}
    for row in long_rows:
        key = (row["target_gene"], row["molecule_chembl_id"], row["canonical_smiles"])
        merged_by_gene_molecule[key].append(row["pchembl_value"])
        details[key] = row

    deduped_long_rows = []
    for key, values in merged_by_gene_molecule.items():
        row = dict(details[key])
        row["pchembl_value"] = median(values)
        row["n_measurements"] = max(row["n_measurements"], len(values))
        deduped_long_rows.append(row)

    # Drop rows whose SMILES would crash RDKit / astartes / chemprop downstream.
    n_before_sanitize = len(deduped_long_rows)
    deduped_long_rows, n_smiles_dropped = sanitize_smiles_rows(deduped_long_rows)
    log(
        f"[smiles-sanitize] dropped {n_smiles_dropped}/{n_before_sanitize} rows "
        "with malformed SMILES (failed RDKit MurckoScaffold check)"
    )

    task_summary = build_task_summary(key_targets, gene_to_target_ids, deduped_long_rows)
    eligible_genes = [
        row["target"]
        for row in task_summary
        if int(row["n_unique_molecules"]) >= args.min_task_samples
    ]

    chemprop_wide = build_wide_table(
        [row for row in deduped_long_rows if row["target_gene"] in eligible_genes],
        eligible_genes,
    )
    inference_template = build_single_inference_template(
        args.query_compound,
        args.query_cas,
        args.query_name,
        args.query_smiles,
        eligible_genes,
    )

    write_csv(
        output_dir / "target_mapping.csv",
        target_mapping_rows,
        [
            "target",
            "consensus_score",
            "platform_vote",
            "pathway_n",
            "sources",
            "target_chembl_ids",
            "n_target_chembl_ids",
        ],
    )
    write_csv(
        output_dir / "chemprop_long.csv",
        deduped_long_rows,
        [
            "molecule_chembl_id",
            "canonical_smiles",
            "target_gene",
            "pchembl_value",
            "n_measurements",
            "standard_types",
            "target_chembl_ids",
        ],
    )
    write_csv(
        output_dir / "task_summary.csv",
        task_summary,
        [
            "target",
            "consensus_score",
            "platform_vote",
            "pathway_n",
            "sources",
            "target_chembl_ids",
            "n_unique_molecules",
            "n_aggregated_rows",
        ],
    )
    if chemprop_wide:
        write_csv(
            output_dir / "chemprop_multitask.csv",
            chemprop_wide,
            ["molecule_chembl_id", "canonical_smiles"] + eligible_genes,
        )
    if inference_template:
        write_csv(
            output_dir / "inference_template.csv",
            inference_template,
            ["compound", "cas", "name", "canonical_smiles"] + eligible_genes,
        )

    summary = {
        "selected_targets": selected_genes,
        "eligible_targets": eligible_genes,
        "n_long_rows": len(deduped_long_rows),
        "n_wide_rows": len(chemprop_wide),
        "n_smiles_dropped": n_smiles_dropped,
        "min_task_samples": args.min_task_samples,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    write_progress(
        progress_path,
        {
            "stage": "completed",
            **summary,
        },
    )

    log(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
