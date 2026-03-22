#!/usr/bin/env python3

import argparse
import csv
import json
import shutil
import subprocess
from pathlib import Path


def copy_if_needed(src: Path, dst: Path) -> Path:
    if src.resolve() == dst.resolve():
        return dst
    shutil.copy2(src, dst)
    return dst


def parse_args():
    parser = argparse.ArgumentParser(
        description="Assemble a case directory from raw Swiss/SEA/ChEMBL/PPB2 inputs."
    )
    parser.add_argument("--case-dir", required=True)
    parser.add_argument("--idmapping", required=True)
    parser.add_argument("--swiss", required=True)
    parser.add_argument("--sea", required=True)
    parser.add_argument("--chembl", required=True)
    parser.add_argument("--ppb2", default=None)
    return parser.parse_args()


def upper_values_from_csv(path, column):
    values = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            value = (row.get(column) or "").strip().upper()
            if value:
                values.append(value)
    return sorted(set(values))


def load_idmapping(path):
    mapping = {}
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            from_id = (row.get("From") or "").strip()
            to_id = (row.get("To") or "").strip().upper()
            if not from_id or not to_id:
                continue
            mapping.setdefault(from_id, set()).add(to_id)
    return mapping


def upper_values_from_chembl_csv(path, accession_to_gene):
    values = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if (row.get("Organism") or "").strip() != "Homo sapiens":
                continue
            accessions = [x.strip() for x in (row.get("Accessions") or "").split("|") if x.strip()]
            for accession in accessions:
                values.extend(accession_to_gene.get(accession, set()))
    return sorted(set(values))


def upper_values_from_sea(path):
    values = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            target_id = (row.get("Target ID") or "").strip()
            if not target_id.endswith("_HUMAN"):
                continue
            name = (row.get("Name") or "").strip().upper()
            if name:
                values.append(name)
    return sorted(set(values))


def write_venn_rds(case_dir, venn_dict):
    venn_json = case_dir / "venn_inputs.json"
    with venn_json.open("w", encoding="utf-8") as f:
        json.dump(venn_dict, f, ensure_ascii=False, indent=2)

    r_code = (
        "library(jsonlite); "
        f"x <- fromJSON('{venn_json.as_posix()}', simplifyVector = FALSE); "
        "x <- lapply(x, unlist); "
        f"saveRDS(x, file='{(case_dir / 'venn.rds').as_posix()}')"
    )
    subprocess.run(["Rscript", "-e", r_code], check=True)


def main():
    args = parse_args()
    case_dir = Path(args.case_dir)
    case_dir.mkdir(parents=True, exist_ok=True)

    idmapping_src = Path(args.idmapping)
    swiss_src = Path(args.swiss)
    sea_src = Path(args.sea)
    chembl_src = Path(args.chembl)
    ppb2_src = Path(args.ppb2) if args.ppb2 else None

    idmapping_dst = case_dir / idmapping_src.name
    swiss_dst = case_dir / swiss_src.name
    sea_dst = case_dir / sea_src.name
    chembl_dst = case_dir / chembl_src.name
    copy_if_needed(idmapping_src, idmapping_dst)
    copy_if_needed(swiss_src, swiss_dst)
    copy_if_needed(sea_src, sea_dst)
    copy_if_needed(chembl_src, chembl_dst)
    if ppb2_src and ppb2_src.exists():
        copy_if_needed(ppb2_src, case_dir / ppb2_src.name)

    accession_to_gene = load_idmapping(idmapping_dst)
    venn = {
        "Swiss": upper_values_from_csv(swiss_dst, "Common name"),
        "SEA": upper_values_from_sea(sea_dst),
        "ChEMBL": upper_values_from_chembl_csv(chembl_dst, accession_to_gene),
        "PPB2": upper_values_from_csv(ppb2_src, "Symbol") if ppb2_src and ppb2_src.exists() else [],
    }
    write_venn_rds(case_dir, venn)

    with (case_dir / "source_summary.json").open("w", encoding="utf-8") as f:
        json.dump(
            {
                "idmapping_file": idmapping_dst.name,
                "swiss_file": swiss_dst.name,
                "sea_file": sea_dst.name,
                "chembl_file": chembl_dst.name,
                "ppb2_file": ppb2_src.name if ppb2_src and ppb2_src.exists() else None,
                "counts": {k: len(v) for k, v in venn.items()},
            },
            f,
            ensure_ascii=False,
            indent=2,
        )


if __name__ == "__main__":
    main()
