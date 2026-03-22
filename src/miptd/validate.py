from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


REQUIRED_DIRS = (
    "a_source_collection",
    "b_venn",
    "c_go_kegg",
    "d_kegg_circlize",
    "e_chemprop",
)

REQUIRED_FILES = (
    "case_manifest.json",
    "a_source_collection/source_summary.json",
    "b_venn/1.p_venn.pdf",
    "c_go_kegg/2.GO_KEGG.pdf",
    "d_kegg_circlize/3.KEGG_circos.pdf",
    "e_chemprop/chemprop_data/inference_template.csv",
    "e_chemprop/chemprop_model/analysis/prediction_summary.json",
    "e_chemprop/Figure1e_filtered.png",
    "e_chemprop/Figure1e_filtered.pdf",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="MIPTD-validate",
        description="Validate a single-CAS MIPTD case directory for structure and single-compound output consistency.",
    )
    parser.add_argument("--case-dir", required=True, help="Path to a CAS_xxxxx_Date case directory.")
    return parser.parse_args()


def read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def fail(message: str) -> None:
    raise SystemExit(message)


def validate_case(case_dir: Path) -> None:
    if not case_dir.exists():
        fail(f"Case directory does not exist: {case_dir}")

    for dirname in REQUIRED_DIRS:
        path = case_dir / dirname
        if not path.is_dir():
            fail(f"Missing required directory: {path}")

    for relpath in REQUIRED_FILES:
        path = case_dir / relpath
        if not path.exists():
            fail(f"Missing required file: {path}")

    manifest = read_json(case_dir / "case_manifest.json")
    compound = manifest.get("compound", {})
    manifest_cas = str(compound.get("cas", "")).strip()
    if not manifest_cas:
        fail("case_manifest.json does not contain compound.cas")

    inference_rows = read_csv(case_dir / "e_chemprop/chemprop_data/inference_template.csv")
    if len(inference_rows) != 1:
        fail(f"inference_template.csv must contain exactly 1 row, found {len(inference_rows)}")
    if inference_rows[0].get("cas", "").strip() != manifest_cas:
        fail("inference_template.csv CAS does not match case_manifest.json")

    summary = read_json(case_dir / "e_chemprop/chemprop_model/analysis/prediction_summary.json")
    compounds = summary.get("compounds", [])
    if len(compounds) != 1:
        fail(f"prediction_summary.json must contain exactly 1 compound, found {len(compounds)}")
    if str(compounds[0].get("cas", "")).strip() != manifest_cas:
        fail("prediction_summary.json CAS does not match case_manifest.json")

    print(json.dumps(
        {
            "status": "ok",
            "case_dir": str(case_dir),
            "cas": manifest_cas,
            "compound": compound.get("name") or compound.get("compound"),
            "validated_dirs": list(REQUIRED_DIRS),
            "validated_files": list(REQUIRED_FILES),
        },
        ensure_ascii=False,
        indent=2,
    ))


def main() -> None:
    args = parse_args()
    validate_case(Path(args.case_dir).resolve())


if __name__ == "__main__":
    main()
