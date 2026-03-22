from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from miptd.validate import validate_case


def build_minimal_case(case_dir: Path, cas: str = "491-71-4") -> Path:
    required_dirs = [
        "a_source_collection",
        "b_venn",
        "c_go_kegg",
        "d_kegg_circlize",
        "e_chemprop/chemprop_data",
        "e_chemprop/chemprop_model/analysis",
    ]
    for rel in required_dirs:
        (case_dir / rel).mkdir(parents=True, exist_ok=True)

    (case_dir / "case_manifest.json").write_text(
        json.dumps({"compound": {"cas": cas, "name": "Chrysoeriol"}}, ensure_ascii=False),
        encoding="utf-8",
    )
    (case_dir / "a_source_collection/source_summary.json").write_text("{}", encoding="utf-8")
    (case_dir / "b_venn/1.p_venn.pdf").write_bytes(b"%PDF-1.4\n")
    (case_dir / "c_go_kegg/2.GO_KEGG.pdf").write_bytes(b"%PDF-1.4\n")
    (case_dir / "d_kegg_circlize/3.KEGG_circos.pdf").write_bytes(b"%PDF-1.4\n")
    (case_dir / "e_chemprop/Figure1e_filtered.png").write_bytes(b"\x89PNG\r\n\x1a\n")
    (case_dir / "e_chemprop/Figure1e_filtered.pdf").write_bytes(b"%PDF-1.4\n")
    (case_dir / "e_chemprop/chemprop_data/inference_template.csv").write_text(
        "cas,smiles\n" f"{cas},CCO\n",
        encoding="utf-8",
    )
    (case_dir / "e_chemprop/chemprop_model/analysis/prediction_summary.json").write_text(
        json.dumps({"compounds": [{"cas": cas, "name": "Chrysoeriol"}]}, ensure_ascii=False),
        encoding="utf-8",
    )
    return case_dir


class ValidateTests(unittest.TestCase):
    def test_validate_case_succeeds_for_minimal_valid_case(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            case_dir = build_minimal_case(Path(tmpdir) / "CAS_491-71-4_2026-03-22")
            validate_case(case_dir)

    def test_validate_case_fails_when_inference_has_multiple_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            case_dir = build_minimal_case(Path(tmpdir) / "CAS_491-71-4_2026-03-22")
            (case_dir / "e_chemprop/chemprop_data/inference_template.csv").write_text(
                "cas,smiles\n491-71-4,CCO\n491-71-4,CCC\n",
                encoding="utf-8",
            )
            with self.assertRaises(SystemExit) as ctx:
                validate_case(case_dir)
            self.assertIn("exactly 1 row", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
