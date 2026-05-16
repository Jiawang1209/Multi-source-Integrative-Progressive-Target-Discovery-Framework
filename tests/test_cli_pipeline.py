from __future__ import annotations

import json
import io
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from miptd.cli import main as cli_main
from miptd.models import CompoundRecord
from miptd.pipeline import PipelineRunner, run_single_case


class CliAndPipelineTests(unittest.TestCase):
    def test_cli_accepts_disase_keywords_alias(self) -> None:
        with patch("miptd.cli.run_single_case") as mock_run:
            mock_run.return_value = type("CasePathsStub", (), {"case_root": Path("/tmp/case")})()
            with patch.object(
                sys,
                "argv",
                [
                    "MIPTD",
                    "--cas",
                    "491-71-4",
                    "--disase-keywords",
                    "NAFLD,liver",
                ],
            ):
                stdout = io.StringIO()
                with redirect_stdout(stdout):
                    cli_main()

            mock_run.assert_called_once()
            kwargs = mock_run.call_args.kwargs
            self.assertEqual(kwargs["cas"], "491-71-4")
            self.assertEqual(kwargs["disease_keywords"], ["NAFLD", "liver"])
            self.assertEqual(stdout.getvalue().strip(), "/tmp/case")

    def test_run_single_case_dry_run_creates_case_layout_and_log(self) -> None:
        fake_compound = CompoundRecord(
            compound="Chrysoeriol",
            cas="491-71-4",
            name="Chrysoeriol",
            formula="C16H12O6",
            canonical_smiles="COC1=CC(=CC(=C1O)O)C2=CC(=O)C3=C(O2)C=C(C=C3)O",
            inchikey="AABBCCDDEEFFGG-UHFFFAOYSA-N",
            identity_source="pubchem",
            pubchem_cid="5280666",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_root = Path(tmpdir)
            with patch("miptd.pipeline.load_compound_record", return_value=fake_compound):
                case_paths = run_single_case(
                    cas="491-71-4",
                    disease_keywords=["NAFLD", "liver"],
                    output_root=output_root,
                    project_root=REPO_ROOT,
                    run_date="2026-03-22",
                    dry_run=True,
                )

            self.assertTrue(case_paths.case_root.is_dir())
            self.assertTrue(case_paths.a_dir.is_dir())
            self.assertTrue(case_paths.b_dir.is_dir())
            self.assertTrue(case_paths.c_dir.is_dir())
            self.assertTrue(case_paths.d_dir.is_dir())
            self.assertTrue(case_paths.e_dir.is_dir())

            manifest_path = case_paths.case_root / "case_manifest.json"
            run_log_path = case_paths.case_root / "run.log"
            self.assertTrue(manifest_path.exists())
            self.assertTrue(run_log_path.exists())

            log_text = run_log_path.read_text(encoding="utf-8")
            self.assertIn("[MIPTD] Step A:", log_text)
            self.assertIn("[MIPTD] Step E.5:", log_text)
            self.assertIn("DRY-RUN", log_text)

    def test_step_a_can_continue_with_single_source_failure(self) -> None:
        fake_compound = CompoundRecord(
            compound="Chrysoeriol",
            cas="491-71-4",
            name="Chrysoeriol",
            formula="C16H12O6",
            canonical_smiles="COC1=CC(=CC(=C1O)O)C2=CC(=O)C3=C(O2)C=C(C=C3)O",
            inchikey="AABBCCDDEEFFGG-UHFFFAOYSA-N",
            identity_source="pubchem",
            pubchem_cid="5280666",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_root = Path(tmpdir)
            with patch("miptd.pipeline.load_compound_record", return_value=fake_compound):
                runner = PipelineRunner(
                    cas="491-71-4",
                    disease_keywords=["NAFLD", "liver"],
                    output_root=output_root,
                    project_root=REPO_ROOT,
                    run_date="2026-03-22",
                    dry_run=False,
                )

            def fake_run(cmd: list[str]) -> None:
                joined = " ".join(cmd)
                if "fetch_swiss_targets.py" in joined:
                    raise RuntimeError("Swiss timeout")
                if "fetch_sea_targets.py" in joined:
                    (runner.case_paths.a_dir / "sea_fetch" / "sea-results.xls").write_text(
                        "Target ID,Name\nABC_HUMAN,AKT1\n",
                        encoding="utf-8",
                    )
                    return
                if "fetch_chembl_targets.py" in joined:
                    (runner.case_paths.a_dir / "chembl_fetch" / "chembl_targets.csv").write_text(
                        "ChEMBL ID,Name,Accessions,Type,Organism,Compounds,Activities,Tax ID,Species Group Flag\n"
                        "CHEMBL1,Target,P31749,SINGLE PROTEIN,Homo sapiens,,,,\n",
                        encoding="utf-8",
                    )
                    return
                if "fetch_ppb2_targets.py" in joined:
                    (runner.case_paths.a_dir / "ppb2_fetch" / "ppb2_targets.csv").write_text(
                        "Symbol\nPPARG\n",
                        encoding="utf-8",
                    )
                    return
                if "build_case_sources.py" in joined:
                    (runner.case_paths.a_dir / "source_summary.json").write_text(
                        json.dumps(
                            {
                                "counts": {
                                    "Swiss": 0,
                                    "SEA": 1,
                                    "ChEMBL": 1,
                                    "PPB2": 1,
                                }
                            }
                        ),
                        encoding="utf-8",
                    )
                    return
                raise AssertionError(f"Unexpected command: {joined}")

            runner._run = fake_run  # type: ignore[method-assign]
            with patch("miptd.pipeline.SOURCE_FETCH_RETRY_DELAYS_SECONDS", (0, 0)):
                with patch("miptd.pipeline.time.sleep"):
                    runner._step_a_source_collection()

            swiss_placeholder = runner.case_paths.a_dir / "swiss_fetch" / "SwissTargetPrediction_491-71-4.csv"
            self.assertTrue(swiss_placeholder.exists())
            self.assertEqual(
                json.loads(runner.status_path.read_text(encoding="utf-8"))["sources"]["Swiss"]["status"],
                "failed",
            )
            self.assertEqual(runner.status["steps"]["step_a"]["status"], "degraded")
            self.assertEqual(runner.status["non_empty_sources"], ["ChEMBL", "PPB2", "SEA"])

    def test_source_fetch_retries_before_marking_ok(self) -> None:
        fake_compound = CompoundRecord(
            compound="Chrysoeriol",
            cas="491-71-4",
            name="Chrysoeriol",
            formula="C16H12O6",
            canonical_smiles="COC1=CC(=CC(=C1O)O)C2=CC(=O)C3=C(O2)C=C(C=C3)O",
            inchikey="AABBCCDDEEFFGG-UHFFFAOYSA-N",
            identity_source="pubchem",
            pubchem_cid="5280666",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_root = Path(tmpdir)
            with patch("miptd.pipeline.load_compound_record", return_value=fake_compound):
                runner = PipelineRunner(
                    cas="491-71-4",
                    disease_keywords=["NAFLD", "liver"],
                    output_root=output_root,
                    project_root=REPO_ROOT,
                    run_date="2026-03-22",
                    dry_run=False,
                )

            output_file = runner.case_paths.a_dir / "sea_fetch" / "sea-results.xls"
            attempts = {"count": 0}

            def flaky_run(cmd: list[str]) -> None:
                attempts["count"] += 1
                if attempts["count"] < 3:
                    raise RuntimeError(f"SEA transient error {attempts['count']}")
                output_file.parent.mkdir(parents=True, exist_ok=True)
                output_file.write_text("Target ID,Name\nABC_HUMAN,AKT1\n", encoding="utf-8")

            runner._run = flaky_run  # type: ignore[method-assign]
            with patch("miptd.pipeline.SOURCE_FETCH_RETRY_DELAYS_SECONDS", (0, 0)):
                with patch("miptd.pipeline.time.sleep") as mock_sleep:
                    runner._fetch_source_with_fallback(
                        "SEA",
                        output_file,
                        ["python3", "fetch_sea_targets.py"],
                        ["Target ID", "Name"],
                    )

            self.assertEqual(attempts["count"], 3)
            self.assertTrue(output_file.exists())
            source_status = json.loads(runner.status_path.read_text(encoding="utf-8"))["sources"]["SEA"]
            self.assertEqual(source_status["status"], "ok")
            self.assertFalse(source_status["used_placeholder"])
            self.assertEqual(source_status["attempt_count"], 3)
            self.assertEqual(len(source_status["attempt_errors"]), 2)
            self.assertEqual(mock_sleep.call_count, 2)

    def test_source_fetch_writes_placeholder_after_retry_exhaustion(self) -> None:
        fake_compound = CompoundRecord(
            compound="Chrysoeriol",
            cas="491-71-4",
            name="Chrysoeriol",
            formula="C16H12O6",
            canonical_smiles="COC1=CC(=CC(=C1O)O)C2=CC(=O)C3=C(O2)C=C(C=C3)O",
            inchikey="AABBCCDDEEFFGG-UHFFFAOYSA-N",
            identity_source="pubchem",
            pubchem_cid="5280666",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_root = Path(tmpdir)
            with patch("miptd.pipeline.load_compound_record", return_value=fake_compound):
                runner = PipelineRunner(
                    cas="491-71-4",
                    disease_keywords=["NAFLD", "liver"],
                    output_root=output_root,
                    project_root=REPO_ROOT,
                    run_date="2026-03-22",
                    dry_run=False,
                )

            output_file = runner.case_paths.a_dir / "ppb2_fetch" / "ppb2_targets.csv"
            attempts = {"count": 0}

            def failing_run(cmd: list[str]) -> None:
                attempts["count"] += 1
                raise RuntimeError(f"PPB2 transient error {attempts['count']}")

            runner._run = failing_run  # type: ignore[method-assign]
            with patch("miptd.pipeline.SOURCE_FETCH_RETRY_DELAYS_SECONDS", (0, 0)):
                with patch("miptd.pipeline.time.sleep"):
                    runner._fetch_source_with_fallback(
                        "PPB2",
                        output_file,
                        ["python3", "fetch_ppb2_targets.py"],
                        ["Symbol"],
                    )

            self.assertEqual(attempts["count"], 3)
            self.assertEqual(output_file.read_text(encoding="utf-8"), "Symbol\n")
            source_status = json.loads(runner.status_path.read_text(encoding="utf-8"))["sources"]["PPB2"]
            self.assertEqual(source_status["status"], "failed")
            self.assertTrue(source_status["used_placeholder"])
            self.assertEqual(source_status["attempt_count"], 3)
            self.assertEqual(len(source_status["attempt_errors"]), 3)


if __name__ == "__main__":
    unittest.main()
