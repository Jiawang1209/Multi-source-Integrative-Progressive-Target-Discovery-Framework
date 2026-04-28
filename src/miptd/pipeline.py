from __future__ import annotations

import csv
import json
import subprocess
from dataclasses import asdict
from pathlib import Path

from .config import (
    DEFAULT_BATCH_SIZE,
    DEFAULT_EPOCHS,
    DEFAULT_MIN_PATHWAY_N,
    DEFAULT_MIN_PLATFORM_VOTE,
    DEFAULT_MIN_TASK_SAMPLES,
    DEFAULT_MIN_TEST_MOLECULES,
    DEFAULT_MIN_TEST_R2,
    DEFAULT_MIN_TRAINING_MOLECULES,
    DEFAULT_PATIENCE,
    DEFAULT_TOP_N_KEGG_TARGETS,
    project_root_from_file,
)
from .discovery import (
    create_case_paths,
    discover_chembl_catalog,
    discover_latest_idmapping,
    load_compound_record,
)
from .models import CasePaths
from .utils import copy_file, ensure_dir, join_keywords_regex, print_stage, python_script_cmd, read_csv_rows, r_script_cmd, run_command, write_text


class PipelineRunner:
    def __init__(
        self,
        cas: str,
        disease_keywords: list[str],
        output_root: Path,
        project_root: Path | None = None,
        run_date: str | None = None,
        chembl_catalog: Path | None = None,
        idmapping_tsv: Path | None = None,
        dry_run: bool = False,
    ) -> None:
        self.project_root = project_root or project_root_from_file()
        self.compound = load_compound_record(cas)
        self.chembl_catalog = chembl_catalog or discover_chembl_catalog(self.project_root)
        self.idmapping_tsv = idmapping_tsv or discover_latest_idmapping(self.project_root)
        self.case_paths = create_case_paths(output_root, cas, run_date)
        self.disease_keywords = disease_keywords
        self.disease_regex = join_keywords_regex(disease_keywords)
        self.dry_run = dry_run
        self.scripts_dir = self.project_root / "scripts"
        self.chemprop_python = self._detect_chemprop_python()
        self.run_log_path = self.case_paths.case_root / "run.log"
        self.status_path = self.case_paths.case_root / "status.json"
        self.status = {
            "pipeline_status": "initialized",
            "case_root": str(self.case_paths.case_root),
            "compound": asdict(self.compound),
            "steps": {},
            "sources": {},
        }

    def run(self) -> CasePaths:
        print_stage("Initialize", f"Resolve input CAS and prepare single-case workspace for {self.compound.cas} ({self.compound.name})", log_path=self.run_log_path)
        print_stage("Output", f"Write all results under {self.case_paths.case_root}", log_path=self.run_log_path)
        self._write_case_manifest()
        self._write_compound_identity()
        self._write_status()
        try:
            self._step_a_source_collection()
            self._step_bcd()
            self._step_e()
        except Exception as exc:
            self.status["pipeline_status"] = "failed"
            self.status["error"] = str(exc)
            self._write_status()
            raise
        self.status["pipeline_status"] = "success"
        self._write_status()
        print_stage("Completed", f"Single-CAS analysis finished: {self.case_paths.case_root}", log_path=self.run_log_path)
        return self.case_paths

    def _write_case_manifest(self) -> None:
        manifest = {
            "compound": asdict(self.compound),
            "case_root": str(self.case_paths.case_root),
            "disease_keywords": self.disease_keywords,
            "disease_regex": self.disease_regex,
            "chembl_catalog": str(self.chembl_catalog),
            "idmapping_tsv": str(self.idmapping_tsv),
        }
        write_text(self.case_paths.case_root / "case_manifest.json", json.dumps(manifest, ensure_ascii=False, indent=2))

    def _write_compound_identity(self) -> None:
        identity_payload = {
            "compound": self.compound.compound,
            "cas": self.compound.cas,
            "name": self.compound.name,
            "formula": self.compound.formula,
            "canonical_smiles": self.compound.canonical_smiles,
            "inchikey": self.compound.inchikey,
            "identity_source": self.compound.identity_source,
            "pubchem_cid": self.compound.pubchem_cid,
        }
        identity_json = json.dumps(identity_payload, ensure_ascii=False, indent=2)
        csv_line = ",".join([
            "compound",
            "cas",
            "name",
            "formula",
            "canonical_smiles",
            "inchikey",
            "identity_source",
            "pubchem_cid",
        ]) + "\n" + ",".join([
            self._csv_escape(identity_payload["compound"]),
            self._csv_escape(identity_payload["cas"]),
            self._csv_escape(identity_payload["name"]),
            self._csv_escape(identity_payload["formula"]),
            self._csv_escape(identity_payload["canonical_smiles"]),
            self._csv_escape(identity_payload["inchikey"]),
            self._csv_escape(identity_payload["identity_source"]),
            self._csv_escape(identity_payload["pubchem_cid"]),
        ]) + "\n"
        for base_dir in (self.case_paths.case_root, self.case_paths.a_dir):
            write_text(base_dir / "compound_identity.json", identity_json)
            write_text(base_dir / "compound_identity.csv", csv_line)

    @staticmethod
    def _csv_escape(value: str) -> str:
        text = value or ""
        if any(ch in text for ch in [",", "\"", "\n"]):
            return "\"" + text.replace("\"", "\"\"") + "\""
        return text

    def _step_a_source_collection(self) -> None:
        self._update_step_status("step_a", "running")
        print_stage("Step A", "Resolve compound identity and collect human targets from Swiss / SEA / ChEMBL / PPB2", log_path=self.run_log_path)
        swiss_fetch_dir = ensure_dir(self.case_paths.a_dir / "swiss_fetch")
        sea_fetch_dir = ensure_dir(self.case_paths.a_dir / "sea_fetch")
        swiss_file = swiss_fetch_dir / f"SwissTargetPrediction_{self.compound.cas}.csv"
        sea_file = sea_fetch_dir / "sea-results.xls"

        chembl_fetch_dir = ensure_dir(self.case_paths.a_dir / "chembl_fetch")
        ppb2_fetch_dir = ensure_dir(self.case_paths.a_dir / "ppb2_fetch")

        print_stage("Step A.1", "Query SwissTargetPrediction with the resolved SMILES", log_path=self.run_log_path)
        self._fetch_source_with_fallback(
            source_name="Swiss",
            output_file=swiss_file,
            cmd=python_script_cmd(
                self.scripts_dir / "fetch_swiss_targets.py",
                "--cas",
                self.compound.cas,
                "--smiles",
                self.compound.canonical_smiles,
                "--result-timeout-seconds",
                "900",
                "--output-dir",
                str(swiss_fetch_dir),
            ),
            placeholder_columns=["Common name"],
        )

        print_stage("Step A.2", "Query SEA with the resolved SMILES and keep human targets", log_path=self.run_log_path)
        self._fetch_source_with_fallback(
            source_name="SEA",
            output_file=sea_file,
            cmd=python_script_cmd(
                self.scripts_dir / "fetch_sea_targets.py",
                "--cas",
                self.compound.cas,
                "--smiles",
                self.compound.canonical_smiles,
                "--compound-name",
                self.compound.name,
                "--output-dir",
                str(sea_fetch_dir),
            ),
            placeholder_columns=["Target ID", "Name"],
        )

        chembl_targets_csv = chembl_fetch_dir / "chembl_targets.csv"
        print_stage("Step A.3", "Query ChEMBL using compound name, InChIKey, and SMILES", log_path=self.run_log_path)
        self._fetch_source_with_fallback(
            source_name="ChEMBL",
            output_file=chembl_targets_csv,
            cmd=python_script_cmd(
                self.scripts_dir / "fetch_chembl_targets.py",
                "--query",
                self.compound.name,
                "--match-inchikey",
                self.compound.inchikey,
                "--match-smiles",
                self.compound.canonical_smiles,
                "--output-dir",
                str(chembl_fetch_dir),
            ),
            placeholder_columns=[
                "ChEMBL ID",
                "Name",
                "Accessions",
                "Type",
                "Organism",
                "Compounds",
                "Activities",
                "Tax ID",
                "Species Group Flag",
            ],
        )

        ppb2_targets_csv = ppb2_fetch_dir / "ppb2_targets.csv"
        print_stage("Step A.4", "Query PPB2 with the resolved SMILES and map targets to human genes", log_path=self.run_log_path)
        self._fetch_source_with_fallback(
            source_name="PPB2",
            output_file=ppb2_targets_csv,
            cmd=python_script_cmd(
                self.scripts_dir / "fetch_ppb2_targets.py",
                "--smiles",
                self.compound.canonical_smiles,
                "--compound-id",
                self.compound.cas,
                "--output-dir",
                str(ppb2_fetch_dir),
                "--chembl-target-catalog",
                str(self.chembl_catalog),
                "--idmapping-tsv",
                str(self.idmapping_tsv),
            ),
            placeholder_columns=["Symbol"],
        )

        if self.dry_run:
            self._update_step_status("step_a", "dry_run")
            self._run(
                python_script_cmd(
                    self.scripts_dir / "build_case_sources.py",
                    "--case-dir",
                    str(self.case_paths.a_dir),
                    "--idmapping",
                    str(self.case_paths.a_dir / self.idmapping_tsv.name),
                    "--swiss",
                    str(self.case_paths.a_dir / swiss_file.name),
                    "--sea",
                    str(self.case_paths.a_dir / sea_file.name),
                    "--chembl",
                    str(self.case_paths.a_dir / chembl_targets_csv.name),
                    "--ppb2",
                    str(self.case_paths.a_dir / ppb2_targets_csv.name),
                )
            )
            return

        copy_file(self.idmapping_tsv, self.case_paths.a_dir / self.idmapping_tsv.name)
        copy_file(swiss_file, self.case_paths.a_dir / swiss_file.name)
        copy_file(sea_file, self.case_paths.a_dir / sea_file.name)
        copy_file(chembl_targets_csv, self.case_paths.a_dir / chembl_targets_csv.name)
        copy_file(ppb2_targets_csv, self.case_paths.a_dir / ppb2_targets_csv.name)

        print_stage("Step A.5", "Build standardized source tables, source summary, and Venn inputs", log_path=self.run_log_path)
        self._run(
            python_script_cmd(
                self.scripts_dir / "build_case_sources.py",
                "--case-dir",
                str(self.case_paths.a_dir),
                "--idmapping",
                str(self.case_paths.a_dir / self.idmapping_tsv.name),
                "--swiss",
                str(self.case_paths.a_dir / swiss_file.name),
                "--sea",
                str(self.case_paths.a_dir / sea_file.name),
                "--chembl",
                str(self.case_paths.a_dir / chembl_targets_csv.name),
                "--ppb2",
                str(self.case_paths.a_dir / ppb2_targets_csv.name),
            )
        )
        source_summary = json.loads((self.case_paths.a_dir / "source_summary.json").read_text(encoding="utf-8"))
        counts = source_summary.get("counts", {})
        non_empty_sources = sorted([name for name, count in counts.items() if int(count) > 0])
        self.status["source_counts"] = counts
        self.status["non_empty_sources"] = non_empty_sources
        step_status = "ok" if all(item.get("status") == "ok" for item in self.status["sources"].values()) else "degraded"
        self._update_step_status("step_a", step_status, {"counts": counts, "non_empty_sources": non_empty_sources})
        if len(non_empty_sources) < 2:
            raise RuntimeError(
                f"At least 2 source databases with non-empty human targets are required to continue; got {len(non_empty_sources)}: {', '.join(non_empty_sources) or 'none'}."
            )

    def _step_bcd(self) -> None:
        self._update_step_status("step_bcd", "running")
        print_stage("Step B/C/D", "Build the Venn diagram, run GO/KEGG enrichment, filter disease pathways, and generate the KEGG circlize plot", log_path=self.run_log_path)
        bcd_work_dir = ensure_dir(self.case_paths.work_dir / "bcd")
        self._run(
            r_script_cmd(
                self.scripts_dir / "run_figure1_bcd.R",
                "--case-dir",
                str(bcd_work_dir),
                "--swiss",
                str(self.case_paths.a_dir / f"SwissTargetPrediction_{self.compound.cas}.csv"),
                "--sea",
                str(self.case_paths.a_dir / "sea-results.xls"),
                "--chembl",
                str(self.case_paths.a_dir / "chembl_targets.csv"),
                "--ppb2",
                str(self.case_paths.a_dir / "ppb2_targets.csv"),
                "--idmapping",
                str(self.case_paths.a_dir / self.idmapping_tsv.name),
                "--liver-regex",
                self.disease_regex,
            )
        )
        if self.dry_run:
            self._update_step_status("step_bcd", "dry_run")
            return

        for filename in ("venn.rds", "venn_inputs.json", "1.p_venn.pdf"):
            copy_file(bcd_work_dir / filename, self.case_paths.b_dir / filename)
        for filename in ("2.GO_result_2.xlsx", "2.KEGG_result_2.xlsx", "GO_KEGG_plot.xlsx", "2.GO_KEGG.pdf"):
            copy_file(bcd_work_dir / filename, self.case_paths.c_dir / filename)
        for filename in ("key_targets_from_kegg.csv", "3.KEGG_circos.pdf"):
            copy_file(bcd_work_dir / filename, self.case_paths.d_dir / filename)
        self._update_step_status("step_bcd", "ok")

    def _step_e(self) -> None:
        self._update_step_status("step_e", "running")
        print_stage("Step E", "Build and apply the Chemprop multi-task model from D-panel key KEGG genes", log_path=self.run_log_path)
        targets_from_d_csv = self.case_paths.e_dir / f"CAS_{self.compound.cas}_targets_from_d.csv"
        self._run(
            python_script_cmd(
                self.scripts_dir / "build_chemprop_targets_from_kegg.py",
                "--kegg-genes-csv",
                str(self.case_paths.d_dir / "key_targets_from_kegg.csv"),
                "--output-csv",
                str(targets_from_d_csv),
            )
        )

        chemprop_data_dir = ensure_dir(self.case_paths.e_dir / "chemprop_data")
        print_stage("Step E.1", "Build Chemprop candidate targets, training tables, and single-compound inference input", log_path=self.run_log_path)
        self._run(
            python_script_cmd(
                self.scripts_dir / "prepare_chemprop_data.py",
                "--key-targets-csv",
                str(targets_from_d_csv),
                "--chembl-target-catalog",
                str(self.chembl_catalog),
                "--idmapping-tsv",
                str(self.idmapping_tsv),
                "--query-compound",
                self.compound.compound,
                "--query-cas",
                self.compound.cas,
                "--query-name",
                self.compound.name,
                "--query-smiles",
                self.compound.canonical_smiles,
                "--output-dir",
                str(chemprop_data_dir),
                "--top-n",
                str(DEFAULT_TOP_N_KEGG_TARGETS),
                "--min-pathway-n",
                str(DEFAULT_MIN_PATHWAY_N),
                "--min-platform-vote",
                str(DEFAULT_MIN_PLATFORM_VOTE),
                "--min-task-samples",
                str(DEFAULT_MIN_TASK_SAMPLES),
                python_bin=self.chemprop_python,
            )
        )
        self._update_step_status("step_e", "ok")

        model_dir = ensure_dir(self.case_paths.e_dir / "chemprop_model")
        print_stage("Step E.2", "Train the Chemprop multi-task model on human ChEMBL activity data", log_path=self.run_log_path)
        self._run(
            python_script_cmd(
                self.scripts_dir / "train_chemprop_multitask.py",
                "--data-csv",
                str(chemprop_data_dir / "chemprop_multitask.csv"),
                "--inference-csv",
                str(chemprop_data_dir / "inference_template.csv"),
                "--output-dir",
                str(model_dir),
                "--epochs",
                str(DEFAULT_EPOCHS),
                "--patience",
                str(DEFAULT_PATIENCE),
                "--batch-size",
                str(DEFAULT_BATCH_SIZE),
                "--num-workers",
                "0",
                "--accelerator",
                "cpu",
                "--devices",
                "auto",
            )
        )

        evaluation_dir = ensure_dir(model_dir / "evaluation")
        print_stage("Step E.3", "Evaluate Chemprop model performance on the held-out test set", log_path=self.run_log_path)
        self._run(
            python_script_cmd(
                self.scripts_dir / "evaluate_chemprop_model.py",
                "--model-dir",
                str(model_dir),
                "--output-dir",
                str(evaluation_dir),
                python_bin=self.chemprop_python,
            )
        )

        analysis_dir = ensure_dir(model_dir / "analysis")
        print_stage("Step E.4", "Rank predicted targets and apply model-quality filtering rules", log_path=self.run_log_path)
        self._run(
            python_script_cmd(
                self.scripts_dir / "analyze_chemprop_predictions.py",
                "--predictions-csv",
                str(model_dir / "predictions.csv"),
                "--task-summary-csv",
                str(chemprop_data_dir / "task_summary.csv"),
                "--metrics-csv",
                str(evaluation_dir / "test_metrics_by_target.csv"),
                "--output-dir",
                str(analysis_dir),
                "--min-test-r2",
                str(DEFAULT_MIN_TEST_R2),
                "--min-test-molecules",
                str(DEFAULT_MIN_TEST_MOLECULES),
                "--min-training-molecules",
                str(DEFAULT_MIN_TRAINING_MOLECULES),
                python_bin=self.chemprop_python,
            )
        )

        print_stage("Step E.5", "Render the filtered Figure 1e target-prioritization panel", log_path=self.run_log_path)
        self._run(
            python_script_cmd(
                self.scripts_dir / "plot_figure1e.py",
                "--ranked-targets-csv",
                str(analysis_dir / "ranked_targets.csv"),
                "--metrics-csv",
                str(evaluation_dir / "test_metrics_by_target.csv"),
                "--cas",
                self.compound.cas,
                "--compound-name",
                self.compound.name,
                "--output-png",
                str(self.case_paths.e_dir / "Figure1e_filtered.png"),
                "--output-pdf",
                str(self.case_paths.e_dir / "Figure1e_filtered.pdf"),
                "--top-n",
                "8",
                "--filter-status",
                "keep",
                python_bin=self.chemprop_python,
            )
        )

    def _detect_chemprop_python(self) -> str:
        candidates = [
            "/Users/liuyue/miniconda3/envs/Chemprop/bin/python",
        ]
        for candidate in candidates:
            if Path(candidate).exists():
                return candidate
        return "python3"

    def _run(self, cmd: list[str]) -> None:
        if self.dry_run:
            with self.run_log_path.open("a", encoding="utf-8") as f:
                f.write(f"DRY-RUN {' '.join(cmd)}\n")
            return
        run_command(cmd, cwd=self.project_root, log_path=self.run_log_path)

    def _write_status(self) -> None:
        write_text(self.status_path, json.dumps(self.status, ensure_ascii=False, indent=2))

    def _update_step_status(self, step_name: str, status: str, extra: dict | None = None) -> None:
        payload = {"status": status}
        if extra:
            payload.update(extra)
        self.status["steps"][step_name] = payload
        self._write_status()

    def _update_source_status(
        self,
        source_name: str,
        status: str,
        output_file: Path,
        message: str | None = None,
        used_placeholder: bool = False,
    ) -> None:
        payload = {
            "status": status,
            "output_file": str(output_file),
            "used_placeholder": used_placeholder,
        }
        if message:
            payload["message"] = message
        self.status["sources"][source_name] = payload
        self._write_status()

    def _fetch_source_with_fallback(
        self,
        source_name: str,
        output_file: Path,
        cmd: list[str],
        placeholder_columns: list[str],
    ) -> None:
        if output_file.exists():
            self._update_source_status(source_name, "ok", output_file)
            return
        if self.dry_run:
            self._update_source_status(source_name, "dry_run", output_file)
            self._run(cmd)
            return
        try:
            self._run(cmd)
            if not output_file.exists():
                raise FileNotFoundError(f"{source_name} did not produce the expected output file: {output_file}")
            self._update_source_status(source_name, "ok", output_file)
        except (subprocess.CalledProcessError, FileNotFoundError, RuntimeError) as exc:
            self._write_placeholder_source_file(output_file, placeholder_columns)
            self._update_source_status(
                source_name,
                "failed",
                output_file,
                message=str(exc),
                used_placeholder=True,
            )

    def _write_placeholder_source_file(self, output_file: Path, columns: list[str]) -> None:
        ensure_dir(output_file.parent)
        with output_file.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(columns)


def run_single_case(
    cas: str,
    disease_keywords: list[str],
    output_root: Path,
    project_root: Path | None = None,
    run_date: str | None = None,
    dry_run: bool = False,
) -> CasePaths:
    runner = PipelineRunner(
        cas=cas,
        disease_keywords=disease_keywords,
        output_root=output_root,
        project_root=project_root,
        run_date=run_date,
        dry_run=dry_run,
    )
    return runner.run()
