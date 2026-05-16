from __future__ import annotations

import importlib.util
import sys
import types
import unittest
from pathlib import Path
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
PREPARE_SCRIPT = REPO_ROOT / "scripts" / "prepare_chemprop_data.py"


def load_prepare_module():
    sys.modules.setdefault("rdkit", types.SimpleNamespace(Chem=types.SimpleNamespace(), RDLogger=types.SimpleNamespace(DisableLog=lambda *_: None)))
    sys.modules.setdefault("rdkit.Chem", types.SimpleNamespace(MolFromSmiles=lambda *_: object()))
    sys.modules.setdefault("rdkit.RDLogger", types.SimpleNamespace(DisableLog=lambda *_: None))
    sys.modules.setdefault(
        "rdkit.Chem.Scaffolds",
        types.SimpleNamespace(MurckoScaffold=types.SimpleNamespace(MurckoScaffoldSmiles=lambda **_: "C")),
    )
    spec = importlib.util.spec_from_file_location("prepare_chemprop_data", PREPARE_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class PrepareChempropDataTests(unittest.TestCase):
    def test_fetch_activities_skips_failed_activity_page(self) -> None:
        prepare = load_prepare_module()

        def fake_fetch_json(url: str, max_retries: int, request_timeout: float) -> dict:
            if "standard_type=IC50" in url:
                raise RuntimeError("ChEMBL temporary 500")
            return {"activities": [], "page_meta": {}}

        with patch.object(prepare, "fetch_json", side_effect=fake_fetch_json):
            rows = prepare.fetch_activities_for_target(
                "CHEMBL123",
                sleep_seconds=0,
                max_retries=1,
                request_timeout=1,
            )

        self.assertEqual(rows, [])


if __name__ == "__main__":
    unittest.main()
