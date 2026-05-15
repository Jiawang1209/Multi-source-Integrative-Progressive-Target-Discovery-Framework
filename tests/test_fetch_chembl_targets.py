from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
CHEMBL_SCRIPT = REPO_ROOT / "scripts" / "fetch_chembl_targets.py"


def load_chembl_module():
    spec = importlib.util.spec_from_file_location("fetch_chembl_targets", CHEMBL_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class ChemblFetchTests(unittest.TestCase):
    def test_choose_molecule_rejects_name_search_hit_when_exact_identifier_mismatches(self) -> None:
        chembl = load_chembl_module()
        molecules = [
            {
                "molecule_chembl_id": "CHEMBL_WRONG",
                "molecule_structures": {
                    "standard_inchi_key": "WRONG-INCHIKEY",
                    "canonical_smiles": "CCC",
                },
            }
        ]

        with self.assertRaisesRegex(ValueError, "No ChEMBL molecule matched the provided exact identifiers"):
            chembl.choose_molecule(
                molecules,
                match_inchikey="RIGHT-INCHIKEY",
                match_smiles="CCO",
            )

    def test_choose_molecule_still_falls_back_without_exact_identifiers(self) -> None:
        chembl = load_chembl_module()
        molecules = [{"molecule_chembl_id": "CHEMBL1", "molecule_structures": {}}]

        self.assertEqual(chembl.choose_molecule(molecules)["molecule_chembl_id"], "CHEMBL1")


if __name__ == "__main__":
    unittest.main()
