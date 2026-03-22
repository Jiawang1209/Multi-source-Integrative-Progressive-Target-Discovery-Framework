from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from miptd.discovery import create_case_paths, discover_chembl_catalog, discover_latest_idmapping


class DiscoveryTests(unittest.TestCase):
    def test_package_resources_are_discoverable(self) -> None:
        idmapping = discover_latest_idmapping(REPO_ROOT)
        chembl_catalog = discover_chembl_catalog(REPO_ROOT)

        self.assertTrue(idmapping.exists())
        self.assertTrue(chembl_catalog.exists())
        self.assertEqual(idmapping.parent.name, "resources")
        self.assertEqual(chembl_catalog.parent.name, "resources")
        self.assertIn("src/miptd/resources", str(idmapping))
        self.assertIn("src/miptd/resources", str(chembl_catalog))

    def test_create_case_paths_creates_expected_structure(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            base_dir = Path(tmpdir)
            case_paths = create_case_paths(base_dir, "491-71-4", "2026-03-22")

            self.assertTrue(case_paths.case_root.is_dir())
            self.assertEqual(case_paths.case_root.name, "CAS_491-71-4_2026-03-22")
            self.assertTrue(case_paths.a_dir.is_dir())
            self.assertTrue(case_paths.b_dir.is_dir())
            self.assertTrue(case_paths.c_dir.is_dir())
            self.assertTrue(case_paths.d_dir.is_dir())
            self.assertTrue(case_paths.e_dir.is_dir())
            self.assertTrue(case_paths.work_dir.is_dir())


if __name__ == "__main__":
    unittest.main()
