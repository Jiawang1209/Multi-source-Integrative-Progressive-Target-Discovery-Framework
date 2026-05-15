from __future__ import annotations

import importlib.util
import sys
import types
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SWISS_SCRIPT = REPO_ROOT / "scripts" / "fetch_swiss_targets.py"


def load_swiss_module():
    sys.modules.setdefault("requests", types.SimpleNamespace(get=lambda *args, **kwargs: None))
    sys.modules.setdefault("bs4", types.SimpleNamespace(BeautifulSoup=lambda *args, **kwargs: None))
    sys.modules.setdefault(
        "playwright.sync_api",
        types.SimpleNamespace(Error=RuntimeError, TimeoutError=RuntimeError, sync_playwright=None),
    )
    spec = importlib.util.spec_from_file_location("fetch_swiss_targets", SWISS_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class SwissFetchTests(unittest.TestCase):
    def test_extract_result_url_uses_http_base_for_relative_result_path(self) -> None:
        swiss = load_swiss_module()

        self.assertEqual(
            swiss.extract_result_url_from_html("location.href='result.php?job=123&organism=Homo_sapiens'"),
            "http://www.swisstargetprediction.ch/result.php?job=123&organism=Homo_sapiens",
        )


if __name__ == "__main__":
    unittest.main()
