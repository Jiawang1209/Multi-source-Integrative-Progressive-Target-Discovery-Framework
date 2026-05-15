from __future__ import annotations

import importlib.util
import sys
import types
import unittest
from pathlib import Path
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
SEA_SCRIPT = REPO_ROOT / "scripts" / "fetch_sea_targets.py"


def load_sea_module():
    sys.modules.setdefault("requests", types.SimpleNamespace(Session=lambda: None))
    sys.modules.setdefault("bs4", types.SimpleNamespace(BeautifulSoup=lambda *args, **kwargs: None))
    spec = importlib.util.spec_from_file_location("fetch_sea_targets", SEA_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class FakeResponse:
    def __init__(self, text: str = "", url: str = "https://sea.compbio.ucsf.edu/jobs/search_test") -> None:
        self.text = text
        self.url = url

    def raise_for_status(self) -> None:
        return None


class FakeSession:
    def __init__(self) -> None:
        self.headers = {}
        self.post_data = None

    def get(self, url: str, timeout: int) -> FakeResponse:
        return FakeResponse(
            '<input name="csrf_token" type="hidden" value="csrf-test">',
            url=url,
        )

    def post(self, url: str, data: dict[str, str], headers: dict[str, str], allow_redirects: bool, timeout: int) -> FakeResponse:
        self.post_data = data
        return FakeResponse()


class SeaFetchTests(unittest.TestCase):
    def test_submit_online_includes_current_fingerprint_field_and_safe_label(self) -> None:
        sea = load_sea_module()
        fake_session = FakeSession()

        with patch.object(sea.requests, "Session", return_value=fake_session):
            sea.submit_online("CCO", "Timosaponin A1")

        assert fake_session.post_data is not None
        self.assertEqual(fake_session.post_data["query_custom_fp_type"], "rdkit_ecfp")
        self.assertEqual(fake_session.post_data["query_custom_targets_paste"], "CCO Timosaponin_A1")


if __name__ == "__main__":
    unittest.main()
