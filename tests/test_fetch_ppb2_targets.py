from __future__ import annotations

import http.client
import importlib.util
import unittest
from pathlib import Path
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
PPB2_SCRIPT = REPO_ROOT / "scripts" / "fetch_ppb2_targets.py"


def load_ppb2_module():
    spec = importlib.util.spec_from_file_location("fetch_ppb2_targets", PPB2_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class FakeUrlopen:
    def __init__(self) -> None:
        self.calls = 0

    def __call__(self, req, timeout: int):
        self.calls += 1
        if self.calls == 1:
            raise http.client.RemoteDisconnected("temporary disconnect")
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def read(self) -> bytes:
        return b"<html>ok</html>"


class Ppb2FetchTests(unittest.TestCase):
    def test_fetch_text_retries_remote_disconnect(self) -> None:
        ppb2 = load_ppb2_module()
        fake_urlopen = FakeUrlopen()

        with patch.object(ppb2.urllib.request, "urlopen", side_effect=fake_urlopen), patch.object(ppb2.time, "sleep", return_value=None):
            text = ppb2.fetch_text("https://ppb2.gdb.tools/predict", max_retries=2)

        self.assertEqual(text, "<html>ok</html>")
        self.assertEqual(fake_urlopen.calls, 2)


if __name__ == "__main__":
    unittest.main()
