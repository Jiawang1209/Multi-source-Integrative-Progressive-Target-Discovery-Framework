from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
TRAIN_SCRIPT = REPO_ROOT / "scripts" / "train_chemprop_multitask.py"


def load_train_module():
    spec = importlib.util.spec_from_file_location("train_chemprop_multitask", TRAIN_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class TrainChempropMultitaskTests(unittest.TestCase):
    def test_default_training_metrics_omit_r2_for_sparse_multitask_splits(self) -> None:
        train = load_train_module()

        with patch.object(
            sys,
            "argv",
            [
                "train_chemprop_multitask.py",
                "--data-csv",
                "data.csv",
                "--inference-csv",
                "inference.csv",
                "--output-dir",
                "model",
            ],
        ):
            args = train.parse_args()

        self.assertEqual(args.metrics, ["rmse", "mae"])


if __name__ == "__main__":
    unittest.main()
