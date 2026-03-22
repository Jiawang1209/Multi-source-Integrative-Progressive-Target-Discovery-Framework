from __future__ import annotations

from importlib.resources import files
from pathlib import Path


STEP_DIRS = {
    "a": "a_source_collection",
    "b": "b_venn",
    "c": "c_go_kegg",
    "d": "d_kegg_circlize",
    "e": "e_chemprop",
}

DEFAULT_TOP_N_KEGG_TARGETS = 40
DEFAULT_MIN_PATHWAY_N = 1
DEFAULT_MIN_PLATFORM_VOTE = 2
DEFAULT_MIN_TASK_SAMPLES = 100
DEFAULT_MIN_TEST_R2 = 0.30
DEFAULT_MIN_TEST_MOLECULES = 30
DEFAULT_MIN_TRAINING_MOLECULES = 100
DEFAULT_EPOCHS = 20
DEFAULT_PATIENCE = 5
DEFAULT_BATCH_SIZE = 64


def project_root_from_file() -> Path:
    return Path(__file__).resolve().parents[2]


def package_resource_path(*parts: str) -> Path:
    return Path(files("miptd").joinpath("resources", *parts))
