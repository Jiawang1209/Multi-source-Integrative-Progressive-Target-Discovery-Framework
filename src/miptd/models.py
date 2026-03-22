from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class CompoundRecord:
    compound: str
    cas: str
    name: str
    formula: str
    canonical_smiles: str
    inchikey: str
    identity_source: str = ""
    pubchem_cid: str = ""


@dataclass(frozen=True)
class CasePaths:
    case_root: Path
    a_dir: Path
    b_dir: Path
    c_dir: Path
    d_dir: Path
    e_dir: Path
    work_dir: Path
