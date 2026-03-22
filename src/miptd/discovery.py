from __future__ import annotations

from datetime import date
from pathlib import Path
import json
import time
import urllib.parse
import urllib.request

from .config import STEP_DIRS, package_resource_path, project_root_from_file
from .models import CasePaths, CompoundRecord
from .utils import ensure_dir


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_USER_AGENT = "miptd/0.1"


def fetch_json(url: str, max_retries: int = 5) -> dict:
    last_error = None
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": PUBCHEM_USER_AGENT})
            with urllib.request.urlopen(req, timeout=60) as response:
                return json.load(response)
        except Exception as exc:
            last_error = exc
            time.sleep(min(10.0, 1.5 * (attempt + 1)))
    raise last_error


def resolve_compound_from_pubchem(cas: str) -> CompoundRecord:
    cas_encoded = urllib.parse.quote(cas)
    cid_url = f"{PUBCHEM_BASE}/compound/name/{cas_encoded}/cids/JSON"
    cid_payload = fetch_json(cid_url)
    cids = ((cid_payload.get("IdentifierList") or {}).get("CID") or [])
    if not cids:
        raise FileNotFoundError(f"CAS {cas} could not be resolved in PubChem.")

    cid = cids[0]
    prop_url = (
        f"{PUBCHEM_BASE}/compound/cid/{cid}/property/"
        "Title,IUPACName,MolecularFormula,CanonicalSMILES,ConnectivitySMILES,InChIKey/JSON"
    )
    prop_payload = fetch_json(prop_url)
    properties = ((prop_payload.get("PropertyTable") or {}).get("Properties") or [])
    if not properties:
        raise FileNotFoundError(f"PubChem returned CID {cid} for CAS {cas}, but no property block was available.")

    row = properties[0]
    title = (row.get("Title") or "").strip()
    iupac_name = (row.get("IUPACName") or "").strip()
    compound_name = title or iupac_name or cas
    return CompoundRecord(
        compound=compound_name,
        cas=cas,
        name=compound_name,
        formula=(row.get("MolecularFormula") or "").strip(),
        canonical_smiles=((row.get("CanonicalSMILES") or row.get("ConnectivitySMILES") or "").strip()),
        inchikey=(row.get("InChIKey") or "").strip(),
        identity_source="pubchem",
        pubchem_cid=str(cid),
    )


def load_compound_record(cas: str) -> CompoundRecord:
    return resolve_compound_from_pubchem(cas)


def create_case_paths(base_dir: Path, cas: str, run_date: str | None = None) -> CasePaths:
    if run_date is None:
        run_date = date.today().isoformat()
    case_root = ensure_dir(base_dir / f"CAS_{cas}_{run_date}")
    a_dir = ensure_dir(case_root / STEP_DIRS["a"])
    b_dir = ensure_dir(case_root / STEP_DIRS["b"])
    c_dir = ensure_dir(case_root / STEP_DIRS["c"])
    d_dir = ensure_dir(case_root / STEP_DIRS["d"])
    e_dir = ensure_dir(case_root / STEP_DIRS["e"])
    work_dir = ensure_dir(case_root / "_work")
    return CasePaths(case_root=case_root, a_dir=a_dir, b_dir=b_dir, c_dir=c_dir, d_dir=d_dir, e_dir=e_dir, work_dir=work_dir)


def discover_latest_idmapping(project_root: Path | None = None) -> Path:
    package_path = package_resource_path("idmapping.tsv")
    if package_path.exists():
        return package_path
    root = project_root or project_root_from_file()
    preferred = root / "resources" / "idmapping.tsv"
    if preferred.exists():
        return preferred
    candidates = sorted(root.glob("**/idmapping*.tsv"))
    if not candidates:
        raise FileNotFoundError("No idmapping TSV file was found in the project tree.")
    return candidates[-1]


def discover_chembl_catalog(project_root: Path | None = None) -> Path:
    package_path = package_resource_path("ChEMBL_target_catalog.csv")
    if package_path.exists():
        return package_path
    root = project_root or project_root_from_file()
    preferred = root / "resources" / "ChEMBL_target_catalog.csv"
    if preferred.exists():
        return preferred
    candidates = sorted(root.glob("**/ChEMBL_*.csv"))
    if not candidates:
        raise FileNotFoundError("No ChEMBL target catalog file was found in the project tree.")
    return candidates[0]
