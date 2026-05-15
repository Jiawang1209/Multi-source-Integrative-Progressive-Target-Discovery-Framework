#!/usr/bin/env python3

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


BROWSER_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/137.0.0.0 Safari/537.36"
    ),
    "Accept": "application/json,text/plain,*/*",
    "Accept-Language": "en-US,en;q=0.9",
    "Cache-Control": "no-cache",
    "Pragma": "no-cache",
    "Referer": "https://www.ebi.ac.uk/chembl/",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fetch human target evidence for a compound from ChEMBL."
    )
    parser.add_argument("--query", required=True, help="Compound name, InChIKey, or other ChEMBL search query.")
    parser.add_argument("--match-inchikey", default=None, help="Exact standard InChIKey to prioritize.")
    parser.add_argument("--match-smiles", default=None, help="Canonical SMILES to prioritize.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sleep-seconds", type=float, default=0.2)
    parser.add_argument("--max-retries", type=int, default=6)
    parser.add_argument("--request-timeout", type=float, default=60.0)
    return parser.parse_args()


def log(message):
    print(message, flush=True)


def fetch_json(url, max_retries=4, request_timeout=60.0):
    last_error = None
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers=BROWSER_HEADERS)
            with urllib.request.urlopen(req, timeout=request_timeout) as response:
                return json.load(response)
        except Exception as exc:
            last_error = exc
            log(f"[retry] {attempt + 1}/{max_retries} failed for {url}: {exc}")
            time.sleep(1.5 * (attempt + 1))
    raise last_error


def search_molecules(query, max_retries, request_timeout):
    encoded = urllib.parse.quote(query)
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q={encoded}"
    payload = fetch_json(url, max_retries=max_retries, request_timeout=request_timeout)
    return payload.get("molecules", [])


def choose_molecule(molecules, match_inchikey=None, match_smiles=None):
    exact_identifier_requested = bool(match_inchikey or match_smiles)
    if match_inchikey:
        for mol in molecules:
            key = ((mol.get("molecule_structures") or {}).get("standard_inchi_key") or "").upper()
            if key == match_inchikey.upper():
                return mol
    if match_smiles:
        for mol in molecules:
            smiles = ((mol.get("molecule_structures") or {}).get("canonical_smiles") or "")
            if smiles == match_smiles:
                return mol
    if not molecules:
        raise ValueError("No ChEMBL molecules matched the query.")
    if exact_identifier_requested:
        raise ValueError("No ChEMBL molecule matched the provided exact identifiers.")
    return molecules[0]


def fetch_all_activities(molecule_chembl_id, sleep_seconds, max_retries, request_timeout):
    rows = []
    offset = 0
    page_i = 0
    while True:
        page_i += 1
        query = urllib.parse.urlencode(
            {
                "molecule_chembl_id": molecule_chembl_id,
                "limit": 1000,
                "offset": offset,
            }
        )
        url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?{query}"
        log(f"[activities] molecule={molecule_chembl_id} page={page_i} offset={offset}")
        payload = fetch_json(url, max_retries=max_retries, request_timeout=request_timeout)
        batch = payload.get("activities", [])
        if not batch:
            break
        rows.extend(batch)
        page_meta = payload.get("page_meta", {})
        if not page_meta.get("next"):
            break
        offset += page_meta.get("limit", 1000)
        time.sleep(sleep_seconds)
    return rows


def fetch_target(target_chembl_id, max_retries, request_timeout):
    url = f"https://www.ebi.ac.uk/chembl/api/data/target/{urllib.parse.quote(target_chembl_id)}.json"
    return fetch_json(url, max_retries=max_retries, request_timeout=request_timeout)


def clean_accessions(target_components):
    values = []
    for component in target_components or []:
        accession = component.get("accession")
        if accession:
            values.append(accession)
    deduped = []
    seen = set()
    for value in values:
        if value not in seen:
            seen.add(value)
            deduped.append(value)
    return "|".join(deduped)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log(f"[search] query={args.query}")
    molecules = search_molecules(args.query, max_retries=args.max_retries, request_timeout=args.request_timeout)
    log(f"[search] matched_molecules={len(molecules)}")
    chosen = choose_molecule(molecules, match_inchikey=args.match_inchikey, match_smiles=args.match_smiles)
    molecule_id = chosen["molecule_chembl_id"]
    log(f"[selected] molecule_chembl_id={molecule_id} pref_name={chosen.get('pref_name')}")

    activities = fetch_all_activities(
        molecule_id,
        args.sleep_seconds,
        max_retries=args.max_retries,
        request_timeout=args.request_timeout,
    )
    log(f"[activities] total_rows={len(activities)}")
    human_target_ids = []
    for row in activities:
        target_id = row.get("target_chembl_id")
        target_organism = row.get("target_organism")
        if not target_id or target_organism != "Homo sapiens":
            continue
        human_target_ids.append(target_id)
    log(f"[activities] human_target_rows={len(human_target_ids)}")

    target_counts = {}
    for target_id in human_target_ids:
        target_counts[target_id] = target_counts.get(target_id, 0) + 1

    target_rows = []
    for target_id, n_activities in sorted(target_counts.items(), key=lambda x: (-x[1], x[0])):
        log(f"[target] fetch {target_id} activities={n_activities}")
        payload = fetch_target(target_id, max_retries=args.max_retries, request_timeout=args.request_timeout)
        target_rows.append(
            {
                "ChEMBL ID": payload.get("target_chembl_id"),
                "Name": payload.get("pref_name"),
                "Accessions": clean_accessions(payload.get("target_components")),
                "Type": payload.get("target_type"),
                "Organism": payload.get("organism"),
                "Compounds": "",
                "Activities": n_activities,
                "Tax ID": payload.get("tax_id"),
                "Species Group Flag": payload.get("species_group_flag"),
            }
        )

    with (output_dir / "selected_molecule.json").open("w", encoding="utf-8") as f:
        json.dump(chosen, f, ensure_ascii=False, indent=2)

    with (output_dir / "raw_activities.json").open("w", encoding="utf-8") as f:
        json.dump(activities, f, ensure_ascii=False)

    with (output_dir / "chembl_targets.csv").open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "ChEMBL ID",
            "Name",
            "Accessions",
            "Type",
            "Organism",
            "Compounds",
            "Activities",
            "Tax ID",
            "Species Group Flag",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(target_rows)

    summary = {
        "query": args.query,
        "molecule_chembl_id": molecule_id,
        "pref_name": chosen.get("pref_name"),
        "standard_inchi_key": (chosen.get("molecule_structures") or {}).get("standard_inchi_key"),
        "canonical_smiles": (chosen.get("molecule_structures") or {}).get("canonical_smiles"),
        "n_total_activities": len(activities),
        "n_human_target_ids": len(human_target_ids),
        "n_unique_human_targets": len(target_rows),
    }
    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    log(f"[done] unique_human_targets={len(target_rows)}")


if __name__ == "__main__":
    main()
