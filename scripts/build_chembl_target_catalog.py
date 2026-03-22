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
    "Referer": "https://www.ebi.ac.uk/chembl/",
}

FIELDNAMES = [
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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a Homo sapiens ChEMBL target catalog resource from the ChEMBL target API."
    )
    parser.add_argument("--output-csv", default="src/miptd/resources/ChEMBL_target_catalog.csv")
    parser.add_argument(
        "--summary-json",
        default="src/miptd/resources/ChEMBL_target_catalog.summary.json",
    )
    parser.add_argument("--limit", type=int, default=1000)
    parser.add_argument("--sleep-seconds", type=float, default=0.2)
    parser.add_argument("--request-timeout", type=float, default=120.0)
    parser.add_argument("--max-retries", type=int, default=6)
    return parser.parse_args()


def log(message: str) -> None:
    print(message, flush=True)


def fetch_json(url: str, request_timeout: float, max_retries: int) -> dict:
    last_error = None
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers=BROWSER_HEADERS)
            with urllib.request.urlopen(req, timeout=request_timeout) as response:
                return json.load(response)
        except Exception as exc:
            last_error = exc
            log(f"[retry] {attempt + 1}/{max_retries} failed for {url}: {exc}")
            time.sleep(min(10.0, 1.5 * (attempt + 1)))
    raise last_error


def dedup_accessions(target_components) -> str:
    values = []
    seen = set()
    for component in target_components or []:
        accession = (component.get("accession") or "").strip()
        if accession and accession not in seen:
            seen.add(accession)
            values.append(accession)
    return "|".join(values)


def build_rows(limit: int, sleep_seconds: float, request_timeout: float, max_retries: int) -> list[dict[str, str]]:
    rows = []
    offset = 0
    page = 0
    while True:
        page += 1
        query = urllib.parse.urlencode({"limit": limit, "offset": offset})
        url = f"https://www.ebi.ac.uk/chembl/api/data/target.json?{query}"
        payload = fetch_json(url, request_timeout=request_timeout, max_retries=max_retries)
        batch = payload.get("targets", [])
        if not batch:
            break
        for target in batch:
            if (target.get("organism") or "").strip() != "Homo sapiens":
                continue
            rows.append(
                {
                    "ChEMBL ID": target.get("target_chembl_id") or "",
                    "Name": target.get("pref_name") or "",
                    "Accessions": dedup_accessions(target.get("target_components")),
                    "Type": target.get("target_type") or "",
                    "Organism": target.get("organism") or "",
                    "Compounds": "",
                    "Activities": "",
                    "Tax ID": target.get("tax_id") or "",
                    "Species Group Flag": (
                        target.get("species_group_flag")
                        if target.get("species_group_flag") is not None
                        else ""
                    ),
                }
            )
        log(f"[chembl-catalog] page={page} offset={offset} total_rows={len(rows)}")
        page_meta = payload.get("page_meta", {})
        if not page_meta.get("next"):
            break
        offset += page_meta.get("limit", limit)
        time.sleep(sleep_seconds)
    rows.sort(key=lambda row: (row["ChEMBL ID"], row["Name"]))
    return rows


def main():
    args = parse_args()
    output_csv = Path(args.output_csv)
    summary_json = Path(args.summary_json)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)

    rows = build_rows(
        limit=args.limit,
        sleep_seconds=args.sleep_seconds,
        request_timeout=args.request_timeout,
        max_retries=args.max_retries,
    )

    with output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter=";")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "file": str(output_csv),
        "n_rows": len(rows),
        "organism": "Homo sapiens",
        "fields": FIELDNAMES,
    }
    summary_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    log(f"[done] wrote {output_csv}")
    log(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
