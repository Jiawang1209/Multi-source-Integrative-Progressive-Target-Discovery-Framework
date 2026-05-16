#!/usr/bin/env python3

import argparse
import csv
import json
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests
from bs4 import BeautifulSoup


BROWSER_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/137.0.0.0 Safari/537.36"
    ),
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.9",
    "Cache-Control": "no-cache",
    "Pragma": "no-cache",
    "Upgrade-Insecure-Requests": "1",
}


SEA_SEARCH_URL = "https://sea.compbio.ucsf.edu/search"
DEFAULT_FINGERPRINT_TYPE = "rdkit_ecfp"
SEA_SUBMIT_RETRY_DELAYS_SECONDS = (5, 15, 30, 60)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Resolve SEA results for a compound through online requests."
    )
    parser.add_argument("--cas", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--smiles", required=True)
    parser.add_argument("--compound-name", default=None)
    parser.add_argument("--organism", default="Homo sapiens")
    return parser.parse_args()


def log(message: str) -> None:
    print(message, flush=True)


def write_progress(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def sea_safe_compound_id(compound_name: Optional[str]) -> str:
    label = compound_name or "compound_1"
    label = re.sub(r"[\s,\t]+", "_", label.strip())
    label = re.sub(r"[^A-Za-z0-9_.:-]+", "_", label)
    return label.strip("_") or "compound_1"


def request_with_retries(action, label: str):
    max_attempts = len(SEA_SUBMIT_RETRY_DELAYS_SECONDS) + 1
    for attempt in range(1, max_attempts + 1):
        try:
            return action()
        except Exception as exc:
            if attempt == max_attempts:
                raise
            delay_seconds = SEA_SUBMIT_RETRY_DELAYS_SECONDS[attempt - 1]
            log(f"[sea] {label} retry {attempt}/{max_attempts} after transient error: {exc}; sleep {delay_seconds}s")
            time.sleep(delay_seconds)
    raise RuntimeError(f"SEA {label} retry loop ended unexpectedly.")


def submit_online(smiles: str, compound_name: Optional[str]) -> Tuple[requests.Session, str]:
    session = requests.Session()
    session.headers.update(BROWSER_HEADERS)
    log("[sea] open search page")
    home = request_with_retries(lambda: session.get(SEA_SEARCH_URL, timeout=120), "open search page")
    home.raise_for_status()
    match = re.search(r'name="csrf_token" type="hidden" value="([^"]+)"', home.text)
    if not match:
        raise RuntimeError("SEA csrf_token not found.")
    label = sea_safe_compound_id(compound_name)
    payload = {
        "csrf_token": match.group(1),
        "query_custom_fp_type": DEFAULT_FINGERPRINT_TYPE,
        "ref_type": "library",
        "query_type": "custom",
        "query_custom_targets_paste": f"{smiles} {label}",
    }
    log("[sea] submit query")
    result = request_with_retries(
        lambda: session.post(
            SEA_SEARCH_URL,
            data=payload,
            headers={
                "Content-Type": "application/x-www-form-urlencoded",
                "Origin": "https://sea.compbio.ucsf.edu",
                "Referer": SEA_SEARCH_URL,
            },
            allow_redirects=True,
            timeout=300,
        ),
        "submit query",
    )
    result.raise_for_status()
    return session, result.url


def wait_for_result_page(session: requests.Session, job_url: str, timeout_seconds: int = 180, poll_seconds: float = 2.0) -> str:
    deadline = time.time() + timeout_seconds
    last_html = ""
    poll_i = 0
    while time.time() < deadline:
        poll_i += 1
        log(f"[sea] poll result page {poll_i}: {job_url}")
        try:
            response = session.get(job_url, timeout=120)
            response.raise_for_status()
        except Exception as exc:
            log(f"[sea] poll retry after transient error: {exc}")
            time.sleep(poll_seconds)
            continue
        html_text = response.text
        last_html = html_text
        if "In Progress" not in html_text and "pending" not in html_text.lower():
            return html_text
        time.sleep(poll_seconds)
    raise RuntimeError("SEA job did not finish before timeout.")


def parse_result_table(html_text: str) -> List[Dict[str, str]]:
    soup = BeautifulSoup(html_text, "html.parser")
    table = None
    for candidate in soup.find_all("table"):
        headers = [th.get_text(" ", strip=True) for th in candidate.find_all("th")]
        if headers[:6] == ["Query", "Target Key", "Target Name", "Description", "P-Value", "MaxTC"]:
            table = candidate
            break
    if table is None:
        raise RuntimeError("SEA result table not found.")

    rows = []
    tbody = table.find("tbody")
    if tbody is None:
        return rows
    for tr in tbody.find_all("tr"):
        cells = tr.find_all("td")
        if len(cells) == 6:
            query_id = cells[0].get_text(" ", strip=True)
            target_id = cells[1].get_text(" ", strip=True)
            name = cells[2].get_text(" ", strip=True)
            description = cells[3].get_text(" ", strip=True)
            pvalue = cells[4].get_text(" ", strip=True)
            max_tc = cells[5].get_text(" ", strip=True)
        elif len(cells) == 5:
            query_id = "compound_1"
            target_id = cells[0].get_text(" ", strip=True)
            name = cells[1].get_text(" ", strip=True)
            description = cells[2].get_text(" ", strip=True)
            pvalue = cells[3].get_text(" ", strip=True)
            max_tc = cells[4].get_text(" ", strip=True)
        else:
            continue
        rows.append(
            {
                "Query ID": query_id,
                "Target ID": target_id,
                "Affinity Threshold (nM)": "5",
                "P-Value": pvalue,
                "Max Tc": max_tc,
                "Cut Sum": "",
                "Z-Score": "",
                "Name": name,
                "Description": description,
                "Query Smiles": "",
            }
        )
    return rows


def filter_rows_by_organism(rows: List[Dict[str, str]], organism: str) -> List[Dict[str, str]]:
    if organism.lower() != "homo sapiens":
        return rows
    return [row for row in rows if row["Target ID"].endswith("_HUMAN")]


def write_csv(path: Path, rows: List[Dict[str, str]], smiles: Optional[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "Query ID",
            "Target ID",
            "Affinity Threshold (nM)",
            "P-Value",
            "Max Tc",
            "Cut Sum",
            "Z-Score",
            "Name",
            "Description",
            "Query Smiles",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            row["Query Smiles"] = smiles or ""
            writer.writerow(row)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / "sea-results.xls"
    progress_path = output_dir / "progress.json"
    try:
        write_progress(progress_path, {"stage": "submit_online"})
        session, result_url = submit_online(args.smiles, args.compound_name)
        write_progress(progress_path, {"stage": "wait_result_page", "result_url": result_url})
        html_text = wait_for_result_page(session, result_url)
        write_progress(progress_path, {"stage": "parse_html", "result_url": result_url})
        raw_rows = parse_result_table(html_text)
        rows = filter_rows_by_organism(raw_rows, args.organism)
        write_csv(destination, rows, args.smiles)
        (output_dir / "raw_result.html").write_text(html_text, encoding="utf-8")
        write_progress(
            progress_path,
            {
                "stage": "completed_online",
                "result_url": result_url,
                "organism": args.organism,
                "n_rows_raw": len(raw_rows),
                "n_rows": len(rows),
            },
        )
        with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
            json.dump(
                {
                    "source": "online",
                    "result_url": result_url,
                    "organism": args.organism,
                    "n_rows_raw": len(raw_rows),
                    "n_rows": len(rows),
                },
                f,
                ensure_ascii=False,
                indent=2,
            )
    except Exception as exc:
        write_progress(progress_path, {"stage": "failed_online", "error": str(exc)})
        raise


if __name__ == "__main__":
    main()
