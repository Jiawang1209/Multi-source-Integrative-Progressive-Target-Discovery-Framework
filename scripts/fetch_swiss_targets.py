#!/usr/bin/env python3

import argparse
import csv
import json
import re
import time
from pathlib import Path
from typing import Dict, List, Optional

import requests
from bs4 import BeautifulSoup

try:
    from playwright.sync_api import Error as PlaywrightError
    from playwright.sync_api import TimeoutError as PlaywrightTimeoutError
    from playwright.sync_api import sync_playwright
except ImportError:  # pragma: no cover - handled at runtime
    PlaywrightError = RuntimeError
    PlaywrightTimeoutError = RuntimeError
    sync_playwright = None


SWISS_HOME_URL = "http://www.swisstargetprediction.ch/"
SWISS_BASE_URL = "http://www.swisstargetprediction.ch/"
RESULT_URL_RE = re.compile(r"https?://www\.swisstargetprediction\.ch/result\.php\?job=\d+&organism=[A-Za-z_]+")
RESULT_PATH_RE = re.compile(r"result\.php\?job=\d+&organism=[A-Za-z_]+")
CHROME_EXECUTABLE = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Resolve SwissTargetPrediction results for a compound through online browser automation."
    )
    parser.add_argument("--cas", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--smiles", required=True)
    parser.add_argument("--organism", default="Homo_sapiens")
    parser.add_argument("--result-timeout-seconds", type=int, default=300)
    return parser.parse_args()


def log(message: str) -> None:
    print(message, flush=True)


def write_progress(path: Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)


def summarize(path: Path, source: str, extra: Optional[Dict] = None) -> Dict:
    with path.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    payload = {
        "source": source,
        "file": path.name,
        "n_rows": len(rows),
        "columns": list(rows[0].keys()) if rows else [],
    }
    if extra:
        payload.update(extra)
    return payload


def require_playwright() -> None:
    if sync_playwright is None:
        raise RuntimeError("Python package `playwright` is required for Swiss online mode.")
    if not CHROME_EXECUTABLE.exists():
        raise RuntimeError(f"Chrome executable not found: {CHROME_EXECUTABLE}")


def accept_swiss_dialog(page) -> None:
    try:
        dialog = page.wait_for_event("dialog", timeout=3000)
        dialog.accept()
    except PlaywrightTimeoutError:
        return


def extract_result_url_from_html(html_text: str) -> Optional[str]:
    direct = RESULT_URL_RE.search(html_text)
    if direct:
        return direct.group(0)
    path_match = RESULT_PATH_RE.search(html_text)
    if path_match:
        return f"{SWISS_BASE_URL}{path_match.group(0)}"
    return None


def submit_online(smiles: str, organism: str, timeout_seconds: int = 180, progress_path: Optional[Path] = None) -> str:
    require_playwright()
    with sync_playwright() as playwright:
        if progress_path:
            write_progress(progress_path, {"stage": "launch_browser"})
        log("[swiss] launch browser")
        browser = playwright.chromium.launch(executable_path=str(CHROME_EXECUTABLE), headless=True)
        context = browser.new_context(ignore_https_errors=True)
        page = context.new_page()
        try:
            if progress_path:
                write_progress(progress_path, {"stage": "goto_home", "url": SWISS_HOME_URL})
            log("[swiss] open homepage")
            page.goto(SWISS_HOME_URL, wait_until="domcontentloaded", timeout=120000)
            if progress_path:
                write_progress(progress_path, {"stage": "handle_dialog"})
            accept_swiss_dialog(page)
            if progress_path:
                write_progress(progress_path, {"stage": "fill_smiles"})
            log("[swiss] fill smiles")
            page.locator("#smilesBox").fill(smiles)
            page.eval_on_selector(
                "#smilesBox",
                """element => {
                    element.dispatchEvent(new Event('input', { bubbles: true }));
                    element.dispatchEvent(new KeyboardEvent('keyup', { bubbles: true }));
                    element.dispatchEvent(new Event('change', { bubbles: true }));
                    if (typeof checkForm === 'function') {
                        checkForm();
                    }
                }""",
            )
            if progress_path:
                write_progress(progress_path, {"stage": "select_organism", "organism": organism})
            page.locator(f"input[name='organism'][value='{organism}']").check()
            if progress_path:
                write_progress(progress_path, {"stage": "submit"})
            log("[swiss] submit prediction")
            page.wait_for_function("() => !document.getElementById('submitButton').disabled", timeout=30000)
            page.locator("#submitButton").click()
            if progress_path:
                write_progress(progress_path, {"stage": "wait_result_url"})
            log("[swiss] wait for result url")
            deadline = time.time() + timeout_seconds
            while time.time() < deadline:
                current_url = page.url
                if RESULT_URL_RE.match(current_url):
                    if progress_path:
                        write_progress(progress_path, {"stage": "result_url_ready", "result_url": current_url})
                    log(f"[swiss] result url={current_url}")
                    return current_url
                try:
                    html_text = page.content()
                except PlaywrightError:
                    page.wait_for_timeout(1000)
                    continue
                discovered = extract_result_url_from_html(html_text)
                if discovered:
                    if progress_path:
                        write_progress(progress_path, {"stage": "result_url_ready", "result_url": discovered})
                    log(f"[swiss] discovered result url from html={discovered}")
                    return discovered
                page.wait_for_timeout(5000)
            raise PlaywrightTimeoutError(f"Timeout {timeout_seconds * 1000}ms exceeded while waiting for Swiss result URL.")
        finally:
            context.close()
            browser.close()


def fetch_result_html(result_url: str) -> str:
    log(f"[swiss] fetch result html: {result_url}")
    response = requests.get(result_url, timeout=120)
    response.raise_for_status()
    return response.text


def parse_result_table(html_text: str) -> List[Dict[str, str]]:
    soup = BeautifulSoup(html_text, "html.parser")
    table = None
    for candidate in soup.find_all("table"):
        headers = [th.get_text(" ", strip=True) for th in candidate.find_all("th")]
        if headers[:7] == [
            "Target",
            "Common name",
            "Uniprot ID",
            "ChEMBL ID",
            "Target Class",
            "Probability*",
            "Known actives (3D/2D)",
        ]:
            table = candidate
            break
    if table is None:
        raise RuntimeError("SwissTargetPrediction result table not found in HTML.")

    rows = []
    for tr in table.find_all("tr"):
        cells = tr.find_all("td")
        if len(cells) < 7:
            continue
        rows.append(
            {
                "Target": cells[0].get_text(" ", strip=True),
                "Common name": cells[1].get_text(" ", strip=True),
                "Uniprot ID": cells[2].get_text(" ", strip=True),
                "ChEMBL ID": cells[3].get_text(" ", strip=True),
                "Target Class": cells[4].get_text(" ", strip=True),
                "Probability*": cells[5].get_text(" ", strip=True),
                "Known actives (3D/2D)": cells[6].get_text(" ", strip=True),
            }
        )
    return rows


def write_csv(path: Path, rows: List[Dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "Target",
            "Common name",
            "Uniprot ID",
            "ChEMBL ID",
            "Target Class",
            "Probability*",
            "Known actives (3D/2D)",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / f"SwissTargetPrediction_{args.cas}.csv"
    progress_path = output_dir / "progress.json"
    try:
        result_url = submit_online(
            args.smiles,
            args.organism,
            timeout_seconds=args.result_timeout_seconds,
            progress_path=progress_path,
        )
        html_text = fetch_result_html(result_url)
        write_progress(progress_path, {"stage": "parse_html", "result_url": result_url})
        rows = parse_result_table(html_text)
        write_csv(destination, rows)
        (output_dir / "raw_result.html").write_text(html_text, encoding="utf-8")
        write_progress(progress_path, {"stage": "completed_online", "result_url": result_url, "n_rows": len(rows)})
        with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
            json.dump(
                summarize(
                    destination,
                    "online",
                    {
                        "result_url": result_url,
                        "job_id": re.search(r"job=(\d+)", result_url).group(1) if "job=" in result_url else None,
                    },
                ),
                f,
                ensure_ascii=False,
                indent=2,
            )
    except Exception as exc:
        write_progress(progress_path, {"stage": "failed_online", "error": str(exc)})
        raise


if __name__ == "__main__":
    main()
