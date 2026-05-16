#!/usr/bin/env python3

import argparse
import csv
import html
import http.client
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from urllib.error import HTTPError, URLError


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
    "Referer": "https://ppb2.gdb.tools/",
}


METHOD_CONFIG = {
    "NN_ECfp4": {"fp": "ECfp4", "method": "Sim", "scoringmethod": "TANIMOTO"},
    "NN_Xfp": {"fp": "Xfp", "method": "Sim", "scoringmethod": "CBD"},
    "NN_MQN": {"fp": "MQN", "method": "Sim", "scoringmethod": "CBD"},
    "NNML_ECfp4": {"fp": "ECfp4", "method": "SimPlusNaiveBayes", "scoringmethod": "TANIMOTO"},
    "NNML_Xfp": {"fp": "Xfp", "method": "SimPlusNaiveBayes", "scoringmethod": "CBD"},
    "NNML_MQN": {"fp": "MQN", "method": "SimPlusNaiveBayes", "scoringmethod": "CBD"},
    "ML_ECfp4": {"fp": "ECfp4", "method": "NaiveBayes", "scoringmethod": "TANIMOTO"},
    "DNN_ECfp4": {"fp": "ECfp4", "method": "DNN", "scoringmethod": "TANIMOTO"},
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fetch target predictions from the PPB2 web interface."
    )
    parser.add_argument("--smiles", required=True)
    parser.add_argument("--compound-id", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--method-id", default="NNML_ECfp4", choices=sorted(METHOD_CONFIG))
    parser.add_argument("--max-ranks", type=int, default=20)
    parser.add_argument("--chembl-target-catalog", default=None)
    parser.add_argument("--idmapping-tsv", default=None)
    return parser.parse_args()


def log(message):
    print(message, flush=True)


def fetch_text(url, max_retries=5):
    last_error = None
    for attempt in range(max_retries):
        try:
            log(f"[ppb2] fetch html attempt={attempt + 1} url={url}")
            req = urllib.request.Request(url, headers=BROWSER_HEADERS)
            with urllib.request.urlopen(req, timeout=180) as response:
                return response.read().decode("utf-8", errors="replace")
        except (HTTPError, URLError, TimeoutError, http.client.RemoteDisconnected) as exc:
            last_error = exc
            log(f"[ppb2] retry {attempt + 1}/{max_retries} failed: {exc}")
            time.sleep(2.0 * (attempt + 1))
    raise last_error


def fetch_json(url, max_retries=5):
    last_error = None
    for attempt in range(max_retries):
        try:
            log(f"[ppb2] fetch json attempt={attempt + 1} url={url}")
            req = urllib.request.Request(url, headers=BROWSER_HEADERS)
            with urllib.request.urlopen(req, timeout=180) as response:
                return json.load(response)
        except (HTTPError, URLError, TimeoutError, http.client.RemoteDisconnected) as exc:
            last_error = exc
            log(f"[ppb2] retry {attempt + 1}/{max_retries} failed: {exc}")
            time.sleep(2.0 * (attempt + 1))
    raise last_error


def build_predict_url(smiles, method_id):
    config = METHOD_CONFIG[method_id]
    query = urllib.parse.urlencode(
        {
            "smi": smiles,
            "fp": config["fp"],
            "method": config["method"],
            "scoringmethod": config["scoringmethod"],
        }
    )
    return f"https://ppb2.gdb.tools/predict?{query}"


def parse_rows(html_text, max_ranks):
    rows = []
    pattern = re.compile(
        r"<tr>\s*<td>\s*(?P<rank>\d+)\s*</td>\s*"
        r"<td>\s*<a [^>]*>\s*(?P<chembl>CHEMBL\d+)</a>\s*</td>\s*"
        r"<td>\s*(?P<name>.*?)\s*</td>\s*"
        r'<td>\s*<button [^>]*showMe\(\'(?P<blob>.*?)\'\)">',
        re.S,
    )
    for match in pattern.finditer(html_text):
        rank = int(match.group("rank"))
        if rank > max_ranks:
            break
        chembl_id = match.group("chembl").strip()
        name = html.unescape(re.sub(r"\s+", " ", match.group("name")).strip())
        blob = html.unescape(match.group("blob"))
        score = None
        pieces = blob.split("\t")
        if len(pieces) >= 3:
            try:
                score = float(pieces[2].split("=NearestNeighbours=")[0])
            except Exception:
                score = None
        rows.append(
            {
                "Rank": rank,
                "ChEMBL ID": chembl_id,
                "Common name": name,
                "PPB2_score": score,
            }
        )
    return rows


def load_idmapping(path):
    mapping = {}
    if not path:
        return mapping
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            from_id = (row.get("From") or "").strip()
            to_id = (row.get("To") or "").strip().upper()
            if not from_id or not to_id:
                continue
            mapping.setdefault(from_id, set()).add(to_id)
    return mapping


def load_target_catalog(path, accession_to_gene):
    target_map = {}
    if not path:
        return target_map
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=";")
        for row in reader:
            chembl_id = (row.get("ChEMBL ID") or "").strip()
            organism = (row.get("Organism") or "").strip()
            if not chembl_id or organism != "Homo sapiens":
                continue
            accessions = [x.strip() for x in (row.get("Accessions") or "").split("|") if x.strip()]
            genes = []
            for accession in accessions:
                genes.extend(sorted(accession_to_gene.get(accession, set())))
            dedup_genes = []
            seen = set()
            for gene in genes:
                if gene not in seen:
                    seen.add(gene)
                    dedup_genes.append(gene)
            target_map[chembl_id] = {
                "UniProt": "|".join(accessions),
                "Symbol": "|".join(dedup_genes),
            }
    return target_map


def fetch_target_mapping_from_chembl(target_chembl_id, accession_to_gene):
    url = f"https://www.ebi.ac.uk/chembl/api/data/target/{urllib.parse.quote(target_chembl_id)}.json"
    payload = fetch_json(url)
    if (payload.get("organism") or "").strip() != "Homo sapiens":
        return {"UniProt": "", "Symbol": ""}
    accessions = []
    genes = []
    for component in payload.get("target_components") or []:
        accession = (component.get("accession") or "").strip()
        if accession:
            accessions.append(accession)
            genes.extend(sorted(accession_to_gene.get(accession, set())))
    dedup_accessions = []
    seen_accessions = set()
    for accession in accessions:
        if accession not in seen_accessions:
            seen_accessions.add(accession)
            dedup_accessions.append(accession)
    dedup_genes = []
    seen_genes = set()
    for gene in genes:
        if gene not in seen_genes:
            seen_genes.add(gene)
            dedup_genes.append(gene)
    return {"UniProt": "|".join(dedup_accessions), "Symbol": "|".join(dedup_genes)}


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    predict_url = build_predict_url(args.smiles, args.method_id)
    log(f"[ppb2] compound={args.compound_id} method={args.method_id}")
    html_text = fetch_text(predict_url)
    rows = parse_rows(html_text, args.max_ranks)
    log(f"[ppb2] parsed_rows={len(rows)}")
    accession_to_gene = load_idmapping(args.idmapping_tsv)
    target_catalog = load_target_catalog(args.chembl_target_catalog, accession_to_gene)

    csv_path = output_dir / "ppb2_targets.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = ["PPB2", "ChEMBL ID", "Common name", "Rank", "PPB2_score", "UniProt", "Symbol"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            log(f"[ppb2] map target chembl={row['ChEMBL ID']} rank={row['Rank']}")
            mapping = target_catalog.get(row["ChEMBL ID"], {})
            if not mapping.get("Symbol"):
                mapping = fetch_target_mapping_from_chembl(row["ChEMBL ID"], accession_to_gene)
            writer.writerow(
                {
                    "PPB2": args.compound_id,
                    "ChEMBL ID": row["ChEMBL ID"],
                    "Common name": row["Common name"],
                    "Rank": row["Rank"],
                    "PPB2_score": row["PPB2_score"],
                    "UniProt": mapping.get("UniProt", ""),
                    "Symbol": mapping.get("Symbol", ""),
                }
            )

    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(
            {
                "compound_id": args.compound_id,
                "method_id": args.method_id,
                "predict_url": predict_url,
                "n_rows": len(rows),
            },
            f,
            ensure_ascii=False,
            indent=2,
        )

    with (output_dir / "raw_result.html").open("w", encoding="utf-8") as f:
        f.write(html_text)
    log(f"[ppb2] done output={csv_path}")


if __name__ == "__main__":
    main()
