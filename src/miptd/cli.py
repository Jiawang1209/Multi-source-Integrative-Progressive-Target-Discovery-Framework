from __future__ import annotations

import argparse
from pathlib import Path

from .pipeline import run_single_case


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="MIPTD",
        description="Run the MIPTD single-CAS target discovery pipeline from CAS number and disease keywords.",
    )
    parser.add_argument("--cas", required=True, help="CAS number, e.g. 491-71-4")
    parser.add_argument(
        "--disease-keywords",
        "--disase-keywords",
        dest="disease_keywords",
        required=True,
        help="Comma-separated disease keywords used in enrichment filtering, e.g. 'NAFLD,liver,hepatic,steatosis'",
    )
    parser.add_argument(
        "--output-root",
        default=".",
        help="Base directory where the CAS_xxxxx_Date project folder will be created. Defaults to current directory.",
    )
    parser.add_argument(
        "--run-date",
        default=None,
        help="Optional date for the case directory name. Defaults to today in YYYY-MM-DD format.",
    )
    parser.add_argument(
        "--project-root",
        default=None,
        help="Optional project root containing scripts/ and required reference files. Defaults to the installed package repository root.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the planned commands without executing them.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    keywords = [item.strip() for item in args.disease_keywords.split(",") if item.strip()]
    case_paths = run_single_case(
        cas=args.cas,
        disease_keywords=keywords,
        output_root=Path(args.output_root).resolve(),
        project_root=Path(args.project_root).resolve() if args.project_root else None,
        run_date=args.run_date,
        dry_run=args.dry_run,
    )
    print(case_paths.case_root)


if __name__ == "__main__":
    main()
