from __future__ import annotations

import csv
from datetime import datetime
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def copy_file(src: Path, dst: Path) -> Path:
    ensure_dir(dst.parent)
    shutil.copy2(src, dst)
    return dst


def run_command(cmd: list[str], cwd: Path, log_path: Path | None = None) -> None:
    env = os.environ.copy()
    mpl_dir = cwd / "tmp" / "mplconfig"
    xdg_dir = cwd / "tmp" / "fontconfig"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    xdg_dir.mkdir(parents=True, exist_ok=True)
    env.setdefault("MPLCONFIGDIR", str(mpl_dir))
    env.setdefault("XDG_CACHE_HOME", str(xdg_dir))
    if log_path:
        ensure_dir(log_path.parent)
        with log_path.open("a", encoding="utf-8") as f:
            f.write(f"{timestamp()} [CMD] {' '.join(cmd)}\n")
            f.flush()
            subprocess.run(cmd, cwd=cwd, check=True, env=env, stdout=f, stderr=subprocess.STDOUT)
            f.write(f"{timestamp()} [CMD] completed\n")
    else:
        subprocess.run(cmd, cwd=cwd, check=True, env=env)


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def print_stage(title: str, detail: str | None = None, log_path: Path | None = None) -> None:
    line = f"[MIPTD] {title}"
    if detail:
        line = f"{line}: {detail}"
    print(line, flush=True)
    if log_path:
        ensure_dir(log_path.parent)
        with log_path.open("a", encoding="utf-8") as f:
            f.write(f"{timestamp()} {line}\n")


def python_script_cmd(script: Path, *args: str, python_bin: str | None = None) -> list[str]:
    return [python_bin or sys.executable, str(script), *args]


def r_script_cmd(script: Path, *args: str) -> list[str]:
    return ["Rscript", str(script), *args]


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_text(path: Path, content: str) -> None:
    ensure_dir(path.parent)
    path.write_text(content, encoding="utf-8")


def join_keywords_regex(keywords: Iterable[str]) -> str:
    clean = [keyword.strip() for keyword in keywords if keyword and keyword.strip()]
    if not clean:
        raise ValueError("At least one disease keyword is required.")
    return "|".join(clean)
