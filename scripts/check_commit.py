#!/usr/bin/env python3
"""Validate the actual Git index, including previously tracked/force-added files.

This is a code-only allowlist guard, not a detector of every possible value
embedded in source text. Keep observations and clinical records out of code.
"""
import re
import subprocess
import sys
from pathlib import Path

ALLOWED = {
    '.gitattributes',
    '.githooks/pre-commit',
    '.gitignore',
    'R/config.R',
    'R/gene_sets.R',
    'README.md',
    'config.example.R',
    'docs/ENVIRONMENT.md',
    'docs/INPUTS.md',
    'docs/METHODS.md',
    'docs/PANEL_MAP.md',
    'docs/PROVENANCE.md',
    'docs/VALIDATION.md',
    'figure1/bulk.R',
    'figure1/copykat.R',
    'figure1/sc_preprocess.R',
    'figure1/single_cell.R',
    'figure4/analyze.R',
    'figure4/plot.R',
    'figure4/preprocess_and_fit.R',
    'scripts/check_commit.py',
    'scripts/check_environment.R',
    'scripts/run.R',
}
MAX_BYTES = 1024 * 1024

def git(*args):
    return subprocess.check_output(["git", *args])

def main():
    root = Path(git("rev-parse", "--show-toplevel").decode().strip())
    if not (root / "scripts/check_commit.py").is_file():
        raise SystemExit("Run this guard in the code-only repository.")
    entries = git("ls-files", "--stage", "-z").split(b"\0")
    failures = []
    checked = 0
    for entry in entries:
        if not entry:
            continue
        meta, raw_path = entry.split(b"\t", 1)
        mode, oid, stage = meta.split()
        path = raw_path.decode("utf-8", errors="strict")
        if stage != b"0":
            failures.append(f"{path}: unresolved merge stage")
            continue
        if path not in ALLOWED:
            failures.append(f"{path}: not on the code/documentation allowlist")
            continue
        if mode not in (b"100644", b"100755"):
            failures.append(f"{path}: symlinks/submodules are not allowed")
            continue
        size = int(git("cat-file", "-s", oid.decode()))
        if size > MAX_BYTES:
            failures.append(f"{path}: larger than the code-file limit")
            continue
        data = git("cat-file", "blob", oid.decode())
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            failures.append(f"{path}: binary or non-UTF-8 content")
            continue
        if b"\0" in data:
            failures.append(f"{path}: binary NUL content")
        # Local user paths belong in ignored runtime configuration.
        if re.search(r"/(?:Users|home)/[A-Za-z0-9_.-]+/", text):
            failures.append(f"{path}: contains a machine-specific home path")
        checked += 1
    if failures:
        print("Commit blocked: this repository permits code/documentation only.", file=sys.stderr)
        for failure in failures:
            print("  " + failure, file=sys.stderr)
        return 1
    print(f"PASS: {checked} indexed code/documentation files; no disallowed paths or binary blobs.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
