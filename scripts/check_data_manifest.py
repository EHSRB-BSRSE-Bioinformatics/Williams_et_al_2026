"""Verify that each Git-distributed data manifest entry matches a tracked file."""

import csv
import hashlib
from pathlib import Path
import re
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parents[1]
    tracked = set(
        subprocess.check_output(["git", "ls-files", "-z"], cwd=root).split(b"\0")
    )
    errors = []
    checked = 0
    with (root / "docs/data_manifest.tsv").open(encoding="utf-8", newline="") as manifest:
        for row in csv.DictReader(manifest, delimiter="\t"):
            if row["distribution"] != "Git":
                continue
            path = row["path"]
            expected = row["sha256"]
            if path.encode("utf-8") not in tracked or not (root / path).is_file():
                errors.append(f"{path}: missing tracked file")
                continue
            if not re.fullmatch("[0-9a-f]{64}", expected):
                errors.append(f"{path}: invalid SHA-256 in data manifest")
                continue
            with (root / path).open("rb") as source:
                observed = hashlib.file_digest(source, "sha256").hexdigest()
            if observed != expected:
                errors.append(f"{path}: expected {expected}, got {observed}")
            checked += 1
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"Validated {checked} tracked data manifest checksums.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
