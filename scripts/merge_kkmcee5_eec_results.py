#!/usr/bin/env python3
"""Add the completed KKMCee 5 pair to the six-sample EEC CSV products."""

import argparse
import csv
import io
from pathlib import Path
import subprocess


NEW_SAMPLE = "KKMCee 5.00.02 + Pythia 8.316"
PRODUCTS = {
    "eec_isr_correction.csv": ("eec_isr_correction_s6.csv", 201),
    "eec_isr_correction_regions.csv": ("eec_isr_correction_regions_s6.csv", 6),
    "eec_verification.csv": ("eec_verification_s6.csv", 1),
    "isr_model_summary.csv": ("isr_model_summary_kkmcee5.csv", 2),
}


def read_csv(text, path):
    lines = text.splitlines()
    reader = csv.DictReader(io.StringIO(text))
    rows = list(reader)
    if not reader.fieldnames or len(lines) != len(rows) + 1:
        raise ValueError(f"expected one CSV record per line: {path}")
    return reader.fieldnames, list(zip(rows, lines[1:])), lines[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--new", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--baseline-git-ref", help="read tracked baseline from this ref")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    repo = None
    if args.baseline_git_ref:
        repo = Path(subprocess.check_output(
            ["git", "-C", str(args.baseline), "rev-parse", "--show-toplevel"],
            text=True).strip())

    for target, (new_name, expected_rows) in PRODUCTS.items():
        baseline_path = args.baseline / target
        if repo:
            relative = baseline_path.resolve().relative_to(repo.resolve())
            baseline_text = subprocess.check_output(
                ["git", "-C", str(repo), "show", f"{args.baseline_git_ref}:{relative.as_posix()}"],
                text=True)
        else:
            baseline_text = baseline_path.read_text()
        fields, baseline, header = read_csv(baseline_text, baseline_path)
        new_path = args.new / new_name
        new_fields, addition, _ = read_csv(new_path.read_text(), new_path)
        if fields != new_fields:
            raise ValueError(f"CSV headers differ for {target}")
        baseline = [(row, line) for row, line in baseline if row["sample"] != NEW_SAMPLE]
        if len({row["sample"] for row, _ in baseline}) != 6:
            raise ValueError(f"expected six original samples in {target}")
        if len(addition) != expected_rows or {row["sample"] for row, _ in addition} != {NEW_SAMPLE}:
            raise ValueError(f"unexpected KKMCee 5 rows in {new_name}")
        with (args.output / target).open("w", newline="\n") as stream:
            stream.write(header + "\n")
            for _, line in baseline + addition:
                stream.write(line + "\n")
        print(f"{target}: {len(baseline)} baseline + {len(addition)} KKMCee 5 rows")


if __name__ == "__main__":
    main()
