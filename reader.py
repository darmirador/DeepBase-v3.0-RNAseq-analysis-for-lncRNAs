#!/usr/bin/env python3
"""
reader.py processes all metafile_NNNN.sampleID.csv files containing lncRNA values
from Deepbase v3.0 in the folder and returns per-cancer and summary CSV files.

Per-cancer file:  ENSGID-output_NNNN.csv
  - column 1: tissue type (normal vs cancer)
  - column 2: TCGA barcode
  - column 3: transcript levels (FPKM)

Summary file:     ENSGID.csv

Sample Input: python3 reader.py -i ENSG00000226496.2
"""

import os
import csv
import argparse
from pathlib import Path

import numpy as np
from scipy import stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize lncRNA expression (normal vs cancer) from Deepbase metafiles."
    )
    parser.add_argument(
        "-i",
        "--ensg",
        required=True,
        help="ENSG accession number of the lncRNA (e.g. ENSG00000226496.2)",
    )
    parser.add_argument(
        "-d",
        "--directory",
        default="/Users/darmirador/Desktop/Python Scripts/normalVsCancer_lncRNA",
        help="Directory containing metafile_*.csv files (default: your previous working directory).",
    )
    return parser.parse_args()


def get_metafiles(path: Path):
    return sorted(
        f for f in path.iterdir()
        if f.is_file() and f.name.startswith("metafile") and f.suffix == ".csv"
    )


def cancer_type_from_filename(fname: Path) -> str:
    # expects metafile_NNNN.sampleID.csv
    parts = fname.name.split("_")
    if len(parts) < 2:
        return "UNKNOWN"
    return parts[1].split(".")[0]


def significance_from_p(p_value: float) -> str:
    if np.isnan(p_value):
        return "n/a"
    if p_value < 0.0001:
        return "****"
    elif p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "ns"


def process_file(
    csv_path: Path,
    accession_number: str,
    cancer_type: str,
    out_dir: Path,
    summary_writer: csv.writer,
):
    # read header for TCGA IDs and find the row for the ENSG of interest
    with csv_path.open("r", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)  # first row: TCGA IDs
        tcga_ids = header[1:]  # skip first column (gene ID)

        expression_row = None
        for row in reader:
            if row and row[0] == accession_number:
                expression_row = row[1:]  # values only
                break

    if expression_row is None:
        print(f"[WARN] {accession_number} not found in {csv_path.name}; skipping.")
        return

    if len(expression_row) != len(tcga_ids):
        print(
            f"[WARN] Length mismatch in {csv_path.name}; "
            f"{len(tcga_ids)} IDs vs {len(expression_row)} expression values."
        )
        return

    normal_ids, normal_vals = [], []
    cancer_ids, cancer_vals = [], []

    for barcode, expr in zip(tcga_ids, expression_row):
        if expr == "" or expr is None:
            continue
        value = float(expr)

        # TCGA-nn-nnnn-xxn-nnn-nnnn-nn → xx (sample type)
        tcga_parts = barcode.split("-")
        if len(tcga_parts) < 4:
            continue
        sample_type_code = tcga_parts[3][:2]

        if sample_type_code in {"01", "02"}:  # primary / recurrent tumor
            cancer_ids.append(barcode)
            cancer_vals.append(value)
        elif sample_type_code == "11":        # solid tissue normal
            normal_ids.append(barcode)
            normal_vals.append(value)

    # write per-cancer output file
    output_filename = f"{accession_number}-output_{cancer_type}.csv"
    output_path = out_dir / output_filename

    with output_path.open("w", newline="") as out_fh:
        writer = csv.writer(out_fh)
        writer.writerow(["tissue", "barcode", "fpkm"])
        for b, v in zip(normal_ids, normal_vals):
            writer.writerow(["normal", b, v])
        for b, v in zip(cancer_ids, cancer_vals):
            writer.writerow(["cancer", b, v])

    normal_vals = np.array(normal_vals, dtype=float)
    cancer_vals = np.array(cancer_vals, dtype=float)

    n_normal = normal_vals.size
    n_cancer = cancer_vals.size

    if n_normal > 1 and n_cancer > 1:
        p_value = stats.ttest_ind(normal_vals, cancer_vals).pvalue
    else:
        p_value = float("nan")

    normal_mean = float(np.mean(normal_vals)) if n_normal > 0 else float("nan")
    cancer_mean = float(np.mean(cancer_vals)) if n_cancer > 0 else float("nan")
    significance = significance_from_p(p_value)

    summary_writer.writerow(
        [
            cancer_type,
            normal_mean,
            cancer_mean,
            n_normal,
            n_cancer,
            p_value,
            significance,
        ]
    )

    print(f"Done processing {csv_path.name}... Saved in new file {output_filename}")


def main():
    args = parse_args()
    accession_number = args.ensg
    base_dir = Path(args.directory).expanduser().resolve()
    meta_files = get_metafiles(base_dir)

    if not meta_files:
        print(f"[ERROR] No metafile_*.csv files found in {base_dir}")
        return

    summary_path = base_dir / f"{accession_number}.csv"
    with summary_path.open("w", newline="") as summary_fh:
        summary_writer = csv.writer(summary_fh)
        summary_writer.writerow(
            [
                "cancer_type",
                "normal-mean",
                "cancer-mean",
                "n_normal",
                "n_cancer",
                "p-value",
                "significance_level",
            ]
        )

        for csv_path in meta_files:
            cancer_type = cancer_type_from_filename(csv_path)
            process_file(
                csv_path,
                accession_number=accession_number,
                cancer_type=cancer_type,
                out_dir=base_dir,
                summary_writer=summary_writer,
            )

    print(f"Done processing {accession_number}.csv file for overall statistics.")


if __name__ == "__main__":
    main()
