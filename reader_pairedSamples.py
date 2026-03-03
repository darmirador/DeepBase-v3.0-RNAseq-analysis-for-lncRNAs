#!/usr/bin/env python3
"""
reader_pairedSamples.py

Creates:
  1) [ENSGID]-pairing_[TCGA].csv  : patient, paired_normal, paired_tumor
  2) [ENSGID]-pairing.csv         : cancer_type, patient, paired_normal, paired_tumor
  3) [ENSGID]-pairing_stats.csv   : cancer_type,normal-mean,cancer-mean,n_normal,n_cancer,p-value,significance_level
"""

import os
import csv
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

try:
    from scipy import stats
    _HAVE_SCIPY = True
except Exception:
    stats = None
    _HAVE_SCIPY = False


def parse_args():
    parser = argparse.ArgumentParser(
        description="Pair normal and tumor samples for an lncRNA across TCGA cancer types."
    )
    parser.add_argument(
        "-p", "--path",
        default=".",
        help="Directory containing [ENSGID]-output_[TCGA].csv files."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Ensembl lncRNA accession (e.g., ENSG00000270195.1)."
    )
    return parser.parse_args()


def extract_cancer_type(filename: str, accession: str) -> str | None:
    base = os.path.basename(filename)
    prefix = f"{accession}-output_"
    if not (base.startswith(prefix) and base.endswith(".csv")):
        return None
    return base[len(prefix):].split(".")[0]


def significance_from_p(p: float) -> str:
    if p is None or np.isnan(p):
        return "n/a"
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def paired_p_value(normal: np.ndarray, tumor: np.ndarray) -> float:
    if normal.size < 2 or tumor.size < 2:
        return float("nan")
    if not _HAVE_SCIPY:
        return float("nan")
    return float(stats.ttest_rel(tumor, normal, nan_policy="omit").pvalue)


def normalize_tissue_label(x: str) -> str:
    """
    Expect upstream to be 'normal' vs 'cancer' (or 'tumor').
    This makes it robust to casing/whitespace; anything else is left as-is.
    """
    s = str(x).strip().lower()
    if s == "cancer":
        return "tumor"
    return s


def main():
    args = parse_args()
    accession = args.input
    base_dir = Path(args.path).expanduser().resolve()

    files = []
    for fname in os.listdir(base_dir):
        cancer = extract_cancer_type(fname, accession)
        if cancer:
            files.append((base_dir / fname, cancer))
    files.sort(key=lambda x: x[1])

    if not files:
        raise SystemExit(f"[ERROR] No files matching {accession}-output_*.csv found in {base_dir}")

    pairing_all_path = base_dir / f"{accession}-pairing.csv"
    with pairing_all_path.open("w", newline="") as fh:
        csv.writer(fh).writerow(["cancer_type", "patient", "paired_normal", "paired_tumor"])

    stats_path = base_dir / f"{accession}-pairing.csv"
    with stats_path.open("w", newline="") as fh:
        csv.writer(fh).writerow([
            "cancer_type",
            "normal-mean",
            "cancer-mean",
            "n_normal",
            "n_cancer",
            "p-value",
            "significance_level"
        ])

    for in_path, cancer_type in files:
        df = pd.read_csv(in_path)

        # Normalize expected columns
        required = {"tissue", "barcode", "fpkm"}
        missing = required - set(df.columns)
        if missing:
            print(f"[WARN] {in_path.name} missing columns {sorted(missing)}; skipping.")
            continue

        df["patient"] = df["barcode"].astype(str).str[:12]
        df["tissue"] = df["tissue"].map(normalize_tissue_label)

        # Pivot to one row per patient
        wide = df.pivot_table(index="patient", columns="tissue", values="fpkm", aggfunc="first")

        # Ensure columns exist even if absent (prevents KeyError later)
        if "normal" not in wide.columns:
            wide["normal"] = np.nan
        if "tumor" not in wide.columns:
            wide["tumor"] = np.nan

        paired = wide.dropna(subset=["normal", "tumor"]).copy()

        # Build per-cancer pairing DF (may be empty, but columns will exist)
        out_df = paired.reset_index()[["patient", "normal", "tumor"]].rename(
            columns={"normal": "paired_normal", "tumor": "paired_tumor"}
        )

        out_pair_path = base_dir / f"{accession}-pairing_{cancer_type}.csv"
        out_df.to_csv(out_pair_path, index=False)

        # Append to combined pairing file
        if not out_df.empty:
            with pairing_all_path.open("a", newline="") as fh:
                w = csv.writer(fh)
                for _, r in out_df.iterrows():
                    w.writerow([cancer_type, r["patient"], r["paired_normal"], r["paired_tumor"]])

        # Stats-only summary
        normal_vals = paired["normal"].to_numpy(dtype=float) if not paired.empty else np.array([], dtype=float)
        tumor_vals = paired["tumor"].to_numpy(dtype=float) if not paired.empty else np.array([], dtype=float)

        n_pairs = int(normal_vals.size)
        normal_mean = float(np.mean(normal_vals)) if n_pairs > 0 else float("nan")
        cancer_mean = float(np.mean(tumor_vals)) if n_pairs > 0 else float("nan")
        p = paired_p_value(normal_vals, tumor_vals)
        sig = significance_from_p(p)

        with stats_path.open("a", newline="") as fh:
            csv.writer(fh).writerow([cancer_type, normal_mean, cancer_mean, n_pairs, n_pairs, p, sig])

        print(f"Generated: {out_pair_path.name}")

    print(f"\nSummary pairing file generated: {pairing_all_path.name}")
    print(f"Stats-only file generated: {stats_path.name}")
    if not _HAVE_SCIPY:
        print("[WARN] SciPy unavailable/incompatible → p-values are NaN in pairing_stats.csv.")


if __name__ == "__main__":
    main()