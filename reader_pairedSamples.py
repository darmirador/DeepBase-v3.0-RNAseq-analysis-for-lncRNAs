#!/usr/bin/env python3
"""
reader_pairedSamples.py
Description:
    This script pairs matched normal and tumor RNA-seq samples for a given
    lncRNA across multiple TCGA cancer types.

    INPUT FORMAT REQUIREMENT:
        The script expects input files in the form:
            [ENSGID]-output_[TCGA].csv
        
        For example:
            ENSG00000270195.1-output_LIHC.csv
            ENSG00000270195.1-output_KIRC.csv

        Each input CSV must contain the following columns:
            tissue  → "normal" or "tumor"
            barcode → TCGA barcode string
            fpkm    → expression value

    OUTPUT:
        1. Per-cancer pairing files:
            [ENSGID]-pairing_[TCGA].csv
           containing:
            patient, paired_normal, paired_tumor

        2. One combined summary file:
            [ENSGID]-pairing.csv
           containing:
            cancer_type, patient, paired_normal, paired_tumor

    PAIRING LOGIC:
        Normal and tumor samples are matched by comparing the first 12
        characters of the TCGA barcode — the unique patient identifier.
"""

import os
import csv
import argparse
import pandas
from pathlib import Path


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Pair normal and tumor samples for an lncRNA across TCGA cancer types."
    )
    parser.add_argument(
        "-p", "--path",
        default=".",     # User’s working directory by default
        help="Directory containing [ENSGID]-output_[TCGA].csv files."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Ensembl lncRNA accession (e.g., ENSG00000270195.1)."
    )
    return parser.parse_args()


def extract_cancer_type(filename: str, accession: str) -> str:
    """
    Extract cancer type from filename format:
        ENSG00000270195.1-output_LIHC.csv  →  LIHC
    """
    base = os.path.basename(filename)
    pattern = accession + "-output_"

    # Skip files that don't match expected format
    if not base.startswith(pattern):
        return None

    # Remove prefix and strip .csv
    rest = base.replace(pattern, "")
    cancerType = rest.split(".")[0]
    return cancerType


def main():
    args = parse_args()
    pathname = args.path
    accessionNumber = args.input

    # Use working directory
    base_dir = Path(pathname).resolve()
    os.chdir(base_dir)

    directoryFiles = os.listdir(base_dir)

    # ----------------------------------------------------------
    # Step 1: Identify input files belonging to this lncRNA
    # ----------------------------------------------------------
    file_pairs = []
    for fname in directoryFiles:
        cancerType = extract_cancer_type(fname, accessionNumber)
        if cancerType:
            file_pairs.append((fname, cancerType))

    # ----------------------------------------------------------
    # Step 2: Create overall summary file
    # ----------------------------------------------------------
    summary_name = f"{accessionNumber}-pairing.csv"
    with open(summary_name, "w", newline="") as summary_fh:
        summary_writer = csv.writer(summary_fh)
        summary_writer.writerow(["cancer_type", "patient", "paired_normal", "paired_tumor"])

    # ----------------------------------------------------------
    # Step 3: Process each cancer-specific input file
    # ----------------------------------------------------------
    for filename, cancerType in sorted(file_pairs):

        # Output format: [ENSGID]-pairing_[TCGA].csv
        out_name = f"{accessionNumber}-pairing_{cancerType}.csv"

        # Create file with header
        with open(out_name, "w", newline="") as csvOut:
            csvOut.write("patient,paired_normal,paired_tumor\n")

        # Read input file
        csvData = pandas.read_csv(filename)

        cancerState = csvData["tissue"].tolist()
        patientMatched = csvData["barcode"].tolist()
        transcriptMatched = csvData["fpkm"].tolist()

        patient = []
        normal = []
        tumor = []

        totalLength = len(cancerState)
        normalCount = cancerState.count("normal")

        # ------------------------------------------------------
        # Pairing logic:
        # Normal samples come first; tumor samples follow.
        # Match by comparing first 12 characters of TCGA barcode.
        # ------------------------------------------------------
        for i in range(normalCount):
            for j in range(totalLength - 1, normalCount, -1):

                # If patient codes match, we have a pair
                if patientMatched[i][:12] == patientMatched[j][:12]:
                    patient.append(patientMatched[i][:12])
                    normal.append(transcriptMatched[i])
                    tumor.append(transcriptMatched[j])
                    break

        # ------------------------------------------------------
        # Step 4: Write per-cancer paired output file
        # ------------------------------------------------------
        with open(out_name, "a", newline="") as csvOut:
            writer = csv.writer(csvOut)
            for p, n, t in zip(patient, normal, tumor):
                writer.writerow([p, n, t])

        # ------------------------------------------------------
        # Step 5: Append to overall summary file
        # ------------------------------------------------------
        with open(summary_name, "a", newline="") as summary_fh:
            summary_writer = csv.writer(summary_fh)
            for p, n, t in zip(patient, normal, tumor):
                summary_writer.writerow([cancerType, p, n, t])

        print(f"Generated: {out_name}")

    print(f"\nSummary file generated: {summary_name}")


if __name__ == "__main__":
    main()
