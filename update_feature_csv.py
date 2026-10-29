#!/usr/bin/env python3
"""Harmonize target gene annotations in a feature CSV against a supplied GTF."""

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path


CONTROL_LABELS = {
    "control", "ctrl", "intergenic", "negativecontrol", "nontargeting",
    "nontargetting", "notarget", "notargeting", "notargetting",
    "safecontrol", "safetargeting",
}
ATTRIBUTE = re.compile(r'(?:^|;)\s*(\w+)\s+"([^"]*)"')
REQUIRED_COLUMNS = {"target_gene_id", "target_gene_name"}


def label_key(value):
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def parse_gtf(gtf_path):
    """Read gene records; retain all IDs for symbols with multiple mappings."""
    id_to_name = {}
    name_to_ids = defaultdict(set)
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GTF record at line {line_number}.")
            if fields[2] != "gene":
                continue
            attributes = dict(ATTRIBUTE.findall(fields[8]))
            gene_id = attributes.get("gene_id", "")
            gene_name = attributes.get("gene_name", "")
            if not gene_id or not gene_name:
                raise ValueError(f"Missing gene_id or gene_name at line {line_number}.")
            if gene_id in id_to_name and id_to_name[gene_id] != gene_name:
                raise ValueError("The GTF assigns conflicting names to one gene ID.")
            id_to_name[gene_id] = gene_name
            name_to_ids[gene_name].add(gene_id)
    if not id_to_name:
        raise ValueError("No gene records with gene_id and gene_name were found.")
    return id_to_name, name_to_ids


def convert_row(row, id_to_name, name_to_ids):
    updated = dict(row)
    gene_id = row["target_gene_id"].strip()
    gene_name = row["target_gene_name"].strip()
    keys = {label_key(gene_id), label_key(gene_name)}
    if keys & CONTROL_LABELS:
        pair, status = ("Non-Targeting", "Non-Targeting"), "control"
    elif keys & {"ignore"}:
        pair, status = ("Ignore", "Ignore"), "ignored"
    elif gene_id in id_to_name:
        pair = (gene_id, id_to_name[gene_id])
        status = "unchanged" if pair == (gene_id, gene_name) else "name_updated"
    elif len(name_to_ids.get(gene_name, set())) == 1:
        pair = (next(iter(name_to_ids[gene_name])), gene_name)
        status = "id_updated"
    else:
        pair = ("NA", "NA")
        status = "ambiguous" if len(name_to_ids.get(gene_name, set())) > 1 else "unresolved"
    updated["target_gene_id"], updated["target_gene_name"] = pair
    return updated, status


def update_csv(csv_path, out_path, id_to_name, name_to_ids, strict=False):
    output = Path(out_path)
    if Path(csv_path).resolve() == output.resolve() or output.exists():
        raise ValueError("Choose a new output path; existing files are not overwritten.")
    with open(csv_path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        if not fields or len(set(fields)) != len(fields):
            raise ValueError("The CSV must have a nonempty, unique header.")
        if not REQUIRED_COLUMNS.issubset(fields):
            raise ValueError("The CSV requires target_gene_id and target_gene_name columns.")
        rows, counts = [], Counter()
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed CSV record near line {reader.line_num}.")
            converted, status = convert_row(row, id_to_name, name_to_ids)
            rows.append(converted)
            counts[status] += 1
    if not rows:
        raise ValueError("The feature CSV contains no data rows.")
    if strict and (counts["unresolved"] or counts["ambiguous"]):
        raise ValueError(
            f"Conversion stopped: {counts['unresolved']} unresolved and "
            f"{counts['ambiguous']} ambiguous rows. No output was written."
        )
    with output.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gtf_path", help="Destination annotation GTF or GTF.gz.")
    parser.add_argument("csv_in", help="Input feature-reference CSV.")
    parser.add_argument("csv_out", help="New output CSV; must not already exist.")
    parser.add_argument(
        "--strict", action="store_true",
        help="Stop without writing output if any target is unresolved or ambiguous.",
    )
    args = parser.parse_args()
    try:
        id_to_name, name_to_ids = parse_gtf(args.gtf_path)
        counts = update_csv(
            args.csv_in, args.csv_out, id_to_name, name_to_ids, strict=args.strict,
        )
    except (OSError, ValueError, csv.Error) as error:
        parser.exit(1, f"Error: {error}\n")
    print(f"Reference genes: {len(id_to_name)}; feature rows: {sum(counts.values())}")
    for status, count in sorted(counts.items()):
        print(f"{status}: {count}")
    if counts["unresolved"] or counts["ambiguous"]:
        print("WARNING: NA annotations require manual review before downstream use.")


if __name__ == "__main__":
    main()
