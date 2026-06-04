#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Add ZS:Z:<spot> tags to SAM records using a qname-to-spot table generated
from R2 barcode decoding.
"""

import argparse
import csv
import gzip
import sys


def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def norm_qname(q: str) -> str:
    q = q.strip().split()[0]

    if q.endswith("/1") or q.endswith("/2"):
        q = q[:-2]

    return q


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Stream SAM from stdin, add spatial spot tag ZS to reads with status=OK, "
            "and discard non-OK reads. The SAM stream must preserve read order. "
            "When using Bowtie2, run with --reorder and do not use --no-unal, -k, or -a."
        )
    )

    ap.add_argument(
        "--map",
        required=True,
        help="qname_to_spot.tsv or qname_to_spot.tsv.gz."
    )

    args = ap.parse_args()

    with open_text(args.map, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        required_cols = {"qname", "spot", "status"}
        if reader.fieldnames is None or not required_cols.issubset(set(reader.fieldnames)):
            raise RuntimeError(
                f"Map file must contain columns: {sorted(required_cols)}. "
                f"Observed columns: {reader.fieldnames}"
            )

        try:
            current = next(reader)
        except StopIteration:
            raise RuntimeError("Map file is empty.")

        n_sam = 0
        n_tagged = 0
        n_drop = 0

        for line in sys.stdin:
            if line.startswith("@"):
                sys.stdout.write(line)
                continue

            n_sam += 1

            if current is None:
                raise RuntimeError(
                    f"SAM has more records than the map file. "
                    f"The map file ended before SAM record {n_sam}."
                )

            fields = line.rstrip("\n").split("\t")
            sam_qname = norm_qname(fields[0])
            map_qname = current["qname"]

            if sam_qname != map_qname:
                raise RuntimeError(
                    f"QNAME mismatch at record {n_sam}: SAM={sam_qname}, MAP={map_qname}. "
                    f"The SAM stream and qname-to-spot table must have the same read order."
                )

            status = current["status"]
            spot = current["spot"]

            if status == "OK":
                fields.append(f"ZS:Z:{spot}")
                sys.stdout.write("\t".join(fields) + "\n")
                n_tagged += 1
            else:
                n_drop += 1

            try:
                current = next(reader)
            except StopIteration:
                current = None

        if current is not None:
            extra = 1 + sum(1 for _ in reader)
            print(
                f"[WARN] Map file has extra records after SAM ended: {extra}",
                file=sys.stderr
            )

    print(
        f"[DONE] sam_records={n_sam:,} tagged_OK={n_tagged:,} dropped_nonOK={n_drop:,}",
        file=sys.stderr
    )


if __name__ == "__main__":
    main()