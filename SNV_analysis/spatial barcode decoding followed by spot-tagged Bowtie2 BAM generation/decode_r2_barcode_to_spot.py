#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Decode DBiT-style spatial barcodes from R2 and write a qname-to-spot table.
"""

import argparse
import gzip
import re
import sys


def iter_fastq(handle):
    while True:
        h = handle.readline()
        if not h:
            return

        s = handle.readline()
        p = handle.readline()
        q = handle.readline()

        if not q:
            raise RuntimeError("Incomplete FASTQ record.")

        yield h, s, p, q


def norm_read_id(header_line: str) -> str:
    h = header_line.strip()

    if h.startswith("@"):
        h = h[1:]

    h = re.split(r"\s+", h, maxsplit=1)[0]

    if h.endswith("/1") or h.endswith("/2"):
        h = h[:-2]

    return h


def read_positions(pos_path: str):
    with open(pos_path, "rt", encoding="utf-8", errors="replace") as f:
        txt = f.read()

    toks = re.split(r"[,\s]+", txt.strip())
    return set(t for t in toks if re.match(r"^\d+x\d+$", t))


def read_barcode_map(map_path: str):

    bc2spot = {}

    with open(map_path, "rt", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()

    if not lines:
        raise RuntimeError(f"Empty barcode map: {map_path}")

    first = lines[0].strip()
    is_header = ("barcode" in first.lower()) and ("x" in first.lower())
    start = 1 if is_header else 0

    for line in lines[start:]:
        line = line.strip()

        if not line:
            continue

        parts = [x.strip() for x in line.split(",")]

        if len(parts) < 2:
            continue

        bc = parts[0].upper()
        spot = parts[1]

        if bc and spot:
            bc2spot[bc] = spot

    return bc2spot


def build_best_match_map(whitelist_barcodes, max_mm=1):

    best = {}

    def add(obs, bc, dist):
        if obs not in best:
            best[obs] = [bc, dist, False]
            return

        cur_bc, cur_dist, _ = best[obs]

        if dist < cur_dist:
            best[obs] = [bc, dist, False]
        elif dist == cur_dist and bc != cur_bc:
            best[obs][2] = True

    for bc in whitelist_barcodes:
        bc = bc.strip().upper()

        if not bc:
            continue

        add(bc, bc, 0)

        if max_mm >= 1:
            for i, ch in enumerate(bc):
                for alt in "ACGT":
                    if alt == ch:
                        continue

                    obs = bc[:i] + alt + bc[i + 1:]
                    add(obs, bc, 1)

    return best


def parse_xy(spot: str):
    m = re.match(r"^(\d+)x(\d+)$", spot)

    if m:
        return m.group(1), m.group(2)

    return "", ""


def get_slice(seq: str, s0: int, s1: int) -> str:
    return seq.strip()[s0:s1].upper()


def main():
    ap = argparse.ArgumentParser(
        description="Assign reads to spatial spots by decoding combined B+A barcodes from R2."
    )

    ap.add_argument("--r2", required=True, help="R2 FASTQ.GZ file.")
    ap.add_argument("--barcode-map", required=True, help="Spatial barcode CSV file.")
    ap.add_argument("--positions", required=True, help="Position file used to retain selected spots.")
    ap.add_argument("--out", required=True, help="Output qname-to-spot TSV.GZ file.")

    ap.add_argument("--x-slice", default="32:40", help="Barcode B slice. Default: 32:40.")
    ap.add_argument("--y-slice", default="70:78", help="Barcode A slice. Default: 70:78.")
    ap.add_argument("--max-mismatch", type=int, default=1, help="Maximum barcode mismatches. Default: 1.")

    args = ap.parse_args()

    x0, x1 = map(int, args.x_slice.split(":"))
    y0, y1 = map(int, args.y_slice.split(":"))

    keep_spots = read_positions(args.positions)
    bc2spot_all = read_barcode_map(args.barcode_map)

    bc2spot = {bc: sp for bc, sp in bc2spot_all.items() if sp in keep_spots}

    if not bc2spot:
        raise RuntimeError("No barcodes were retained after filtering by position file.")

    print(
        f"[INFO] keep_spots={len(keep_spots)} kept_barcodes={len(bc2spot)}",
        file=sys.stderr
    )

    best_map = build_best_match_map(sorted(bc2spot.keys()), max_mm=args.max_mismatch)

    n_total = 0
    n_ok = 0
    n_short = 0
    n_nomap = 0
    n_ambig = 0

    with gzip.open(args.r2, "rt", encoding="utf-8", errors="replace") as f, \
            gzip.open(args.out, "wt", encoding="utf-8") as out:

        out.write("\t".join(["qname", "raw_bc", "corr_bc", "spot", "x", "y", "status"]) + "\n")

        for h, s, _, _ in iter_fastq(f):
            n_total += 1
            qname = norm_read_id(h)
            seq = s.strip()

            if len(seq) < max(x1, y1):
                out.write("\t".join([qname, "", "", "", "", "", "FAIL_SHORT"]) + "\n")
                n_short += 1
                continue

            bc_b = get_slice(seq, x0, x1)
            bc_a = get_slice(seq, y0, y1)
            raw_bc = bc_b + bc_a

            hit = best_map.get(raw_bc)

            if hit is None:
                out.write("\t".join([qname, raw_bc, "", "", "", "", "FAIL_NOMAP"]) + "\n")
                n_nomap += 1
                continue

            corr_bc, _, amb = hit

            if amb:
                out.write("\t".join([qname, raw_bc, "", "", "", "", "FAIL_AMBIG"]) + "\n")
                n_ambig += 1
                continue

            spot = bc2spot.get(corr_bc, "")

            if not spot:
                out.write("\t".join([qname, raw_bc, corr_bc, "", "", "", "FAIL_SPOT"]) + "\n")
                n_nomap += 1
                continue

            x, y = parse_xy(spot)
            out.write("\t".join([qname, raw_bc, corr_bc, spot, x, y, "OK"]) + "\n")
            n_ok += 1

            if n_total % 1000000 == 0:
                print(f"[INFO] processed={n_total:,} ok={n_ok:,}", file=sys.stderr)

    print(
        f"[DONE] total={n_total:,} ok={n_ok:,} short={n_short:,} "
        f"no_map={n_nomap:,} ambiguous={n_ambig:,}",
        file=sys.stderr
    )


if __name__ == "__main__":
    main()