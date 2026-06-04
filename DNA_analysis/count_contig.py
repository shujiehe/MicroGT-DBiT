#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Count contig-by-spot alignments from mapped spatial DNA reads.

The input BAM is expected to be the sorted and indexed BAM from the mapping step.
Read names should contain spatial coordinates in the format "...|:_:|X_Y".

Outputs:
  expmat_contig.tsv
  countcontig.log
"""

import argparse
import logging
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict


def run_capture(cmd, log=None):
    if log:
        log.info("RUN: %s", " ".join(cmd))
    return subprocess.check_output(cmd, text=True)


def file_exists(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


_RE_MAIN = re.compile(r".+\|:_:\|([0-9]+_[0-9]+)$")
_RE_FALLBACK = re.compile(r".*?([0-9]+[_xX][0-9]+)$")

PARSE_STATS = {"fallback": 0}


def parse_spot(qname):
    """
    Parse the spatial spot key from a read name.

    Main expected format:
      ...|:_:|X_Y

    Fallback accepted format:
      ...X_Y or ...XxY
    """
    name = qname[1:] if qname.startswith("@") else qname

    match = _RE_MAIN.match(name)
    if match:
        return match.group(1)

    match = _RE_FALLBACK.match(name)
    if match:
        PARSE_STATS["fallback"] += 1
        spot_raw = match.group(1)
        return spot_raw.replace("x", "_").replace("X", "_")

    return None


def load_barcodes(path):
    """
    Load barcode coordinates from a whitelist-like table.

    Expected format:
      <B+A>\\t<x>\\t<y>

    Internal spot keys are stored as X_Y.
    """
    spots = []

    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < 3:
                continue

            x = fields[1].strip()
            y = fields[2].strip()
            spots.append(f"{x}_{y}")

    return spots


def spot_display_name(spot):
    return spot.replace("_", "x", 1) if "_" in spot else spot


def get_bam_contigs(bam, log):
    """
    Get contig names from samtools idxstats.

    The returned order follows the BAM index order.
    """
    try:
        output = run_capture(["samtools", "idxstats", bam], log=log)
    except subprocess.CalledProcessError:
        sys.exit("[ERROR] samtools idxstats failed. Please index the BAM first.")

    lines = output.strip().splitlines()
    contigs = [
        line.split("\t", 1)[0]
        for line in lines
        if line and not line.startswith("*")
    ]

    seen = set()
    ordered = []

    for contig in contigs:
        if contig and contig != "*" and contig not in seen:
            seen.add(contig)
            ordered.append(contig)

    return ordered


def primary_selfcheck(bam, log):
    """
    Report whether secondary or supplementary records are present.

    The counting step filters them with -F 2308.
    """
    def count_flag(flag):
        try:
            return int(run_capture(["samtools", "view", "-c", "-f", flag, bam]))
        except Exception:
            return -1

    secondary = count_flag("256")
    supplementary = count_flag("2048")

    log.info(
        "Primary-alignment self-check: secondary=%s, supplementary=%s",
        secondary,
        supplementary
    )

    detected = []

    if isinstance(secondary, int) and secondary > 0:
        detected.append("secondary")

    if isinstance(supplementary, int) and supplementary > 0:
        detected.append("supplementary")

    if detected:
        log.warning(
            "Detected non-primary records in BAM: %s. "
            "They will be removed during counting with -F 2308.",
            ", ".join(detected)
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Count primary alignments in a contig-by-spot matrix."
    )

    parser.add_argument(
        "--bam",
        required=True,
        help="Input sorted and indexed BAM file, e.g. align.sorted.bam"
    )

    parser.add_argument(
        "--barcodes",
        required=True,
        help="Barcode coordinate table: <B+A>\\t<x>\\t<y>"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of threads for samtools view. Default: 16"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    log_fp = os.path.join(args.outdir, "countcontig.log")
    logging.basicConfig(
        filename=log_fp,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )
    log = logging.getLogger("count_contig")

    log.info("=== Contig-by-spot counting started ===")
    log.info("BAM: %s", args.bam)
    log.info("Barcode table: %s", args.barcodes)
    log.info("Threads: %d", args.threads)

    if not file_exists(args.bam):
        sys.exit(f"[ERROR] BAM not found or empty: {args.bam}")

    spots = load_barcodes(args.barcodes)
    if not spots:
        sys.exit("[ERROR] Barcode table is empty or incorrectly formatted.")

    spot_headers = [spot_display_name(spot) for spot in spots]
    spot_set = set(spots)

    contigs = get_bam_contigs(args.bam, log)
    if not contigs:
        sys.exit("[ERROR] No contigs found in BAM idxstats output.")

    primary_selfcheck(args.bam, log)

    counts = defaultdict(int)

    log.info("Streaming primary alignments with samtools view -F 2308")

    proc = subprocess.Popen(
        ["samtools", "view", "-@", str(args.threads), "-F", "2308", args.bam],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    n_in = 0
    n_badname = 0
    n_offspot = 0
    bad_examples = Counter()

    assert proc.stdout is not None

    with proc.stdout:
        for line in proc.stdout:
            if not line or line.startswith("@"):
                continue

            n_in += 1
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 11:
                continue

            contig = fields[2]
            if contig == "*" or not contig:
                continue

            qname = fields[0]
            spot = parse_spot(qname)

            if not spot:
                n_badname += 1
                if qname in bad_examples or len(bad_examples) < 10:
                    bad_examples[qname] += 1
                continue

            if spot not in spot_set:
                n_offspot += 1
                continue

            counts[(contig, spot)] += 1

    return_code = proc.wait()
    if return_code != 0:
        sys.exit(f"[ERROR] samtools view exited with return code {return_code}")

    log.info(
        "BAM scanned: total_alignments=%d, bad_read_names=%d, spots_not_in_barcode_table=%d",
        n_in,
        n_badname,
        n_offspot
    )

    if PARSE_STATS.get("fallback", 0) > 0:
        log.warning(
            "Fallback spot parser was used for %d reads. "
            "Please confirm that read names follow the expected format: ...|:_:|X_Y",
            PARSE_STATS["fallback"]
        )

    if bad_examples:
        log.warning("Examples of read names that could not be parsed:")
        for qname, count in bad_examples.most_common(10):
            log.warning("  %s  count=%d", qname, count)

    out_tsv = os.path.join(args.outdir, "expmat_contig.tsv")

    with open(out_tsv, "w") as out:
        out.write("contig\t" + "\t".join(spot_headers) + "\n")

        for contig in contigs:
            row = [contig]

            for spot in spots:
                row.append(str(counts.get((contig, spot), 0)))

            out.write("\t".join(row) + "\n")

    total_nonzero = sum(1 for value in counts.values() if value > 0)
    total_alignments = sum(counts.values())

    log.info("Matrix written: %s", out_tsv)
    log.info("Nonzero cells: %d", total_nonzero)
    log.info("Total primary alignments counted: %d", total_alignments)
    log.info("=== Contig-by-spot counting completed ===")


if __name__ == "__main__":
    main()